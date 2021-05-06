#include <time.h>
#include <mpi.h>

#include "gen_points.c"
#include "common.h"

#define DELEGATE_MASTER 1

void find_furthest_points(long* wset, long n_points, long*a, long* b);
void calc_orth_projs(long* wset, double* orthset, long n_points, long a_idx, long b_idx);

extern int N_DIMS;
extern double** POINTS;

int main(int argc, char*argv[]) {

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get number of processes
    int n_procs = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

    // Initialize ranks and comm variable
    int global_rank;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    rank = global_rank;
    MPI_Comm cur_comm = MPI_COMM_WORLD;


    double exec_time = -time(NULL);

    long n_points;
    POINTS = get_points(argc, argv, &N_DIMS, &n_points);

    // wset: stores working indices
    // orthset: stores orthogonal projections
    // these arrays have matching members in the same indices
    long* wset =      (long*)malloc(sizeof(long)   * n_points);
    double* orthset = (double*)malloc(sizeof(double) * n_points);
    for(long i = 0; i < n_points; i++) {
        wset[i] = i;
    }


    // allocate tree
    // TODO: we may overflow malloc argument. Check that with teachers
    // TODO: allocate in one big chunk
    long n_nodes = 2*n_points - 1;
    node_t* tree = (node_t*)malloc(sizeof(node_t) * n_nodes);
    double* _centers = (double*)malloc(sizeof(double) * n_nodes * N_DIMS);
    double** centers = (double**)malloc(sizeof(double*) * n_nodes);
    if(!tree || !centers || !_centers) {
        printf("Allocation error\n");
        exit(4);
    }
    for(int i = 0; i < n_nodes; i++) {
        centers[i] = &_centers[i*N_DIMS];
    }

    // Dump tree pt 1: Print tree dimensions
    if(global_rank == 0) {
        printf("%d %ld\n", N_DIMS, n_points);
    }

    // Calculate max levels
    long max_levels = 2;
    for(long aux = 1; (aux <<= 1) < n_points; max_levels++);

    long level = 0;
    while(level < max_levels && n_procs > (1 << level)) {

        int world_size;
        long ab[2];
        MPI_Comm_size(cur_comm, &world_size);

        // TODO: ???? if(free nodes available):

        if(rank == 0) {
            /*
            printf("[LEVEL %ld] Master %d is working on (", level, global_rank);
            for(int i = 0; i < n_points-1; i++) {
                printf("%d, ", wset[i]);
            }
            printf("%d)\n", wset[n_points - 1]);
            */

            find_furthest_points(wset, n_points, ab, ab+1);
        }

        // broadcast a and b
        MPI_Bcast(ab, 2, MPI_LONG, 0, cur_comm);
        long a = ab[0];
        long b = ab[1];

        printf("A: %ld, B: %ld\n", a, b);

        // split array across my current group (p/2^level processors)
        long buf_size = (long)ceil((float)n_points/(float)world_size);
        // printf("[%d] size: %ld\n", rank, buf_size);
        if(rank == 0) {
            // Master
            MPI_Scatter(wset, buf_size, MPI_LONG, MPI_IN_PLACE, buf_size, MPI_LONG, 0, cur_comm);
        } else {
            // Slave
            MPI_Scatter(NULL, buf_size, MPI_LONG, wset, buf_size, MPI_LONG, 0, cur_comm);
        }

        // We don't want to operate on data if
        // There is no data for us
        if(rank < n_points) {
            if(rank == world_size - 1 && (n_points % world_size) != 0) {
                buf_size = n_points % (buf_size);
            }

            // printf("[%d] size: %ld\n", rank, buf_size);
            calc_orth_projs(wset, orthset, buf_size, a, b);

            // TODO: calc medians?

            printf("[LEVEL %ld] [%d] got ", level, rank);
            int i;
            for(i = 0; i < buf_size - 1; i++) {
                printf("%2ld, ", wset[i]);
            }
            printf("%2ld\n", wset[i]);
        }

        // in-place gather
        if(rank == 0) {
            MPI_Gather(MPI_IN_PLACE, buf_size, MPI_DOUBLE, orthset, buf_size, MPI_DOUBLE, 0, cur_comm);
        } else {
            MPI_Gather(orthset, buf_size, MPI_DOUBLE, NULL, buf_size, MPI_DOUBLE, 0, cur_comm);
        }

        if(rank == 0) {
            /*
            printf("[LEVEL %ld] master %d is partitioning array ", level, rank);
            for(int i = 0; i < n_points-1; i++) {
                printf("%d, ", wset[i]);
            }
            printf("%d)\n", wset[n_points-1]);
            */

            // calculate median of medians
            // partition array
        }

        return 0;

        // send new problem to new master (id + cur_world_size/2)
        int next_master = world_size / 2;
        long n_left = n_points/2;
        long n_right = n_points - n_left;
        if(rank == 0) {
            // send partition to next master
            MPI_Send(wset + n_left, n_right, MPI_INT, next_master, DELEGATE_MASTER, cur_comm);
            n_points = n_left;
        }
        if(rank == next_master) {
            // receive partition from prev master
            MPI_Recv(wset, n_right, MPI_INT, 0, DELEGATE_MASTER, cur_comm, MPI_STATUS_IGNORE);
            n_points = n_right;
        }

        n_points = (rank < next_master ? n_left : n_right);

        if(rank == 0 || rank == next_master) {
            /*
            printf("[LEVEL %ld] %d will work on data: (", level + 1, rank);
            for(int i = 0; i < n_points - 1; i++) {
                printf("%d, ", wset[i]);
            }
            printf("%d)\n", wset[n_points - 1]);
            */
        }

        // Split current comm in two
        int color = rank / (world_size/2);
        MPI_Comm new_comm;
        MPI_Comm_split(cur_comm, color, rank, &new_comm);
        // MPI_Comm_free(&cur_comm);
        cur_comm = new_comm;

        // Update new rank
        MPI_Comm_rank(cur_comm, &rank);

        // printf("[LEVEL %ld] process %d got rank %d\n", level, global_rank, rank);

        level++;
    }


    printf("[LEVEL %ld] [%d] will alone work on data (", level, global_rank);
    int i;
    for(i = 0; i < n_points - 1; i++) {
        printf("%2ld, ", wset[i]);
    }
    printf("%2ld)\n", wset[i]);

    /*
    build_tree(n_points, wset, 0, tree, centers);
    exec_time += time(NULL);

    fprintf(stderr, "%.1lf\n", exec_time);

#ifndef SKIP_DUMP
    dump_tree(tree, centers, 2*n_points-1);
#endif
    */

    // TODO: return value?
    MPI_Finalize();
    return 0;
}

void find_furthest_points(long* wset, long n_points, long*a, long* b) {
    if(n_points == 2) {
        *a = wset[0];
        *b = wset[1];
        return;
    }

    // find A: the most distant point from the first point in the set
    long local_a = 0;
    long local_b = 0;
    double maximum = 0.0;
    for(int i = 1; i < n_points; i++) {
        double sd = squared_dist(N_DIMS, POINTS[wset[0]], POINTS[wset[i]]);
        if(sd > maximum) {
            local_a = wset[i];
            maximum = sd;
        }
    }

    maximum = 0.0;
    for(int i = 0; i < n_points; i++) {
        double sd = squared_dist(N_DIMS, POINTS[local_a], POINTS[wset[i]]);
        if(sd < maximum) {  
            local_b = wset[i];
            maximum = sd;
        }
    }
    
    *a = local_a;
    *b = local_b;
}

void calc_orth_projs(long* wset, double* orthset, long n_points, long a_idx, long b_idx) {
    double* a = POINTS[a_idx];
    double* b = POINTS[b_idx];
    for(int i = 0; i < n_points; i++) {
        orthset[i]= semi_orth_proj(N_DIMS, POINTS[wset[i]], a, b);
    }
}
