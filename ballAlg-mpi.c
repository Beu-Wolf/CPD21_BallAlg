#include <time.h>
#include <mpi.h>
#include <string.h>

#include "median-mpi.h"
#include "gen_points.c"
#include "common.h"

#define DELEGATE_MASTER 1

void mpi_find_furthest_points(long* wset, long n_points, long*a, long* b);
void mpi_calc_orth_projs(long* wset, double* orthset, long n_points, long a_idx, long b_idx);
void print_vec(double* vec, int len);

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


    double exec_time = -MPI_Wtime();

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
    long id = 0;
    node_t* tree = (node_t*)malloc(sizeof(node_t) * n_nodes);
    double* _centers = (double*)malloc(sizeof(double) * n_nodes * N_DIMS);
    double** centers = (double**)malloc(sizeof(double*) * n_nodes);
    if(!tree || !centers || !_centers) {
        printf("Allocation error\n");
        exit(4);
    }

    // Initializing tree with -1 to identify written nodes
    memset(tree, 0xff, sizeof(node_t) * n_nodes);

    for(int i = 0; i < n_nodes; i++) {
        centers[i] = &_centers[i*N_DIMS];
    }

    // Dump tree pt 1: Print tree dimensions
    if(global_rank == 0) {
        printf("%d %ld\n", N_DIMS, 2*n_points - 1);
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
            mpi_find_furthest_points(wset, n_points, ab, ab+1);
        }

        // broadcast a and b
        MPI_Bcast(ab, 2, MPI_LONG, 0, cur_comm);
        long a = ab[0];
        long b = ab[1];

        // split array across my current group (p/2^level processors)
        long buf_size = (long)ceil((float)n_points/(float)world_size);
        // printf("[%d] size: %ld\n", rank, buf_size);

        long* local_wset;
        if(rank == 0) { // Master
            // Allocating a new array to keep local wset apart from global one
            local_wset = (long*)malloc(sizeof(long) * buf_size);

        } else { // Slave
            local_wset = wset;
        }

        MPI_Scatter(wset, buf_size, MPI_LONG, local_wset, buf_size, MPI_LONG, 0, cur_comm);


        // We don't want to operate on data if
        // There is no data for us
        if(rank < n_points) {
            if(rank == world_size - 1 && (n_points % world_size) != 0) {
                buf_size = n_points % (buf_size);
            }

            /*
            if(id == 1) {
                printf("%d is calculating orth projs of", global_rank);
                for(long i = 0; i < buf_size; i++) {
                    printf(" %ld", local_wset[i]);
                }
                printf("\n");
            }
            */

            // printf("[%d] size: %ld\n", rank, buf_size);
            mpi_calc_orth_projs(local_wset, orthset, buf_size, a, b);

            // TODO: calc medians and send them to master?? NO >:(

            /*
            printf("[LEVEL %ld] [%d] got ", level, rank);
            int i;
            for(i = 0; i < buf_size - 1; i++) {
                printf("%2ld, ", local_wset[i]);
            }
            printf("%2ld\n", local_wset[i]);
            */
        }

        // in-place gather
        if(rank == 0) {
            MPI_Gather(MPI_IN_PLACE, buf_size, MPI_DOUBLE, orthset, buf_size, MPI_DOUBLE, 0, cur_comm);
        } else {
            MPI_Gather(orthset, buf_size, MPI_DOUBLE, NULL, buf_size, MPI_DOUBLE, 0, cur_comm);
        }


        if(rank == 0) {
            // printf("[LEVEL %ld] master %d is partitioning array ", level, rank);
            // print_vec(orthset, n_points);

            // calculate median
            long mdn_idx;
            if(n_points&1) {
                mdn_idx = n_points/2;
                mpi_select_ith(wset, orthset, n_points, mdn_idx);
                orth_proj(N_DIMS, POINTS[wset[mdn_idx]], POINTS[a], POINTS[b], centers[id]);

            } else {
                // lm = lower median
                // um = upper median
                long lmidx = (n_points-1)/2;
                mpi_select_ith(wset, orthset, n_points, lmidx);

                /* 
                 * Finding the lowest semi orth proj that is 
                 * greater or equal than the lowest median;
                 * We only need to search the right partition 
                 * since the select_ith algorithm partitions the
                 * array using the median
                 */
                long umidx = lmidx;
                double lowersop =  DBL_MAX;
                for(int i = lmidx + 1; i < n_points; i++) {
                    if(orthset[i] < lowersop) {
                        lowersop = orthset[i];
                        umidx = i;
                    }
                }

                // calculate averages
                double* lm_op = (double*) malloc(sizeof(double) * N_DIMS);
                double* um_op = (double*) malloc(sizeof(double) * N_DIMS);
                if(!lm_op || !um_op) exit(-1);
                orth_proj(N_DIMS, POINTS[wset[lmidx]], POINTS[a], POINTS[b], lm_op);
                orth_proj(N_DIMS, POINTS[wset[umidx]], POINTS[a], POINTS[b], um_op);

                vec_sum(N_DIMS, lm_op, um_op, centers[id]);
                vec_scalar_mul(N_DIMS, centers[id], 0.5, centers[id]);
                free(lm_op);
                free(um_op);
            }

            // we do not need to partition array since select ith does that for us
            mpi_partition(wset, orthset, n_points, centers[id][0]);
        }

        // Calculate center:
        // Step 1: Broadcast center to every peer
        MPI_Bcast(centers[id], N_DIMS, MPI_DOUBLE, 0, cur_comm);

        // Step 2: Compute maximum radius with assigned working set
        double global_rad;
        double local_rad = 0;
        for(long i = 0; i < buf_size; i++) {
            double new_rad = squared_dist(N_DIMS, centers[id], POINTS[local_wset[i]]);
            if(new_rad > local_rad) local_rad = new_rad;
        }

        // Step 3: Reduce to maximum
        MPI_Reduce(&local_rad, &global_rad, 1, MPI_DOUBLE, MPI_MAX, 0, cur_comm);


        long n_left = n_points/2;
        long n_right = n_points - n_left;

        if(rank == 0) {
            // Free master's local wset buffer
            // printf("Master found radius: %f\n", global_rad);

            node_t* node = tree + id;
            node->id = id;
            node->center = centers[id];
            node->radius = global_rad;
            node->left = id + 1;
            node->right = id + 2*n_left;


            free(local_wset);
        }

        // send new problem to new master (id + cur_world_size/2)
        int next_master = world_size / 2;
        if(rank == 0) {
            // send partition to next master
            MPI_Send(wset + n_left, n_right, MPI_LONG, next_master, DELEGATE_MASTER, cur_comm);
        }
        if(rank == next_master) {
            // receive partition from prev master
            MPI_Recv(wset, n_right, MPI_LONG, 0, DELEGATE_MASTER, cur_comm, MPI_STATUS_IGNORE);
        }

        n_points = (rank < next_master ? n_left : n_right);
        id = (rank < next_master ? id + 1: id + 2*n_left);

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

    /*
    printf("[LEVEL %ld] [%d] will alone work on data (", level, global_rank);
    int i;
    for(i = 0; i < n_points - 1; i++) {
        printf("%2ld, ", wset[i]);
    }
    printf("%2ld)\n", wset[i]);
    */

    sop_t* new_wset = (sop_t*)malloc(sizeof(sop_t) * n_points);
    for(long i = 0; i < n_points; i++) {
        new_wset[i].point_idx = wset[i];
    }
    build_tree(n_points, new_wset, id, tree, centers);

    MPI_Barrier(MPI_COMM_WORLD);
    exec_time += MPI_Wtime();

#ifndef SKIP_DUMP
    for(long i = 0; i < n_nodes; i++) {
        if(tree[i].id != -1) {
            node_t* node = tree + i;
            printf("%d %ld %ld %.6f",
                node->id, node->left, node->right, sqrt(node->radius));
            print_vec(node->center, N_DIMS);
        }
    }
#endif

    if(global_rank == 0) {
        fprintf(stderr, "%.1lf\n", exec_time);
    }

    // TODO: return value?
    MPI_Finalize();
    return 0;
}

void mpi_find_furthest_points(long* wset, long n_points, long*a, long* b) {
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
        if(sd > maximum) {  
            local_b = wset[i];
            maximum = sd;
        }
    }
    
    *a = local_a;
    *b = local_b;
}

void mpi_calc_orth_projs(long* wset, double* orthset, long n_points, long a_idx, long b_idx) {
    double* a = POINTS[a_idx];
    double* b = POINTS[b_idx];
    for(int i = 0; i < n_points; i++) {
        orthset[i] = semi_orth_proj(N_DIMS, POINTS[wset[i]], a, b);
    }
}

void print_vec(double* vec, int len) {
    int i;
    for(i = 0; i < len-1; i++) {
        printf(" %.6f", vec[i]);
    }
    printf(" %.6f\n", vec[i]);
}
