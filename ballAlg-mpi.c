#include <time.h>
#include <mpi.h>

#include "gen_points.c"
#include "common.h"

#define DELEGATE_MASTER 1


extern int N_DIMS;
extern double** POINTS;

int main(int argc, char*argv[]){
    
    /*
    long n_points;

    double exec_time = -time(NULL);
    POINTS = get_points(argc, argv, &N_DIMS, &n_points);

    sop_t* wset = (sop_t*)malloc(sizeof(sop_t) * n_points);
    for(long i = 0; i < n_points; i++) {
        wset[i].point_idx = i;
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
    */


    // Initialize MPI
    MPI_Init(&argc, &argv);

    int total_procs = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &total_procs);

    int global_rank;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    rank = global_rank;


    long n_points = 16;
    int wset[n_points];
    // TODO: (remove) populate array for testing
    if(rank == 0) {
        for(int i = 0; i < n_points; i++) {
            wset[i] = i;
        }
    }

    // Calculate max levels
    long max_levels = 2;
    for(long aux = 1; (aux <<= 1) < n_points; max_levels++);

    long level = 0;



    MPI_Comm cur_comm = MPI_COMM_WORLD;

    while(level < max_levels && total_procs > 1 << level) {

        int world_size = 0;
        MPI_Comm_size(cur_comm, &world_size);

        // TODO: ???? if(free nodes available):

        if(rank == 0) {
            printf("[LEVEL %ld] Master %d is working on (", level, global_rank);
            for(int i = 0; i < n_points-1; i++) {
                printf("%d, ", wset[i]);
            }
            printf("%d)\n", wset[n_points - 1]);
        }

        // find furthest points
        // TODO

        // split array across my current group (p/2^level processors)
        long buf_size = (long)ceil((float)n_points/(float)world_size);
        int recv_buf[buf_size];
        // printf("[%d] size: %ld\n", rank, buf_size);
        MPI_Scatter(wset, buf_size, MPI_INT, recv_buf, buf_size, MPI_INT, 0, cur_comm);

        // We don't want to operate on data if
        // There is no data for us
        if(rank < n_points) {
            if(rank == world_size - 1 && (n_points % world_size) != 0) {
                buf_size = n_points % (buf_size);
            }
            // printf("[%d] size: %ld\n", rank, buf_size);

            // calc orth projs/medians

            /*
            printf("[LEVEL %ld] [%d] got ", level, rank);
            int i;
            for(i = 0; i < buf_size - 1; i++) {
                printf("%2d, ", recv_buf[i]);
            }
            printf("%2d\n", recv_buf[i]);
            */

            // TODO: remove
            for(int i = 0; i < buf_size; i++) {
                recv_buf[i] = recv_buf[i] * 2;
            }
        }

        MPI_Gather(recv_buf, buf_size, MPI_INT, wset, buf_size, MPI_INT, 0, cur_comm);

        // TODO: remove
        MPI_Barrier(cur_comm);

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
        printf("%2d, ", wset[i]);
    }
    printf("%2d)\n", wset[i]);

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

