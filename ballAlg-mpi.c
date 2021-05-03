#include <time.h>
#include <mpi.h>

#include "gen_points.c"
#include "common.h"


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


    long n_points = 16;
    int wset[n_points];

    // Calculate max levels
    long max_levels = 2;
    for(long aux = 1; (aux <<= 1) < n_points; max_levels++);

    long level = 0;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    while(level < max_levels && total_procs > 1 << level) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int world_size = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        // TODO: ???? if(free nodes available):

        if(rank == 0 && level == 0) {
            for(int i = 0; i < n_points; i++) {
                wset[i] = i;
            }

            printf("[LEVEL %ld] Master %d is working on (", level, rank);
            for(int i = 0; i < n_points-1; i++) {
                printf("%d, ", wset[i]);
            }
            printf("%d)\n", wset[n_points - 1]);
        }

        // TODO: remove
        MPI_Barrier(MPI_COMM_WORLD);


        // find furthest points
        // TODO

        // split array across my current group (p/2^level processors)
        long buf_size = (long)ceil((float)n_points/(float)world_size);
        int recv_buf[buf_size];
        // printf("[%d] size: %ld\n", rank, buf_size);
        MPI_Scatter(wset, buf_size, MPI_INT, recv_buf, buf_size, MPI_INT, 0, MPI_COMM_WORLD);

        // We don't want to operate on data if
        // There is no data for us
        if(rank < n_points) {
            if(rank == world_size - 1 && (n_points % buf_size) != 0) {
                buf_size = n_points % (buf_size);
            }
            // printf("[%d] size: %ld\n", rank, buf_size);

            // calc orth projs/medians

            printf("[LEVEL %ld] [%d] got ", level, rank);
            int i;
            for(i = 0; i < buf_size - 1; i++) {
                printf("%2d, ", recv_buf[i]);
            }
            printf("%2d\n", recv_buf[i]);

            // TODO: remove
            for(int i = 0, j = buf_size - 1; i < j; i++, j--) {
                int tmp = recv_buf[i];
                recv_buf[i] = recv_buf[j];
                recv_buf[j] = tmp;
            }

        }


        MPI_Gather(recv_buf, buf_size, MPI_INT, wset, buf_size, MPI_INT, 0, MPI_COMM_WORLD);

        // TODO: remove
        MPI_Barrier(MPI_COMM_WORLD);

        if(rank == 0) {
            printf("[LEVEL %ld] master %d got ", level, rank);
            for(int i = 0; i < n_points; i++) {
                printf("%d, ", wset[i]);
            }
            printf("\n");

            // calculate median of medians
            // partition array
        }


        int color = rank / (world_size/2);
        MPI_Comm row_comm;
        MPI_Comm_split(MPI_COMM_WORLD, color, rank, &row_comm);

        int row_rank, row_size;
        MPI_Comm_rank(row_comm, &row_rank);
        MPI_Comm_size(row_comm, &row_size);

        MPI_Comm_free(&row_comm);

        // send new problem to new master (id + cur_world_size/2)


        // new_group_id = 2*my_id / world_size (my_id / (world_size/2))
        // am_i_master = (my_id % world_size/2) == 0

        level++;
    }

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

