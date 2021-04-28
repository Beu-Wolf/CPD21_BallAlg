#include <time.h>
#include <mpi.h>

#include "gen_points.c"
#include "common.h"


extern int N_DIMS;
extern double** POINTS;

int main(int argc, char*argv[]){
    
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


//    if master of current level:
//        while not finished:
//            if(free nodes available):
//                find furthest points
//                split array across my current group (p/2^level processors)
//                wait for orth projs/medians
//                calculate median of medians
//                partition array
//                send new problem to new master (id + cur_world_size/2)
//            else:
//                receive (probe?) array -> calculate orht projs
//                calculate median
//                send orth projs (median points at the end of the array)
//                new_group_id = 2*my_id / world_size (my_id / (world_size/2))
//                am_i_master = (my_id % world_size/2) == 0

    /*
    build_tree(n_points, wset, 0, tree, centers);
    exec_time += time(NULL);

    fprintf(stderr, "%.1lf\n", exec_time);

#ifndef SKIP_DUMP
    dump_tree(tree, centers, 2*n_points-1);
#endif
    */
    return 0;
}

