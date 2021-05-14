#include <time.h>

#include "gen_points.c"
#include "common.h"

extern int N_DIMS;
extern double** POINTS;

void find_furthest_points(sop_t* wset, long n_points, long* a, long* b);

int main(int argc, char*argv[]){
    
    long n_points;

    double exec_time = -time(NULL);
    char useless;
    POINTS = get_points(argc, argv, 0, 1, &N_DIMS, &n_points, &useless);

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

    build_tree(n_points, wset, 0, tree, centers);
    exec_time += time(NULL);

    fprintf(stderr, "%.1lf\n", exec_time);

#ifndef SKIP_DUMP
    dump_tree(tree, centers, 2*n_points-1);
#endif
    return 0;
}


