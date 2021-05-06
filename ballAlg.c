#include <time.h>

#include "gen_points.c"
#include "common.h"

int N_DIMS;
double** POINTS;

void find_furthest_points(sop_t* wset, long n_points, long* a, long* b);
void calc_orth_projs(sop_t* wset, long n_points, long a_idx, long b_idx);

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

    build_tree(n_points, wset, 0, tree, centers);
    exec_time += time(NULL);

    fprintf(stderr, "%.1lf\n", exec_time);

#ifndef SKIP_DUMP
    dump_tree(tree, centers, 2*n_points-1);
#endif
    return 0;
}


void find_furthest_points(sop_t* wset, long n_points, long* a, long* b) {
    if(n_points == 2) {
        *a = wset[0].point_idx;
        *b = wset[1].point_idx;
        return;
    }

    // find A: the most distant point from the first point in the set
    long local_a = 0;
    long local_b = 0;
    double maximum = 0.0;

    max_struct_t priv;
    priv.index = 0;
    priv.max=0.0;
    // #pragma omp declare reduction(test:max_struct_t:omp_out=max_with_index(omp_out, omp_in)) initializer(omp_priv={0.0, 0})

    // #pragma omp taskloop shared(wset) reduction(test:priv) if(n_points >= 128)
    for(int i = 1; i < n_points; i++) {
        max_struct_t t;
        t.max = squared_dist(N_DIMS, POINTS[wset[0].point_idx], POINTS[wset[i].point_idx]);
        t.index = wset[i].point_idx;    
        if (priv.max < t.max) {
            priv=t;
        }
    }

    local_a = priv.index;

    // find B: the most distant point from a
    priv.index = 0;
    priv.max=0.0;

    //#pragma omp taskloop shared(wset) reduction(test:priv) if(n_points >= 128)
    for(int i = 0; i < n_points; i++) {
        max_struct_t t;
        t.max = squared_dist(N_DIMS, POINTS[local_a], POINTS[wset[i].point_idx]);
        t.index = wset[i].point_idx;
        if(priv.max < t.max) {  
            priv = t;
        }
    }

    local_b = priv.index;
    
    *a = local_a;
    *b = local_b;
}

void calc_orth_projs(sop_t* wset, long n_points, long a_idx, long b_idx) {
    double* a = POINTS[a_idx];
    double* b = POINTS[b_idx];
    //#pragma omp taskloop if (n_points >= 128)
    for(int i = 0; i < n_points; i++) {
        wset[i].sop = semi_orth_proj(N_DIMS, POINTS[wset[i].point_idx], a, b);
    }
}
