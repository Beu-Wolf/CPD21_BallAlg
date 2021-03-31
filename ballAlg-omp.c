#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// #define DEBUG
#include "vectors.h"
#include "median.h"
#include "gen_points.c"
#include "sop.h"

typedef struct _node {
    int id;
    double radius;
    double* center;
    long left;
    long right;
} node_t;


// core functions
void build_tree(int n_points, sop_t* wset, long id, node_t* tree, double** centers);
void calc_orth_projs(sop_t* wset, long n_points, long a_idx, long b_idx);
void find_furthest_points(sop_t* wset, long n_points, long* a, long* b);
void dump_tree(node_t* tree, double** centers, long len);

int N_DIMS;
double** POINTS;

int main(int argc, char*argv[]){
    
    long n_points;

    double exec_time = -omp_get_wtime();
    POINTS = get_points(argc, argv, &N_DIMS, &n_points);

    sop_t* wset = (sop_t*)malloc(sizeof(sop_t) * n_points);
    

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

    for(long i = 0; i < n_points; i++) {
        wset[i].point_idx = i;
    }
    for(int i = 0; i < n_nodes; i++) {
        centers[i] = &_centers[i*N_DIMS];
    }

    build_tree(n_points, wset, 0, tree, centers);
    
    exec_time += omp_get_wtime();

    fprintf(stderr, "%.1lf\n", exec_time);
    dump_tree(tree, centers, 2*n_points-1);
    return 0;
}

void build_tree(int n_points, sop_t* wset, long id, node_t* tree, double** centers) {

    if(n_points == 1) {
        // create leaf node
        node_t *leaf = &(tree[id]);
        leaf->id = id;
        leaf->center = POINTS[wset[0].point_idx];
        leaf->radius = 0.0;
        leaf->left = -1;
        leaf->right = -1;
        return;
    }

    // find furthest points
    long a_idx, b_idx;
    find_furthest_points(wset, n_points, &a_idx, &b_idx);

    // orthogonal projection
    calc_orth_projs(wset, n_points, a_idx, b_idx);

    // partitions the array into two subsets according to median
    double mdn_sop = 0.0;
    if(n_points&1) { // odd
        sop_t mdn = select_ith(wset, n_points, n_points/2);
        mdn_sop = mdn.sop;
        orth_proj(N_DIMS, POINTS[mdn.point_idx], POINTS[a_idx], POINTS[b_idx], centers[id]);
    } else {
        // lm = lower median
        // um = upper median
        sop_t lm = select_ith(wset, n_points, (n_points-1)/2);
        sop_t um;
        int num_equals = 0; // number of times we found lm. um may be equal to lm
        um.sop =  DBL_MAX;
        for(int i = 0; i < n_points; i++) {
            if(wset[i].sop == lm.sop) num_equals++;
            else if(wset[i].sop > lm.sop && wset[i].sop < um.sop) um = wset[i];
        }

        if(num_equals > 1) {
            // if um is equal to lm, then the center is given by lm
            orth_proj(N_DIMS, POINTS[lm.point_idx], POINTS[a_idx], POINTS[b_idx], centers[id]);
            mdn_sop = lm.sop;
        } else {
            // calculate averages
            double* lm_op = (double*) malloc(sizeof(double) * N_DIMS);
            double* um_op = (double*) malloc(sizeof(double) * N_DIMS);
            if(!lm_op || !um_op) exit(-1);
            orth_proj(N_DIMS, POINTS[lm.point_idx], POINTS[a_idx], POINTS[b_idx], lm_op);
            orth_proj(N_DIMS, POINTS[um.point_idx], POINTS[a_idx], POINTS[b_idx], um_op);

            vec_sum(N_DIMS, lm_op, um_op, centers[id]);
            vec_scalar_mul(N_DIMS, centers[id], 0.5, centers[id]);
            mdn_sop = lm.sop + um.sop;
            free(lm_op);
            free(um_op);
        }
    }

    partition(wset, n_points, mdn_sop);


    double sq_radius = 0.0;
    double* center = centers[id];
    for(int i = 0; i < n_points; i++) {
        double new_rad = squared_dist(N_DIMS, center, POINTS[wset[i].point_idx]);
        if(sq_radius < new_rad) sq_radius = new_rad;
    }

    // left/right point indices (partitions)
    long n_left = n_points/2;
    long n_right = n_points - n_left;

    // create node
    node_t* node = tree + id;
    node->id = id;
    node->center = centers[id];
    node->radius = sq_radius;
    node->left = id + 1;
    node->right = id + 2*n_left;

    if (n_points == 2) {
        // if nr_points is 2, there is no need to create new tasks, the recursion is way to fast
        build_tree(n_left, wset, node->left, tree, centers);
        build_tree(n_right, wset + n_left, node->right, tree, centers);
    } else {
        // left partition
        build_tree(n_left, wset, node->left, tree, centers);

        // right partition
        build_tree(n_right, wset + n_left, node->right, tree, centers);
    }
}

typedef struct {
    double max;
    int index;
} max_struct_t;

max_struct_t max_with_index(max_struct_t r, max_struct_t n) {
    max_struct_t a;
    if (r.max > n.max) {
        a = r;
    } else {
        a = n;
    }
    return a;
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
    for(int i = 0; i < n_points; i++) {
        wset[i].sop = semi_orth_proj(N_DIMS, POINTS[wset[i].point_idx], a, b);
    }
}

void dump_tree(node_t* tree, double** centers, long len) {
    printf("%d %ld\n", N_DIMS, len);
    for(int i = 0; i < len; i++) {
        printf("%d %ld %ld %.6f",
            i, tree[i].left, tree[i].right,
            sqrt(tree[i].radius));
        for(int j = 0; j < N_DIMS; j++) {
            printf(" %.6f", tree[i].center[j]);
        }
        printf("\n");
    }
}


void print_point(double* point, int dims) {
    fprintf(stderr, "[point]: (");
    for(int i = 0; i < dims-1; i++) {
        fprintf(stderr, "%f, ", point[i]);
    }
    fprintf(stderr, "%f)\n", point[dims-1]);
}
