#include <omp.h>
#include <math.h>
#include <stdlib.h>

#include "common.h"
#include "vectors.h"
#include "median.h"
#include "sop.h"

void build_tree_aux(int n_points, sop_t* wset, long id, node_t* tree, double** centers);
void build_tree(int n_points, sop_t* wset, long id, node_t* tree, double** centers);

// Used to solve problems bigger than the capacity of a single processor
void build_tree_BIGGUS(int n_points, sop_t* wset, long id, node_t* tree, double** centers, long write_idx);
void build_tree_aux_BIGGUS(int n_points, sop_t* wset, long id, node_t* tree, double** centers, long write_idx);

void calc_orth_projs(sop_t* wset, long n_points, long a_idx, long b_idx);
void find_furthest_points(sop_t* wset, long n_points, long* a, long* b);

extern int N_DIMS;
extern double** POINTS;


void build_tree(int n_points, sop_t* wset, long id, node_t* tree, double** centers) {
#pragma omp parallel
    {
        #pragma omp master
        build_tree_aux(n_points, wset, id, tree, centers);
    }
}
void build_tree_aux(int n_points, sop_t* wset, long id, node_t* tree, double** centers) {

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

    int trash;

    partition(wset, n_points, mdn_sop, &trash);


    double sq_radius = 0.0;
    double* center = centers[id];
    //#pragma omp taskloop shared(wset, center) reduction(max:sq_radius) if(n_points >= 128)
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

    
    // left partition
    #pragma omp task untied if (n_points >= 128)
    build_tree_aux(n_left, wset, node->left, tree, centers);

    // right partition
    #pragma omp task untied if (n_points >= 128)
    build_tree_aux(n_right, wset + n_left, node->right, tree, centers);
}

void build_tree_BIGGUS(int n_points, sop_t* wset, long id, node_t* tree, double** centers, long write_idx) {
#pragma omp parallel
    {
        #pragma omp master
        build_tree_aux_BIGGUS(n_points, wset, id, tree, centers, write_idx);
    }
}
void build_tree_aux_BIGGUS(int n_points, sop_t* wset, long id, node_t* tree, double** centers, long write_idx) {

    if(n_points == 1) {
        // create leaf node
        node_t *leaf = &(tree[write_idx]);
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
        orth_proj(N_DIMS, POINTS[mdn.point_idx], POINTS[a_idx], POINTS[b_idx], centers[write_idx]);
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
            orth_proj(N_DIMS, POINTS[lm.point_idx], POINTS[a_idx], POINTS[b_idx], centers[write_idx]);
            mdn_sop = lm.sop;
        } else {
            // calculate averages
            double* lm_op = (double*) malloc(sizeof(double) * N_DIMS);
            double* um_op = (double*) malloc(sizeof(double) * N_DIMS);
            if(!lm_op || !um_op) exit(-1);
            orth_proj(N_DIMS, POINTS[lm.point_idx], POINTS[a_idx], POINTS[b_idx], lm_op);
            orth_proj(N_DIMS, POINTS[um.point_idx], POINTS[a_idx], POINTS[b_idx], um_op);

            vec_sum(N_DIMS, lm_op, um_op, centers[write_idx]);
            vec_scalar_mul(N_DIMS, centers[write_idx], 0.5, centers[write_idx]);
            mdn_sop = lm.sop + um.sop;
            free(lm_op);
            free(um_op);
        }
    }

    int trash;
    partition(wset, n_points, mdn_sop, &trash);


    double sq_radius = 0.0;
    double* center = centers[write_idx];
    //#pragma omp taskloop shared(wset, center) reduction(max:sq_radius) if(n_points >= 128)
    for(int i = 0; i < n_points; i++) {
        double new_rad = squared_dist(N_DIMS, center, POINTS[wset[i].point_idx]);
        if(sq_radius < new_rad) sq_radius = new_rad;
    }

    // left/right point indices (partitions)
    long n_left = n_points/2;
    long n_right = n_points - n_left;

    // create node
    node_t* node = tree + write_idx;
    node->id = id;
    node->center = centers[write_idx];
    node->radius = sq_radius;
    node->left = id + 1;
    node->right = id + 2*n_left;

    long left_widx = write_idx + 1;
    long right_widx = write_idx + 2*n_left;
    
    // left partition
    #pragma omp task untied if (n_points >= 128)
    build_tree_aux_BIGGUS(n_left, wset, node->left, tree, centers, left_widx);

    // right partition
    #pragma omp task untied if (n_points >= 128)
    build_tree_aux_BIGGUS(n_right, wset + n_left, node->right, tree, centers, right_widx);
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
