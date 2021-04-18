#include <omp.h>
#include <math.h>
#include <stdlib.h>

#include "common.h"
#include "vectors.h"
#include "median.h"
#include "sop.h"

void build_tree_aux(int n_points, sop_t* wset, long id, node_t* tree, double** centers);
void build_tree(int n_points, sop_t* wset, long id, node_t* tree, double** centers);
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

    partition(wset, n_points, mdn_sop);


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
    #pragma omp task untied if (n_points >= 8)
    build_tree_aux(n_left, wset, node->left, tree, centers);

    #pragma omp task untied if (n_points >= 8)
    // right partition
    build_tree_aux(n_right, wset + n_left, node->right, tree, centers);
    
}
