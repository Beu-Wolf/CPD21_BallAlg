#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "vectors.h"
#include "median.h"
#include "sop.h"


typedef struct {
    long size;
    long tree_id;
    long start_idx;
} aux_t; // parent node info

extern int N_DIMS;
extern double** POINTS;

void build_tree(int n_points, sop_t* wset, long id, node_t* tree, double** centers);
void build_tree_aux(int n_points, sop_t* wset, long id, node_t* tree, double** centers);
void calc_orth_projs(sop_t* wset, long n_points, long a_idx, long b_idx, char is_parallel);
void find_furthest_points(sop_t* wset, long n_points, long* a, long* b, char is_parallel);

void build_tree(int n_points, sop_t* wset, long id, node_t* tree, double** centers) {
    // number of nodes in the tree
    long n_nodes = 2*n_points - 1;

    // calculate number of tree levels
    long n_levels = 2;
    for(long aux = 1; (aux <<= 1) < n_points; n_levels++);


    // TODO: localidade: METER TUDO NUMA ESTRUTURA PARA FICAR TUDO NA CACHE
    // TODO (really minor): memory is still being used after processing leaves (compared
    //       to recursive version)
    aux_t* part_sizes = (aux_t*) malloc((1 << n_levels) * sizeof(aux_t));

    // Create basic nodes for every thread
    // (currently assuming a power of 2 number of threads)
    #pragma omp parallel
    {
        // These variables are declared here to be used outside master scope
        long level;
        long n_thread_levels;
        long n_threads;
        long n_level_nodes;

        #pragma omp single copyprivate(level, n_thread_levels, n_threads, n_level_nodes)
        {
            // TODO: initialize this value in for loop
            n_thread_levels = 2;

#ifndef SERIAL
            n_threads = omp_get_num_threads();
#else
            n_threads = 1;
#endif
            if(n_threads == 1) n_thread_levels = 1;
            else for(long aux = 1; (aux<<=1) < n_threads; n_thread_levels++);

            for(level = 0; level < n_thread_levels - 1; level++) {
                n_level_nodes = 1 << level;
                for(long level_node_idx = 0; level_node_idx < n_level_nodes; level_node_idx++) {
                    // calculate cenas
                    long cur_node_idx = (1 << level) + level_node_idx - 1;
                    aux_t cur_info;
                    if(level == 0) {
                        cur_info.size = n_points;
                        cur_info.tree_id = 0;
                        cur_info.start_idx = 0;
                        // x = n_points / 2;
                    } else {
                        // parent's index relative to the previous level
                        long prev_level_parent_idx = level_node_idx / 2; 
                        long parent_idx = (1 << (level-1)) + prev_level_parent_idx - 1;

                        aux_t parent = part_sizes[parent_idx];

                        if(parent.size == 1) {
                            cur_info.size = 1;
                            continue;
                        }

                        long x = parent.size / 2;
                        long y = parent.size % 2;
                        if(level_node_idx % 2 == 0) { // left child
                            cur_info.size = x;
                            cur_info.tree_id = parent.tree_id + 1;
                            cur_info.start_idx = parent.start_idx;

                        } else { // right child
                            cur_info.size = x + y;
                            cur_info.tree_id = parent.tree_id + 2*x; // 2*x - 1 + 1
                            cur_info.start_idx = parent.start_idx + x;
                        }
                    }
                    part_sizes[cur_node_idx] = cur_info;

                    sop_t* my_wset = (wset + cur_info.start_idx);
                    long my_n_points = cur_info.size;

                    if(cur_info.size == 1) {
                        // create leaf node
                        node_t *leaf = &(tree[cur_info.tree_id]);
                        leaf->id = cur_info.tree_id;
                        leaf->center = POINTS[my_wset[0].point_idx];
                        leaf->radius = 0.0;
                        leaf->left = -1;
                        leaf->right = -1;
                        continue;
                    }

                    // find furthest points
                    long a_idx, b_idx;
                    find_furthest_points(my_wset, my_n_points, &a_idx, &b_idx, 1);

                    // orthogonal projection
                    calc_orth_projs(my_wset, my_n_points, a_idx, b_idx, 1);

                    // partitions the array into two subsets according to median
                    double mdn_sop = 0.0;
                    if(my_n_points&1) { // odd
                        sop_t mdn = select_ith(my_wset, my_n_points, my_n_points/2);
                        mdn_sop = mdn.sop;
                        orth_proj(N_DIMS, POINTS[mdn.point_idx], POINTS[a_idx], POINTS[b_idx], centers[cur_info.tree_id]);
                    } else {
                        // lm = lower median
                        // um = upper median
                        sop_t lm = select_ith(my_wset, my_n_points, (my_n_points-1)/2);
                        sop_t um;
                        int num_equals = 0; // number of times we found lm. um may be equal to lm
                        um.sop =  DBL_MAX;

                        for(int i = 0; i < my_n_points; i++) {
                            if(my_wset[i].sop == lm.sop) {
                                num_equals++;

                            } else if(my_wset[i].sop > lm.sop && my_wset[i].sop < um.sop) {
                                um = my_wset[i];
                            }
                        }

                        if(num_equals > 1) {
                            // if um is equal to lm, then the center is given by lm
                            orth_proj(N_DIMS, POINTS[lm.point_idx], POINTS[a_idx], POINTS[b_idx], centers[cur_info.tree_id]);
                            mdn_sop = lm.sop;
                        } else {
                            // calculate averages
                            double* lm_op = (double*) malloc(sizeof(double) * N_DIMS);
                            double* um_op = (double*) malloc(sizeof(double) * N_DIMS);
                            if(!lm_op || !um_op) exit(-1);
                            orth_proj(N_DIMS, POINTS[lm.point_idx], POINTS[a_idx], POINTS[b_idx], lm_op);
                            orth_proj(N_DIMS, POINTS[um.point_idx], POINTS[a_idx], POINTS[b_idx], um_op);

                            vec_sum(N_DIMS, lm_op, um_op, centers[cur_info.tree_id]);
                            vec_scalar_mul(N_DIMS, centers[cur_info.tree_id], 0.5, centers[cur_info.tree_id]);
                            mdn_sop = lm.sop + um.sop;
                            free(lm_op);
                            free(um_op);
                        }
                    }

                    partition(my_wset, my_n_points, mdn_sop);

                    double sq_radius = 0.0;
                    double* center = centers[cur_info.tree_id];
                    for(int i = 0; i < my_n_points; i++) {
                        double new_rad = squared_dist(N_DIMS, center, POINTS[my_wset[i].point_idx]);
                        if(sq_radius < new_rad) sq_radius = new_rad;
                    }

                    // create node
                    node_t* node = tree + cur_info.tree_id;
                    node->id = cur_info.tree_id;
                    node->center = centers[cur_info.tree_id];
                    node->radius = sq_radius;
                    node->left = cur_info.tree_id + 1;
                    // (size / 2) * 2 will make the LSB 0, so we just & with 0b1111...110
                    node->right = cur_info.tree_id + (cur_info.size & ~1);

                } // END OF INNER FOR
            }

            // update number of nodes in next level
            n_level_nodes = 1 << level;
        } // END OF SINGLE (implicit barrier wait for every thread)
        #pragma omp barrier

        // now we reached a level where we have 1 thread per node
        #pragma omp for
        for(long level_node_idx = 0; level_node_idx < n_level_nodes; level_node_idx++) {

            long cur_node_idx = (1 << level) + level_node_idx - 1;
            aux_t cur_info;
            if(level == 0) {
                cur_info.size = n_points;
                cur_info.tree_id = 0;
                cur_info.start_idx = 0;
            } else {
                // parent's index relative to the previous level
                long prev_level_parent_idx = level_node_idx / 2; 
                long parent_idx = (1 << (level-1)) + prev_level_parent_idx - 1;

                aux_t parent = part_sizes[parent_idx];

                if(parent.size == 1) {
                    cur_info.size = 1;
                    continue;
                }

                long x = parent.size / 2;
                long y = parent.size % 2;
                if(level_node_idx % 2 == 0) { // left child
                    cur_info.size = x;
                    cur_info.tree_id = parent.tree_id + 1;
                    cur_info.start_idx = parent.start_idx;

                } else { // right child
                    cur_info.size = x + y;
                    cur_info.tree_id = parent.tree_id + 2*x; // 2*x - 1 + 1
                    cur_info.start_idx = parent.start_idx + x;
                }
            }

            sop_t* my_wset = (wset + cur_info.start_idx);
            long my_n_points = cur_info.size;

            // each thread will serially operate on this subtree
            build_tree_aux(my_n_points, my_wset, cur_info.tree_id, tree, centers);
        }
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
    find_furthest_points(wset, n_points, &a_idx, &b_idx, 0);

    // orthogonal projection
    calc_orth_projs(wset, n_points, a_idx, b_idx, 0);

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

    build_tree_aux(n_left, wset, node->left, tree, centers);

    build_tree_aux(n_right, wset + n_left, node->right, tree, centers);
    
}
