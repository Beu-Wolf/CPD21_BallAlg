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

typedef struct {
    long size;
    long tree_id;
    long start_idx;
} aux_t; // parent node info

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
    // number of nodes in the tree
    long n_nodes = 2*n_points - 1;

    // calculate number of tree levels
    long n_levels = 2;
    for(long aux = 1; (aux <<= 1) < n_points; n_levels++);

    // TODO: localidade: METER TUDO NUMA ESTRUTURA PARA FICAR TUDO NA CACHE
    // TODO (really minor): memory is still being used after processing leaves (compared
    //       to recursive version)
    aux_t* part_sizes = (aux_t*) malloc((1 << n_levels) * sizeof(aux_t));

    // printf("[%ld POINTS] Number of levels: %ld\n", n_points, n_levels);
    for(long level = 0; level < n_levels; level++) {
        long n_level_nodes = 1 << level;
        // printf("level %ld has %ld nodes\n", level, n_level_nodes);

#pragma omp parallel for
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
                long prev_level_parent_idx = level_node_idx / 2; // parent's index relative to the previous level
                long parent_idx = (1 << (level-1)) + prev_level_parent_idx - 1;
                aux_t parent = part_sizes[parent_idx];

                if(parent.size == 1) {
                    continue;
                }

                // declared in the scope above to be used later
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

            /*
            printf("%ld, %ld: [%ld, %ld[ (%ld), id=%ld\n",
                    level, level_node_idx, 
                    cur_info.start_idx, cur_info.start_idx + cur_info.size, cur_info.size,
                    cur_info.tree_id);
                    */

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
            find_furthest_points(my_wset, my_n_points, &a_idx, &b_idx);

            // orthogonal projection
            calc_orth_projs(my_wset, my_n_points, a_idx, b_idx);

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
            node->right = cur_info.tree_id + (cur_info.size / 2) * 2; // size / 2 * 2
        }
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
