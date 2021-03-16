#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define DEBUG
#include "vectors.h"
#include "gen_points.c"

typedef struct _node {
    int id;
    double radius;
    double* center;
    long point_idx;
    int left;
    int right;
} node_t;


typedef struct {
    long point_idx; // points to POINTS
    double* orth_proj;
} proj_point_t;


// core functions
void build_tree(int n_points, long* point_indices, int id, node_t* tree, double* centers);
void find_furthest_points(long* point_indices, long n_points, long* a, long* b);
void dump_tree(node_t* tree);

int cmp_orth(const void* _a, const void* _b);

int N_DIMS;
double** POINTS;

int main(int argc, char*argv[]){
    
    long n_points;

    double exec_time = -omp_get_wtime();
    POINTS = get_points(argc, argv, &N_DIMS, &n_points);

    long point_indices[n_points];
    for(long i = 0; i < n_points; i++) {
        point_indices[i] = i;
    }

    // allocate tree
    // TODO: we may overflow malloc argumtn. Check that with teachers
    node_t* tree = (node_t*)malloc(sizeof(node_t) * 2*n_points-1);
    double* centers = (double*)malloc(sizeof(double) * (2*n_points - 1) * N_DIMS);
    if(!tree) {
        printf("Allocation error\n");
        exit(4);
    }

    build_tree(n_points, point_indices, 0, tree, centers);
    exec_time += omp_get_wtime();

    fprintf(stderr, "%.1lf\n", exec_time);
    dump_tree(tree);
    return 0;
}

void build_tree(int n_points, long* point_indices, int id, node_t* tree, double* centers) {

    printf("I'm in...\n");

    if(n_points == 1) {
        // create leaf node
        node_t *leaf = &(tree[id]);
        leaf->id = id;
        leaf->point_idx = point_indices[0];
        leaf->radius = 0.0;
        leaf->left = -1;
        leaf->right = -1;
        return;
    }

    // find furthest points
    long a_idx, b_idx;
    find_furthest_points(point_indices, n_points, &a_idx, &b_idx);

    // orthogonal projection
    double orth_projs[n_points][N_DIMS];
    proj_point_t aux[n_points];
    for(long i = 0; i < n_points; i++) {
        // TODO: check if it's cheaper to orth proj
        if(point_indices[i] == a_idx || point_indices[i] == b_idx) {
            vec_copy(N_DIMS, POINTS[point_indices[i]], orth_projs[i]);
        } else {
            orth_proj(N_DIMS, POINTS[point_indices[i]], POINTS[a_idx], POINTS[b_idx], orth_projs[i]);
        }

        aux[i].point_idx = point_indices[i];
        aux[i].orth_proj = orth_projs[i];
    }

    // compute center
    // TODO: find median in O(n): https://en.wikipedia.org/wiki/Selection_algorithm
    qsort(aux, n_points, sizeof(proj_point_t), cmp_orth);

    // calculate center
    double* center = centers + (id * N_DIMS);
    double* sm = aux[n_points/2 - 1].orth_proj;
    double* bm = aux[n_points/2].orth_proj;
    if(n_points % 2 == 0) {
        for(int i = 0; i < N_DIMS; i++) {
            center[i] = (sm[i] + bm[i]) / 2.0;
        }
    } else {
        vec_copy(N_DIMS, aux[n_points / 2].orth_proj, center);
    }

    double radius = 0.0;
    // TODO: calculate radius (from center)


    // create node
    node_t* node = &(tree[id]);
    node->id = id;
    node->radius = radius;
    node->center = center;
    node->point_idx = -1;
    node->left = id + 1;
    node->right = id + n_points/2;
    
    // left/right point indices (partitions)
    long n_left = n_points/2;
    long n_right = n_points - n_left;
    long lpi[n_left];
    long rpi[n_right];

    int i = 0;
    int j = 0;
    while(i < n_points/2) {
        lpi[j++] = aux[i].point_idx;
        i++;
    }
    j = 0;
    while(i < n_points) {
        rpi[j++] = aux[i].point_idx;
        i++;
    }

    // left partition
    build_tree(n_left, lpi, node->left, tree, centers);

    // right partition
    build_tree(n_right, rpi, node->right, tree, centers);
}

void find_furthest_points(long* point_indices, long n_points, long* a, long* b) {
    double max = 0.0;
    if(n_points == 2) {
        *a = 0;
        *b = 1;
    }

    // TODO: ...
}

int cmp_orth(const void* _a, const void* _b) {
    proj_point_t a = *(proj_point_t*)_a;
    proj_point_t b = *(proj_point_t*)_b;
    if(a.orth_proj[0] > b.orth_proj[0]) return 1;
    if(a.orth_proj[0] < b.orth_proj[0]) return -1;
    return 0;
}

void dump_tree(node_t* tree) {}


void print_point(double* point, int dims) {
    printf("[point]: (");
    for(int i = 0; i < dims-1; i++) {
        printf("%f, ", point[i]);
    }
    printf("%f)\n", point[dims-1]);
}
