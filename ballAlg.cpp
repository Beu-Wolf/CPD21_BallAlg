#include <omp.h>
#include <math.h>
#include <stdio.h>

#define DEBUG
#include "vectors.h"
#include "gen_points.c"

// core functions
int build_tree(double** pts, int len);
void dump_tree(int tree);

int main(int argc, char*argv[]){
    
    int n_dims;
    long n_points;

    double exec_time = -omp_get_wtime();
    double **pts = get_points(argc, argv, &n_dims, &n_points);
    int root = build_tree(pts, 1);
    exec_time += omp_get_wtime();

    fprintf(stderr, "%.1lf\n", exec_time);
    dump_tree(root);
    return 0;
}

int build_tree(double** points, int len) {

    if(len == 1) {
        // create leaf node
        return 0;
    }

    // find furthest points

    // orthogonal projection

    // compute center

    // create L and R
    // node.left  = build_tree(L)
    // node.right = build_tree(R)
    return 0;
}

// to the stdout!
void dump_tree(int tree) {}


void print_point(double* point, int dims) {
    printf("[point]: (");
    for(int i = 0; i < dims-1; i++) {
        printf("%f, ", point[i]);
    }
    printf("%f)\n", point[dims-1]);
}
