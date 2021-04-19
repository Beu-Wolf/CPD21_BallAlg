#ifndef __COMMON_H__
#define __COMMON_H__

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "vectors.h"
#include "sop.h"

typedef struct _node {
    int id;
    double radius;
    double* center;
    long left;
    long right;
} node_t;

void build_tree(int n_points, sop_t* wset, long id, node_t* tree, double** centers);
void dump_tree(node_t* tree, double** centers, long len);
void print_point(double* point, int dims);

void find_furthest_points(sop_t* wset, long n_points, long* a, long* b);
void calc_orth_projs(sop_t* wset, long n_points, long a_idx, long b_idx);

#endif
