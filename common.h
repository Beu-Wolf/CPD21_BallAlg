#ifndef __COMMON_H__
#define __COMMON_H__

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "vectors.h"
#include "sop.h"

#define MIN(x, y) (x > y ? y : x)
#define MAX(x, y) (x > y ? x : y)

typedef struct _node {
    int id;
    double radius;
    double* center;
    long left;
    long right;
} node_t;

typedef struct {
    double max;
    int index;
} max_struct_t;

void build_tree(int n_points, sop_t* wset, long id, node_t* tree, double** centers);
void dump_tree(node_t* tree, double** centers, long len);
void print_point(double* point, int dims);

#endif
