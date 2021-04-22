#include "common.h"

int N_DIMS;
double** POINTS;

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

void find_furthest_points(sop_t* wset, long n_points, long* a, long* b, char is_parallel) {
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

    #pragma omp declare reduction(test:max_struct_t:omp_out=max_with_index(omp_out, omp_in)) initializer(omp_priv={0.0, 0})
    #pragma omp taskloop shared(wset) reduction(test:priv) if(n_points >= 128 && is_parallel)
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


void calc_orth_projs(sop_t* wset, long n_points, long a_idx, long b_idx, char is_parallel) {
    double* a = POINTS[a_idx];
    double* b = POINTS[b_idx];
    #pragma omp taskloop if (n_points >= 128 && is_parallel)
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
    return;
}

void print_point(double* point, int dims) {
    printf("(");
    int i;
    for(i = 0; i < dims-1; i++) {
        printf("%.6f, ", point[i]);
    }
    printf("%.6f)", point[i]);
}
