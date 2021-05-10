#include "common.h"

int N_DIMS;
double** POINTS;

max_struct_t max_with_index(max_struct_t r, max_struct_t n) {
    max_struct_t a;
    if (r.max > n.max) {
        a = r;
    } else {
        a = n;
    }
    return a;
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
