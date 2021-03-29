#include <stdio.h>
#include <assert.h>
#include "vectors.h"

// test vec sub
void test_vec_sub(int dim, double* a, double* b, double* e) {
    double res[dim];
    vec_sub(dim, a, b, res);
    for(int i = 0; i < dim; i++) {
        assert(res[i] == e[i]);
    }
}

// test dot product
void test_dot_prod(int dim, double* a, double* b, double e) {
    assert(dot_prod(dim, a, b) ==  e);
}

// test orth_proj
void test_orth_proj(int dim, double* p, double* a, double* b, double e) {
    double res = orth_proj_marosca(dim, p, a, b);
    printf("%06f\n", res);
    assert(res >= e - 0.050000000001);
    assert(res <= e + 0.050000000001);
}

int main() {
    double a[3] = { 3.0, 6.0, 7.0 };
    double b[3] = { 2.0, 4.0, 4.0 };
    double c[3] = { 1.0, 2.0, 3.0 };

    double d[2] = { 3.4, 7.7 };
    double e[2] = { 9.1, 2.0 };
    double f[2] = { 7.8, 8.0 };

    test_vec_sub(3, a, b, c);

    test_dot_prod(3, a, b, 58.0);

    test_orth_proj(2, f, d, e, 5.5);


    return 0;
}

/*
void test_median(int np, int seed) {
    item_t* vec;
    
    vec = gen_points(np, seed);
    // print_vec(vec, np);
    item_t m = median(vec, np);
    item_t correct = nlogn_median(vec, np);
    assert(m.sop == correct.sop);
    free(vec);
}

void test_select_ith(int np, int seed) {
    item_t* vec;
    
    vec = gen_points(np, seed);
    print_vec(vec, np);
    item_t ith = select_ith(vec, np, (np-1)/2);
    printf("SUCCESS! Got %f\n", ith.sop);
    free(vec);
}

void test_partition(int np, int seed) {
    // generating new array, to have array values
    item_t* vec = gen_points(np, seed);

    printf("Using vec:\n");
    print_vec(vec, np);
    print_separator();
    
    // test_part_aux(np, seed, 5.0);
    // test_part_aux(np, seed, -1);
    // test_part_aux(np, seed, RANGE+1);
    for(int i = 0; i < np; i++) {
        // test_part_aux(np, seed, vec[i]);
    }
}
void test_part_aux(int np, int seed, item_t ref) {
    item_t* vec;
    int i;
    vec = gen_points(np, seed);
    i = partition(vec, np, ref);
    printf("Partitioning with %f => %d\n", ref, i);
    print_vec(vec, np);
    print_separator();
    free(vec);
}
*/
