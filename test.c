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
