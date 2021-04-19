#include "vectors.h"

#include <stdlib.h>

/*
 * performs a·b
 */
double dot_prod(int dim, double* a, double* b) {
    double res = 0;
    #pragma omp simd reduction(+:res)
    for(int i = 0; i < dim; i++) {
        res += a[i] * b[i];
    }
    return res;
}

/*
 * performs c = a - b
 */
void vec_sub(int dim, double* a, double* b, double* c) {
    #pragma omp simd
    for(int i = 0; i < dim; i++) {
        c[i] = a[i]-b[i];
    }
}

void vec_sum(int dim, double* a, double* b, double* c) {
    #pragma omp simd
    for(int i = 0; i < dim; i++) {
        c[i] = a[i]+b[i];
    }
}

void vec_scalar_mul(int dim, double* a, double scalar, double* c) {
    #pragma omp simd
    for(int i = 0; i < dim; i++) {
        c[i] = scalar*a[i];
    }
}

void vec_copy(int dim, double* src, double* dst) {
    #pragma omp simd
    for(int i = 0; i < dim; i++) {
        dst[i] = src[i];
    }
}

/*
 * returns the the orhogonal projection of p on line ab
 *      po = a + [(p−a)·(b−a) / (b−a)·(b−a) ]*(b−a)
 */
void orth_proj(int dim, double* p, double* a, double* b, double* ret) {
    double* bma = (double*) malloc(sizeof(double) * dim);
    double* pma = (double*) malloc(sizeof(double) * dim);
    if(!bma || !pma) exit(-1);
    vec_sub(dim, b, a, bma);
    vec_sub(dim, p, a, pma);

     // [(p−a)·(b−a) / (b−a)·(b−a) ]
    double scalar = (dot_prod(dim, pma, bma)/dot_prod(dim, bma, bma));

    // scalar * (b-a)
    vec_scalar_mul(dim, bma, scalar, bma);

    // scalar*(b-a) + a
    vec_sum(dim, a, bma, ret);
    free(bma);
    free(pma);
}

double semi_orth_proj(int dim, double* p, double* a, double* b) {
    double* bma = (double*) malloc(sizeof(double) * dim);
    double* pma = (double*) malloc(sizeof(double) * dim);
    if(!bma || !pma) exit(-1);
    vec_sub(dim, b, a, bma);
    vec_sub(dim, p, a, pma);

     // [(p−a)·(b−a) / (b−a)·(b−a) ]
    double scalar = (dot_prod(dim, pma, bma)/dot_prod(dim, bma, bma));

    free(bma);
    free(pma);

    return a[0] + (scalar * (b[0] - a[0]));
}

/*
 * returns the square of the distance between a and b
 */
double squared_dist(int dim, double* a, double* b) {
    double square_sum = 0.0;
    // #pragma omp simd
    for(int i = 0; i < dim; i++) {
        double diff = a[i] - b[i];
        square_sum += diff*diff;
    }

    return square_sum;
}
