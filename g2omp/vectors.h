#ifndef __VECTORS_H__
#define __VECTORS_H__

void vec_sub(int dim, double* a, double* b, double* c);
void vec_sum(int dim, double* a, double* b, double* c);
void vec_scalar_mul(int dim, double* a, double scalar, double* c);
void vec_copy(int dim, double* src, double* dst);
void orth_proj(int dim, double* p, double* a, double* b, double* ret);
double semi_orth_proj(int dim, double* p, double* a, double* b);
double squared_dist(int dim, double* a, double* b);

#endif
