/*
 * performs a·b
 */
double dot_prod(int dim, double* a, double* b) {
    double res = 0;
    for(int i = 0; i < dim; i++) {
        res += a[i] * b[i];
    }
    return res;
}

/*
 * performs c = a - b
 */
void vec_sub(int dim, double* a, double* b, double* c) {
    for(int i = 0; i < dim; i++) {
        c[i] = a[i]-b[i];
    }
}

void vec_sum(int dim, double* a, double* b, double* c) {
    for(int i = 0; i < dim; i++) {
        c[i] = a[i]+b[i];
    }
}

void vec_scalar_prod(int dim, double* a, double scalar, double* c) {
    for(int i = 0; i < dim; i++) {
        c[i] = scalar*a[i];
    }
}

void vec_copy(int dim, double* src, double* dst) {
    for(int i = 0; i < dim; i++) {
        dst[i] = src[i];
    }
}

/*
 * returns the the orhogonal projection of p on line ab
 *      po = a + [(p−a)·(b−a) / (b−a)·(b−a) ]*(b−a)
 */
void orth_proj(int dim, double* p, double* a, double* b, double* ret) {
    double bma[dim];
    double pma[dim];
    vec_sub(dim, b, a, bma);
    vec_sub(dim, p, a, pma);

     // [(p−a)·(b−a) / (b−a)·(b−a) ]
    double scalar = (dot_prod(dim, pma, bma)/dot_prod(dim, bma, bma));

    // scalar * (b-a)
    vec_scalar_prod(dim, bma, scalar, bma);

    // scalar*(b-a) + a
    vec_sum(dim, a, bma, ret);
}


/*
 * returns the square of the distance between a and b
 */
double squared_dist(int dim, double* a, double* b) {
    double square_sum = 0;
    for(int i = 0; i < dim; i++) {
        double diff = a[i] - b[i];
        square_sum += diff*diff;
    }

    return square_sum;
}
