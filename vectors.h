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


/*
 * returns the firt coordinate of 
 * the orhogonal projection of p on line ab
 *      po = a + [(p−a)·(b−a) / (b−a)·(b−a) ]*(b−a)
 */
double orth_proj_marosca(int dim, double* p, double* a, double* b) {
    double bma[dim];
    double pma[dim];
    vec_sub(dim, b, a, bma);
    vec_sub(dim, p, a, pma);

    return a[0] + (b[0] - a[0])*(dot_prod(dim, pma, bma)/dot_prod(dim, bma, bma));
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
