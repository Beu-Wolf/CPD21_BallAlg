#include <stdio.h>
#include <float.h>
#include <assert.h>
#include <stdlib.h>

#define RANDOM(len) (random() % len)

#define SWAP(type, a, b) { \
    type tmp = a; \
    a = b; \
    b = tmp; \
}

double mpi_pick_pivot(long* wset, double* orthset, long len);
void mpi_select_ith(long* wset, double* orthset, long len, long ith);
int mpi_partition(long* wset, double* orthset, long len, double ref);
/*
item_t median(item_t* vec, int len);
int cmp_item(const void* _a, const void* _b);
item_t nlogn_median(item_t* vec, int len);
void print_vec(item_t* vec, int len);
*/

// returns idx of ith smallest element of vec (after scrambled)
void mpi_select_ith(long* wset, double* orthset, long len, long ith) {
    // printf("Looking for %dth smallest number\n", ith+1); // ith starts in 0
    if(len == 1) return;

    if(len == 2) {
        if(orthset[0] > orthset[1]) { // sort if unsorted
            SWAP(long, wset[0], wset[1]);
            SWAP(double, orthset[0], orthset[1]);
        }
        return;
    }

    int idx = 0;

    double pivot = mpi_pick_pivot(wset, orthset, len);
    idx = mpi_partition(wset, orthset, len, pivot);

    if(ith == idx) {
        return;
    } else if(ith < idx) {
        // printf("Searching in left\n");
        mpi_select_ith(wset, orthset, idx, ith);
        return;
    } else {
        // printf("Searching in right\n");
        mpi_select_ith(wset+idx, orthset+idx, len-idx, ith-idx);
        return;
    }
}

double mpi_pick_pivot(long* wset, double* orthset, long len) {
    // randomized
    return orthset[(long)(RANDOM(len))];
    /*
    // pick median of medians
    if(len < 5) return nlogn_median(vec, len);

    int n_sub_arrays = len/5;
    item_t* medians = (item_t *) malloc(n_sub_arrays*sizeof(item_t));
    // n/5 * O(1) => O(n)
    for(int i = 0; i < n_sub_arrays; i++) {
        // O(1), constant values
        medians[i] = nlogn_median(&(vec[i*5]), 5);
    }

    item_t t = select_ith(medians, n_sub_arrays, n_sub_arrays/2);
    free(medians);
    return t;
    */
}

// returns first index of second partition
int mpi_partition(long* wset, double* orthset, long len, double ref) {
    int i = -1;
    int j = len;

    while(i < j) {
        while(i < j && orthset[++i] < ref) {
            if (i >= len) break;
            // printf("accessing i %d\n", i);
        }
        while(j > i && orthset[--j] >= ref) {
            // printf("accessing j %d\n", j);
        }
        if(i >= j) break;
        // printf("swap %d %d\n", i, j);
        // does this copy the entire structure?
        SWAP(long, wset[i], wset[j]);
        SWAP(double, orthset[i], orthset[j]);
    }
    return i;
}

// returns index of pivot
int mpi_qs_partition(double* points, double* orthset, int n_dims, long low, long high) {
    double ref = orthset[high];

    // printf("partitioning with pivot %.6f\n", ref);

    double aux[n_dims];

    int i = low-1;

    for (int j = low; j <= high-1; j++) {
        if(orthset[j] < ref) {
            i++;
            memcpy(aux, &points[i * n_dims], sizeof(double)*n_dims);
            memcpy(&points[i * n_dims], &points[j * n_dims], sizeof(double)*n_dims);
            memcpy(&points[j * n_dims], aux, sizeof(double)*n_dims);
            SWAP(double, orthset[i], orthset[j]);
        }
    }

    i++;
    SWAP(double, orthset[i], orthset[high]);
    
    memcpy(aux, &points[i * n_dims], sizeof(double)*n_dims);
    memcpy(&points[i * n_dims], &points[high * n_dims], sizeof(double)*n_dims);
    memcpy(&points[high * n_dims], aux, sizeof(double)*n_dims);
   
    return i;
}

void mpi_quicksort_aux(double* points, double* orthset, long low, long high, int n_dims) {
    if(low >= high) return;

    long i = mpi_qs_partition(points, orthset, n_dims, low, high);
    mpi_quicksort_aux(points, orthset, low, i-1, n_dims);
    mpi_quicksort_aux(points, orthset, i+1, high, n_dims);
}

void mpi_quicksort(double* points, double* orthset, int n_dims, long n_points) {
    mpi_quicksort_aux(points, orthset, 0, n_points-1, n_dims);
}

/*
int cmp_item(const void* _a, const void* _b) {
    item_t a = *(item_t*)_a;
    item_t b = *(item_t*)_b;

    if(a.sop > b.sop) return  1;
    if(a.sop < b.sop) return -1;
    return 0;
}

void insertion_sort(item_t* vec, int len) {
    int i, j;
    for (i = 1; i < len; i++) {
        item_t eli = vec[i];
        j = i-1;
        while(j >= 0 && eli.sop < vec[j].sop) {
            vec[j+1] = vec[j];
            j--;
        }
        vec[j+1] = eli;
    }
}

item_t nlogn_median(item_t* vec, int len) {
    // qsort(vec, len, sizeof(item_t), cmp_item);
    insertion_sort(vec, len);
    item_t res;
    res.sop = (vec[len/2].sop + vec[(len-1)/2].sop)/2;
    return res;
}


void print_vec(item_t* vec, int len) {
    printf("[");
    int i;
    for(i = 0; i < len-1; i++) {
        printf("%f, ", vec[i].sop);
    }
    printf("%f]\n", vec[i].sop);
}
*/
