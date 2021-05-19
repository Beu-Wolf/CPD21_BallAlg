#include <stdio.h>
#include <float.h>
#include <assert.h>
#include <stdlib.h>

#include "sop.h"

typedef sop_t item_t;

item_t pick_pivot(item_t* vec, int len);
item_t nlogn_median(item_t* vec, int len);
int cmp_item(const void* _a, const void* _b);
item_t median(item_t* vec, int len);
item_t select_ith(item_t* vec, int len, int i);
int partition(item_t* vec, int len, double ref, int* different_vals);


#define RANDOM(len) (random() % len)

// void print_vec(item_t* vec, int len);

// returns ith smallest element of vec
item_t select_ith(item_t* vec, int len, int ith) {
    // printf("Looking for %dth smallest number\n", ith+1); // ith starts in 0
    if(len == 1) return vec[0];

    if(len == 2) {
        if(vec[0].sop > vec[1].sop) { // sort if unsorted
            item_t t = vec[0];
            vec[0] = vec[1];
            vec[1] = t;
        }
        return vec[ith]; // return ith
    }

    int idx = 0;

    int different_vals = 0;
    while(idx == 0 || idx == len) { // wait to find a good partition
        item_t pivot = pick_pivot(vec, len);

        // printf("random idx: %d (%f)\n", r, vec[r]);
        idx = partition(vec, len, pivot.sop, &different_vals);
        // print_vec(vec, len);
        if (!different_vals) break; // If we see all sops are the same, we can just return the first
        if((idx == 0 || idx == len) && different_vals) continue;

        // printf("partitioned at %d\n", idx);

        // printf("left:  ");
        // print_vec(vec, idx);
        // printf("right: ");
        // print_vec(vec+idx, len-idx);
    }

    if(ith < idx) {
        // printf("Searching in left\n");
        return select_ith(vec, idx, ith);
    } else if (!different_vals){
        return vec[0];
    } else {
        // printf("Searching in right\n");
        return select_ith(vec+idx, len-idx, ith-idx);
    }
}

item_t pick_pivot(item_t* vec, int len) {
    // randomized
    return vec[(int)(RANDOM(len))];

    // pick median of medians
    /*if(len < 5) return nlogn_median(vec, len);

    int n_sub_arrays = len/5;
    item_t* medians = (item_t *) malloc(n_sub_arrays*sizeof(item_t));
    // n/5 * O(1) => O(n)
    for(int i = 0; i < n_sub_arrays; i++) {
        // O(1), constant values
        medians[i] = nlogn_median(&(vec[i*5]), 5);
    }

    item_t t = select_ith(medians, n_sub_arrays, n_sub_arrays/2);
    free(medians);
    return t;*/
}

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

// returns first index of second partition
int partition(item_t* vec, int len, double ref, int* different_vals) {
    int i = -1;
    int j = len;

    double sop = vec[0].sop;
    if (sop != ref) *different_vals = 1;

    while(i < j) {
        while(i < j && vec[++i].sop < ref) {
            if (vec[i].sop != sop && different_vals == 0) *different_vals = 1;
            if (i >= len) break;
            // printf("accessing i %d\n", i);
        }
        while(j > i && vec[--j].sop >= ref) {
            if (vec[j].sop != sop && different_vals == 0) *different_vals = 1;
            // printf("accessing j %d\n", j);
        }
        if (vec[i].sop != sop || vec[j].sop != sop || vec[i].sop != vec[j].sop) *different_vals = 1;
        if(i >= j) break;
        // printf("swap %d %d\n", i, j);
        // does this copy the entire structure?
        item_t t = vec[i];
        vec[i] = vec[j];
        vec[j] = t;
    }

    return i;
}

/*
void print_vec(item_t* vec, int len) {
    printf("[");
    int i;
    for(i = 0; i < len-1; i++) {
        printf("%f, ", vec[i].sop);
    }
    printf("%f]\n", vec[i].sop);
}
*/

