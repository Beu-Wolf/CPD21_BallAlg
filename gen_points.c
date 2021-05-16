#include <stdio.h>
#include <stdlib.h>

#define RANGE 10
#define _1GB 134217728

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n) - 1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n) - BLOCK_LOW(id,p,n) + 1)
#define BLOCK_OWNER(index,p,n) (((p)*((index)+1)-1)/(n)) 


extern void print_point(double *, int);

double **create_array_pts(int n_dims, long np)
{
    double *_p_arr;
    double **p_arr;

    _p_arr = (double *) malloc(n_dims * np * sizeof(double));
    p_arr = (double **) malloc(np * sizeof(double *));
    if((_p_arr == NULL) || (p_arr == NULL)){
        printf("Error allocating array of points, exiting.\n");
        exit(4);
    }

    for(long i = 0; i < np; i++)
        p_arr[i] = &_p_arr[i * n_dims];

    return p_arr;
}


double **get_points(int argc, char *argv[], int world_rank, int world_size, int *n_dims, long *np, char* is_big)
{
    double **pt_arr;
    unsigned seed;
    long i;
    int j;

    long requested_np;

    if(argc != 4){
        printf("Usage: %s <n_dims> <n_points> <seed>\n", argv[0]);
        exit(1);
    }

    *n_dims = atoi(argv[1]);
    if(*n_dims < 2){
        printf("Illegal number of dimensions (%d), must be above 1.\n", *n_dims);
        exit(2);
    }

    *np = atol(argv[2]);
    if(*np < 1){
        printf("Illegal number of points (%ld), must be above 0.\n", *np);
        exit(3);
    }

    seed = atoi(argv[3]);
    srandom(seed);

    *is_big = (*np) * (*n_dims) > _1GB;
    //*is_big = 1;

    if(*is_big) {
        // Big version: processes have different points
        long blockSize = BLOCK_SIZE(world_rank, world_size, *np);
        // printf("process %d: bsize: %ld, blow: %ld\n", world_rank, blockSize, BLOCK_LOW(world_rank, world_size, *np));
        pt_arr = (double **) create_array_pts(*n_dims, blockSize);
        for(i = 0; i < BLOCK_LOW(world_rank, world_size, *np); i++)
            for(j = 0; j < *n_dims; j++)
                // discard other processes' points
                random();

        for(i = 0; i < BLOCK_SIZE(world_rank, world_size, *np); i++)
            for(j = 0; j < *n_dims; j++)
                pt_arr[i][j] = RANGE * ((double) random()) / RAND_MAX;
    } else {
        // Small version: every process has every point
        pt_arr = (double **) create_array_pts(*n_dims, *np);
        for(i = 0; i < *np; i++)
            for(j = 0; j < *n_dims; j++)
                pt_arr[i][j] = RANGE * ((double) random()) / RAND_MAX;

    }

#ifdef DEBUG
    for(i = 0; i < *np; i++)
        print_point(pt_arr[i], *n_dims);
#endif

    return pt_arr;
}
