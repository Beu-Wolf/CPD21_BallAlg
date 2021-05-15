#include <time.h>
#include <mpi.h>
#include <string.h>
#include <unistd.h>

#include "median-mpi.h"
#include "common.h"
#include "gen_points.c"
#include "sop.h"

#define DELEGATE_MASTER 1
#define OET_POINTS_1 2
#define OET_POINTS_2 3
#define OET_ORTHS_1 4
#define OET_ORTHS_2 5
#define EVEN_MEDIAN 6

int solve_small(int n_procs, int global_rank, double exec_time, long n_points);

void build_tree_BIGGUS(int n_points, sop_t* wset, long id, node_t* tree, double** centers, long write_idx);

void reduce_max_distance(void *in, void *inout, int *len, MPI_Datatype *dptr);
int solve_big(int n_procs, int global_rank, double exec_time, long n_points);

void mpi_find_furthest_points(long* wset, long n_points, long*a, long* b);
void mpi_calc_orth_projs(long* wset, double* orthset, long n_points, long a_idx, long b_idx);
void print_vec(double* vec, int len);

MPI_Op REDUCE_MAX_DISTANCE;

extern int N_DIMS;
extern double** POINTS;
double* SWAPPIE_SWAPPIE;
double* SWAPPIE_SWAPPIE_ORTH;

double* global_aux;

int main(int argc, char*argv[]) {

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Create custom reduce function
    MPI_Op_create(reduce_max_distance, 1, &REDUCE_MAX_DISTANCE);


    // Get number of processes
    int n_procs = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

    // Initialize ranks and comm variable
    int global_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);

    double exec_time = -MPI_Wtime();

    long n_points;
    char is_big;
    POINTS = get_points(argc, argv, global_rank, n_procs, &N_DIMS, &n_points, &is_big);

    if(is_big) {
        solve_big(n_procs, global_rank, exec_time, n_points);
    } else {
        solve_small(n_procs, global_rank, exec_time, n_points);
    }

    //free(global_aux);
    return 0;
}

void reduce_max_distance(void *in, void *inout, int *len, MPI_Datatype *dptr) {
    double sqdst1 = squared_dist(*len, in, global_aux);
    double sqdst2 = squared_dist(*len, inout, global_aux);
    if(sqdst1 > sqdst2) {
        vec_copy(*len, in, inout);
    }
    
}

int solve_big(int n_procs, int global_rank, double exec_time, long _n_points) {

    int rank = global_rank;
    MPI_Comm cur_comm = MPI_COMM_WORLD;

    long _n_nodes = 2*_n_points - 1;

    long local_n_points = BLOCK_SIZE(rank, n_procs, _n_points);
    long local_n_nodes = 2*local_n_points + n_procs;

    SWAPPIE_SWAPPIE = (double*)malloc(sizeof(double) * N_DIMS * local_n_points);
    SWAPPIE_SWAPPIE_ORTH = (double*)malloc(sizeof(double) * local_n_points);

    // orthset: stores orthogonal projections
    double* orthset = (double*)malloc(sizeof(double) * local_n_points);

    global_aux = (double*)malloc(sizeof(double) * N_DIMS);

    // allocate tree
    node_t* tree = (node_t*)malloc(sizeof(node_t) * local_n_nodes);
    double* _centers = (double*)malloc(sizeof(double) * local_n_nodes * N_DIMS);
    double** centers = (double**)malloc(sizeof(double*) * local_n_nodes);
    if(!tree || !centers || !_centers) {
        printf("Allocation error\n");
        exit(4);
    }
    for(int i = 0; i < local_n_nodes; i++) {
        centers[i] = &_centers[i*N_DIMS];
    }

    // Initializing local tree with -1 to identify written nodes
    memset(tree, 0xff, sizeof(node_t) * local_n_nodes);


    // Dump tree pt 1: Print tree dimensions
    if(global_rank == 0) {
        printf("%d %ld\n", N_DIMS, 2*_n_points - 1);
    }

    // Calculate max levels
    long max_levels = 2;
    for(long aux = 1; (aux <<= 1) < _n_points; max_levels++);

    long id = 0;
    long level = 0;
    long write_idx = 0;
    while(level < max_levels && n_procs > (1 << level)) {

        int comm_size;
        MPI_Comm_size(cur_comm, &comm_size);

        // Step 1: Find furthest points

        // 1.1: Send first point to every process
        double a[N_DIMS];
        double b[N_DIMS];

        if(rank == 0) {
            vec_copy(N_DIMS, POINTS[0], global_aux);
        }
        MPI_Bcast(global_aux, N_DIMS, MPI_DOUBLE, 0, cur_comm);

        // 1.2 Each process finds its local a
        double max_sq_dist = 0;
        long max_idx = 0;
        for(long i = 0; i < local_n_points; i++) {
            double tmp_dist = squared_dist(N_DIMS, global_aux, POINTS[i]);
            if(tmp_dist > max_sq_dist) {
                max_idx = i;
                max_sq_dist = tmp_dist;
            }
        }
        // printf("[%d] aux point: ", global_rank);
        // print_vec(global_aux, N_DIMS);

        // 1.3 Find global a

        MPI_Allreduce(POINTS[max_idx], a, N_DIMS, MPI_DOUBLE, REDUCE_MAX_DISTANCE, cur_comm);

        // 1.4 Find local b
        vec_copy(N_DIMS, a, global_aux);
        max_sq_dist = 0;
        max_idx = 0;
        for(long i = 0; i < local_n_points; i++) {
            double tmp_dist = squared_dist(N_DIMS, global_aux, POINTS[i]);
            if(tmp_dist > max_sq_dist) {
                max_idx = i;
                max_sq_dist = tmp_dist;
            }
        }

        // 1.5 Find global b
        MPI_Allreduce(POINTS[max_idx], b, N_DIMS, MPI_DOUBLE, REDUCE_MAX_DISTANCE, cur_comm);

        // Step 2: Calculate semi orth projs
        for(int i = 0; i < local_n_points; i++) {
            orthset[i] = semi_orth_proj(N_DIMS, POINTS[i], a, b);
        }

        MPI_Barrier(cur_comm);

        // Step 3: Sort global array
        // 3.1: Sort local array
        mpi_quicksort(POINTS[0], orthset, N_DIMS, local_n_points);

        // 3.2: OET
        for(int phase = 0; phase < comm_size; phase++) {
            int peer = rank;

            //              phase%2 = 0 | phase%2 = 1
            // rank%2 = 0       ++      |     --
            // rank%2 = 1       --      |     ++
            //
            if(phase % 2 == 0) { // even phase
                if(rank % 2 == 0) peer++;
                else peer--;
            } else { // odd phase
                if(rank % 2 == 0) peer--;
                else peer++;
            }

            if(peer >= 0 && peer < comm_size) {
                // 3.2.1: Send my data to my peer
                //printf("BEFORE MALOC: [%d] phase %d\n", rank, phase);
                long peer_block_size = BLOCK_SIZE(peer, comm_size, _n_points);

                double* peer_points = (double*)malloc(sizeof(double) * N_DIMS * peer_block_size);
                double* peer_orth = (double*)malloc(sizeof(double) * peer_block_size);

                
                if(rank % 2 == 0) { // these nodes send first
                    // HACK: POINTS[0] points to the beginning of the point buffer
                    MPI_Send(POINTS[0], local_n_points * N_DIMS, MPI_DOUBLE, peer, OET_POINTS_1, cur_comm);
                    MPI_Recv(peer_points, peer_block_size * N_DIMS, MPI_DOUBLE, peer, OET_POINTS_2, cur_comm, NULL);

                    MPI_Send(orthset, local_n_points, MPI_DOUBLE, peer, OET_ORTHS_1, cur_comm);
                    MPI_Recv(peer_orth, peer_block_size, MPI_DOUBLE, peer, OET_ORTHS_2, cur_comm, NULL);

                } else { // these nodes receive first
                    MPI_Recv(peer_points, peer_block_size * N_DIMS, MPI_DOUBLE, peer, OET_POINTS_1, cur_comm, NULL);
                    MPI_Send(POINTS[0], local_n_points * N_DIMS, MPI_DOUBLE, peer, OET_POINTS_2, cur_comm);

                    MPI_Recv(peer_orth, peer_block_size, MPI_DOUBLE, peer, OET_ORTHS_1, cur_comm, NULL);
                    MPI_Send(orthset, local_n_points, MPI_DOUBLE, peer, OET_ORTHS_2, cur_comm);
                }



                // 3.2.2: Merge-sort arrays
                // Copy local points to swappie swappie
                memcpy(SWAPPIE_SWAPPIE, POINTS[0], sizeof(double) * N_DIMS * local_n_points);
                if(rank < peer) {
                    long my_idx = 0, peer_idx = 0;
                    for(long k = 0; k < local_n_points; k++) {
                        if(peer_orth[peer_idx] < orthset[my_idx] && peer_idx < peer_block_size) {
                            memcpy(POINTS[k], peer_points + (peer_idx * N_DIMS), sizeof(double) * N_DIMS);
                            SWAPPIE_SWAPPIE_ORTH[k] = peer_orth[peer_idx];
                            peer_idx++;
                        } else {
                            memcpy(POINTS[k], SWAPPIE_SWAPPIE + (my_idx * N_DIMS), sizeof(double) * N_DIMS);
                            SWAPPIE_SWAPPIE_ORTH[k] = orthset[my_idx];
                            my_idx++;
                        }
                    }
                } else {
                    long my_idx = local_n_points - 1, peer_idx = peer_block_size - 1;
                    for(long k = local_n_points - 1; k >= 0; k--) {
                        if(peer_orth[peer_idx] > orthset[my_idx] && peer_idx >= 0) {
                            memcpy(POINTS[k], peer_points + (peer_idx * N_DIMS), sizeof(double) * N_DIMS);
                            SWAPPIE_SWAPPIE_ORTH[k] = peer_orth[peer_idx];
                            peer_idx--;
                        } else {
                            memcpy(POINTS[k], SWAPPIE_SWAPPIE + (my_idx * N_DIMS), sizeof(double) * N_DIMS);
                            SWAPPIE_SWAPPIE_ORTH[k] = orthset[my_idx];
                            my_idx--;
                        }
                    }
                }
                // do the swappie swappie
                double* tmp = orthset;
                orthset = SWAPPIE_SWAPPIE_ORTH;
                vec_copy(local_n_points, SWAPPIE_SWAPPIE_ORTH, tmp);
                SWAPPIE_SWAPPIE_ORTH = tmp;

                free(peer_points);
                free(peer_orth);
            }


            MPI_Barrier(cur_comm);
        }

        /*
        for(int p = 0; p < comm_size; p++) {
            if(p == rank) {
                for(long i = 0; i < local_n_points; i++) {
                    printf("[%d] %.6f | ", rank, orthset[i]);
                    print_vec(POINTS[i], N_DIMS);
                }
            }
            fflush(stdout);
            MPI_Barrier(cur_comm);
        }
        */

        // Step 4: Calculate center
        long mdnidx = (_n_points - 1) / 2;
        int median_owner = 0;
        if(_n_points % 2 == 1) {
            median_owner = BLOCK_OWNER(mdnidx, comm_size, _n_points);
            if(rank == median_owner) {
                orth_proj(N_DIMS, POINTS[0], a, b, centers[write_idx]);
            }
        } else {
            long lmid = mdnidx;
            long umid = lmid + 1;
            int lower_owner = BLOCK_OWNER(lmid, comm_size, _n_points);
            int upper_owner = BLOCK_OWNER(umid, comm_size, _n_points);
            median_owner = upper_owner;

            assert(BLOCK_OWNER(umid, comm_size, _n_points) - BLOCK_OWNER(lmid, comm_size, _n_points) == 1);
            assert(BLOCK_HIGH(lower_owner, comm_size, _n_points) == lmid);
            assert(BLOCK_LOW(upper_owner, comm_size, _n_points) == umid);


            if(rank == lower_owner) {
                // Send lower point to the dude that has the upper point
                MPI_Send(POINTS[local_n_points - 1], N_DIMS, MPI_DOUBLE, upper_owner, EVEN_MEDIAN, cur_comm);
            }
            if(rank == upper_owner) {
                double lower_point[N_DIMS];
                double lm_op[N_DIMS];
                double um_op[N_DIMS];

                MPI_Recv(lower_point, N_DIMS, MPI_DOUBLE, lower_owner, EVEN_MEDIAN, cur_comm, NULL);

                orth_proj(N_DIMS, lower_point, a, b, lm_op);
                orth_proj(N_DIMS, POINTS[0], a, b, um_op);
                vec_sum(N_DIMS, lm_op, um_op, centers[write_idx]);
                vec_scalar_mul(N_DIMS, centers[write_idx], 0.5, centers[write_idx]);
            }
        }

        // Step 5: Calculate greatest radius
        // Step 5.1: Broadcast center
        double center[N_DIMS];
        if(rank == median_owner) {
            memcpy(center, centers[write_idx], sizeof(double) * N_DIMS);
        }
        MPI_Bcast(center, N_DIMS, MPI_DOUBLE, median_owner, cur_comm);

        // Step 5.2: Calculate local radius
        double global_rad;
        double local_rad = 0;
        for(long i = 0; i < local_n_points; i++) {
            double new_rad = squared_dist(N_DIMS, center, POINTS[i]);
            if(new_rad > local_rad) local_rad = new_rad;
        }

        // Step 5.3: Reduce to maximum
        MPI_Reduce(&local_rad, &global_rad, 1, MPI_DOUBLE, MPI_MAX, median_owner, cur_comm);

        MPI_Barrier(MPI_COMM_WORLD);

        // Step 6: Store tree node in single node
        long n_left = _n_points/2;
        long n_right = _n_points - n_left;
        if(rank == median_owner) {
            node_t* node = tree + write_idx;
            node->id = id;
            node->center = centers[write_idx];
            node->radius = global_rad;
            node->left = id + 1;
            node->right = id + 2*n_left;

            write_idx++;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Step 7: Divide problem to conquer professor's heart
        _n_points = (rank < median_owner ? n_left : n_right);
        id = (rank < median_owner ? id + 1: id + 2*n_left);

        // Split current comm in two
        int color = rank / (comm_size/2);
        MPI_Comm new_comm;
        MPI_Comm_split(cur_comm, color, rank, &new_comm);
        // MPI_Comm_free(&cur_comm);
        cur_comm = new_comm;

        // Update new rank
        MPI_Comm_rank(cur_comm, &rank);

        level++;
    }
    sop_t* wset = (sop_t*)malloc(sizeof(sop_t) * local_n_points);
    for(long i = 0; i < local_n_points; i++) {
        wset[i].point_idx = i;
    }

    build_tree_BIGGUS(local_n_points, wset, id, tree + write_idx, centers + write_idx, 0);

    MPI_Barrier(MPI_COMM_WORLD);
    exec_time += MPI_Wtime();

    /*
    free(SWAPPIE_SWAPPIE);
    free(SWAPPIE_SWAPPIE_ORTH);
    free(orthset);
    free(global_aux);
    */
#ifndef SKIP_DUMP
    for (long r = 0; r < n_procs; r++) {
        if (global_rank == r) {
            for (long i = 0; i < local_n_nodes; i++) {
                if(tree[i].id != -1) {
                    node_t* node = tree + i;
                    printf("%d %ld %ld %.6f",
                        node->id, node->left, node->right, sqrt(node->radius));
                    print_vec(node->center, N_DIMS);
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    if(global_rank == 0) {
        fprintf(stderr, "%.1lf\n", exec_time);
    }

    MPI_Finalize();
    return 0;
}

int solve_small(int n_procs, int global_rank, double exec_time, long n_points) {

    int rank = global_rank;
    MPI_Comm cur_comm = MPI_COMM_WORLD;

    // wset: stores working indices
    // orthset: stores orthogonal projections
    // these arrays have matching members in the same indices
    long* wset =      (long*)malloc(sizeof(long)   * n_points);
    double* orthset = (double*)malloc(sizeof(double) * n_points);
    for(long i = 0; i < n_points; i++) {
        wset[i] = i;
    }


    // allocate tree
    long n_nodes = 2*n_points - 1;
    long id = 0;
    node_t* tree = (node_t*)malloc(sizeof(node_t) * n_nodes);
    double* _centers = (double*)malloc(sizeof(double) * n_nodes * N_DIMS);
    double** centers = (double**)malloc(sizeof(double*) * n_nodes);
    if(!tree || !centers || !_centers) {
        printf("Allocation error\n");
        exit(4);
    }

    // Initializing tree with -1 to identify written nodes
    memset(tree, 0xff, sizeof(node_t) * n_nodes);

    for(int i = 0; i < n_nodes; i++) {
        centers[i] = &_centers[i*N_DIMS];
    }

    // Dump tree pt 1: Print tree dimensions
    if(global_rank == 0) {
        printf("%d %ld\n", N_DIMS, 2*n_points - 1);
    }

    // Calculate max levels
    long max_levels = 2;
    for(long aux = 1; (aux <<= 1) < n_points; max_levels++);

    long level = 0;
    while(level < max_levels && n_procs > (1 << level)) {

        int comm_size;
        MPI_Comm_size(cur_comm, &comm_size);

        long ab[2];

        if(rank == 0) {
            mpi_find_furthest_points(wset, n_points, ab, ab+1);
        }

        // broadcast a and b
        MPI_Bcast(ab, 2, MPI_LONG, 0, cur_comm);
        long a = ab[0];
        long b = ab[1];

        // split array across my current group (p/2^level processors)
        long buf_size = (long)ceil((float)n_points/(float)comm_size);
        // printf("[%d] size: %ld\n", rank, buf_size);

        long* local_wset;
        if(rank == 0) { // Master
            // Allocating a new array to keep local wset apart from global one
            local_wset = (long*)malloc(sizeof(long) * buf_size);

        } else { // Slave
            local_wset = wset;
        }
        // printf("hello from %d\n", global_rank);
        MPI_Scatter(wset, buf_size, MPI_LONG, local_wset, buf_size, MPI_LONG, 0, cur_comm);


        // We don't want to operate on data if
        // There is no data for us
        if(n_points - (rank*buf_size) > 0) {
            // rank*buf_size = the number of points attributed to the nodes
            // with a smaller rank
            buf_size = MIN(buf_size, n_points - rank*buf_size);

            mpi_calc_orth_projs(local_wset, orthset, buf_size, a, b);
            // TODO: calc medians and send them to master?? NO >:(

        } else {
            buf_size = 0;
        }

        int recvcounts[comm_size];
        int displs[comm_size];
        int count = 0;
        int remainder;
        if (rank == 0) {
            for (int i = 0; i < comm_size; i++) {
                remainder = n_points - (i*buf_size);
                if (remainder > 0) {
                    recvcounts[i] = MIN(buf_size, remainder);
                } else {
                    recvcounts[i] = 0;
                }
                displs[i] = count;
                count += recvcounts[i];
            }
        }
        // in-place gather
        if(rank == 0) {
            MPI_Gatherv(MPI_IN_PLACE, buf_size, MPI_DOUBLE, orthset, recvcounts, displs, MPI_DOUBLE, 0, cur_comm);
        } else {
            MPI_Gatherv(orthset, buf_size, MPI_DOUBLE, orthset, NULL, NULL, MPI_DOUBLE, 0, cur_comm);
        }
        // printf("After gather: %d", global_rank);

        if(rank == 0) {
            // printf("[LEVEL %ld] master %d is partitioning array ", level, rank);
            // print_vec(orthset, n_points);

            // calculate median
            long mdn_idx;
            if(n_points&1) {
                mdn_idx = n_points/2;
                mpi_select_ith(wset, orthset, n_points, mdn_idx);
                orth_proj(N_DIMS, POINTS[wset[mdn_idx]], POINTS[a], POINTS[b], centers[id]);

            } else {
                // lm = lower median
                // um = upper median
                long lmidx = (n_points-1)/2;
                mpi_select_ith(wset, orthset, n_points, lmidx);

                /* 
                 * Finding the lowest semi orth proj that is 
                 * greater or equal than the lowest median;
                 * We only need to search the right partition 
                 * since the select_ith algorithm partitions the
                 * array using the median
                 */
                long umidx = lmidx;
                double lowersop =  DBL_MAX;
                for(int i = lmidx + 1; i < n_points; i++) {
                    if(orthset[i] < lowersop) {
                        lowersop = orthset[i];
                        umidx = i;
                    }
                }

                // calculate averages
                double* lm_op = (double*) malloc(sizeof(double) * N_DIMS);
                double* um_op = (double*) malloc(sizeof(double) * N_DIMS);
                if(!lm_op || !um_op) exit(-1);
                orth_proj(N_DIMS, POINTS[wset[lmidx]], POINTS[a], POINTS[b], lm_op);
                orth_proj(N_DIMS, POINTS[wset[umidx]], POINTS[a], POINTS[b], um_op);

                vec_sum(N_DIMS, lm_op, um_op, centers[id]);
                vec_scalar_mul(N_DIMS, centers[id], 0.5, centers[id]);
                free(lm_op);
                free(um_op);
            }

            // we do not need to partition array since select ith does that for us
            mpi_partition(wset, orthset, n_points, centers[id][0]);
        }

        // Calculate center:
        // Step 1: Broadcast center to every peer
        MPI_Bcast(centers[id], N_DIMS, MPI_DOUBLE, 0, cur_comm);

        // Step 2: Compute maximum radius with assigned working set
        double global_rad;
        double local_rad = 0;
        for(long i = 0; i < buf_size; i++) {
            double new_rad = squared_dist(N_DIMS, centers[id], POINTS[local_wset[i]]);
            if(new_rad > local_rad) local_rad = new_rad;
        }

        // Step 3: Reduce to maximum
        MPI_Reduce(&local_rad, &global_rad, 1, MPI_DOUBLE, MPI_MAX, 0, cur_comm);


        long n_left = n_points/2;
        long n_right = n_points - n_left;

        if(rank == 0) {
            // Free master's local wset buffer
            // printf("Master found radius: %f\n", global_rad);

            node_t* node = tree + id;
            node->id = id;
            node->center = centers[id];
            node->radius = global_rad;
            node->left = id + 1;
            node->right = id + 2*n_left;


            free(local_wset);
        }

        // send new problem to new master (id + comm_size/2)
        int next_master = comm_size / 2;
        if(rank == 0) {
            // send partition to next master
            MPI_Send(wset + n_left, n_right, MPI_LONG, next_master, DELEGATE_MASTER, cur_comm);
        }
        if(rank == next_master) {
            // receive partition from prev master
            MPI_Recv(wset, n_right, MPI_LONG, 0, DELEGATE_MASTER, cur_comm, MPI_STATUS_IGNORE);
        }

        n_points = (rank < next_master ? n_left : n_right);
        id = (rank < next_master ? id + 1: id + 2*n_left);

        // Split current comm in two
        int color = rank / (comm_size/2);
        MPI_Comm new_comm;
        MPI_Comm_split(cur_comm, color, rank, &new_comm);
        // MPI_Comm_free(&cur_comm);
        cur_comm = new_comm;

        // Update new rank
        MPI_Comm_rank(cur_comm, &rank);

        level++;
    }

    sop_t* new_wset = (sop_t*)malloc(sizeof(sop_t) * n_points);
    for(long i = 0; i < n_points; i++) {
        new_wset[i].point_idx = wset[i];
    }
    build_tree(n_points, new_wset, id, tree, centers);

    MPI_Barrier(MPI_COMM_WORLD);
    exec_time += MPI_Wtime();

#ifndef SKIP_DUMP
    for (long r = 0; r < n_procs; r++) {
        if (global_rank == r) {
            for (long i = 0; i < n_nodes; i++) {
                if(tree[i].id != -1) {
                    node_t* node = tree + i;
                    printf("%d %ld %ld %.6f",
                        node->id, node->left, node->right, sqrt(node->radius));
                    print_vec(node->center, N_DIMS);
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    if(global_rank == 0) {
        fprintf(stderr, "%.1lf\n", exec_time);
    }

    // TODO: return value?

    MPI_Finalize();
    return 0;
}

void mpi_find_furthest_points(long* wset, long n_points, long*a, long* b) {
    if(n_points == 2) {
        *a = wset[0];
        *b = wset[1];
        return;
    }

    // find A: the most distant point from the first point in the set
    long local_a = 0;
    long local_b = 0;
    double maximum = 0.0;
    for(int i = 1; i < n_points; i++) {
        double sd = squared_dist(N_DIMS, POINTS[wset[0]], POINTS[wset[i]]);
        if(sd > maximum) {
            local_a = wset[i];
            maximum = sd;
        }
    }

    maximum = 0.0;
    for(int i = 0; i < n_points; i++) {
        double sd = squared_dist(N_DIMS, POINTS[local_a], POINTS[wset[i]]);
        if(sd > maximum) {  
            local_b = wset[i];
            maximum = sd;
        }
    }
    
    *a = local_a;
    *b = local_b;
}

void mpi_calc_orth_projs(long* wset, double* orthset, long n_points, long a_idx, long b_idx) {
    double* a = POINTS[a_idx];
    double* b = POINTS[b_idx];
    for(int i = 0; i < n_points; i++) {
        orthset[i] = semi_orth_proj(N_DIMS, POINTS[wset[i]], a, b);
    }
}

void print_vec(double* vec, int len) {
    int i;
    for(i = 0; i < len-1; i++) {
        printf(" %.6f", vec[i]);
    }
    printf(" %.6f\n", vec[i]);
}
