# !/bin/bash

for i in 1 2 4 8 16 32 64; do
    for j in "20 1000000 0" "3 5000000 0" "4 10000000 0" "3 20000000 0" "4 20000000 0"; do
        srun -o "./mpi_outputs/$i.$j.txt" -n $i ./ballAlg-mpi $j
    done
done