#!/bin/bash

# Usage: ./test/stress.sh

N_ITER=1
tests=(
    ## Professor tests
    # "20 1000000 0" # 7.3
    # "3 5000000 0"  # 23.2
    # "4 10000000 0" # 57.2
    # "3 20000000 0" # 122.5
    "4 20000000 0" # 131.6

    ## Stress test diagonal
    "5000000 10 0"
    "500000 100 0"
    "50000 1000 0"
    "5000 10000 0"
    "500 100000 0"
    "50 1000000 0"
    "5 10000000 0"
)

make bench
if [[ $? -eq 2 ]]; then
    exit
fi

si="serial_itr"
sr="serial_rec"
pi="parall_itr"
pr="parall_rec"
running=(
    ${sr}
    # ${si}
    # ${pr}
    ${pi}
)

threads=(
    "1"
    "2"
    "4"
    # "8"
)

echo "Running tests... (${N_ITER} iterations each)"
printf "%10s %10s %5s %5s |  " "#dims" "#points" "seed" "#thrds"
for bin in $running; do
    printf "%12s" ${bin}
done
printf "\n"

for arg in "${tests[@]}"; do
    for thread in "${threads[@]}"; do
        split=($arg)
        dim=${split[0]}
        num_points=${split[1]}
        seed=${split[2]}

        printf "%10s %10s %5s  %5s |  " "${dim}" "${num_points}" "${seed}" "${thread}"

        for suf in $running; do
            bin="ballAlg_${suf}"
            sum=0
            for i in $(seq 1 ${N_ITER}); do
                time="$(OMP_NUM_THREADS=${thread} ./${bin} ${arg} 2>&1> /dev/null)"
                sum=$(bc -l <<<"${sum}+${time}")
            done
            avg=$(bc -l <<<"${sum}/${N_ITER}")
            printf "%12.4f" ${avg}
        done
        printf "\n"
    done
done
