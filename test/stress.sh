#!/bin/bash

# Usage: ./test/go.sh [serial]
# if "serial" given, test serial version
# if "tasks" given, test tasks version

N_ITER=1
tests=(
    ## Professor tests
    "20 1000000 0" # 7.3
    "3 5000000 0"  # 23.2
    "4 10000000 0" # 57.2
    "3 20000000 0" # 122.5
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

serial="serial"
tasks="recur"
iter="iter"
running="${serial} ${tasks} ${iter}"

echo "Running tests... (${N_ITER} iterations each)"
printf "%10s %10s %10s  |  " "#dims" "#points" "seed"
for bin in $running; do
    printf "%10s" ${bin}
done
printf "\n"

for arg in "${tests[@]}"; do
    split=($arg)
    dim=${split[0]}
    num_points=${split[1]}
    seed=${split[2]}

    printf "%10s %10s %10s  |  " "${dim}" "${num_points}" "${seed}"

    for suf in $running; do
        bin="ballAlg_${suf}"
        sum=0
        for i in $(seq 1 ${N_ITER}); do
            time="$(./${bin} ${arg} 2>&1> /dev/null)"
            sum=$(bc -l <<<"${sum}+${time}")
        done
        avg=$(bc -l <<<"${sum}/${N_ITER}")
        printf "%10.4f" ${avg}
    done

    printf "\n"
done
