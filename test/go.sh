#!/usr/bin/bash

# Usage: ./test/go.sh

N_ITER=1
tests=(
    ## Test Correctness
    "2 5 0"
    "2 5 394893"
    "2 5 932852"
    "2 5 490358"

    ## Test increasing number of points
    "5 1 0"
    "5 10 0"
    "5 100 0"
    "5 1000 0"
    "5 10000 0"
    "5 100000 0"
    "5 1000000 0"

    ## Test increasing number of dimensions
    "5 100 0"
    "50 100 0"
    "500 100 0"
    "5000 100 0"
    "50000 100 0"
    "500000 100 0" # seg fault :(

    ## Stress test diagonal
    "50 1000000 0"
    "500 100000 0"
    "5000 10000 0"
    "50000 1000 0"
    "500000 100 0"
    "5000000 10 0"
)

make serial
if [[ $? -eq 2 ]]; then
    exit
fi

make
if [[ $? -eq 2 ]]; then
    exit
fi


# Generate serial output if needed
for arg in "${tests[@]}"; do
    split=($arg)
    dim=${split[0]}
    num_points=${split[1]}
    seed=${split[2]}

    out="test/serial_out/${dim}_${num_points}_${seed}.out" 
    if [[ ! -f ${out} ]]; then
        echo "Generating output for ${arg}"
        ./ballAlg-serial ${arg} > ${out} 2> /dev/null
    fi

done


echo "Running tests... (${N_ITER} iterations each)"
echo -e "#dims\t\t#points\t\tseed\t\ttime (s)"
for arg in "${tests[@]}"; do
    split=($arg)
    dim=${split[0]}
    num_points=${split[1]}
    seed=${split[2]}

    out="test/serial_out/${dim}_${num_points}_${seed}.out" 

    echo -ne "${dim}\t\t${num_points}\t\t${seed}\t\t"
    sum=0
    for i in $(seq 1 ${N_ITER}); do
        tmp=$(mktemp)
        time="$(./ballAlg-omp ${arg} 2>&1> ${tmp})"
        diff -q ${tmp} ${out} > /dev/null
        if [[ ! $? -eq 0 ]]; then
            echo "[ERROR] Unexpected output from test [${arg}]"
            exit
        fi
        sum=$(bc -l <<<"${sum}+${time}")
    done
    avg=$(bc -l <<<"${sum}/${N_ITER}")
    echo ${avg}
done
