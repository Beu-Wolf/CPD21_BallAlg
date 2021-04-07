#!/bin/bash

# Usage: ./test/go.sh [serial]
# if "serial" given, test serial version

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
    "500000 100 0"

    ## Stress test diagonal
    "50 1000000 0"
    "500 100000 0"
    "5000 10000 0"
    "50000 1000 0"
    "500000 100 0"
    "5000000 10 0"
)

if [[ $# -eq 1 && $1 == "serial" ]]; then
    serial=1
fi

make serial
if [[ $? -eq 2 ]]; then
    exit
fi

# if running serial no need to compile parallel version
if [[ ! ${serial} ]]; then
    make
    if [[ $? -eq 2 ]]; then
        exit
    fi
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

[[ $serial ]] && bin="ballAlg-serial" || bin="ballAlg-omp"

echo "Running tests... (${N_ITER} iterations each)"
echo -e "#dims\t\t#points\t\tseed\t\ttime (s)"
tmp_out="tmp.txt"
for arg in "${tests[@]}"; do
    split=($arg)
    dim=${split[0]}
    num_points=${split[1]}
    seed=${split[2]}

    out="test/serial_out/${dim}_${num_points}_${seed}.out" 

    echo -ne "${dim}\t\t${num_points}\t\t${seed}\t\t"
    sum=0
    for i in $(seq 1 ${N_ITER}); do
        time="$(./${bin} ${arg} 2>&1> ${tmp_out})"
        diff -q ${tmp_out} ${out} > /dev/null
        if [[ ! $? -eq 0 ]]; then
            echo "[ERROR] Unexpected output from test [${arg}]"
            exit
        fi
        sum=$(bc -l <<<"${sum}+${time}")
    done
    avg=$(bc -l <<<"${sum}/${N_ITER}")
    echo ${avg}
done
rm ${tmp_out}
