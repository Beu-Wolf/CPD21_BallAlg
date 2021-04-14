#!/bin/bash

# Usage: ./test/check.sh

tests=(
    ## Test Correctness
    "2 5 0"
    "2 5 394893"
    "2 5 932852"
    "2 5 490358"
    "2 8 0"
    "2 8 394893"
    "2 8 932852"
    "2 8 490358"
)
                                                                                                        
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

versions="ballAlg-omp ballAlg-tasks"

tmp_out="tmp.txt"
rm -f test/err/*
for bin in ${versions}; do
    printf "\n\n[ ============== %15s ============== ]\n" ${bin}
    printf "%10s %10s %10s  |  RESULT\n" "#dims" "#points" "seed"
    for arg in "${tests[@]}"; do
        split=($arg)
        dim=${split[0]}
        num_points=${split[1]}
        seed=${split[2]}
        
        name="${dim}_${num_points}_${seed}.out"
        out="test/serial_out/${name}" 
        printf "%10s %10s %10s  |  " "${dim}" "${num_points}" "${seed}"

        ./${bin} ${arg} 2> /dev/null 1> ${tmp_out}
        diff -q ${tmp_out} ${out} > /dev/null
        if [[ ! $? -eq 0 ]]; then
            echo "FAILED"
            cp ${tmp_out} "test/err/${bin}_${name}.err"
        else
            echo "OK"
        fi
    done
done
rm -f ${tmp_out}
