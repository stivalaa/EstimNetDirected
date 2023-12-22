#!/bin/sh

## Run testBipartiteAlphaBetaChangeStats and compare against
## b1nodematch / b2nodematch using statnet
## as gold standard to verify implementations changeBipartiteNodematchBetaA,
## changeBipartiteNodematchBetaB, changeBipartiteNodematchAlphaA, and
## changeBipartiteNodematchAlphaB in EstimNetDirected.

BASELINE=statnet_b1nodematch_baseline.out
OUTPUT=b1nodematch_test.out

echo "Bipartite nodematch (alpha and beta)"
for implementation in ./testBipartiteAlphaBetaChangeStats ./testBipartiteAlphaBetaChangeStats_array ./testBipartiteAlphaBetaChangeStats_hash
do
    echo ${implementation}
    cat /dev/null > ${OUTPUT}
    for i in `seq 0 10`
    do
        exponent=`echo "$i / 10" | bc -l`
        ${implementation} ../../examples/bipartite/simulated/bpnet_A12000_B4000_attrs_sim830000000.net ../../examples/bipartite/simulation/catattr_all.txt $exponent >> ${OUTPUT}
        if [ $? -ne 0 ]; then
            echo "**** FAILED ****"
            exit 1
        fi
    done
    ./compareStatsValues.py ${BASELINE} ${OUTPUT}
    if [ $? -ne 0 ]; then
        echo "**** FAILED ****"
        exit 1
    fi
done
echo PASSED
exit 0
