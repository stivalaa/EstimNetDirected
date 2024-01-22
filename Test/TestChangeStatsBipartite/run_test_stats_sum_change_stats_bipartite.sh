#!/bin/sh

## Run testStatsSumChangeStatsBipartite with different implementations
## and different networks and different lambda values to verify that
## summing change statistics over all edges gives same value as
## computing statistic explicitly according to definition.


echo "Testing summing change statistics for bipartite networks"
for implementation in ./testStatsSumChangeStatsBipartite ./testStatsSumChangeStatsBipartite_array ./testStatsSumChangeStatsBipartite_hash
do
    echo ${implementation}
    for i in `seq 1 50`
    do
        netfile=../../examples/bipartite/simulated/bpnet_A12000_B4000_attrs_sim830000000.net
        lambda=`echo "1 + $i / 10" | bc -l`
        echo lambda = ${lambda}
        ${implementation} ${netfile}  ${lambda}
        if [ $? -ne 0 ]; then
            echo "**** FAILED ****"
            exit 1
        fi
    done
done
echo PASSED
exit 0
