#!/bin/sh

## Run testStatsSumChangeStatsUndirected with different implementations
## and different networks and different lambda values to verify that
## summing change statistics over all edges gives same value as
## computing statistic explicitly according to definition.

. ./netfiles.sh

echo "Testing summing change statistics for undirected networks"
for implementation in ./testStatsSumChangeStatsUndirected ./testStatsSumChangeStatsUndirected_array ./testStatsSumChangeStatsUndirected_hash
do
    echo ${implementation}
    for netfile in ${NETFILES}
    do
        echo ${netfile}
        for i in `seq 1 50`
        do
            lambda=`echo "1 + $i / 10" | bc -l`
            #echo lambda = ${lambda}
            ${implementation} ${flags} ${netfile}  ${lambda} > /dev/null
            if [ $? -ne 0 ]; then
                echo "**** FAILED ****"
                exit 1
            fi
        done
    done
done
echo PASSED
exit 0
