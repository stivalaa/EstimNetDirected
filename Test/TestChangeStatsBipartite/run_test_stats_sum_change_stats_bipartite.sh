#!/bin/sh

## Run testStatsSumChangeStatsBipartite with different implementations
## and different networks and different lambda values to verify that
## summing change statistics over all edges gives same value as
## computing statistic explicitly according to definition.

. ./netfiles.sh

echo "Testing summing change statistics for bipartite networks"
for implementation in ./testStatsSumChangeStatsBipartite ./testStatsSumChangeStatsBipartite_array ./testStatsSumChangeStatsBipartite_hash
do
    echo ${implementation}
    for netfile in ${NETFILES}
    do
        echo ${netfile}
        # Only test also the slow implementations of the statistics funcions
        # on the network that is small enough that they are not too slow
        if [ `basename ${netfile}` == "inouye_pyke_pollinators_bipartite.net" ]
        then
            flags="-s"
        else
            flags=""
        fi
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
