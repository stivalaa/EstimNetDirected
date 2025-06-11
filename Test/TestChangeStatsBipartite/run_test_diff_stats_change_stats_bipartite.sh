#!/bin/sh

## Run testDiffStatsChangeStatsBipartite with different implementations
## and different networks and different lambda values to verify that
## summing change statistics over all edges gives same value as
## computing statistic explicitly according to definition.

. ./netfiles.sh

echo "Testing change statistics against diff in stats for bipartite networks"
for implementation in ./testDiffStatsChangeStatsBipartite ./testDiffStatsChangeStatsBipartite_array ./testDiffStatsChangeStatsBipartite_hash
do
    echo ${implementation}
    for netfile in ${FAST_NETFILES}
    do
        echo ${netfile}
        # Only test also the slow implementations of the statistics funcions
        # on the network that is small enough that they are not too slow
        if [ `basename ${netfile}` = "inouye_pyke_pollinators_bipartite.net" ]
        then
            flags="-s"
        else
            flags=""
        fi
        for lambda in 2 5 10
        do
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
