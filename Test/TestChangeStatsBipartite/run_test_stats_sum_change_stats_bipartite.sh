#!/bin/sh

## Run testStatsSumChangeStatsBipartite with different implementations
## and different networks and different lambda values to verify that
## summing change statistics over all edges gives same value as
## computing statistic explicitly according to definition.


echo "Testing summing change statistics for bipartite networks"
for implementation in ./testStatsSumChangeStatsBipartite ./testStatsSumChangeStatsBipartite_array ./testStatsSumChangeStatsBipartite_hash
do
    echo ${implementation}
    for netfile in chain_bipartite.net fourcycle3_bipartite.net fourcycle4_fourcycle_components_bipartite.net fourcycle6_bipartite.net fourcycle_bipartite.net opsahl_bipartite.net ring_bipartite.net  star_bipartite.net twopath_bipartite.net grid_open_bipartite.net grid_bipartite.net fourcycle3_revmode_bipartite.net fourcycle3_open_bipartite.net fourcycle3_leaf_bipartite.net  fourcycle3_leaf_revmode_bipartite.net inouye_pyke_pollinators_snowball_sample_2_waves_seed_3.net  ../../examples/bipartite/inouye_pyke_pollinators/inouye_pyke_pollinators_bipartite.net   ../../examples/bipartite/simulated/bpnet_A6000_B750_sparse_sim100000000.net robertson_pollinators_bipartite.net #  ../../examples/bipartite/simulated/bpnet_A12000_B4000_attrs_sim830000000.net
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
