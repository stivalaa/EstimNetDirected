#!/bin/sh

## Run testStatsSumChangeStatsAttrBipartite with different
## implementations to verify that summing change statistics over all
## dges gives same value as computing statistic explicitly according to
## definition.

echo "Testing summing change statistics for bipartite attribute node networks"
for implementation in ./testStatsSumChangeStatsAttrBipartite ./testStatsSumChangeStatsAttrBipartite_array ./testStatsSumChangeStatsAttrBipartite_hash
do
    echo ${implementation}
    
    stats=`${implementation} ../TestChangeStatsBipartite/twopath_bipartite.net  twopath_binattr.txt`
    if [ $? -ne 0 ]; then
        echo "**** FAILED ****"
        exit 1
    fi
    if [ "${stats}" != "0 2 0 1 " ]; then
        echo "**** FAILED ****"
        exit 1
    fi
    
    ${implementation}  ../../examples/bipartite/simulated/bpnet_A12000_B4000_attrs_sim830000000.net  ../../examples/bipartite/simulation/binattr_all.txt
    if [ $? -ne 0 ]; then
        echo "**** FAILED ****"
        exit 1
    fi
done
echo PASSED
exit 0
