#!/bin/sh

## Run testDiffStatsChangeStatsAttrBipartite_array ith different
## implementations to verify that change statistic gives same result
## as subtracting old from new statistic explicitly according to definition.

echo "Testing diff and change statistics for bipartite attribute node networks"
for implementation in ./testDiffStatsChangeStatsAttrBipartite ./testDiffStatsChangeStatsAttrBipartite_array ./testDiffStatsChangeStatsAttrBipartite_hash
do
    echo ${implementation}
    
    ${implementation}  ../../examples/bipartite/simulated/bpnet_A12000_B4000_attrs_sim830000000.net  ../../examples/bipartite/simulation/binattr_all.txt
    if [ $? -ne 0 ]; then
        echo "**** FAILED ****"
        exit 1
    fi
done
echo PASSED
exit 0
