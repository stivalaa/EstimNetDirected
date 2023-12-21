#!/bin/sh

## Run tstBipartiteAlphaBetaChangeStats and compare against
## b1nodematch / b2nodematch using statnet
## as gold standard to verify implementations changeBipartiteNodematchBetaA,
## changeBipartiteNodematchBetaB, changeBipartiteNodematchAlphaA, and
## changeBipartiteNodematchAlphaB in EstimNetDirected.

rc=0

for i in `seq 0 100`
do
  alpha=`echo "$i / 100" | bc -l`
  echo $alpha
  ./testBipartiteAlphaBetaChangeStats ../../examples/bipartite/simulated/bpnet_A12000_B4000_attrs_sim830000000.net ../../examples/bipartite/simulation/catattr_all.txt $alpha
done

exit $rc
