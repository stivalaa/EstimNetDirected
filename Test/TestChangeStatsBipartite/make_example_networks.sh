#!/bin/sh

Rscript makeExampleNetworks.R  | tee makeExampleNetworks.out
Rscript makeExampleNetworkTwoComponents.R  | tee makeExampleNetworkTwoComponents.out
for i in *.eps
do
  ./convert_eps_to_png_and_crop.sh $i
done
