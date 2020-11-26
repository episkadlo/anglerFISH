#!/bin/bash

lowerLimit=$1
upperLimit=$2
genomeName=$3
refGenome=$4

for i in $(seq $lowerLimit $upperLimit); do
  indexName="${genomeName}_${i}"
  jellyfish count -m $i --out-counter-len 1 -L 2 $refGenome -s 100M -t 10 -o $indexName.jf
done
