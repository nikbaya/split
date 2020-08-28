#! /usr/bin/env bash

set -e

# hard-parallelize jobs across clusters
# - loops over specified number of iterations (maxi)
# - spins up cluster
# - NOTE: Only use for last part of meta_split_parallel.py (permutation/meta-analysis step)

maxi=$((10))


for i in `seq 1 ${maxi}`; do
# for i in "${array[@]}"; do
	# cluster start ukbb-nb-1-$i --num-workers 2 --num-preemptible-workers 0 --max-idle 10m &
	cluster submit ukbb-nb-1-$i meta_split_parallel.py --args "--phen 50 --batch 1 --iters 100 --parsplit ${maxi} --paridx ${i}" &
done


# eof
