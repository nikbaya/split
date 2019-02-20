#! /usr/bin/env bash

set -e

# hard-parallelize jobs across clusters
# - loops specified number of batches
# - spins up cluster
# - NOTE: Only use for last part of meta_split_parallel.py (permutation/meta-analysis step)

maxi=$((10))
# array=( 8 )

batch=1

for i in `seq 6 ${maxi}`; do
# for i in "${array[@]}"; do
	# cluster start ukbb-nb${batch}-$i -m n1-standard-16 --num-workers 2 --num-preemptible-workers 0 --max-idle 10m &
	# cluster submit ukbb-nb${batch}-$i meta_downsample.py --args "--phen 20160 --batch ${batch} --subbatch $i --start_idx 0 --stop_idx 34 --paridx ${i} --parsplit ${maxi}" &
	cluster submit ukbb-nb${batch}-$i meta_downsample.py --args "--phen 50 --batch ${batch} --subbatch $i --start_idx 0 --stop_idx 34 --paridx 1 --parsplit 1" &
	# cluster submit ukbb-nb1-5 meta_downsample.py --args "--phen 20160 --batch 1 --start_idx 0 --stop_idx 34 --paridx 5 --parsplit 5" &
done


# eof
