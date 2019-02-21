# split
Used for running GWAS on split datasets.
<br>

## Chunking/Meta-analysis instructions

### Preprocessing:
Submit a job with the three preprocess function calls at the bottom of the script uncommented. 
<br>
If the job gets stuck, try submitting jobs with only one preprocess function call uncommented at a time. 
<br>
If jobs with only one uncommented preprocess function call still don't work, try increasing the number of workers (may need to go up to 100 workers). 
<br>
Below I show how I would submit jobs with only one preprocess function call uncommented.

##### preprocess1

	cd ${directory containing meta_split.py}
	cluster start ${clustername} --num-workers 20 --num-preemptible-workers 100 --max-idle 10m
	cluster submit ${clustername} meta_split.py --args "--phen ${phen} --n_chunks ${n_chunks} --batch ${batch} --variant_set ${variant_set}"

read code for descriptions of each argument passed to the script 
<br>
if the job takes too long or seems to be stuck, try restarting the cluster and setting num-preemptibleworkers to 0 and increasing num-workers to a higher number (up to 100)


##### preprocess2 

	cluster start ${clustername} --num-workers 20 --num-preemptible-workers 100 --max-idle 10m
	cluster submit ${clustername} meta_split.py --args "--phen ${phen} --n_chunks ${n_chunks} --batch ${batch} --variant_set ${variant_set}"

if the job takes too long or seems to be stuck, try restarting the cluster and setting num-preemptibleworkers to 0 and increasing num-workers to a higher number (up to 100)

##### preprocess3 

	cluster start ${clustername} --num-workers 50
	cluster submit ${clustername} meta_split.py --args "--phen ${phen} --n_chunks ${n_chunks} --batch ${batch} --variant_set ${variant_set}"

*NOTE: DO NOT USE PREEMPTIBLE WORKERS FOR THIS STEP*
<br>
You may find that you will need to increase num-workers to 100 if the job gets stuck

<br>
### Meta-split: 
Be careful when running a job with both metasplit1 and metasplit2 uncommented in the same script. I typically run them separately.



##### metasplit1 

	cluster start ${clustername} --num-workers 20 --max-idle 10m
	cluster submit ${clustername} meta_split.py --args "--phen ${phen} --n_chunks ${n_chunks} --batch ${batch} --variant_set ${variant_set}"

Takes about 1 hour to finish using 20 n1-standard-8 workers


##### metasplit2

	cluster start ${clustername} --num-workers 20 --max-idle 10m
	cluster submit ${clustername} meta_split.py --args "--phen ${phen} --n_chunks ${n_chunks} --batch ${batch} --variant_set ${variant_set}"

Takes 40 n1-standard-8 workers about 5 hours. Takes 20 n1-standard-8 workers slightly more than 8 hours. 
<br>
You can try playing around with the number of workers to suit your needs. Typically, fewer workers takes longer but are more cost-effective. 
<br>
However, I generally wouldn't recommend using fewer than 20 workers because this may cause the job to become stuck.


##### metasplit3 
	cluster start ${clustername} --num-workers 2 --max-idle 10m
	cluster submit ${clustername} meta_split.py --args "--phen ${phen} --n_chunks ${n_chunks} --batch ${batch} --variant_set ${variant_set}"

*NOTE: I would recommend using meta_split_parallel.py instead of this last function. It will be much faster because the meta-analysis of the permutations can be hard parallelized over clusters without an increase in cost.*



<br>
### meta_split_parallel.py:
Substitute this script for the last step in the meta-analysis/split pipeline. It is much faster because it hard parallelizes tasks over clusters.

###### (Step 1): Start up your clusters using a parallel for-loop

	maxi=5
	for i in `seq 1 ${maxi}`; do
		cluster start ${clusterprefix}-$i --max-idle 10m &
	done

*NOTE: Wait for the clusters to all be running before proceeding to the next step*
<br>
*NOTE: However, don't let the clusters sit idle for more than 10 minutes because they will automatically be deleted*

###### (Step 2): Submit tasks using parallel for-loop

	cd ${directory containing meta_split_parallel.py}
	maxi=5
	for i in `seq 1 ${maxi}`; do
		cluster submit ${clusterprefix}-$i meta_split_parallel.py --args "--phen ${phen} --n_chunks ${300} --batch ${batch} --reps 100 --parsplit ${maxi} --paridx $i" &
	done

