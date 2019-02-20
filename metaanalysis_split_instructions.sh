# Meta-analysis/split pipeline cloud submission instructions

#–––––––––––– Preprocessing: ––––––––––––
# Submit a job with the three preprocess function calls at the bottom of the script uncommented. 
# If the job gets stuck, try submitting jobs with only one preprocess function call uncommented at a time. 
# If jobs with only one uncommented preprocess function call still don't work, try increasing the number of workers (may need to go up to 100 workers). 
# Below I show how I would submit jobs with only one preprocess function call uncommented.

#–––––– preprocess1 (only uncomment preprocess1 at the bottom of the script meta_split.py) ––––––

	cluster start ${clustername} --num-workers 20 --num-preemptible-workers 100 --max-idle 10m
	cluster submit ${clustername} meta_split.py --args "--phen ${phen} --n_chunks ${n_chunks} --batch ${batch} --variant_set ${variant_set}"

	# read code for descriptions of each argument passed to the script
	# if the job takes too long or seems to be stuck, try restarting the cluster and setting num-preemptibleworkers to 0 and increasing num-workers to a higher number (up to 100)


#–––––– preprocess2 (only uncomment preprocess2 at the bottom of the script meta_split.py) ––––––

	cluster start ${clustername} --num-workers 20 --num-preemptible-workers 100 --max-idle 10m
	cluster submit ${clustername} meta_split.py --args "--phen ${phen} --n_chunks ${n_chunks} --batch ${batch} --variant_set ${variant_set}"

	# if the job takes too long or seems to be stuck, try restarting the cluster and setting num-preemptibleworkers to 0 and increasing num-workers to a higher number (up to 100)

#–––––– preprocess3 (only uncomment preprocess3 at the bottom of the script meta_split.py) ––––––

	cluster start ${clustername} --num-workers 50
	cluster submit ${clustername} meta_split.py --args "--phen ${phen} --n_chunks ${n_chunks} --batch ${batch} --variant_set ${variant_set}"

	# NOTE: DO NOT USE PREEMPTIBLE WORKERS FOR THIS STEP
	# You may find that you will need to increase num-workers to 100 if the job gets stuck










#–––––––––––– Meta-split: ––––––––––––
# Be careful when running a job with both metasplit1 and metasplit2 uncommented in the same script. I typically run them separately.



#–––––– metasplit1 (only uncomment metasplit1 at the bottom of the script meta_split.py) ––––––

	cluster start ${clustername} --num-workers 20 --max-idle 10m
	cluster submit ${clustername} meta_split.py --args "--phen ${phen} --n_chunks ${n_chunks} --batch ${batch} --variant_set ${variant_set}"

	# Takes about 1 hour to finish using 20 n1-standard-8 workers


#–––––– metasplit2 (only uncomment metasplit2 at the bottom of the script meta_split.py) ––––––

	cluster start ${clustername} --num-workers 20 --max-idle 10m
	cluster submit ${clustername} meta_split.py --args "--phen ${phen} --n_chunks ${n_chunks} --batch ${batch} --variant_set ${variant_set}"

	# Takes 40 n1-standard-8 workers about 5 hours. Takes 20 n1-standard-8 workers slightly more than 8 hours. 
	# You can try playing around with the number of workers to suit your needs. Typically, fewer workers takes longer but are more cost-effective. 
	# However, I generally wouldn't recommend using fewer than 20 workers because this may cause the job to become stuck.


#–––––– metasplit3 (only uncomment metasplit3 at the bottom of the script meta_split.py) ––––––
	cluster start ${clustername} --num-workers 2 --max-idle 10m
	cluster submit ${clustername} meta_split.py --args "--phen ${phen} --n_chunks ${n_chunks} --batch ${batch} --variant_set ${variant_set}"

	# NOTE: I would recommend using meta_split_parallel.py instead of this last function. It will be much faster because the meta-analysis of the permutations can be hard parallelized over clusters without an increase in cost.




#–––––––––––– meta_split_parallel.py: ––––––––––––
	# Substitute this script for the last step in the meta-analysis/split pipeline. It is much faster because it hard parallelizes tasks over clusters.

	# (Step 1): Start up your clusters using a parallel for-loop
	# The higher maxi is, the fewer "reps" (replicates) each cluster will have to process, making the meta-split faster.
	# However, this only takes between 1.3 and 1.6 hours if maxi=5.

	maxi=5 # Make sure this is a factor of the "reps" argument passed to meta_split_parallel.py
	for i in `seq 1 ${maxi}`; do
		cluster start ${clusterprefix}-$i --num-workers 2 --max-idle 10m &
	done

	# NOTE: Wait for the clusters to all be running before proceeding to the next step
	# NOTE: However, don't let the clusters sit idle for more than 10 minutes because they will automatically be deleted

	# (Step 2): Submit tasks using parallel for-loop

	maxi=5
	for i in `seq 1 ${maxi}`; do
		cluster submit ${clusterprefix}-$i meta_split_parallel.py --args "--phen ${phen} --n_chunks ${n_chunks} --batch ${batch} --reps ${reps} --parsplit ${maxi} --paridx $i" &
	done




#–––––––––––– rg_split_parallel.py: ––––––––––––
# This is the script used to calculate genetic correlation using the ldsc from MTAG.
# It is a bit confusing to read through the script, but basically all you need to know is that you only need to check a few variables: phenfile, out_bucket, rg_outname, gs_phenfile_bucket, and gs_sumstat_dir

# DESCRIPTION OF VARIABLES
# phenfile: the name of your "h2part" file. This is a tsv file you must create for each phenotype you use. Below is my python code for generating the file. 

		import pandas as pd
		import os
		phen = '50_raw'
		variant_set = 'hm3'
		batch = 1
		n_chunks = 300

		phenotypes = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/imputed-v3-results/phenotypes.both_sexes.tsv',sep='\t')

		#vals = phenotypes.loc[phenotypes['phenotype'] == str(phen)]

		#os.system('gsutil cp gs://nbaya/rg_sex/template.h2part.tsv ~/Documents/lab/ukbb-sexdiff/rg_sex/')
		tb1 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/rg_sex/template.h2part.tsv',sep = '\t')
		tb1 = tb1.iloc[:,8:42]
		tb1.at[0,'female_file'] = variant_set+'_'+str(phen)+'_meta_A_n'+str(n_chunks)+'_batch_'+str(batch)+'_s0.tsv.bgz'
		tb1.at[0,'male_file'] = variant_set+'_'+str(phen)+'_meta_B_n'+str(n_chunks)+'_batch_'+str(batch)+'_s0.tsv.bgz'
		tb1.at[0,'phen'] = str(phen)+'_s0'
		tb1.at[0,'female_n'] = int(360388/2) #vals['n_non_missing']/2 #{'50_raw': 360338, '30100': 350423}
		tb1.at[0,'male_n'] = int(360388/2) #vals['n_non_missing']/2 #int(350423/2) #
		tb1.at[0,'female_n_cas'] = float('NaN') #vals['n_cases']/2 #float('NaN') #
		tb1.at[0,'male_n_cas'] = float('NaN') #vals['n_cases']/2 #float('NaN') #
		tb1.at[0,'female_n_con'] = float('NaN') #vals['n_controls']/2 #float('NaN') #
		tb1.at[0,'male_n_con'] = float('NaN') #vals['n_controls']/2 #float('NaN') #
		tb = tb1.copy()

    
		for i in range(1,100):
		    tb = tb.append(tb1, ignore_index=True)
		    tb.at[i,'female_file'] = variant_set+'_'+str(phen)+'_meta_A_n'+str(n_chunks)+'_batch_'+str(batch)+'_s'+str(i)+'.tsv.bgz'
		    tb.at[i,'male_file'] = variant_set+'_'+str(phen)+'_meta_B_n'+str(n_chunks)+'_batch_'+str(batch)+'_s'+str(i)+'.tsv.bgz'
		    tb.at[i,'phen'] = str(phen)+'_s'+str(i)
    

		filename = str(phen)+'.h2part.n'+str(n_chunks)+'_batch_'+str(batch)+'.tsv'
		local_wd = '~/Documents/lab/ukbb-sexdiff/rg_sex/'
		cloud_wd = 'gs://nbaya/rg_sex/'
		tb.to_csv(local_wd+filename, sep='\t',index=False)

		os.system('gsutil cp '+local_wd+filename+' '+cloud_wd)

		# NOTE: The "phenoytpes" variable above is just a table containing statistics about all phenotypes in the UKB. You will need to change the values in columns like female_n and male_n to match the phenotypes you are using. You may also need to set the n_cas and n_con values to something other than NaN if your phenotype is binary.

# out_bucket: the cloud bucket that you want your output file containing genetic correlation estimates to be saved
# rg_outname: the name of the output file
# gs_phenfile_bucket: the cloud bucket containing your phenfile, aka your "h2part" file.
# gs_sumstat_dir: the cloud bucket that contains all meta-split summary stats that you generated from meta_split_parallel.py

# submission command
cluster start ${clustername} --num-workers 2 --max-idle 10m --version 0.1 --spark 2.0.2 
cluster submit ${clustername} rg_split_parallel.py --args "--phen ${phen} --n_chunks ${n_chunks} --batch ${batch} " 

# NOTE: There is no need to hard parallelize this step because it is already parallelized. 
# This typically takes less than 15 minutes, even for 100 replicates. 


