#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 08:37:36 2018

Used for the downsampling (or upsampling) simulation, to test the number of 
individuals needed to have a stable heritability estimate.

Generates summary statistics for sets of chunks, starting with a set of 1 chunk
and progressively adding more chunks to the set. Note that a chunk is what we
call a group of individuals who were part of the same linear regression.

This file uses meta_split.py (and possibly meta_split_parallel.py) to create the
gmt matrix table.

@author: nbaya
"""

import hail as hl
import numpy as np
import datetime
import argparse

phen_dict = {
    '50':'height irnt',
    '20160':'smoking',
    '2443':'diabetes',
    '23107':'legimp_r',
    '6138_100':'qual_none',
    '50_raw':'height_raw',
    '30100':'platelet_vol',
    '50_sim_inf':'height_sim_infinitesimal',
    '50_sim_inf_h2_0.1':'height_sim_infinitesimal, h2=0.1',
    '50_raw_res':'Residuals from 50_raw linreg'
}


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--phen', type=str, required=True, help="phenotype code")
parser.add_argument('--n_chunks', type=int, required=True, help="number of subgroups (or 'chunks')")
parser.add_argument('--isnested', type=int, required=True, help="whether to do nested sets, default: 1")
parser.add_argument('--batch', type=str, required=True, help="batch number")
parser.add_argument('--subbatch', type=str, required=True, help="sub-batch number")

args = parser.parse_args()

#Globals
phen = args.phen #code for phenotype
desc = phen_dict[phen]
n_chunks = args.n_chunks #number of subgroups (or "chunks")
isnested = bool(args.isnested)
batch = args.batch #batch number
subbatch = args.subbatch #batch number


# Note: start_idx and stop_idx are indexing the subsets list (actually an numpy array)
if n_chunks == 300:
    subsets = np.asarray(np.concatenate((np.linspace(1,20,20),np.linspace(21,42,8),
                                     np.linspace(48,84,7)),axis=0).astype(int).tolist()) #for n=300
    idx = range(0,35)
elif n_chunks == 150:
    subsets = np.asarray(np.concatenate((np.linspace(1,20,20),np.linspace(21,42,10)),axis=0).astype(int).tolist()) #for n=150
    idx = range(0,30)

print('\n####################')
print('Phenotype: '+phen)
print('Description: '+desc)
print('n chunks: '+str(n_chunks))
print('Batch: '+batch)
print('Sub-Batch: '+subbatch)
print('Subsets: '+str(subsets[list(idx)]))
print('####################')

"""
╔═════════════════════════╗
║ Part 1: Generate subset ║
╚═════════════════════════╝
This step determines the order in which group IDs are added to the sample subset (nested downsampling)
"""
if isnested:
    groups = list(range(n_chunks))
    #randstate = np.random.RandomState(int(batch)) #seed with batch number (used for batches 1-5 for phen 50, 20160)
    seed_id = int(batch+subbatch.zfill(4)) #OPTION 2: create a seed_id unique to every split
    randstate = np.random.RandomState(seed_id) #OPTION 2: seed with seed_id
    randstate.shuffle(groups)

"""
╔═══════════════════════════════════════════════════╗
║ Part 2: Calculate meta-analyzed summary statistics║
╚═══════════════════════════════════════════════════╝
Uses inverse-variance weighted meta-analysis
"""
print('Starting Part 5: Calculate meta-analyzed sumstats')
print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))

gmt = hl.read_matrix_table('gs://nbaya/split/meta_split/ukb31063.hm3_'+phen+'_gmt'+str(n_chunks)+'_batch_'+batch+'.mt')
gmt = gmt.add_col_index()
gmt = gmt.rename({'rsid': 'SNP'})

for subset_i in subsets[list(idx)].tolist():

    starttime = datetime.datetime.now()
    print('\n####################')
    print('Running meta-analysis for subset '+ str(subset_i))
    print('####################')
    
    if not isnested:
        groups = list(range(n_chunks))
        seed_id = int(batch+subbatch.zfill(4)+str(subset_i).zfill(4)) #OPTION 2: create a seed_id unique to every split
        randstate = np.random.RandomState(seed_id) #OPTION 2: seed with seed_id
        randstate.shuffle(groups)
    subset = hl.literal(set(groups[0:subset_i])) #set of group ids to use in sample
    gmt_subset = gmt.filter_cols(subset.contains(gmt['group_id']))
    pi = np.array(['A']*n_chunks).tolist()
    gmt_lab = gmt_subset.annotate_cols(label = hl.literal(pi)[hl.int32(gmt_subset.col_idx)])
    ht = gmt_lab.group_cols_by(gmt_lab.label).aggregate(
            unnorm_meta_beta=hl.agg.sum(gmt_lab.beta / gmt_lab.standard_error ** 2),
            inv_se2 = hl.agg.sum(1 / gmt_lab.standard_error ** 2),
        group_id = hl.agg.sum(gmt_lab.group_id)).make_table()
    ht = ht.annotate(A_Z = ht['A.unnorm_meta_beta'] / hl.sqrt(ht['A.inv_se2']))
    ht = ht.drop('A.unnorm_meta_beta','A.inv_se2').key_by('SNP')
    
    variants = hl.import_table('gs://nbaya/rg_sex/50_snps_alleles_N.tsv.gz',types={'N': hl.tint32}) #can use variants table from height phenotype because all phenotypes have the same hm3 variants
    variants = variants.key_by('SNP')
#    mt_all = hl.read_matrix_table('gs://nbaya/split/ukb31063.hm3_variants.gwas_samples_'+phen+'_grouped_batch_'+batch+'.mt') #matrix table containing individual samples (used for phen=50,20160)
    mt_all = hl.read_matrix_table('gs://nbaya/split/ukb31063.hm3_variants.gwas_samples_'+phen+'_grouped'+str(n_chunks)+'_batch_'+batch+'.mt') #matrix table containing individual samples
    n_samples = int(mt_all.count_cols()*subset_i/n_chunks)
    variants = variants.annotate(N = hl.int32(n_samples))
    print('\nN samples for phen '+phen+', subset '+str(subset_i)+' = '+str(n_samples))
    
    metaA = variants.annotate(Z = ht[variants.SNP].A_Z)
                              
#    metaA_path = 'gs://nbaya/rg_sex/'+phen+'_sample_A_batch_'+batch+'.'+subbatch+'_set'+str(subset_i)+'_nested'+str(isnested)+'.tsv.bgz' # used for n = 300 version

    metaA_path = 'gs://nbaya/rg_sex/'+phen+'_sample_n'+str(n_chunks)+'_batch_'+batch+'.'+subbatch+'_set'+str(subset_i)+'_nested'+str(isnested)+'.tsv.bgz' # used for n = 300 version
    metaA.export(metaA_path)
    
    endtime = datetime.datetime.now()
    elapsed = endtime-starttime
    print('####################')
    print('Completed sample subset '+str(subset_i))
    print('File written to:')
    print(metaA_path)
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    print('Iteration time: '+str(round(elapsed.seconds/60, 2))+' minutes')
    print('####################')    


print('####################')
print('Finished Part 5')
print('Meta-sampling complete')
print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
print('####################')