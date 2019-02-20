#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 16:29:43 2019

Iterates through list of phenotypes (default: heritable phenotypes in ukb_phenos_for_sex_rg.tsv),
generating genetic correlation between randomly split halves of the population.
Iteration is hard parallelized across clusters if parsplit is more than 1.

@author: nbaya
"""


import hail as hl
import numpy as np
import datetime
import argparse
import os

parser = argparse.ArgumentParser(add_help=False)

parser.add_argument('--n_chunks', type=int, required=True, help="number of partitions. Default: 2")
parser.add_argument('--batch', type=int, required=True, help="batch number for reproducibility of stochastic steps (used as seed number). Default: 1")
parser.add_argument('--phen_set', type=str, required=True, help="set of phenotypes to work with. Options: phesant, finngen, icd10, more_phesant, stillmore1, stillmore2, stillmore3")
parser.add_argument('--parsplit', type=int, required=True, help="number of parallel batches, suggested: 50")
parser.add_argument('--paridx', type=int, required=True, help="which of the parallel batches to run")

args = parser.parse_args()

n_chunks = args.n_chunks #number of partitions of the dataset
batch = args.batch #This is used for the seed to replicate the random partitioning. default: 1
phen_set = args.phen_set 

output_bucket = 'gs://nbaya/rg_sex/allphens/'

heritable_phens_tb = hl.import_table('gs://nbaya/split/ukb_phenos_for_sex_rg.tsv')
#heritable_phens_tb = hl.import_table('/Users/nbaya/Downloads/ukb_phenos_for_sex_rg.tsv')
heritable_phens_ls = [x.phenotype.strip('_irnt') for x in heritable_phens_tb.select('phenotype').collect()] #list of heritable phenotypes for UKBB round 2
 
if phen_set == 'phesant':
    phen_tb_all = hl.import_table('gs://phenotype_31063/ukb31063.phesant_phenotypes.both_sexes.tsv.bgz',
                                          missing='',impute=True,types={'"userId"': hl.tstr}).rename({ '"userId"': 's'})
elif phen_set == 'finngen':
    phen_tb_all = hl.import_table('gs://ukb31063-mega-gwas/phenotype-files/curated-phenotypes/2018-04-06_ukb-finngen-pheno-for-analysis.tsv',
                                      missing='',impute=True, types={'eid':hl.tstr}).rename({'eid':'s'})
elif phen_set == 'icd10':
    phen_tb_all = hl.import_table('gs://nbaya/sex_linreg/icd10_phenotypes.both_sexes.tsv.bgz',
                                      missing='',impute=True,types={'s': hl.tstr})
elif phen_set == 'more_phesant':
    phen_tb_all = hl.import_table('gs://ukb31063-mega-gwas/phenotype-files/still-more-phesant/cts_irnt.both_sexes.tsv',
                                      missing='',impute=True, types = {'userId':hl.tstr}).rename({ 'userId"': 's'})
elif phen_set == 'stillmore1':
    phen_tb_all = hl.import_table('gs://ukb31063-mega-gwas/phenotype-files/still-more-phesant/neale_lab_parsed_and_restricted_to_QCed_samples_cat_variables_both_sexes.1.tsv',
                                  impute=True, missing="",types={'all_sexes$userId': hl.tstr}).rename({'all_sexes$userId':'s'})
elif phen_set == 'stillmore2':
    phen_tb_all = hl.import_table('gs://ukb31063-mega-gwas/phenotype-files/still-more-phesant/neale_lab_parsed_and_restricted_to_QCed_samples_cat_variables_both_sexes.2.tsv',
                                  impute=True, missing="",types={'all_sexes$userId': hl.tstr}).rename({'all_sexes$userId':'s'})
elif phen_set == 'stillmore3':
    phen_tb_all = hl.import_table('gs://ukb31063-mega-gwas/phenotype-files/still-more-phesant/neale_lab_parsed_and_restricted_to_QCed_samples_cat_variables_both_sexes.3.tsv',
                                  impute=True, missing="",types={'all_sexes$userId': hl.tstr}).rename({'all_sexes$userId':'s'})
    
withdrawn = hl.import_table('gs://nbaya/w31063_20181016.csv',missing='',no_header=True)
withdrawn_set = set(withdrawn.f0.take(withdrawn.count()))
phen_tb_all = phen_tb_all.filter(hl.literal(withdrawn_set).contains(phen_tb_all['s']),keep=False) 
phen_tb_all = phen_tb_all.key_by('s')

phen_tb_all_cols = [x.strip('"') for x in list(phen_tb_all.row)]

completed_ls = os.popen('gsutil ls gs://nbaya/rg_sex/allphens/*').read().split('\n')
completed_ls = completed_ls[1:len(completed_ls)-1]
completed_phens = [x.split('_meta')[0].split('/')[5] for x in completed_ls if 'meta_B' in x]

phens_ls = list(set(heritable_phens_ls).intersection(set(phen_tb_all_cols)).difference(set(completed_phens)))

if len(phens_ls) == 0:
    print('\n####################')
    print(phen_set+' is complete')
    print('\n####################')

idx_ls = range(args.paridx-1,len(phens_ls), args.parsplit) #hard-parallelized

variants = hl.read_matrix_table('gs://nbaya/split/ukb31063.hm3_variants.gwas_samples_repart.mt') #hm3 variants matrix table

print('\n####################')
print('phen set: '+phen_set)
print('parsplit: '+str(args.parsplit))
print('paridx: '+str(args.paridx))
print('n chunks: '+str(n_chunks))
print('batch: '+str(batch))
print('running '+str(len(phens_ls))+' phenotypes: '+'\t'.join(phens_ls[args.paridx-1::args.parsplit]))
print('####################')
          
for phen_i in idx_ls:
    phen = phens_ls[phen_i]
        
    print('####################')
    print('Starting phen '+phen+' ('+str(phen_i+1)+' of '+str(len(phens_ls))+' '+phen_set+' phens)')
    print('####################')
    phen_starttime = datetime.datetime.now()

    """
    ╔═══════════════╗
    ║ Preprocessing ║
    ╚═══════════════╝
    """
    
    if phen_set == 'phesant':
        phen_tb_all1 = phen_tb_all.rename({'"'+phen+'"': 'phen'})
    else:
        phen_tb_all1 = phen_tb_all.rename({phen: 'phen'})
    phen_tb = phen_tb_all1.select(phen_tb_all1['phen'])
             
    mt1 = variants.annotate_cols(phen_str = hl.str(phen_tb[variants.s]['phen']).replace('\"',''))
    mt1 = mt1.filter_cols(mt1.phen_str == '',keep=False)
    
    if phen_tb.phen.dtype == hl.dtype('bool'):
        mt1 = mt1.annotate_cols(phen = hl.bool(mt1.phen_str)).drop('phen_str')
    else:
        mt1 = mt1.annotate_cols(phen = hl.float64(mt1.phen_str)).drop('phen_str')
    
    n_samples = mt1.count_cols()
    print('\n>>> phen '+phen+': N samples = '+str(n_samples)+' <<<') #expect n samples to match n_non_missing from phenotypes.both_sexes.tsv
    
    mt2 = mt1.add_col_index()
    group_size = int(n_samples/n_chunks)+1     #the ideal number of samples in each group
    group_ids = np.ndarray.tolist(np.ndarray.flatten(np.asarray([range(n_chunks)]*group_size))) #list of group ids to be paired to each sample (Note: length of group_ids > # of cols in mt, but it doesn't affect the result)
    group_ids = group_ids[0:n_samples]
    randstate = np.random.RandomState(int(batch)) #seed with batch number
    randstate.shuffle(group_ids)
    mt3 = mt2.annotate_cols(group_id = hl.literal(group_ids)[hl.int32(mt2.col_idx)]) #assign group ids
    
    mt = mt3.rename({'dosage': 'x', 'phen': 'y'})
    
    """
    ╔═══════════════════════╗
    ║ Run linear regression ║
    ╚═══════════════════════╝
    """
        
    if n_chunks == 2: #traditional split of population into even halves (no meta-analysis)
        mt_A = mt.filter_cols(mt.group_id == 0)
        mt_B = mt.filter_cols(mt.group_id == 1)
        
        cov_list_A = [ mt_A['isFemale'], mt_A['age'], mt_A['age_squared'], mt_A['age_isFemale'],
                mt_A['age_squared_isFemale'] ]+ [mt_A['PC{:}'.format(i)] for i in range(1, 21)] 
    
        cov_list_B = [ mt_B['isFemale'], mt_B['age'], mt_B['age_squared'], mt_B['age_isFemale'],
                mt_B['age_squared_isFemale'] ]+ [mt_B['PC{:}'.format(i)] for i in range(1, 21)] 
             
        ht_A = hl.linear_regression_rows(
                y=mt_A.y,
                x=mt_A.x,
                covariates=[1]+cov_list_A,
                pass_through = ['rsid'])

        ht_B = hl.linear_regression_rows(
                y=mt_B.y,
                x=mt_B.x,
                covariates=[1]+cov_list_B,
                pass_through = ['rsid'])
        
        ht_A = ht_A.rename({'rsid':'SNP'}).key_by('SNP')
        ht_B = ht_B.rename({'rsid':'SNP'}).key_by('SNP')
        
        ht_A = ht_A.select(Z = ht_A.beta/ht_A.standard_error)
        ht_B = ht_B.select(Z = ht_B.beta/ht_B.standard_error)
        
        sumstats_template = hl.import_table('gs://nbaya/rg_sex/50_snps_alleles_N.tsv.gz',types={'N': hl.tint64})
        sumstats_template = sumstats_template.key_by('SNP')
        sumstats_template = sumstats_template.annotate(N = int(n_samples/2))
#            sumstats_template.show()
        
        metaA = sumstats_template.annotate(Z = ht_A[sumstats_template.SNP]['Z'])
        metaB = sumstats_template.annotate(Z = ht_B[sumstats_template.SNP]['Z'])
        
        metaA_path = output_bucket+phen+'_meta_A_nchunks'+str(n_chunks)+'_batch'+str(batch)+'_split.tsv.bgz' 
        metaB_path = output_bucket+phen+'_meta_B_nchunks'+str(n_chunks)+'_batch'+str(batch)+'_split.tsv.bgz' 
        metaA.export(metaA_path)
        metaB.export(metaB_path)
        
        print('\n####################')
        print('Files written to:')
        print(metaA_path+'\t'+metaB_path)
        print('####################')
        
    elif n_chunks > 2:  #jackknife split        
        for split_i in range(n_chunks):
            print('####################')
            print('Starting split '+str(split_i+1)+' of '+str(n_chunks))
            print('####################')
            starttime = datetime.datetime.now()
            
            mt_A = mt.filter_cols(mt.group_id != split_i) #main cohort of n_chunk-1 partitions
            mt_B = mt.filter_cols(mt.group_id == split_i) #single cohort of 1 partition
            
            cov_list_A = [ mt_A['isFemale'], mt_A['age'], mt_A['age_squared'], mt_A['age_isFemale'],
                    mt_A['age_squared_isFemale'] ]+ [mt_A['PC{:}'.format(i)] for i in range(1, 21)] 
        
            cov_list_B = [ mt_B['isFemale'], mt_B['age'], mt_B['age_squared'], mt_B['age_isFemale'],
                    mt_B['age_squared_isFemale'] ]+ [mt_B['PC{:}'.format(i)] for i in range(1, 21)] 
                 
            ht_A = hl.linear_regression_rows(
                    y=mt_A.y,
                    x=mt_A.x,
                    covariates=[1]+cov_list_A,
                    pass_through = ['rsid'])
    
            ht_B = hl.linear_regression_rows(
                    y=mt_B.y,
                    x=mt_B.x,
                    covariates=[1]+cov_list_B,
                    pass_through = ['rsid'])
            
            ht_A = ht_A.rename({'rsid':'SNP'}).key_by('SNP')
            ht_B = ht_B.rename({'rsid':'SNP'}).key_by('SNP')
            
            ht_A = ht_A.select(Z = ht_A.beta/ht_A.standard_error)
            ht_B = ht_B.select(Z = ht_B.beta/ht_B.standard_error)
            
            sumstats_template = hl.import_table('gs://nbaya/rg_sex/50_snps_alleles_N.tsv.gz',types={'N': hl.tint64})
            sumstats_template = sumstats_template.key_by('SNP')
            sumstats_template = sumstats_template.annotate(N = int(n_samples/2))
            sumstats_template.show()
            
            metaA = sumstats_template.annotate(Z = ht_A[sumstats_template.SNP]['Z'])
            metaB = sumstats_template.annotate(Z = ht_B[sumstats_template.SNP]['Z'])
            
            metaA_path = output_bucket+phen+'_meta_A_nchunks'+str(n_chunks)+'_batch'+str(batch)+'_split'+str(split_i)+'of'+str(n_chunks)+'.tsv.bgz' 
            metaB_path = output_bucket+phen+'_meta_B_nchunks'+str(n_chunks)+'_batch'+str(batch)+'_split'+str(split_i)+'of'+str(n_chunks)+'.tsv.bgz' 
            metaA.export(metaA_path)
            metaB.export(metaB_path)
            
            endtime = datetime.datetime.now()
            elapsed = endtime-starttime
            print('\n####################')
            print('Completed split '+str(split_i+1)+' of '+str(n_chunks))
            print('Files written to:')
            print(metaA_path+'\t'+metaB_path)
            print('Split time: '+str(round(elapsed.seconds/60, 2))+' minutes')
            print('####################')  
        
    phen_endtime = datetime.datetime.now()
    phen_elapsed = phen_endtime -phen_starttime
    print('\n####################')
    print('Finished phen '+phen+' ('+str(phen_i+1)+' of '+str(len(phens_ls))+')')
    print('Phen iter time: '+str(round(phen_elapsed.seconds/60, 2))+' minutes')
    print('####################')  