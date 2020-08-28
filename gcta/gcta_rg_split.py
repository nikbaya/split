#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 14:43:46 2019

Script to generate a replicate of random split GCTA rg estimate

@author: nbaya
"""

import argparse
import subprocess
import os
import pandas as pd
from itertools import chain
import hail as hl
import numpy as np
#import matplotlib.pyplot as plt
import datetime


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--min_id', type=int, required=False, help="minimum replicate ID for batch")
parser.add_argument('--max_id', type=int, required=False, help="maximum replicate ID for batch")
parser.add_argument('--parsplit', type=int, required=False, help="number of processes to run in parallel for the batch, should be factor of (max_id-min_d)")
parser.add_argument('--paridx', type=int, required=False, help="index for which parallel process to run")

args = parser.parse_args()
min_id = args.min_id
max_id = args.max_id
parsplit = args.parsplit
paridx = args.paridx

def make_grm(nsamples):
    nsamples=str(int(nsamples/1000))
    wd = '/home/gcta/'
    print(f'\r#########\nCalculating GRM for {nsamples}k samples\n#########')
    if not os.path.isfile(wd+'gcta64'):
        print('###### Installing gcta64')
        subprocess.call(['gsutil','cp','gs://nbaya/split/gcta/gcta64', wd+'gcta64'])
    if not (os.path.isfile(wd+f'gcta_{nsamples}k.bed') & os.path.isfile(wd+f'gcta_{nsamples}k.bim') & os.path.isfile(wd+f'gcta_{nsamples}k.fam')):
        print(f'###### Installing grm files for {nsamples}k individuals from UKB')
        subprocess.call(['gsutil','-m','cp',f'gs://nbaya/split/gcta/gcta_{nsamples}k.b*', wd])
        subprocess.call(['gsutil','-m','cp',f'gs://nbaya/split/gcta/gcta_{nsamples}k.fam', wd])
    try:
        is_complete = subprocess.check_output([f'gsutil','ls',f'gs://nbaya/split/gcta/gcta_{nsamples}k.grm.bin']) != None
    except:
        is_complete = False
    if not is_complete:            
        os.chdir(wd)
        subprocess.call(['chmod', '+x','gcta64'])
        subprocess.call(['./gcta64',
                         '--bfile', f'gcta_{nsamples}k',
                         '--autosome',
                         '--make-grm',
                         '--out',f'gcta_{nsamples}k',
                         '--thread-num', '10'])
        subprocess.call(['gsutil','-m','cp',f'gcta_{nsamples}k.grm*','gs://nbaya/split/gcta/'])
    else:
        print(f'###### Already completed GRM for {nsamples}k samples')


def get_phen_files(nsamples,min_id,max_id,parsplit,paridx):
    nsamples=str(int(nsamples/1000))
    print(f'\r#########\nGetting phen files for {nsamples}k samples\n#########')
    mt0 = hl.read_matrix_table('gs://nbaya/ldscsim/hm3.50_sim_h2_0.08.mt/')
    ht0 = mt0.select_cols(mt0.nonsim_phen).cols()
    ht1 = ht0.rename({'s':'IID','nonsim_phen':'y'})
    ht1 = ht1.annotate(FID = '0')
    ht1 = ht1.key_by(ht1.FID)
    ht1 = ht1.select(ht1.IID,ht1.y)
    ht1 = ht1.key_by(ht1.IID)
    ids = hl.import_table(f'gs://nbaya/split/gcta/gcta_{nsamples}k.grm.id',no_header=True) #GRM ids
    ids = ids.rename({'f0':'FID','f1':'IID'})
    ids = set(ids.IID.take(ids.count()))
    ht2 = ht1.filter(hl.literal(ids).contains(ht1['IID']))
    n = ht2.count()
    rep_ids = range(min_id+paridx-1,max_id+1,parsplit) #replicate "IDs", which were used as seeds to generate the random split
    for rep_id in rep_ids:
        try:
            is_complete = subprocess.check_output(['gsutil','ls',f'gs://nbaya/split/gcta/gcta_{nsamples}k.s{rep_id}.phen']) != None
        except:
            is_complete = False
        if not is_complete:
            start=datetime.now()
            pi = [1]*int(n/2) + [0]*int(n/2)
            randstate = np.random.RandomState(rep_id)
            randstate.shuffle(pi)
            ht = ht2.add_index()
            ht = ht.annotate(label = hl.literal(pi)[hl.int32(ht.idx)])
            ht = ht.annotate(y1 = hl.cond(ht.label==1, ht.y, hl.null('float')))
            ht = ht.annotate(y2 = hl.cond(ht.label==0, ht.y, hl.null('float')))
            ht = ht.drop(ht.idx,ht.label,ht.y)
            ht = ht.order_by(ht.y1)
            ht.show()
            ht.export(f'gs://nbaya/split/gcta/gcta_{nsamples}k.s{rep_id}.phen')
            runtime=datetime.now()-start
            print(f'######\nRuntime for generating phenfile of rep {rep_id}: {round((runtime.total_seconds())/60, 4)} min')
        else:
            print(f'###### Already completed phenfile for replicate #{rep_id}')

    
def run_gcta_rg(nsamples,min_id,max_id,parsplit,paridx):
    '''Calculate rg of random split using GCTA'''
    nsamples=str(int(nsamples/1000))
    rep_ids = range(min_id+paridx-1,max_id+1,parsplit) #replicate "IDs", which were used as seeds to generate the random split
    wd = '/home/gcta/'
    print(f'\r#########\nRunning gcta rg for replicates: {list(rep_ids)}\n#########')
    if not os.path.isfile(wd+'gcta64'):
        print('###### Installing gcta64')
        subprocess.call(['gsutil','cp','gs://nbaya/split/gcta/gcta64', wd+'gcta64'])
    if not (os.path.isfile(wd+f'gcta_{nsamples}k.grm.N.bin') & os.path.isfile(wd+f'gcta_{nsamples}k.grm.bin') & os.path.isfile(wd+f'gcta_{nsamples}k.grm.id')):
        print(f'###### Downloading GRM files for {nsamples}k individuals from UKB')
        subprocess.call(['gsutil','-m','cp',f'gs://nbaya/split/gcta/gcta_{nsamples}k.grm*', wd])
    for rep_id in rep_ids:
        try:
            is_complete = subprocess.check_output(['gsutil','ls',f'gs://nbaya/split/gcta/gcta_{nsamples}k.s{rep_id}.hsq']) != None
        except:
            is_complete = False
        if not is_complete:            
            if not os.path.isfile(wd+f'gcta_{nsamples}k.s{rep_id}.phen'):        
                print(f'###### Downloading phen file for {nsamples}k replicate #{rep_id}')
                subprocess.call(['gsutil','cp',f'gs://nbaya/split/gcta/gcta_{nsamples}k.s{rep_id}.phen',wd])
            os.chdir(wd)
            subprocess.call(['chmod', '+x','gcta64'])
            out=f'gcta_{nsamples}k.s{rep_id}'
            subprocess.call(['./gcta64',
                             '--reml-bivar',
                             '--reml-bivar-no-constrain',
                             '--grm',f'gcta_{nsamples}k',
                             '--pheno',f'gcta_{nsamples}k.s{rep_id}.phen', 
                             '--out',out,
                             '--thread-num', '96'])
            subprocess.call(['gsutil','cp',out+'.hsq','gs://nbaya/split/gcta/'])
        else:
            print(f'###### Already completed GCTA rg calculation for {nsamples}k replicate #{rep_id}')

def run_gcta_he_rg(nsamples,min_id,max_id,parsplit,paridx):
    '''Calculate rg of random split using GCTA Haseman-Elston regression'''
    nsamples=str(int(nsamples/1000))
    rep_ids = range(min_id+paridx-1,max_id+1,parsplit) #replicate "IDs", which were used as seeds to generate the random split
    wd = '/home/gcta/'
    print(f'\r#########\nRunning gcta HEreg rg for replicates: {list(rep_ids)}\n#########')
    if not os.path.isfile(wd+'gcta64'):
        print('###### Installing gcta64')
        subprocess.call(['gsutil','cp','gs://nbaya/split/gcta/gcta64', wd+'gcta64'])
    if not (os.path.isfile(wd+f'gcta_{nsamples}k.grm.N.bin') & os.path.isfile(wd+f'gcta_{nsamples}k.grm.bin') & os.path.isfile(wd+f'gcta_{nsamples}k.grm.id')):
        print(f'###### Downloading GRM files for {nsamples}k individuals from UKB')
        subprocess.call(['gsutil','-m','cp',f'gs://nbaya/split/gcta/gcta_{nsamples}k.grm*', wd])
    for rep_id in rep_ids:
        try:
            is_complete = subprocess.check_output(['gsutil','ls',f'gs://nbaya/split/gcta/gcta_{nsamples}k.s{rep_id}.HEreg']) != None
        except:
            is_complete = False
        if not is_complete:            
            if not os.path.isfile(wd+f'gcta_{nsamples}k.s{rep_id}.phen'):        
                print(f'###### Downloading phen file for {nsamples}k replicate #{rep_id}')
                subprocess.call(['gsutil','cp',f'gs://nbaya/split/gcta/gcta_{nsamples}k.s{rep_id}.phen',wd])
            os.chdir(wd)
            subprocess.call(['chmod', '+x','gcta64'])
            out=f'gcta_{nsamples}k.s{rep_id}'
            subprocess.call(['./gcta64',
                             '--HEreg-bivar',
                             '--grm',f'gcta_{nsamples}k',
                             '--pheno',f'gcta_{nsamples}k.s{rep_id}.phen', 
                             '--out',out,
                             '--thread-num', '10']) #change based on number of cores available in VM
            subprocess.call(['gsutil','cp',out+'.HEreg','gs://nbaya/split/gcta/'])
        else:
            print(f'###### Already completed GCTA rg calculation for {nsamples}k replicate #{rep_id}')

def gwas(mt, x, y, cov_list=[], with_intercept=True, pass_through=[], path_to_save=None, 
         normalize_x=False, is_std_cov_list=False):
    '''Runs GWAS in Hail'''
    
    mt = mt._annotate_all(col_exprs={'__y':y},
                           entry_exprs={'__x':x})
    if normalize_x:
        mt = mt.annotate_rows(__gt_stats = hl.agg.stats(mt.__x))
        mt = mt.annotate_entries(__x= (mt.__x-mt.__gt_stats.mean)/mt.__gt_stats.stdev) 
        mt = mt.drop('__gt_stats')
    
    if is_std_cov_list:
        cov_list = ['isFemale','age','age_squared','age_isFemale',
                    'age_squared_isFemale']+['PC{:}'.format(i) for i in range(1, 21)]
        
    if str in list(map(lambda x: type(x),cov_list)):
        cov_list = list(map(lambda x: mt[x] if type(x) is str else x,cov_list))
        
    cov_list = ([1] if with_intercept else [])+cov_list
    
    pass_through = list(set(['rsid']+pass_through))
    print(f'variables to pass through: {pass_through}')

    gwas_ht = hl.linear_regression_rows(y=mt.__y,
                                        x=mt.__x,
                                        covariates=cov_list,
                                        pass_through = pass_through)
    
    gwas_ht = gwas_ht.annotate_globals(with_intercept = with_intercept)
        
    gwas_ht = gwas_ht.rename({'rsid':'SNP'}).key_by('SNP')
        
    gwas_ht = gwas_ht.select(Z = gwas_ht.t_stat,
                             N = gwas_ht.n)
    
    ss_template = hl.read_table('gs://nbaya/rg_sex/hm3.sumstats_template.ht') # sumstats template as a hail table
    ss_template = ss_template.key_by('SNP')
        
    ss = ss_template.annotate(Z = gwas_ht[ss_template.SNP].Z,
                              N = gwas_ht[ss_template.SNP].N)
    
    if path_to_save is not None:
        ss.export(path_to_save)
        
    return ss

def gwas_on_gcta_splits(min_id,max_id,parsplit,paridx):
    rep_ids = range(min_id+paridx-1,max_id+1,parsplit) #replicate "IDs", which were used as seeds to generate the random split
    ids = hl.import_table('gs://nbaya/split/gcta/gcta_20k.grm.id',no_header=True) #GRM ids
    ids = ids.rename({'f0':'FID','f1':'IID'})
    ids = set(ids.IID.take(ids.count()))
    mt0 = hl.read_matrix_table('gs://nbaya/ldscsim/hm3.50_sim_h2_0.08.mt/')
    mt1 = mt0.filter_cols(hl.literal(ids).contains(mt0.s))
    for i in rep_ids:
        try:
            y1_complete = subprocess.check_output(['gsutil','ls',f'gs://nbaya/split/gcta/20k_sumstats.y1.s{i}.tsv.bgz']) != None
        except:
            y1_complete = False
        try:
            y2_complete = subprocess.check_output(['gsutil','ls',f'gs://nbaya/split/gcta/20k_sumstats.y2.s{i}.tsv.bgz']) != None
        except:
            y2_complete = False
        if not (y1_complete and y2_complete):           
            phen = hl.import_table(f'gs://nbaya/split/gcta/gcta_20k.s{i}.phen',types={'y1':hl.tfloat64,'y2':hl.tfloat64},key='IID')
            mt = mt1.annotate_cols(y1 = phen[mt1.s].y1,
                                   y2 = phen[mt1.s].y2)
            if not y1_complete:
                print(f'\r##########\nRunning GWAS for y1 replicate {i} \n##########')
                gwas(mt=mt,x=mt.dosage,y=mt.y1,is_std_cov_list=True,path_to_save=f'gs://nbaya/split/gcta/20k_sumstats.y1.s{i}.tsv.bgz')
            if not y2_complete:
                print(f'\r##########\nRunning GWAS for y2 replicate {i} \n##########')
                gwas(mt=mt,x=mt.dosage,y=mt.y2,is_std_cov_list=True,path_to_save=f'gs://nbaya/split/gcta/20k_sumstats.y2.s{i}.tsv.bgz')
        else:
            print(f'\r##########\n GWAS already complete for y1 and y2 for replicate {i} \n##########')
                

    
def combine_rgs():
    '''
    Combine all GCTA rg results for splits of the 20k individual dataset
    NOTE: Go to meta_split_supp.py for comparisons of these results to previous 
    random splits for other phenotypes
    '''
    ids = subprocess.check_output(['gsutil','ls',f'gs://nbaya/split/gcta/gcta.s*.hsq']).decode('utf-8')
    ids = ids.split('\n')[:-1]
    ids = [x.split('.')[-2].strip('s') for x in ids]
    hsq_ls = []
    for i in ids: # ids is the list of replicate IDs which have corresponding hsq files (i.e. they have already had rg calculated)
        hsq = subprocess.check_output(['gsutil','cat',f'gs://nbaya/split/gcta/gcta.s{i}.hsq'])
        hsq = hsq.decode('utf-8').split('\n')[1:-1]
        hsq = list(chain.from_iterable([[str(i)]]+list(map(lambda x: x.split('\t')[1:],hsq))))
        hsq_ls.append(hsq)
    df = pd.DataFrame(hsq_ls,columns=['rep_id','V(G)_tr1',  'V(G)_tr1_SE', 'V(G)_tr2','V(G)_tr2_SE',
                               'C(G)_tr12', 'C(G)_tr12_SE','V(e)_tr1','V(e)_tr1_SE',
                               'V(e)_tr2', 'V(e)_tr2_SE', 'Vp_tr1', 'Vp_tr1_SE',
                               'Vp_tr2', 'Vp_tr2_SE', 'V(G)/Vp_tr1','V(G)/Vp_tr1_SE',
                               'V(G)/Vp_tr2', 'V(G)/Vp_tr2_SE', 'rg', 'rg_SE', 'logL', 'n'],
                        dtype=np.float32)    
    df['rep_id'] = df.rep_id.astype(np.int32)
    df['n'] = df.n.astype(np.int32)
    hl.Table.from_pandas(df).export('gs://nbaya/split/gcta/gcta_20k.rgs.tsv.bgz')
    
def combine_he_rgs():
    '''
    Combine all GCTA HEreg rg results for splits of the 20k individual dataset
    NOTE: Go to meta_split_supp.py for comparisons of these results to previous 
    random splits for other phenotypes
    '''
    ids = subprocess.check_output(['gsutil','ls',f'gs://nbaya/split/gcta/gcta_20k.s*.HEreg']).decode('utf-8')
    ids = ids.split('\n')[:-1]
    ids = [x.split('.')[-2].strip('s') for x in ids]
    file_ls = []
    for i in ids: # ids is the list of replicate IDs which have corresponding hsq files (i.e. they have already had rg calculated)
        hereg0 = subprocess.check_output(['gsutil','cat',f'gs://nbaya/split/gcta/gcta_20k.s{i}.HEreg'])
        hereg = hereg0.decode('utf-8').split('\n')[3:-1]
        hereg = list(chain.from_iterable([[str(i)]]+list(map(lambda x: list(filter(None,x.split(' ')))[1:],hereg))))
        file_ls.append(hereg)
    df = pd.DataFrame(file_ls,columns=['rep_id','Intercept_tr1',  'Intercept_tr1_SE_OLS', 
                                       'Intercept_tr1_SE_Jackknife','Intercept_tr1_P_OLS',
                                       'Intercept_tr1_P_Jackknife', 'Intercept_tr2',  'Intercept_tr2_SE_OLS', 
                                       'Intercept_tr2_SE_Jackknife','Intercept_tr2_P_OLS',
                                       'Intercept_tr2_P_Jackknife', 'Intercept_tr12',  'Intercept_tr12_SE_OLS', 
                                       'Intercept_tr12_SE_Jackknife','Intercept_tr12_P_OLS',
                                       'Intercept_tr12_P_Jackknife', 'V(G)/Vp_tr1',  'V(G)/Vp_tr1_SE_OLS', 
                                       'V(G)/Vp_tr1_SE_Jackknife','V(G)/Vp_tr1_P_OLS',
                                       'V(G)/Vp_tr1_P_Jackknife', 'V(G)/Vp_tr2',  'V(G)/Vp_tr2_SE_OLS', 
                                       'V(G)/Vp_tr2_SE_Jackknife','V(G)/Vp_tr2_P_OLS',
                                       'V(G)/Vp_tr2_P_Jackknife','V(G)/Vp_tr12',  'V(G)/Vp_tr12_SE_OLS', 
                                       'V(G)/Vp_tr12_SE_Jackknife','V(G)/Vp_tr12_P_OLS',
                                       'V(G)/Vp_tr12_P_Jackknife','rg','rg_SE_OLS','rg_SE_Jackknife',
                                       'N_tr1','N_tr2'],dtype=np.float32)    
    df['rep_id'] = df.rep_id.astype(np.int32)
    df['N_tr1'] = df.N_tr1.astype(np.int32)
    df['N_tr2'] = df.N_tr2.astype(np.int32)
    hl.Table.from_pandas(df).export('gs://nbaya/split/gcta/gcta_he_20k.rgs.tsv.bgz')
    
    
if __name__ == "__main__":
#    make_grm(nsamples=100e3)
#    get_phen_files(nsamples=100e3,min_id=min_id,max_id=max_id,parsplit=parsplit,paridx=paridx)
#    run_gcta_rg(nsamples=20e3,min_id=min_id,max_id=max_id,parsplit=parsplit,paridx=paridx)
#    run_gcta_he_rg(nsamples=20e3,min_id=min_id,max_id=max_id,parsplit=parsplit,paridx=paridx)
    combine_he_rgs()
#    combine_rgs()
#    gwas_on_gcta_splits(min_id,max_id,parsplit,paridx)