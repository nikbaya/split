#!/usr/bin/env python3
"""
This file is the main pipeline for the chunking/meta-analysis method.

@author: nbaya
"""


import hail as hl
import numpy as np
import datetime
import argparse


phen_dict = {
    '50':'height',
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
parser.add_argument('--phen', type=str, required=True, help="phenotype code (e.g. for height, phen = 50")
parser.add_argument('--n_chunks', type=int, required=True, help="number of subgroups (or 'chunks'). Default: 300")
parser.add_argument('--batch', type=str, required=True, help="batch number for reproducibility of stochastic steps (used as seed number). Default: 1")
parser.add_argument('--variant_set', type=str, required=True, help="set of variants to use")

args = parser.parse_args()

      
"""
╔═══════════════════════╗
║ Part 1: Preprocessing ║
╚═══════════════════════╝
Preprocessing steps are split up to ensure that we can use preemptible workers
in preprocess1 and preprocess2.
"""

def preprocess1(variant_set):
    print('\n##################')
    print('Starting Pre-processing 1: Creating variants table (variant_set: '+variant_set+')')
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    print('\n##################')
    
    if variant_set == 'hm3':
        variants = hl.import_table('gs://nbaya/split/hapmap3_variants.tsv')
    elif variant_set == 'qc_pos':
        variants = hl.import_table('gs://ukb31063-mega-gwas/qc/ukb31063.gwas_variants.autosomes.tsv')
        
    variants = variants.annotate(**hl.parse_variant(variants.v))
    variants = variants.key_by('locus','alleles') 
    
    variants.write('gs://nbaya/split/'+variant_set+'_variants.ht')  

    print('\n##################')
    print('Finished Pre-processing 1: Creating variants table using variant_set: '+variant_set)
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    print('\n##################')
          
def preprocess2(variant_set):
    print('\n##################')
    print('Starting Pre-processing 2: Filtering variants table (variant_set: '+variant_set+')')
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    print('\n##################')
          
    variants = hl.read_table('gs://nbaya/split/'+variant_set+'_variants.ht') # for hm3: import table hapmap3_variants.tsv'
        
    mt = hl.read_matrix_table('gs://phenotype_31063/hail/imputed/ukb31063.dosage.autosomes.mt')
    mt = mt.filter_rows(hl.is_defined(variants[mt.locus,mt.alleles])) #filter to variant_set variants
    
    covs = hl.import_table('gs://phenotype_31063/ukb31063.gwas_covariates.both_sexes.tsv',
                                 key='s', impute=True, types={'s': hl.tstr})
    
    mt = mt.annotate_cols(**covs[mt.s])
    mt = mt.filter_cols(hl.is_defined(mt.PC1), keep=True)
    
    mt.write('gs://nbaya/split/ukb31063.'+variant_set+'_variants.gwas_samples_prerepart.mt')

    print('\n##################')
    print('Finished Pre-processing 2: Filtering variants table (variant_set: '+variant_set+')')
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    print('\n##################')
    
def preprocess3(variant_set):
    print('\n##################')
    print('Starting Pre-processing 3: Repartitioning (variant_set: '+variant_set+')')
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    print('\n##################')
          
    mt = hl.read_matrix_table('gs://nbaya/split/ukb31063.'+variant_set+'_variants.gwas_samples_prerepart.mt')
    if variant_set == 'hm3':
        mt = mt.repartition(1000)        
    elif variant_set == 'qc_pos':
        mt = mt.repartition(5000,shuffle=False) #qc_pos has more variants, thus requires more partitions
        
    
    #mt.write('gs://nbaya/split/ukb31063.hm3_variants.gwas_samples_v1.mt') #old version
    mt.write('gs://nbaya/split/ukb31063.'+variant_set+'_variants.gwas_samples_repart.mt',overwrite=True)
    
    print('\n##################')
    print('Finished Pre-processing 3: Repartitioning (variant_set: '+variant_set+')')
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    print('\n##################')
    


"""
╔══════════════════════════════════════════════════════════════════════╗
║ Part 2: Annotate with phenotype and label with n number of group IDs ║
╚══════════════════════════════════════════════════════════════════════╝
IMPORTANT: First differentiation by phenotype, batch, and n_chunks
Takes ~1h with 20 workers if writing matrix table. Otherwise it takes 2 min with
20 workers to write out the table.
"""
def get_phen_mt(variant_set,phen,batch,n_chunks,write):
    print('Starting Part 2: Splitting into n groups')
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    
    mt0 = hl.read_matrix_table('gs://nbaya/split/ukb31063.'+variant_set+'_variants.gwas_samples_repart.mt')

    if 'sim' in phen:
        print('\nReading simulated phenotype...')
        if variant_set == 'qc_pos':
            mt1 = hl.read_matrix_table('gs://nbaya/rg_sex/qc_pos.50_sim_inf_h2_0.485223.mt') #outdated
        elif variant_set == 'hm3':
#            mt1 = hl.read_matrix_table('gs://nbaya/rg_sex/50_sim_inf_h2_0.485223.mt')      
            sim_phen = phen.split('_sim')[0]
            if sim_phen == '50':
                phen_tb = hl.read_table('gs://nbaya/ldscsim/'+variant_set+'.phen_'+sim_phen+'.sim_h2_'+str(0.485223)+'.ht')
        
    else:
        print('\nReading UKB phenotype...')
#        mt0 = hl.read_matrix_table('gs://nbaya/split/ukb31063.hm3_variants.gwas_samples_v2.mt') #old version
        
        if phen == '50_raw':
            phen_tb0 = hl.import_table('gs://nbaya/ukb31063.50_raw.tsv.bgz',missing='',impute=True,types={'s': hl.tstr}).rename({phen: 'phen'})
        elif phen == '50_raw_res':
            phen_tb0 = hl.read_table('gs://nbaya/split/50_raw_linreg.ht').rename({'res': 'phen'})
        else:
            phen_tb0 = hl.import_table('gs://phenotype_31063/ukb31063.phesant_phenotypes.both_sexes.tsv.bgz',
                                          missing='',impute=True,types={'"userId"': hl.tstr}).rename({ '"userId"': 's', '"'+phen+'"': 'phen'})
    
        phen_tb = phen_tb0.select(phen_tb0['phen'])    

    mt1 = mt0.annotate_cols(phen_str = hl.str(phen_tb[mt0.s]['phen']).replace('\"',''))
    mt1 = mt1.filter_cols(mt1.phen_str == '',keep=False)
    
    if phen_tb.phen.dtype == hl.dtype('bool'):
        mt1 = mt1.annotate_cols(phen = hl.bool(mt1.phen_str)).drop('phen_str')
    else:
        mt1 = mt1.annotate_cols(phen = hl.float64(mt1.phen_str)).drop('phen_str')            
    
    #Remove withdrawn samples
    withdrawn = hl.import_table('gs://nbaya/w31063_20181016.csv',missing='',no_header=True)
    withdrawn_set = set(withdrawn.f0.take(withdrawn.count()))
    mt1 = mt1.filter_cols(hl.literal(withdrawn_set).contains(mt1['s']),keep=False) 
    mt1 = mt1.key_cols_by('s')

    n_samples = mt1.count_cols()
    print('\n>>> N samples = '+str(n_samples)+' <<<') #expect n samples to match n_non_missing from phenotypes.both_sexes.tsv, minus withdrawn samples.
    
    mt2 = mt1.add_col_index()
    group_size = int(n_samples/n_chunks)+1     #the ideal number of samples in each group
    #list of group ids to be paired to each sample (Note: length of group_ids > # of cols in mt, but it doesn't affect the result)
    group_ids = np.ndarray.tolist(np.ndarray.flatten(np.asarray([range(n_chunks)]*group_size))) 
    group_ids = group_ids[0:n_samples]
    randstate = np.random.RandomState(int(batch)) #seed with batch number
    randstate.shuffle(group_ids)
    mt3 = mt2.annotate_cols(group_id = hl.literal(group_ids)[hl.int32(mt2.col_idx)]) #assign group ids # OLD VERSION
    ht_group_ids = mt3.select_cols(mt3.group_id).cols() #assign group ids
    
    if write:
        print('Writing HailTable with group ids...')
    #    mt3.write('gs://nbaya/split/ukb31063.'+variant_set+'_variants.gwas_samples_'+phen+'_grouped'+str(n_chunks)+'_batch_'+batch+'.mt',overwrite=True) #Takes ~30 min with 50 workers OLD VERSION
        ht_group_ids.write('gs://nbaya/split/ukb31063.'+variant_set+'_variants.gwas_samples_'+phen+'_grouped'+str(n_chunks)+'_batch_'+batch+'.ht',overwrite=True) 

    print('Finished Part 2: Splitting into n groups')
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now())) #takes ~1h with 20 workers, 42 min with 30 workers
    
    return mt3

"""
╔══════════════════════════════════════════════╗
║ Part 3: Run linear regression for each group ║
╚══════════════════════════════════════════════╝
"""
def metasplit1(variant_set,phen,batch,n_chunks):
    print('Starting Part 3: Linear regression')
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    
    mt = get_phen_mt(variant_set,phen,batch,n_chunks,write=False)
    
    mt = mt.rename({'dosage': 'x', 'phen': 'y'})
    
    cov_list = [ mt['isFemale'], mt['age'], mt['age_squared'], mt['age_isFemale'],
                mt['age_squared_isFemale'] ]+ [mt['PC{:}'.format(i)] for i in range(1, 21)] 
               
    gmt = (mt.group_cols_by(mt.group_id)
            .aggregate(linreg = hl.agg.linreg(y=mt.y, x = [mt.x, 1] + cov_list)))
    
    gmt.select_entries(beta=gmt.linreg.beta[0],standard_error=gmt.linreg.standard_error[0],
                       beta_int=gmt.linreg.beta[1],standard_error_int=gmt.linreg.standard_error[1]).write('gs://nbaya/split/meta_split/ukb31063.'+variant_set+'_'+phen+'_gmt'+str(n_chunks)+'_batch_'+batch+'.mt',overwrite=True)
    
    print('Finished Part 4') #Takes 1.7 hours with 10 workers + 200 pre-emptible workers, around 4.5 hours with 50 workers, 4.55 hours with 10 workers + 40 pre-emptible
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))

"""
╔═════════════════════════════════════════════════════════╗
║ Part 4: Split and calculate meta-analyzed summary stats ║
╚═════════════════════════════════════════════════════════╝
Split n groups into two populations (A, B) and calculate meta summary statistics (via inverse-variance weighting meta-analysis)
"""
def metasplit2(variant_set, phen,batch,n_chunks):    
    print('Starting Part 4: Splitting into populations A and B, calculatng summary statistics')
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    
    #gmt = hl.read_matrix_table('gs://nbaya/split/meta_split/ukb31063.hm3_'+phen+'_gmt.mt') #old version
    #gmt = hl.read_matrix_table('gs://nbaya/split/meta_split/ukb31063.hm3_'+phen+'_gmt_batch_'+batch+'.mt') #used for n=300, phenotypes 50 and 20160
    gmt = hl.read_matrix_table('gs://nbaya/split/meta_split/ukb31063.'+variant_set+'_'+phen+'_gmt'+str(n_chunks)+'_batch_'+batch+'.mt') #used for n=150 and other phenotypes
    
    gmt = gmt.add_col_index()
    gmt = gmt.rename({'rsid': 'SNP'})
    
    def run_meta_split(i):
        print('####################')
        print('Starting split '+str(i))
        print('####################')
        starttime = datetime.datetime.now()
        pi = ['A']*int(n_chunks/2) + ['B']*int(n_chunks/2)
        #    randstate = np.random.RandomState(i) #OPTION 1: seed with split number (used for 20160 batch1)
        seed_id = int(batch+str(i).zfill(4)) #OPTION 2: create a seed_id unique to every split
        randstate = np.random.RandomState(seed_id) #OPTION 2: seed with seed_id
        
        randstate.shuffle(pi)
        gmt_shuf = gmt.annotate_cols(label = hl.literal(pi)[hl.int32(gmt.col_idx)])
        
        mt = gmt_shuf.group_cols_by(gmt_shuf.label).aggregate(unnorm_meta_beta=hl.agg.sum(gmt_shuf.beta / gmt_shuf.standard_error ** 2),
                                    inv_se2 = hl.agg.sum(1 / gmt_shuf.standard_error ** 2))
               
        ht = mt.make_table()

        ht = ht.annotate(A_Z = ht['A.unnorm_meta_beta'] / hl.sqrt(ht['A.inv_se2']),
                         B_Z = ht['B.unnorm_meta_beta'] / hl.sqrt(ht['B.inv_se2']))
        
        ht = ht.drop('A.unnorm_meta_beta','B.unnorm_meta_beta','A.inv_se2','B.inv_se2').key_by('SNP')
        
        variants = hl.import_table('gs://nbaya/rg_sex/50_snps_alleles_N.tsv.gz',types={'N': hl.tint64})
        variants = variants.key_by('SNP')
        mt_all = hl.read_matrix_table('gs://nbaya/split/ukb31063.hm3_variants.gwas_samples_'+phen+'_grouped'+str(n_chunks)+'_batch_'+batch+'.mt') #matrix table containing individual samples
        variants = variants.annotate(N = hl.int32(mt_all.count_cols()/2))
        variants.show()
        
        metaA = variants.annotate(Z = ht[variants.SNP].A_Z)
        metaB = variants.annotate(Z = ht[variants.SNP].B_Z)
        
        metaA_path = 'gs://nbaya/rg_sex/'+phen+'_meta_A_n'+str(n_chunks)+'_batch_'+batch+'_s'+str(i)+'.tsv.bgz' #used for n_chunks=150 and all other phenotypes
        metaB_path = 'gs://nbaya/rg_sex/'+phen+'_meta_B_n'+str(n_chunks)+'_batch_'+batch+'_s'+str(i)+'.tsv.bgz' #used for n_chunks=150 and all other phenotypes
        metaA.export(metaA_path)
        metaB.export(metaB_path)
        
        endtime = datetime.datetime.now()
        elapsed = endtime-starttime
        print('####################')
        print('Completed iteration '+str(i))
        print('Files written to:')
        print(metaA_path+'\t'+metaB_path)
        print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
        print('Iteration time: '+str(round(elapsed.seconds/60, 2))+' minutes')
        print('####################')    
    
    for i in range(100):
        run_meta_split(i)
    
    print('####################')
    print('Finished Part 5')
    print('Meta-split complete')
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
    print('####################')


"""
╔════════════════╗
║ Run meta_split ║
╚════════════════╝
Note: Be careful when running metasplit1 and metasplit2 on the same cluster. Sometimes the 
last step of metasplit2 will fail to finish.
"""
if __name__ == "__main__":    
    phen = args.phen
    desc = phen_dict[phen]
    batch = args.batch
    n_chunks = args.n_chunks #number of subgroups (or "chunks")
    variant_set = str(args.variant_set)
    
    print('####################')
    print('Phenotype: '+phen)
    print('Description: '+desc)
    print('Batch: '+batch)
    print('n chunks: '+str(n_chunks))
    print('variant set: '+variant_set)
    print('####################')
    
    #preprocess1(variant_set) #create variant table
    #preprocess2(variant_set) #filter variants
    #preprocess3(variant_set) #repartition
#    get_phen_mt(variant_set,phen,batch,n_chunks,write=True)  #Only used to save HailTable of phens. Can skip this step and start with metasplit1, which uses this function to get the matrix table of phenotypes. Takes 2 min with 20 n1-standard-8 workers.
    metasplit1(variant_set,phen,batch,n_chunks) #typically takes 40 n1-standard-8 workers about 5 hours
    #metasplit2(variant_set, phen,batch,n_chunks) #use meta_split_parallel.py to run in hard parallel across clusters
