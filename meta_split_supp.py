# -*- coding: utf-8 -*-
"""
Contains supplementary code to the chunking/meta-analysis pipelines for rg and 
h2 calculations.

Use for reading in data, generating plots and calculating statistics.
"""


# import packages
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import  os


###############################################################################
"""
Non-meta-analysis version of the split
"""

#Read h2part results file for height iterations 0-6 (NON-META-ANALYSIS VERSION)


for v in list(map(str,range(7))):
    if not os.path.isfile('/Users/nbaya/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.phesant.batch_1_v'+v+'.tsv.gz'):
        os.system('gsutil cp gs://nbaya/split/test/ukbb31063.h2part_results.phesant.batch_1_v'+v+'.tsv.gz ~/Documents/lab/ukbb-sexdiff/h2part/')

h2_split_50 = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.phesant.batch_1_v0.tsv.gz',compression='gzip',sep='\t').iloc[:,0:20]

for v in list(map(str,range(1,7))):
    temp = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.phesant.batch_1_v'+v+'.tsv.gz',compression='gzip',sep='\t').iloc[:,0:20]
    h2_split_50 = h2_split_50.append(temp,ignore_index=True)
    
np.mean(h2_split_50.h2_observed)
np.std(h2_split_50.h2_observed)

#Read rg results file for height iterations 0-6 (NON-META-ANALYSIS VERSION)
for v in list(map(str,range(7))): #read in replicates
    if not os.path.isfile('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/ukbb31063.phesant.rg_sex.batch_1_v'+v+'.tsv.gz'):
        os.system('gsutil cp gs://nbaya/rg_sex/batches/ukbb31063.phesant.rg_sex.batch_1_v'+v+'.tsv.gz ~/Documents/lab/ukbb-sexdiff/rg_sex/')
        
rg_split_50 = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/ukbb31063.phesant.rg_sex.batch_1_v0.tsv.gz',compression='gzip',sep='\t').iloc[:,0:20]

for v in list(map(str,range(1,7))): #read in replicates
    temp = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/ukbb31063.phesant.rg_sex.batch_1_v'+v+'.tsv.gz',compression='gzip',sep='\t').iloc[:,0:20]
    rg_split_50 = rg_split_50.append(temp,ignore_index=True)

np.mean(rg_split_50.rg)
np.std(rg_split_50.rg)

#Read h2part results file for smoking iterations 0-4 (NON-META-ANALYSIS VERSION)

for v in list(map(str,range(5))):
    if not os.path.isfile('/Users/nbaya/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.20160_smoking.batch_1_v'+v+'.tsv.gz'):
        os.system('gsutil cp gs://nbaya/split/test/ukbb31063.h2part_results.20160_smoking.batch_1_v'+v+'.tsv.gz ~/Documents/lab/ukbb-sexdiff/h2part/')

h2_split_20160 = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.20160_smoking.batch_1_v'+v+'.tsv.gz',compression='gzip',sep='\t').iloc[:,0:20]

for v in list(map(str,range(1,5))):
    temp = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.20160_smoking.batch_1_v'+v+'.tsv.gz',compression='gzip',sep='\t').iloc[:,0:20]
    h2_split_20160 = h2_split_20160.append(temp,ignore_index=True)
    
np.mean(h2_split_20160.h2_observed)
np.std(h2_split_20160.h2_observed)

#Read rg results file for smoking iterations 0-4 (NON-META-ANALYSIS VERSION)
for v in list(map(str,range(5))):
    if not os.path.isfile('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/ukbb31063.rg_sex.20160_smoking.batch_1_s'+v+'.tsv.gz'):
        os.system('gsutil cp gs://nbaya/rg_sex/batches/ukbb31063.rg_sex.20160_smoking.batch_1_s'+v+'.tsv.gz ~/Documents/lab/ukbb-sexdiff/rg_sex/')
        
rg_split_20160 = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/ukbb31063.rg_sex.20160_smoking.batch_1_s0.tsv.gz',compression='gzip',sep='\t').iloc[:,0:20]

for v in list(map(str,range(1,5))):
    temp = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/ukbb31063.rg_sex.20160_smoking.batch_1_s'+v+'.tsv.gz',compression='gzip',sep='\t').iloc[:,0:20]
    rg_split_20160 = rg_split_20160.append(temp,ignore_index=True)

print(np.mean(rg_split_20160.rg))
print(np.std(rg_split_20160.rg))
print(stats.ttest_1samp(rg_split_20160.rg, 1).pvalue/2)

#get h2_ref
phen='50'
h2 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/rg_sex/ukbb31063.both_sexes.h2part_results.phesant.tsv.gz', 
                              sep='\t',compression='gzip').rename(index=str,columns={'phenotype':'phen'}).iloc[:,0:20]
h2_ref = h2[h2['phen']==phen+'_irnt'].h2_observed[0] # full data ref
#Stats

stats.ttest_1samp(h2_split_50.h2_observed,0.485)
stats.ttest_1samp(h2_50.iloc[34].filter(regex=('observed')),0.485)

stats.ttest_1samp(rg_split_20160.rg,1)

stats.ttest_1samp(rg_split_50.rg,1)
stats.ttest_1samp(rg_split_20160.rg,1)


#Plots

sns.kdeplot(rg_split_20160.rg)
sns.kdeplot(rg_split_50.rg)
plt.legend(['20160 rg','50 rg'])
plt.title('Comparison of rg distr. \n 20160 vs. 50')
fig = plt.gcf()
fig.set_size_inches(6,4)
fig.savefig('/Users/nbaya/Desktop/kde_50vs20160_non_meta_version.png',dpi=300)

rg_split_50.loc[:,'description'] = '50_irnt'
rg_split_20160.loc[:,'description'] = '20160'
rg_split_temp = rg_split_50.append(rg_split_20160)
#rg_split_temp['description'].astype(str)

sns.violinplot(x = 'rg', y= 'description', data = rg_split_temp,palette='Set3')
plt.plot([1, 1], [-12, 12], '--k')
plt.ylabel('phenotype')
fig = plt.gcf()
fig.set_size_inches(6*1.5, 4*1.5)
fig.savefig('/Users/nbaya/Desktop/rg_dist.png',dpi=600)

###############################################################################
"""
Initial versions of meta-analysis
"""

import pandas as pd
import os

meta = pd.read_csv('~/Desktop/ukb31063.hm3_20160_metasplit_final_s1.tsv.bgz', 
                   sep = '\t', compression = 'gzip')
meta = meta.sort_values(by='SNP')

old_metaA = pd.read_csv('~/Documents/lab/ukbb-sexdiff/split/20160_vdsA_s1.tsv.gz'
                       , sep = '\t', compression = 'gzip')
old_metaA = old_metaA.sort_values(by='SNP')
old_metaA = old_metaA.reset_index(drop = True)

old_metaB = pd.read_csv('~/Documents/lab/ukbb-sexdiff/split/20160_vdsB_s1.tsv.gz'
                       , sep = '\t', compression = 'gzip')
old_metaB = old_metaB.sort_values(by='SNP')
old_metaB = old_metaB.reset_index(drop = True)


#os.system('gsutil cp gs://nbaya/split/meta_split/ukb31063.hm3_20160_metasplit_final_s1_test.tsv.bgz ~/Desktop/')

test_meta = pd.read_csv('~/Desktop/ukb31063.hm3_20160_metasplit_final_s1_test2.tsv.bgz',
                        sep = '\t', compression = 'gzip')
test_meta = test_meta.sort_values(by='SNP')
test_meta['A_Z_check'] = test_meta['A_meta_beta']/test_meta['A_meta_se']
test_meta['B_Z_check'] = test_meta['B_meta_beta']/test_meta['B_meta_se']
test_meta = test_meta.reset_index(drop=True)


test_metaA = old_metaA[['SNP','A1','A2','N']].copy()
test_metaA['Z'] = test_meta['A_Z'].copy()

test_metaB = old_metaA[['SNP','A1','A2','N']].copy()
test_metaB['Z'] = test_meta['B_Z'].copy()


test_metaA.to_csv('~/Documents/lab/ukbb-sexdiff/split/20160_meta_A_s1.tsv.gz', sep = '\t', compression = 'gzip', index = False)
test_metaB.to_csv('~/Documents/lab/ukbb-sexdiff/split/20160_meta_B_s1.tsv.gz', sep = '\t', compression = 'gzip', index = False)

os.system('gsutil cp ~/Documents/lab/ukbb-sexdiff/split/20160_meta_**_s1.tsv.gz gs://nbaya/rg_sex/')

###############################################################################
import pandas as pd

iters = 100
batch = '2'
phen = 50

phen_summary = pd.read_csv('~/Documents/lab/ukbb-sexdiff/rg_sex/'+str(phen)+'.h2part.batch_'+batch+'.tsv',sep='\t')

###############################################################################
"""
Compare cost of conventional to meta-analysis method
"""
x = np.linspace(0,200,201)
y_c = (39.368+36.82133)/2*x
y_m = (656.853+681.7584)/2+ 7.98*x #intercept is mean between cost of `50_raw` and `30100`, both with 300 chunks. slope comes from `30100`.
plt.plot(x, y_c)
plt.plot(x, y_m)
plt.legend(['conventional', 'meta-analysis'])
plt.title('Conventional vs. Meta-analysis')
plt.xlim([0, 200])
plt.ylim([0, np.max(y_c)])
plt.xlabel('replicates')
plt.ylabel('cost ($)')
fig = plt.gcf()
fig.set_size_inches(6,4)
fig.savefig('/Users/nbaya/Desktop/cost_comparison.png',dpi=1000)


###############################################################################
"""
Create snps_alleles_N file for meta_split.py method
"""
snps_alleles = old_metaA[['SNP','A1','A2','N']].copy()
snps_alleles.to_csv('~/Documents/lab/ukbb-sexdiff/split/50_snps_alleles_N.tsv.gz',
                    sep = '\t', compression = 'gzip', index = False)
os.system('gsutil cp ~/Documents/lab/ukbb-sexdiff/split/50_snps_alleles_N.tsv.gz gs://nbaya/rg_sex/')

###############################################################################
"""
Create h2part file for rg calculation (aka "phenfile" in the rg ldsc file)
"""

import pandas as pd
import os
phen = '50_sim_inf'
variant_set = 'hm3'
batch = 1
n_chunks = 300

phenotypes = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/imputed-v3-results/phenotypes.both_sexes.tsv',sep='\t')
heritable_phens = pd.read_csv('/Users/nbaya/Downloads/ukb_phenos_for_sex_rg.tsv',sep='\t')

#vals = phenotypes.loc[phenotypes['phenotype'] == str(phen)]

tb1 = pd.DataFrame(np.zeros(shape=[100,10]),columns=['phen','female_file','male_file','desc','female_n','male_n','female_n_cas','male_n_cas','female_n_con','male_n_con'])

for row in range(100):
    tb1.loc[row,'female_file'] = variant_set+'_'+phen+'_meta_A_n'+str(n_chunks)+'_batch_'+str(batch)+'_s'+str(row)+'.tsv.bgz'
    tb1.loc[row,'male_file'] = variant_set+'_'+phen+'_meta_B_n'+str(n_chunks)+'_batch_'+str(batch)+'_s'+str(row)+'.tsv.bgz'
    tb1.loc[:,'desc'] = variant_set+'_'+phen+'_s'+str(row)
    tb1.loc[row,'female_n'] = int(360338/2)
    tb1.loc[row,'male_n'] = int(360338/2)
    tb1.loc[row,'female_n_cas'] = float('NaN')
    tb1.loc[row,'male_n_cas'] = float('NaN')
    tb1.loc[row,'female_n_con'] = float('NaN')
    tb1.loc[row,'male_n_con'] = float('NaN')

tb1.loc[:,'phen'] = variant_set+'_'+phen

tb1 = tb1[['phen','female_file','male_file','desc','female_n','male_n','female_n_cas',
           'male_n_cas','female_n_con','male_n_con']]

filename = variant_set+'.'+phen+'.h2part.nchunks'+str(n_chunks)+'_batch_'+str(batch)+'.tsv'
local_wd = '~/Documents/lab/ukbb-sexdiff/rg_sex/'
cloud_wd = 'gs://nbaya/rg_sex/'
tb1.to_csv(local_wd+filename, sep='\t',index=False)

os.system('gsutil cp '+local_wd+filename+' '+cloud_wd)

###############################################################################
"""
Create h2part file for rg calculation (aka "phenfile" in the rg ldsc file)
"""

import pandas as pd
import os
phen = 'heritable_phens'
variant_set = 'hm3'
batch = 1
n_chunks = 2

phenotypes = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/imputed-v3-results/phenotypes.both_sexes.tsv',sep='\t')
heritable_phens = pd.read_csv('/Users/nbaya/Downloads/ukb_phenos_for_sex_rg.tsv',sep='\t')

#vals = phenotypes.loc[phenotypes['phenotype'] == str(phen)]

tb1 = heritable_phens
tb1['phen'] = heritable_phens['phenotype'].str.strip('_irnt')
tb1['female_file'] = tb1['phen']+'_meta_A_nchunks2_batch1_split.tsv.bgz'
tb1['male_file'] = tb1['phen']+'_meta_B_nchunks2_batch1_split.tsv.bgz'
tb1.loc[:,'desc'] = heritable_phens.loc[:,'description']
tb1['female_n'] = heritable_phens['n']/2 
tb1['male_n'] = heritable_phens['n']/2 
tb1['female_n_cas'] = heritable_phens['n_cases']/2 
tb1['male_n_cas'] = heritable_phens['n_cases']/2 
tb1['female_n_con'] = heritable_phens['n_controls']/2 
tb1['male_n_con'] = heritable_phens['n_controls']/2 

tb1 = tb1[['female_file','male_file','phen','desc','female_n','male_n','female_n_cas',
           'male_n_cas','female_n_con','male_n_con']]

filename = variant_set+'.'+phen+'.h2part.nchunks'+str(n_chunks)+'_batch_'+str(batch)+'.tsv'
local_wd = '~/Documents/lab/ukbb-sexdiff/rg_sex/'
cloud_wd = 'gs://nbaya/rg_sex/allphens/'
tb1.to_csv(local_wd+filename, sep='\t',index=False)

os.system('gsutil cp '+local_wd+filename+' '+cloud_wd)

###############################################################################
"""
Create h2part file for rg calculation (SPECIFICALLY FOR SEX STRAT SUMSTATS)
"""

import pandas as pd
import os
import numpy as np
phen_set = 'sex_strat'
variant_set = 'hm3'
batch = 1
n_chunks = 2

phenotypes = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/imputed-v3-results/phenotypes.both_sexes.tsv',sep='\t')

phens = os.popen('gsutil ls gs://nbaya/rg_sex/sex_strat/*nchunks*').read().split('\n')[:-1]
phens_files = [x.split('/')[5] for x in phens]
phens_set = sorted(list(set(['_'.join(x.split('_')[:-5]) for x in phens_files])))

df = pd.DataFrame(np.zeros(shape=[6*6,9]),columns=['female_file','male_file','desc','female_n','male_n','female_n_cas','male_n_cas','female_n_con','male_n_con'])

for i,phen in enumerate(phens_set):
    df.loc[i*6:i*6+6,'phen'] = phen

    df.loc[i*6,'female_file'] = phen+'_female_A_nchunks2_batch1_split.tsv.bgz' #phen+'_female_B_nchunks2_batch1_split.tsv.bgz' #
    df.loc[i*6,'male_file'] = phen+'_female_B_nchunks2_batch1_split.tsv.bgz'
    
    df.loc[i*6+1,'female_file'] = phen+'_male_A_nchunks2_batch1_split.tsv.bgz'
    df.loc[i*6+1,'male_file'] = phen+'_male_B_nchunks2_batch1_split.tsv.bgz'
    
    df.loc[i*6+2,'female_file'] = phen+'_female_A_nchunks2_batch1_split.tsv.bgz' #phen+'_female_B_nchunks2_batch1_split.tsv.bgz' #phen+'_female_A_nchunks2_batch1_split.tsv.bgz'
    df.loc[i*6+2,'male_file'] = phen+'_male_A_nchunks2_batch1_split.tsv.bgz'
    
    df.loc[i*6+3,'female_file'] = phen+'_female_A_nchunks2_batch1_split.tsv.bgz' #phen+'_female_B_nchunks2_batch1_split.tsv.bgz'#phen+'_female_B_nchunks2_batch1_split.tsv.bgz'
    df.loc[i*6+3,'male_file'] = phen+'_male_B_nchunks2_batch1_split.tsv.bgz'
    
    df.loc[i*6+4,'female_file'] = phen+'_female_A_nchunks2_batch1_split.tsv.bgz'
    df.loc[i*6+4,'male_file'] = phen+'_male_B_nchunks2_batch1_split.tsv.bgz'
    
    df.loc[i*6+5,'female_file'] = phen+'_female_B_nchunks2_batch1_split.tsv.bgz'
    df.loc[i*6+5,'male_file'] = phen+'_male_A_nchunks2_batch1_split.tsv.bgz'
    
    df.loc[i*6:i*6+6,'female_n'] = int(heritable_phens.loc[heritable_phens.phenotype.str.strip('_irnt') == phen]['n'].tolist()[0]/4)
    df.loc[i*6:i*6+6,'male_n'] = int(heritable_phens.loc[heritable_phens.phenotype.str.strip('_irnt') == phen]['n'].tolist()[0]/4 )
    df.loc[i*6:i*6+6,'female_n_cas'] = heritable_phens.loc[heritable_phens.phenotype.str.strip('_irnt') == phen]['n_cases'].tolist()[0]/4 
    df.loc[i*6:i*6+6,'male_n_cas'] = heritable_phens.loc[heritable_phens.phenotype.str.strip('_irnt') == phen]['n_cases'].tolist()[0]/4 
    df.loc[i*6:i*6+6,'female_n_con'] = heritable_phens.loc[heritable_phens.phenotype.str.strip('_irnt') == phen]['n_controls'].tolist()[0]/4 
    df.loc[i*6:i*6+6,'male_n_con'] = heritable_phens.loc[heritable_phens.phenotype.str.strip('_irnt') == phen]['n_controls'].tolist()[0]/4 

df.loc[:,'desc'] = 'description'

tb1 = df

filename = variant_set+'.'+phen_set+'.h2part.nchunks'+str(n_chunks)+'_batch_'+str(batch)+'.tsv'
local_wd = '~/Documents/lab/ukbb-sexdiff/rg_sex/'
cloud_wd = 'gs://nbaya/rg_sex/allphens/'
tb1.to_csv(local_wd+filename, sep='\t',index=False)

os.system('gsutil cp '+local_wd+filename+' '+cloud_wd)

###############################################################################
"""
Create h2part file for rg calculation (specifically for blood pressure phenotypes)
"""

import pandas as pd
import os
import numpy as np
variant_set = 'hm3'

#phenotypes = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/imputed-v3-results/phenotypes.both_sexes.tsv',sep='\t')
#phens = os.popen('gsutil ls gs://nbaya/rg_sex/sex_strat/*nchunks*').read().split('\n')[:-1]
#phens_files = [x.split('/')[5] for x in phens]
#phens_set = sorted(list(set(['_'.join(x.split('_')[:-5]) for x in phens_files])))

phens_ls = ['4079','4080']

df = pd.DataFrame(np.zeros(shape=[len(phens_ls),10]),columns=['phen','female_file','male_file','desc','female_n','male_n','female_n_cas','male_n_cas','female_n_con','male_n_con'])



for i,phen in enumerate(phens_ls):
    df.loc[i,'phen'] = phen

    df.loc[i,'female_file'] = phen+'_female.tsv.bgz' 
    df.loc[i,'male_file'] = phen+'_male.tsv.bgz' 

    df.loc[i,'female_n_cas'] = float('NaN')
    df.loc[i,'male_n_cas'] = float('NaN')
    df.loc[i,'female_n_con'] = float('NaN')
    df.loc[i,'male_n_con'] = float('NaN')

    df.loc[i,'desc'] = phen+' filtered by BP meds'

df.loc[0,'female_n'] = 142710
df.loc[0,'male_n'] = 107511
df.loc[1,'female_n'] = 142710
df.loc[1,'male_n'] = 107511

tb1 = df

filename = 'blood_pressure_phens.sex_strat.h2part.tsv'
local_wd = '~/Documents/lab/ukbb-sexdiff/rg_sex/'
cloud_wd = 'gs://nbaya/rg_sex/sex_strat/'
tb1.to_csv(local_wd+filename, sep='\t',index=False)

os.system('gsutil cp '+local_wd+filename+' '+cloud_wd)


###############################################################################
"""
╔═════════════════╗
║ Read rg results ║
╚═════════════════╝
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# possible phenotypes: 50 (irnt), 20160, 50 (n_chunks = 150), 50_raw, 30100, 2443, 23107, 6138_100
no_intercept = True
no_gcov_int = False
all_filenames = ['ukbb31063.rg_sex.50.batch_1.tsv.gz', 'ukbb31063.rg_sex.20160.batch_1.tsv.gz',
             'ukbb31063.rg_sex.50.n150_batch_1.tsv.gz','ukbb31063.rg_sex.50_raw.batch_1.tsv.gz',
             'ukbb31063.rg_sex.30100.batch_1.tsv.gz', 'ukbb31063.rg_sex.2443.batch_1.tsv.gz',
             'ukbb31063.rg_sex.23107.batch_1.tsv.gz', 'ukbb31063.rg_sex.6138_100.batch_1.tsv.gz',
             'ukbb31063.rg_sex.hm3.50_sim_inf.nchunks300.batch_1'+('.no_intercept' if no_intercept else '')+('.int_gcov_0' if no_gcov_int else '')+'.tsv.gz', 'ukbb31063.rg_sex.50_sim_inf.n150.batch_1.tsv.gz',
             'ukbb31063.rg_sex.50_sim_inf_h2_0.1.n300.batch_1.tsv.gz','ukbb31063.rg_sex.50_raw_res.n300.batch_1.tsv.gz',
             'ukbb31063.rg_sex.qc_pos.50_sim_inf.n300.batch_1.tsv.gz','gcta_20k.rg.tsv.bgz']
#origpheno_filenames = ['ukbb31063.rg_sex.50.batch_1.tsv.gz', 'ukbb31063.rg_sex.20160.batch_1.tsv.gz',
#             'ukbb31063.rg_sex.30100.batch_1.tsv.gz', 'ukbb31063.rg_sex.2443.batch_1.tsv.gz',
#             'ukbb31063.rg_sex.23107.batch_1.tsv.gz', 'ukbb31063.rg_sex.6138_100.batch_1.tsv.gz']
filenames = all_filenames
cloud_wd = 'gs://nbaya/rg_sex/batches/'
wd = '/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/'
for filename in filenames:
    if not os.path.isfile(wd+filename):
        print('Importing '+filename)
        os.system('gsutil cp '+cloud_wd+filename+' '+wd)
rg_df =  pd.read_csv(wd+filenames[0],sep='\t',compression='gzip').sort_values(by='p1').reset_index(drop=True)
for i in range(1,len(filenames)):
    if '20k' in filenames[i]:
        temp = pd.read_csv(wd+filenames[i],sep='\t',compression='gzip').rename(
                index=str,columns={'rg_SE':'se','rep_id':'phenotype','n':'ph1_n'}).reset_index(drop=True)
        temp['p1'] = 'gcta_20k_'+temp['phenotype'].astype(str)+'A'
        temp['p2'] = 'gcta_20k_'+temp['phenotype'].astype(str)+'B'
        temp['z'] = temp['rg']/temp['se']
        temp['ph2_n'] = temp['ph1_n']
        temp = temp.drop(columns=temp.columns.values[1:-8])
    else:
        temp = pd.read_csv(wd+filenames[i],sep='\t',compression='gzip').sort_values(by='p1').reset_index(drop=True)
    rg_df = rg_df.append(temp)
rg_df.loc[:,'z_not_1'] = (rg_df.loc[:,'rg']-1)/rg_df.loc[:,'se']

phen_ids = {#'20160':'20160',
            '50_meta_A_batch':'50_irnt',
#            '50_meta_A_n150':'50 n150',
#            '50_raw_meta':'50_raw',
#            '30100':'30100',
#            '2443':'2443',
#            '23107':'23107',
#            '6138_100':'6138_100',
#            'meta':'All phens',
            'hm3_50_sim_inf_meta_A_n300':'Inf. model, n300',
#            '50_sim_inf_meta_A_n150':'Inf. model, n150',
#            'inf_h2_0.1':'Inf. model, h2=0.1',
#            '50_raw_res':'50_raw_res',
#            'qc_pos_50':'Inf. model with QC variants',
            'gcta_20k':'gcta_20k'
            }

phen_ids_labels = {'20160':'20160',
            '50_meta_A_batch':'Height (50_irnt)',
            '50_meta_A_n150':'50 n150',
            '50_raw_meta':'50_raw',
            '30100':'30100',
            '2443':'2443',
            '23107':'23107',
            '6138_100':'6138_100',
#            'meta':'All phens',
            'hm3_50_sim_inf_meta_A_n300':'Inf. model, n300',
            '50_sim_inf_meta_A_n150':'Inf. model, n150',
            'inf_h2_0.1':'Inf. model, h2=0.1',
            '50_raw_res':'50_raw_res',
            'qc_pos_50':'Inf. model with QC+ variants',
            'gcta_20k':'gcta_20k'
            }

col='rg'
for phen_id, desc in phen_ids.items():
    rg_df.loc[rg_df.p1.str.contains(phen_id),'description'] = desc

###### Code for generating violin plot for poster ######
rg_df_temp = rg_df[rg_df['description'].isin(list(phen_ids.values()))].iloc[:,0:20].append(rg_split_temp)
rg_df_temp['method'] = 'conventional'
rg_df_temp.loc[rg_df_temp.p1.str.contains('meta'),'method'] = 'meta-analysis'
rg_df_temp.loc[rg_df_temp.description.str.contains('50_irnt'),'description'] = '50_irnt\n"Standing height"'
rg_df_temp.loc[rg_df_temp.description.str.contains('20160'),'description'] = '20160\n"Ever smoked"'

sns.violinplot(y = 'description',x='rg', data=rg_df_temp.iloc[::-1], hue='method', split=True, scale='width')
plt.plot([1, 1], [-12, 12], '--k')
plt.legend(loc=4)
plt.ylabel('phenotype')
plt.title('Conventional vs. Meta-analysis')#\nComparison of rg distributions')
plt.tight_layout()
fig = plt.gcf()
fig.set_size_inches(6, 4)
fig.savefig('/Users/nbaya/Desktop/rg_dist.png',dpi=2000)

sns.violinplot(x='description', y='rg', data=rg_df[rg_df['description'].isin(list(phen_ids.values()))], palette='Set3')
plt.ylabel('phenotype')
plt.title('Conventional vs. Meta-analysis\nComparison of rg distributions')
fig = plt.gcf()
fig.set_size_inches(6, 4)
####

###### Code for generating kde plots for all selected phenotypes ######
col='rg'
for phen_id, desc in phen_ids.items():
#    print(rg_df[(rg_df.p1.str.contains(phen_id))][col].shape)
    print(rg_df[(rg_df.description.str.contains(desc))][col].shape)
#    sns.kdeplot(rg_df[rg_df.p1.str.contains(phen_id)][col])
    sns.kdeplot(rg_df[rg_df.description.str.contains(desc)][col])
#    plt.hist(rg_df[rg_df.p1.str.contains(phen_id)][col],alpha=0.5)

plt.legend([value for key, value in phen_ids_labels.items()])
plt.title(col+' distributions')
plt.xlabel(col)
plt.ylabel('density')
plt.ylim([0,90])
fig = plt.gcf()
fig.set_size_inches(6*1.2, 4*1.2)
fig.savefig('/Users/nbaya/Desktop/50_sim_phens_'+col+'.png',dpi=600)

col='h2_int'
col1='ph1_'+col
for phen_id, desc in phen_ids.items():
    sns.kdeplot(rg_df[rg_df.p1.str.contains(phen_id)][col1])
col2='ph2_'+col
for i, phen_id in enumerate(list(phen_ids.keys())):
    sns.kdeplot(rg_df[rg_df.p1.str.contains(phen_id)][col2], c=plt.rcParams['axes.prop_cycle'].by_key()['color'][i],linestyle='--')
#for i in range(4):
#    phen_id = list(phen_ids.keys())[i]
#    mean = np.mean(rg_df[rg_df.p1.str.contains(phen_id)].rg)
#    plt.plot([mean,mean],[0.1,65],c=plt.rcParams['axes.prop_cycle'].by_key()['color'][i],linestyle='--',linewidth=1)
plt.legend([value for key, value in phen_ids.items()],loc=9)
plt.title(col+' distributions')
plt.xlabel(col)
plt.ylabel('density')
fig = plt.gcf()
fig.set_size_inches(6, 4)
fig.savefig('/Users/nbaya/Desktop/50_'+col+'.png',dpi=300)


sns.kdeplot(rg_df[rg_df.p1.str.contains('50_meta_A_batch')]['rg'],color='#ff7f0e',linestyle='-')
sns.kdeplot(rg_df[rg_df.p1.str.contains('30100')]['rg'],color='#2ca02c',linestyle='-')
plt.legend(['50 rg', '30100 rg'],loc=0)    
plt.title('Comparison of stock ldsc rg distr.\n50 vs. 30100')
plt.xlabel('rg')
plt.ylabel('density')
fig = plt.gcf()
fig.set_size_inches(6, 4)
fig.savefig('/Users/nbaya/Desktop/meta_split_50_30100_rg.png',dpi=300)


plt.title('h2_int'+' distributions')
plt.xlabel('h2_int')
plt.ylabel('density')
fig = plt.gcf()
fig.set_size_inches(6, 4)
fig.savefig('/Users/nbaya/Desktop/meta_split_height_h2_int.png',dpi=300)


###### Print statistics ######
col='rg'
for phen_id, desc in phen_ids.items():
    print(desc+' Mean: '+str(np.mean(rg_df[rg_df.p1.str.contains(phen_id)][col])))
    print(desc+' Std: '+str(np.std(rg_df[rg_df.p1.str.contains(phen_id)][col])))
    print(desc+' t-test 1samp=1 pvalue: '+str(stats.ttest_1samp(rg_df[rg_df.p1.str.contains(phen_id)][col],1).pvalue/2))
    print('\r')
#    print(desc+' pearson r: '+str(stats.pearsonr(rg_df[rg_df.p1.str.contains(phen_id)]['rg'],rg_df[rg_df.p1.str.contains(phen_id)]['gcov_int'])))


###### proportion of replicates that have a 95% confidence interval that contains 1 ######
for phen_id, desc in phen_ids.items():    
    s = rg_df[rg_df.p1.str.contains(phen_id)][['rg','se']]
    print(phen_id)
    print(np.mean(abs(s['rg']-1) < 2*s['se']))
    print('\r')
#    meta_rg = (np.sum(s['rg']/s['se']**2)/np.sum(1/s['se']**2))
#    meta_se = (1/np.sum(1/s['se']**2))
#    sns.kdeplot(np.random.normal(meta_rg, meta_se, 1000000))

s = rg_df[rg_df.p1.str.contains('hm3_50_sim_inf_meta_A_n300')]
s = rg_df[rg_df.p1.str.contains('50_sim_inf_meta_A_n150')]
stats.ttest_1samp(s['rg'],1)
sns.kdeplot(s['rg'])

sns.kdeplot(s.ph1_h2_obs)
sns.kdeplot(s.ph2_h2_obs)


###### plot rg vs. rg_se ######
fig,ax = plt.subplots(figsize=(8,6))
for phen_id, desc in phen_ids.items():
    rg_temp = rg_df[rg_df.p1.str.contains(phen_id)]['rg']
    se_temp = rg_df[rg_df.p1.str.contains(phen_id)]['se']
    print(f'{desc}\ncorr: {stats.pearsonr(rg_temp,se_temp)[0]}')
    print(f'p-val: {stats.ttest_1samp(rg_temp,popmean=1)[1]}\n')
    ax.plot(rg_temp,se_temp,'.')
plt.legend([value for key, value in phen_ids.items()])
plt.ylabel('se')
plt.xlabel('rg')
plt.title('rg vs rg_se')
fig.set_size_inches(8,6)
fig.savefig('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/plots/rg_vs_rg_se_realphens.png',dpi=600)
    
###### plot rg's w/ error bars ####
rg_df_temp = rg_df
rg_df_temp['is_used'] = False
for phen_id, desc in phen_ids.items():
    rg_df_temp.loc[rg_df_temp.p1.str.contains(phen_id),'is_used'] = True
rg_df_temp = rg_df_temp[rg_df_temp.is_used==True]
rg_df_temp = rg_df_temp.reset_index()
fig,ax = plt.subplots(figsize=(12,12))
#ax.plot([1,1],[0,rg_df_temp.shape[0]],'k--',alpha=0.5)
ax.plot([0,0],[0,rg_df_temp.shape[0]],'k--',alpha=0.5)
for phen_id, desc in phen_ids.items():
    df_temp = rg_df_temp[rg_df_temp.p1.str.contains(phen_id)]
#    df_temp = df_temp.sort_values(by='rg')
    df_temp = df_temp.sort_values(by='z_not_1')
    min_idx = np.min(df_temp.index.values)
    df_temp = df_temp.reset_index()
#    ax.errorbar(df_temp.rg,df_temp.index+min_idx,xerr=2*df_temp.se,fmt='.')
    ax.errorbar(df_temp.z_not_1,df_temp.index+min_idx,xerr=2,fmt='.')
#ax.legend(['rg=1']+[value for key, value in phen_ids.items()])
ax.legend(['rg=1 (z_not_1=0)']+[value for key, value in phen_ids.items()])
plt.xlabel('rg estimate')
plt.xlabel('z_not_1')
plt.ylabel('replicates')
#plt.xlim([0.9,1.1])
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft=False) # labels along the bottom edge are off
plt.title('rg estimate')
plt.title('z_not_1')
#plt.tight_layout()
fig = plt.gcf()
#fig.set_size_inches(12,9)
fig.savefig('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/plots/test.png',dpi=400)

################################################################################
'''
Code for reading GCTA 20k results
'''
df = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/gcta_20k.rg.tsv.bgz',delimiter='\t',compression='gzip')
ct = 0
for i in range(df.shape[0]):
    if not abs(df.rg[i]-1) < 2*df.rg_SE[i]:
        print(f'rg: {df.rg[i]}, se: {df.rg_SE[i]}')
        ct += 1
print(f'Proportion of rg estimates w/ 95% CI not containing 1: {ct/df.shape[0]}')

df = df.sort_values(by='rg')
fig,ax = plt.subplots(figsize=(8,8))
#    plt.plot(range(df.shape[0]),df.rg,'.')
#    plt.fill_between(range(df.shape[0]),df.rg-2*df.rg_SE,df.rg+2*df.rg_SE,alpha=0.5)
ax.plot([1,1],[0,df.shape[0]],'k--',alpha=0.5)
ax.errorbar(df.rg,range(df.shape[0]),xerr=2*df.rg_SE,fmt='.')
ax.legend(['rg=1 reference','rg estimates (95% CI)'])
plt.xlabel('rg estimate')
plt.ylabel('replicates')
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft=False) # labels along the bottom edge are off
plt.title('gcta rg estimates')
fig = plt.gcf()
fig.set_size_inches(8,8)
fig.savefig('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/plots/gcta_rg_estimates.png',dpi=600)

###############################################################################
"""
╔═════════════════════════════════════════════╗
║ Read meta-analysis rg results (old version) ║
╚═════════════════════════════════════════════╝
"""
rg_50 = rg_df[rg_df.p1.str.contains('50_meta_A_batch')]
rg_20160 = rg_df[rg_df.p1.str.contains('20160')]
rg_n150_50 = rg_df[rg_df.p1.str.contains('n150')]
rg_50_raw = rg_df[rg_df.p1.str.contains('50_raw')]

rg_temp = rg_50.append(rg_20160)

stats.ttest_rel(rg_50['rg'], rg_n150_50['rg'])
stats.ttest_1samp(rg_50['rg'],1)
stats.ttest_1samp(rg_n150_50['rg'],1)

np.mean(rg_20160['rg'])
np.std(rg_20160['rg'])

np.mean(rg_n150_50['rg'])
np.std(rg_n150_50['rg'])

np.mean(rg_20160['z_not_1'])
np.std(rg_20160['z_not_1'])

plt.subplot(1,2,1)
plt.hist(rg_20160['rg'])
plt.title('Histogram of rg values')
plt.xlabel('rg')
plt.ylabel('count')
fig = plt.gcf()
fig.set_size_inches(12, 8*.5)

plt.subplot(1,2,2)
plt.hist(rg_20160['z_not_1'])
plt.title('Histogram of z_not_1 values')
plt.xlabel('z_not_1')
plt.ylabel('count')
fig = plt.gcf()
fig.set_size_inches(12, 8*.5)
fig.savefig('/Users/nbaya/Desktop/meta_split_20160_hist.png',dpi=120)


np.mean(rg_50['rg'])
np.std(rg_50['rg'])
np.mean(rg_50['z_not_1'])
np.std(rg_50['z_not_1'])

plt.subplot(1,2,1)
plt.hist(rg_50['rg'])
plt.title('Histogram of rg values')
plt.xlabel('rg')
plt.ylabel('count')
fig = plt.gcf()
fig.set_size_inches(12, 8*.5)

plt.subplot(1,2,2)
plt.hist(rg_50['z_not_1'])
plt.title('Histogram of z_not_1 values')
plt.xlabel('z_not_1')
plt.ylabel('count')
fig = plt.gcf()
fig.set_size_inches(12, 8*.5)
fig.savefig('/Users/nbaya/Desktop/meta_split_50_hist.png',dpi=300)

np.mean(rg_n150_50['rg'])
np.std(rg_n150_50['rg'])
np.mean(rg_n150_50['z_not_1'])
np.std(rg_n150_50['z_not_1'])

plt.subplot(1,2,1)
plt.hist(rg_n150_50['rg'])
plt.title('Histogram of rg values')
plt.xlabel('rg')
plt.ylabel('count')
fig = plt.gcf()
fig.set_size_inches(12, 8*.5)

plt.subplot(1,2,2)
plt.hist(rg_n150_50['z_not_1'])
plt.title('Histogram of z_not_1 values')
plt.xlabel('z_not_1')
plt.ylabel('count')
fig = plt.gcf()
fig.set_size_inches(12, 8*.5)
fig.savefig('/Users/nbaya/Desktop/meta_split_50_n150_hist.png',dpi=300)




import seaborn as sns
plt.subplot(1,2,1)
#plot kde plot of 50 vs 20160 rg
sns.kdeplot(rg_20160['rg'])
sns.kdeplot(rg_50['rg'])
sns.kdeplot(rg_n150_50['rg'])
plt.legend(['20160 rg','50 rg (n=300)','50 rg (n=150)'])
#plt.legend(['20160 rg','50 rg'])
plt.xlabel('rg')
plt.ylabel('density')
plt.title('Comparison of rg distr. \n 20160 vs. 50 (n=300) vs. 50 (n=150)')
#plt.title('Comparison of rg distr. \n 20160 vs. 50')
fig = plt.gcf()
fig.set_size_inches(6, 4)
fig.savefig('/Users/nbaya/Desktop/kde_50vs20160vs50n150.png',dpi=300)

plt.subplot(1,2,2)
sns.kdeplot(rg_20160['z_not_1'])
sns.kdeplot(rg_50['z_not_1'])
sns.kdeplot(rg_n150_50['z_not_1'])
plt.legend(['20160 z_not_1','50 z_not_1 (n=300)','50 z_not_1 (n=150)'])
plt.title('Comparison of z_not_1 distr. \n 20160 vs. 50 (n=300) vs. 50 (n=150)')
plt.xlabel('z_not_1')
plt.ylabel('density')
fig = plt.gcf()
fig.set_size_inches(12, 8*.5)
fig.savefig('/Users/nbaya/Desktop/kde_50vs20160.png',dpi=300)


sns.kdeplot(rg_20160['z_not_1'])
sns.kdeplot(rg_50['z_not_1'])
plt.legend(['20160','50'])




# Calculate correlation between rg and gcov_int

from scipy import stats

stats.pearsonr(rg1['rg'],rg1['gcov_int'])
stats.pearsonr(rg2['rg'],rg2['gcov_int'])
stats.pearsonr(rg3['rg'],rg3['gcov_int'])

stats.spearmanr(rg1['rg'],rg1['gcov_int'])
stats.spearmanr(rg2['rg'],rg2['gcov_int'])
stats.spearmanr(rg3['rg'],rg3['gcov_int'])

stats.pearsonr(rg1['rg'],rg1['ph1_h2_int'])
stats.pearsonr(rg2['rg'],rg2['ph1_h2_int'])
stats.pearsonr(rg3['rg'],rg3['ph1_h2_int'])

stats.pearsonr(rg1['rg'],rg1['ph2_h2_int'])
stats.pearsonr(rg2['rg'],rg2['ph2_h2_int'])
stats.pearsonr(rg3['rg'],rg3['ph2_h2_int'])


###############################################################################
# Compare distributions of Z scores for halves of a meta split
metaA = pd.read_csv('~/Desktop/50_meta_A_s1.tsv.bgz', sep='\t',compression='gzip')
metaB = pd.read_csv('~/Desktop/50_meta_B_s1.tsv.bgz', sep='\t',compression='gzip')

oldsplitA = pd.read_csv('~/Documents/lab/ukbb-sexdiff/split/vds1_1.tsv.gz',compression='gzip',
                        sep='\t')
oldsplitB = pd.read_csv('~/Documents/lab/ukbb-sexdiff/split/vds2_1.tsv.gz',compression='gzip',
                        sep='\t')
np.mean(metaB['Z'])
np.mean(metaA['Z'])

np.mean(oldsplitA['Z'])
np.mean(oldsplitB['Z'])

metaA = pd.read_csv('~/Desktop/20160_meta_A_s0.tsv.bgz', sep='\t',compression='gzip')
metaB = pd.read_csv('~/Desktop/20160_meta_B_s0.tsv.bgz', sep='\t',compression='gzip')

###############################################################################
# Compare distributions of meta betas for halves of a meta split
import matplotlib
matplotlib.rcParams['agg.path.chunksize'] = 1000000

meta_df = pd.read_csv('~/Desktop/meta_beta_df.tsv.gz', sep='\t',compression='gzip')
meta_df['Z_A'] = meta_df.meta_beta_A/meta_df.meta_se_A
meta_df['Z_B'] = meta_df.meta_beta_B/meta_df.meta_se_B
meta_df['meta_diff'] = (meta_df.meta_beta_A-meta_df.meta_beta_B)
meta_df['meta_se'] = (((meta_df.meta_se_A**2+meta_df.meta_se_A**2)/2)**(1/2))
for i in range(10):
    meta_df_sub = meta_df.iloc[100000*(i):100000*(i+1)]

    sortby = 'meta_diff'
    meta_diff = meta_df_sub.sort_values(by=sortby).meta_diff
    meta_se = meta_df_sub.sort_values(by=sortby).meta_se
#    plt.plot(range(len(meta_diff)),meta_diff,'.-',color=[0, 0, 1])
#    plt.fill_between(range(len(meta_diff)),meta_diff-2*meta_se, meta_diff+2*meta_se,
#                     color=[0, 0.5, 1])
#    plt.plot([0, len(meta_diff)],[0, 0],'r--')
#    plt.xlim([0, len(meta_diff)])
#    fig=plt.gcf()
#    scale=1
#    fig.set_size_inches(12*scale, 4*scale)
    print('\nFor subset '+str(i))
    print('meta beta r: '+ str(stats.pearsonr(meta_df_sub.meta_beta_A, meta_df_sub.meta_beta_B)[1]))
    print('meta se r: '+str(stats.pearsonr(meta_df_sub.meta_se_A, meta_df_sub.meta_se_B)[1]))
    print('ttest meta beta: '+str(stats.ttest_ind(meta_df_sub.meta_se_A, meta_df_sub.meta_se_B)[1]))
    print('ttest meta beta paired: '+str(stats.ttest_rel(meta_df_sub.meta_se_A, meta_df_sub.meta_se_B)[1]))

plt.plot(meta_df.meta_beta_A,meta_df.meta_beta_B,'.',alpha=0.5)
#plt.errorbar(meta_df.meta_beta_A,meta_df.meta_beta_B,xerr=meta_df.meta_se_A,yerr=meta_df.meta_se_B,linestyle='')
plt.plot([-0.1,0.1],[-0.1,0.1])
plt.xlabel('meta_beta_A')
plt.ylabel('meta_beta_B')
plt.title('Comparison of meta betas for split halves A and B')
fig=plt.gcf()
scale=1.5
fig.set_size_inches(6*scale, 4*scale)
fig.savefig('/Users/nbaya/Desktop/meta_beta_comparison.png',dpi=300)

stats.pearsonr(meta_df.meta_beta_A, meta_df.meta_beta_B)[1]
stats.ttest_ind(meta_df.meta_beta_A, meta_df.meta_beta_B)[1]

plt.plot(range(len(meta_df)),meta_df.meta_diff,'.',alpha=0.5,markersize=1)
plt.xlim([0,len(meta_df)])
plt.xlabel('snp index')
plt.ylabel('difference between meta betas (A-B)')
plt.title('Difference between meta betas')
fig=plt.gcf()
scale=1.5
fig.set_size_inches(6*scale, 4*scale)
fig.savefig('/Users/nbaya/Desktop/meta_beta_diff.png',dpi=300)

plt.plot(range(len(meta_df)),abs(meta_df.meta_se_A-meta_df.meta_se_B),'.',alpha=0.5,markersize=1)
plt.xlim([0,len(meta_df)])
plt.xlabel('snp index')
plt.ylabel('difference between meta se (A-B)')
plt.title('Absolute difference between meta se')
fig=plt.gcf()
scale=1.5
fig.set_size_inches(6*scale, 4*scale)
fig.savefig('/Users/nbaya/Desktop/meta_se_absdiff.png',dpi=300)

j=0
eps = 50
plt.plot([-eps,eps],[-eps,eps],'black')
for i in np.logspace(0,1,50):
    meta_df_temp = meta_df[(abs(meta_df.Z_A-meta_df.Z_B) < (i)/10*np.max(abs(meta_df.Z_A-meta_df.Z_B))) &
                           (abs(meta_df.Z_A-meta_df.Z_B) > (j)/10*np.max(abs(meta_df.Z_A-meta_df.Z_B)))]
#    print((j)/10*np.max(abs(meta_df.Z_A-meta_df.Z_B)))
#    print((i)/10*np.max(abs(meta_df.Z_A-meta_df.Z_B)))
#    print(len(meta_df_temp))
    j = i
    plt.plot(meta_df_temp.Z_A,meta_df_temp.Z_B,linestyle='',marker='.',alpha=(i)/10*1,color=[1-i/10,0 , i/10])
plt.xlim([-eps, eps])
plt.ylim([-eps, eps])
plt.xlabel('Z_A')
plt.ylabel('Z_B')
plt.title('Comparison of z scores for split halves A and B')
fig=plt.gcf()
scale=1.5
fig.set_size_inches(6*scale, 6*scale)
fig.savefig('/Users/nbaya/Desktop/zscore_comparison.png',dpi=300)

plt.hist(meta_df.meta_diff,50)

sns.kdeplot(meta_df.meta_diff)
plt.xlabel('meta diff (A-B)')
plt.ylabel('density')
plt.title('kde of meta diff')
fig=plt.gcf()
scale=1.5
fig.set_size_inches(6*scale, 4*scale)
fig.savefig('/Users/nbaya/Desktop/meta_diff_kde.png',dpi=300)




plt.plot(range(len(meta_df.iloc[0:100000])), meta_df.iloc[0:100000].meta_diff)
for i in range(10):
    plt.plot([10000*i,10000*i],[-0.10,0.10],'r--',alpha=0.5)
    
    
    
"""
gmt = mt.add_col_index()
n_chunks = 300
batch='1'
i=1
pi = ['A']*int(n_chunks/2) + ['B']*int(n_chunks/2)
seed_id = int(batch+str(i).zfill(4)) #OPTION 2: create a seed_id unique to every split
randstate = np.random.RandomState(seed_id) #OPTION 2: seed with seed_id

randstate.shuffle(pi)
gmt_shuf = gmt.annotate_cols(label = hl.literal(pi)[hl.int32(gmt.col_idx)])

mt = gmt_shuf.group_cols_by(gmt_shuf.label).aggregate(unnorm_meta_beta=hl.agg.sum(gmt_shuf.beta / gmt_shuf.standard_error ** 2),
                            inv_se2 = hl.agg.sum(1 / gmt_shuf.standard_error ** 2)) 

mt = mt.annotate_entries(meta_beta = mt.unnorm_meta_beta/mt.inv_se2)
mt = mt.annotate_entries(meta_se = hl.sqrt(1/mt.inv_se2))

meta_beta_A = np.asarray(mt.aggregate_entries(hl.agg.filter(mt.label=='A', hl.agg.collect(mt.meta_beta))))
meta_beta_B = np.asarray(mt.aggregate_entries(hl.agg.filter(mt.label=='B', hl.agg.collect(mt.meta_beta))))
meta_se_A = np.asarray(mt.aggregate_entries(hl.agg.filter(mt.label=='A', hl.agg.collect(mt.meta_se))))
meta_se_B = np.asarray(mt.aggregate_entries(hl.agg.filter(mt.label=='B', hl.agg.collect(mt.meta_se))))

meta_beta_df = pd.DataFrame({'meta_beta_A':meta_beta_A,'meta_beta_B':meta_beta_B,'meta_se_A':meta_se_A, 'meta_se_B':meta_se_B})
"""

###############################################################################
"""
Cochrans Q
"""
import pandas as pd
from scipy import stats
Q = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/split/cochrans_q.50_batch_1.tsv.bgz',sep='\t',compression='gzip')
fig, ax = plt.subplots(1,1)
df=299
x = np.linspace(stats.chi2.ppf(0.0001, df),stats.chi2.ppf(0.9999, df), 100)
ax.plot(x, stats.chi2.pdf(x, df),
        'r-', lw=5, alpha=0.6, label='chi2 pdf')
hist = ax.hist(Q.Q_j,100,density=True)
ax.legend(['chi^2 dist for df=299','hist of Cochrans Qs'])
fig.set_size_inches(6*1.5,4*1.5)
plt.title("Cochran's Qs from 300 chunks compared with chi^2 dist for df=299")
plt.xlabel('Q')
plt.ylabel('density')
fig.savefig('/Users/nbaya/Desktop/cochrans_q_300_chunks.png',dpi=600)
        
Q['pval']=(Q.Q_j>299)*(1-stats.chi2.cdf(Q.Q_j,df=299) + (stats.chi2.cdf(299-(Q.Q_j-299),df=299)))+(
        Q.Q_j<=299)*(stats.chi2.cdf(Q.Q_j,df=299)+1-stats.chi2.cdf(299+(299-Q.Q_j),df=299))
Q['pval']=stats.chi2.cdf(Q.Q_j,df=299)
Q['pval']=1-stats.chi2.cdf(Q.Q_j,df=299)

np.isnan(Q.pval)

obs = np.sort(-np.log10(Q.pval.values))
exp = np.sort(-np.log10(np.linspace(1,1/len(obs),len(obs))))
plt.plot(exp, obs,'o',alpha=0.5)
plt.plot([0,6],[0,6],'r-',alpha=0.5)
plt.title("Q-Q plot of 300 chunks Cochran's Qs p-values")
plt.xlabel('expected')
plt.ylabel('observed')
fig=plt.gcf()
fig.set_size_inches(8,6)
fig.savefig('/Users/nbaya/Desktop/cochrans_q_300_chunks_qqplot_HA_ne.png',dpi=600)



Q_split = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/split/cochrans_q.50_batch_1_split.tsv.bgz',sep='\t',compression='gzip')
fig, ax = plt.subplots(1,1)
df=1
x = np.linspace(stats.chi2.ppf(0.00001, df),stats.chi2.ppf(0.99, df), 1000)
ax.plot(x, stats.chi2.pdf(x, df)/400000,
        'r-', lw=5, alpha=0.6, label='chi2 pdf')
#ax.hist(Q_split[Q_split.Q_j<np.max(x)].Q_j,10,density=True)
ax.legend(['chi^2 dist for df=1','hist of Cochrans Qs'])
plt.xlim([0,np.max(x)])
fig.set_size_inches(6*1.5,4*1.5)
plt.title("Cochran's Qs from 2 halves compared with chi^2 dist for df=1")
fig.savefig('/Users/nbaya/Desktop/cochrans_q_2_halves.png')


Q_split['pval']=1-stats.chi2.cdf(Q_split.Q_j,df=1)
obs = np.sort(-np.log10(Q_split.pval.values))
exp = np.sort(-np.log10(np.linspace(1,1/len(obs),len(obs))))
plt.plot(exp, obs,'o',alpha=0.5)
plt.plot([0,6],[0,6],'r-',alpha=0.5)
plt.title("Q-Q plot of 2 halves Cochran's Qs p-values")
plt.xlabel('expected')
plt.ylabel('observed')
fig=plt.gcf()
fig.set_size_inches(8,6)
fig.savefig('/Users/nbaya/Desktop/cochrans_q_300_chunks_qqplot_HA_ne.png',dpi=600)



###############################################################################
"""
╔══════════════════════════════════╗
║ Create split.tsv file for h2part ║
╚══════════════════════════════════╝
"""
import pandas as pd
import os
import math

phen = '50_raw_res'
batch = '1'
n = 300
isnested = True
for subbatch in range(1,51):
    subbatch = str(subbatch)
    
    df = pd.read_csv('~/Documents/lab/ukbb-sexdiff/split/phenotypes.split.tsv',sep='\t').drop(list(range(24)))
    
    phenotypes = pd.read_csv('~/Documents/lab/ukbb-sexdiff/imputed-v3-results/phenotypes.both_sexes.tsv',sep='\t')
    
    vals = phenotypes.loc[phenotypes['phenotype'] == '50']
    
    #os.system('gsutil cp gs://nbaya/rg_sex/smoking.h2part.tsv ~/Documents/lab/ukbb-sexdiff/rg_sex/')
         
    #NOTE: Hardcoded values are from dataset after withdrawn samples are removed
    for i in range(1,n+1):
        df.at[i-1,'phenotype'] = str(phen)+'_sample_n'+str(n)+'_batch_'+batch+'.'+subbatch+'_set'+str(i)+'_nested'+str(isnested)
#        df.at[i-1,'phenotype'] = str(phen)+'_sample_A_batch_'+batch+'.'+subbatch+'_set'+str(i)
#        df.at[i-1,'phenotype'] = str(phen)+'_sample_n'+str(n)+'_A_batch_'+batch+'_set'+str(i) #Used ONLY for the set containing all subsets
        df.at[i-1,'description'] = str(phen)+'_set'+str(i)
        df.at[i-1,'n_non_missing'] = int(360338*i/n) #int(vals['n_non_missing']*i/n) #{'6138_100': 357500, '50_raw': 360338, '30100': 350423}
        df.at[i-1,'n_missing'] = int(806*i/n) #int(vals['n_missing']*i/n) #{'6138_100': 3644, '50_raw': 806, '30100': 10721}
        if math.isnan(vals['n_cases']):
            df.at[i-1,'n_cases'] = float('NaN')
            df.at[i-1,'n_controls'] = float('NaN')
        else:
            df.at[i-1,'n_cases'] = int(vals['n_cases']*i/n) #{6138_100: 61083}
            df.at[i-1,'n_controls'] = int(vals['n_controls']*i/n) #{6138_100: 296417}
            
    df['source'] = "'"+phen+"'"
    
    filename = phen+'.downsample_n'+str(n)+'_nested'+str(isnested)+'.batch_'+batch+'.'+subbatch+'.tsv'
    local_wd = '~/Documents/lab/ukbb-sexdiff/split/'
    cloud_wd = 'gs://nbaya/h2part'
    
    df.to_csv(local_wd+filename,sep='\t', index=False)
    
    os.system('gsutil cp '+local_wd+filename+' '+cloud_wd)

###############################################################################

phen = 'blood_pressure_phens'
    
df = pd.read_csv('~/Documents/lab/ukbb-sexdiff/split/phenotypes.split.tsv',sep='\t').drop(list(range(24)))

#phenotypes = pd.read_csv('~/Documents/lab/ukbb-sexdiff/imputed-v3-results/phenotypes.both_sexes.tsv',sep='\t')

#vals = phenotypes.loc[phenotypes['phenotype'] == '50']

#os.system('gsutil cp gs://nbaya/rg_sex/smoking.h2part.tsv ~/Documents/lab/ukbb-sexdiff/rg_sex/')
     
#NOTE: Hardcoded values are from dataset after withdrawn samples are removed

df.at[0,'phenotype'] = '4079_female'
df.at[1,'phenotype'] = '4079_male'
df.at[2,'phenotype'] = '4080_female'
df.at[3,'phenotype'] = '4080_male'
df.at[:,'description'] = str(phen)
df.at[0,'n_non_missing'] =  142710 #4079 female
df.at[1,'n_non_missing'] =  107511 #4079 male
df.at[2,'n_non_missing'] =  142710 #4080 female
df.at[3,'n_non_missing'] =  107511 #4080 male
df.loc[:,'n_missing'] = 361144-df.loc[:,'n_non_missing'] 
df.at[:,'n_cases'] = float('NaN')
df.at[:,'n_controls'] = float('NaN')

df = df.astype({'n_non_missing': int, 'n_missing': int})

#df.at[:,'n_missing'] = int(806) #int(vals['n_missing']*i/n) #{'6138_100': 3644, '50_raw': 806, '30100': 10721}
#if math.isnan(vals['n_cases']):
#    df.at[:,'n_cases'] = float('NaN')
#    df.at[:,'n_controls'] = float('NaN')
#else:
#    df.at[:,'n_cases'] = int(vals['n_cases']) #{6138_100: 61083}
#    df.at[:,'n_controls'] = int(vals['n_controls']) #{6138_100: 296417}
#
#df.loc[0,'female_n'] = 142710   #4079 female
#df.loc[1,'male_n'] = 107511     #4079 male
#df.loc[2,'female_n'] = 142710   #4080 female
#df.loc[3,'male_n'] = 107511     #4080 male
        
df['source'] = "'"+phen+"'"


filename = phen+'.split.tsv'
local_wd = '~/Documents/lab/ukbb-sexdiff/split/'
cloud_wd = 'gs://nbaya/rg_sex/sex_strat/'

df.to_csv(local_wd+filename,sep='\t', index=False)

os.system('gsutil cp '+local_wd+filename+' '+cloud_wd)

#df = pd.read_csv('~/Downloads/ukbb31063.h2part_results.'+phen+'split.tsv.gz',sep='\t',compression='gzip').iloc[:,0:20]
df.h2_observed

###############################################################################

phen = 'simulated_phens'
    
df = pd.read_csv('~/Documents/lab/ukbb-sexdiff/split/phenotypes.split.tsv',sep='\t').drop(list(range(24)))

df.at[0,'phenotype'] = 'hm3.sim_h2_0.5.no_popstrat.sumstats' #'normal_phenotype_h2_0.08_sumstats'
df.at[1,'phenotype'] = 'hm3.sim_h2_0.5.w_popstrat_PC1.popstrat_s2_3.sumstats' #'ordinal_gamma_h2_0.08_sumstats'
#df.at[2,'phenotype'] = 'prevalence_0.06_h2_0.2_sumstats'
#df.at[3,'phenotype'] = 'prevalence_0.01_h2_0.3_sumstats'
df.at[:,'description'] = str(phen)
df.at[0,'n_non_missing'] =  360338 
df.at[1,'n_non_missing'] =  360338  
#df.at[2,'n_non_missing'] =  360338 
#df.at[3,'n_non_missing'] =  360338 
df.loc[:,'n_missing'] = 361144-df.loc[:,'n_non_missing'] 
df.at[0,'n_cases'] = float('NaN')
df.at[1,'n_controls'] = float('NaN')
#df.at[2,'n_cases'] = int(df.at[2,'n_non_missing']*0.06)
#df.at[3,'n_cases'] = int(df.at[3,'n_non_missing']*0.01)
#df.at[2,'n_controls'] = int(df.at[2,'n_non_missing']*(1-0.06))
#df.at[3,'n_controls'] = int(df.at[3,'n_non_missing']*(1-0.01))

        
df['source'] = "'"+phen+"'"


filename = phen+'.split.tsv'
local_wd = '~/Documents/lab/ukbb-sexdiff/split/'
cloud_wd = 'gs://nbaya/ldscsim/'

df.to_csv(local_wd+filename,sep='\t', index=False)

os.system('gsutil cp '+local_wd+filename+' '+cloud_wd)


###############################################################################
"""
╔═══════════════════════════╗
║ h2part results for height ║
╚═══════════════════════════╝
Batches 1-5, 1.1-1.50 
There were 50 "proper" subbatches (1.1-1.50) for batch 1, which is important for 
the seed given to the subbatch when permuting chunks before meta-analyzing.
I initially ran this for 5 batches before realizing I should be using the subbatch
framework. Those initial five batches did not have the same method of generating
the random seed.
"""
phen = '50'

#for i in range(11,51):
#    os.system('gsutil cp gs://nbaya/split/meta_split/h2part/ukbb31063.h2part_results.'+phen+'.downsample.batch_'+batch+'.'+str(i)+'.tsv.gz ~/Documents/lab/ukbb-sexdiff/h2part/')
#os.system('gsutil cp gs://nbaya/split/meta_split/h2part/ukbb31063.h2part_results.50.downsample.batch_1.1_set300.tsv.gz ~/Documents/lab/ukbb-sexdiff/h2part/')

h2_50_batch1 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.50.downsample.batch_1.tsv.gz', sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20]
h2_50_batch1_set300 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.50.downsample.batch_1.1_set300.tsv.gz', sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20]
h2_50_batch1_n150_set150 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.50.downsample_n150.batch_1.1.tsv.gz', sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20]

h2_50 = h2_50_batch1.copy()[['n','h2_observed','h2_observed_se']]
h2_50 = h2_50.append(h2_50_batch1_set300.iloc[1,:].copy()[['n','h2_observed','h2_observed_se']],ignore_index=True) #only use when showing only batch 1 subbatch results
h2_50 = h2_50.rename(axis='columns',mapper={"h2_observed":"h2_observed_1","h2_observed_se":"h2_se_1"})


#for i in range(2,6): #adds batches 2-5 results
#    other = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.50.downsample.batch_'+str(i)+'.tsv.gz', sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20]
#    h2_50[['h2_observed_'+str(i),'h2_se_'+str(i)]] = other[['h2_observed','h2_observed_se']]


for i in range(1,51): #adds batch 1 subbatch results
    other = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.50.downsample.batch_1.'+str(i)+'.tsv.gz', sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20]
    other = other.append(h2_50_batch1_set300.iloc[1,:].copy()[['n','h2_observed','h2_observed_se']],ignore_index=True) #only use when showing only batch 1 subbatch results
    h2_50[['h2_observed_1_'+str(i),'h2_se_1_'+str(i)]] = other[['h2_observed','h2_observed_se']]

h2 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/rg_sex/ukbb31063.both_sexes.h2part_results.phesant.tsv.gz', 
                              sep='\t',compression='gzip').rename(index=str,columns={'phenotype':'phen'}).iloc[:,0:20]
h2_ref = h2[h2['phen']==phen+'_irnt'].h2_observed[0] # full data ref

n_batches = int((h2_50.shape[1]-1)/2)

combined_mean_h2 = np.mean(h2_50.filter(regex=('observed')),axis=1)
combined_h2_se = np.mean(h2_50.filter(regex=('_se'))**2,axis=1)**(1/2)
combined_h2_se_mean = (np.mean(h2_50.filter(regex=('_se'))**2,axis=1)/h2_50.filter(regex=('_se')).shape[1])**(1/2)

#h2_ref = 0.474646 # set 300 reference
plt.axhline(y=h2_ref,color='k',alpha=1, linewidth=0.5) # h2part reference
#plt.axhline(y=0.474646,color='k',alpha=1, linewidth=0.5) # set 300 reference
plt.plot(h2_50['n'],combined_mean_h2,'.-',color=[0, 0, 1])
plt.fill_between(h2_50['n'],list(map(float,combined_mean_h2-2*combined_h2_se)),
                 list(map(float,combined_mean_h2+2*combined_h2_se)),alpha = 0.1,color=[0, 0.5, 1])
plt.fill_between(h2_50['n'],list(map(float,combined_mean_h2-2*combined_h2_se_mean)),
                 list(map(float,combined_mean_h2+2*combined_h2_se_mean)),alpha = 0.2,color=[0, 0, 1])
plt.ylabel('h2_observed')
plt.xlabel('n_non_missing')
plt.ylim([h2_ref-0.1,h2_ref+0.1])
plt.xlim([0,np.max(h2_50['n'])*1.01])
plt.xlim([0,1e5*1.03])
plt.text(15e3,(h2_ref-0.1*0.8),'n = 25k')
plt.text(40e3,(h2_ref-0.1*0.8),'n = 50k')
plt.text(88.5e3,(h2_ref-0.1*0.8),'n = 100k')
n_batches = int((h2_50.shape[1]-1)/2)
plt.title('Combined h2_observed of 50 ("Standing height")\n(%i batches)' % n_batches)
#plt.legend(['subset 300 h2 ref \n(%f)' % h2_ref,'combined h2','prediction interval CI','CI for mean'],loc=1)
plt.legend(['full data h2 ref \n(%f)' % h2_ref,'combined h2','prediction interval CI','CI for mean'],loc=1)
plt.axvline(x=25e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=50e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=100e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
fig = plt.gcf()
fig.set_size_inches(12, 8*.7)
fig.savefig('/Users/nbaya/Desktop/upsampling_'+phen+'_h2_observed_height_'+str(n_batches)+'_batches.png',dpi=300)


#plt.axhline(y=h2_ref,color='k',alpha=1, linewidth=0.5) # h2part reference
#plt.axhline(y=0.474646,color='k',alpha=1, linewidth=0.5) # set 300 reference
plt.plot(h2_50['n'],h2_50.filter(regex='observed'))
for i in range(1,6):
    plt.fill_between(h2_50['n'],h2_50['h2_observed_'+str(i)]-2*h2_50['h2_se_'+str(i)],
                     h2_50['h2_observed_'+str(i)]+2*h2_50['h2_se_'+str(i)],alpha = 0.2)
for i in range(1,51):
    plt.fill_between(h2_50['n'],h2_50['h2_observed_1_'+str(i)]-2*h2_50['h2_se_1_'+str(i)],
                     h2_50['h2_observed_1_'+str(i)]+2*h2_50['h2_se_1_'+str(i)],alpha = 0.2)
plt.axvline(x=25e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=50e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=100e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.ylabel('h2_observed')
plt.xlabel('n_non_missing')
plt.ylim([h2_ref-0.1,h2_ref+0.1])
plt.xlim([0,np.max(h2_50['n'])*1.01])
plt.text(15e3,(h2_ref-0.1*0.8),'n = 25k')
plt.text(40e3,(h2_ref-0.1*0.8),'n = 50k')
plt.text(88.5e3,(h2_ref-0.1*0.8),'n = 100k')
plt.title('h2_observed of height\n(%i batches)' % n_batches)
plt.legend(['full data h2 ref ref\n(0.485)'],loc=1)
fig = plt.gcf()
fig.set_size_inches(12, 8*.7)
fig.savefig('/Users/nbaya/Desktop/upsampling_h2_observed_height_'+str(n_batches)+'_batches_not_combined.png',dpi=300)

#Read in results of using n_chunks=150
h2_50_n150 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.50.downsample_n150.batch_1.1.tsv.gz', sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20]
h2_50_n150[['phenotype','h2_observed']].iloc[1] #expect to be closer to h2_ref than the n_chunks=300 version (h2_50_batch1_set300)
h2_50_batch1_set300[['phenotype','h2_observed']].iloc[1]
h2_ref

###############################################################################
"""
╔════════════════════════════╗
║ h2part results for smoking ║
╚════════════════════════════╝
Batches: 1-3, 5, 1.1-1.50
See above for an explanation of how the batches differ.

NOTE: There is no batch 4 (probably due to an error while running the task). 
      If you are only interested in using the 50 subbatches derived from batch 1
      this will not affect you.
"""
phen = '20160'
batch='1'

#for i in range(1,51):
#    os.system('gsutil cp gs://nbaya/split/meta_split/h2part/ukbb31063.h2part_results.'+phen+'.downsample.batch_'+batch+'.'+str(i)+'.tsv.gz ~/Documents/lab/ukbb-sexdiff/h2part/')
os.system('gsutil cp gs://nbaya/split/meta_split/h2part/ukbb31063.h2part_results.20160.downsample.batch_1.1_set300.tsv.gz ~/Documents/lab/ukbb-sexdiff/h2part/')

h2_20160_batch1 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.'+phen+'.downsample.batch_1.tsv.gz', 
                              sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20].copy()[['n','h2_observed','h2_observed_se']]
#h2_20160_batch1_set300 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.20160.downsample.batch_1.1_set300.tsv.gz', 
#                              sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20].copy()[['n','h2_observed','h2_observed_se']]

h2_20160 = h2_20160_batch1.rename(axis='columns',mapper={"h2_observed":"h2_observed_1","h2_observed_se":"h2_se_1"})
    
for i in range(2,6):
    if i != 4:
        other = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.'+phen+'.downsample.batch_'+str(i)+'.tsv.gz', sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20]
        h2_20160[['h2_observed_'+str(i),'h2_se_'+str(i)]] = other[['h2_observed','h2_observed_se']]
for i in range(1,51):
    other = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.'+phen+'.downsample.batch_1.'+str(i)+'.tsv.gz', sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20]
    h2_20160[['h2_observed_1_'+str(i),'h2_se_1_'+str(i)]] = other[['h2_observed','h2_observed_se']]
    
#for i in range(300,301):
#    other = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.'+phen+'.downsample.batch_'+batch+'.'+str(i)+'.tsv.gz', sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20]
#    h2_20160[['h2_observed_'+str(i),'h2_se_'+str(i)]] = other[['h2_observed','h2_observed_se']]

h2 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/rg_sex/ukbb31063.both_sexes.h2part_results.phesant.tsv.gz', 
                              sep='\t',compression='gzip').rename(index=str,columns={'phenotype':'phen'}).iloc[:,0:20]
h2_ref = h2[h2['phen']==phen].h2_observed[0] # full data ref

n_batches = int((h2_20160.shape[1]-1)/2)

plt.axhline(y=h2_ref,color='k',alpha=1, linewidth=0.5) # h2part reference
#plt.axhline(y=0.073778,color='k',alpha=1, linewidth=0.5) # set 300 reference
plt.plot(h2_20160['n'],h2_20160.filter(regex=('observed')))
plt.ylim([-0.4, 0.4])
for i in range(1,6):
    plt.fill_between(h2_20160['n'],h2_20160['h2_observed_'+str(i)]-2*h2_20160['h2_se_'+str(i)],
                     h2_20160['h2_observed_'+str(i)]+2*h2_20160['h2_se_'+str(i)],alpha = 0.2)
for i in range(1,51):
    plt.fill_between(h2_20160['n'],h2_20160['h2_observed_1_'+str(i)]-2*h2_20160['h2_se_1_'+str(i)],
                 h2_20160['h2_observed_1_'+str(i)]+2*h2_20160['h2_se_1_'+str(i)],alpha = 0.2)
plt.ylabel('h2_observed')
plt.xlabel('n_non_missing')
plt.legend(['full data h2 ref ref\n(0.074)'],loc=1)
plt.axvline(x=25e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=50e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=100e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.ylim([h2_ref-0.1,h2_ref+0.1])
plt.xlim([0,np.max(h2_20160['n'])*1.01])
plt.text(15e3,(h2_ref-0.1*0.8),'n = 25k')
plt.text(40e3,(h2_ref-0.1*0.8),'n = 50k')
plt.text(88.5e3,(h2_ref-0.1*0.8),'n = 100k')
plt.title('comparison of h2_observed between batches\nfor phen 20160 (smoking)')
fig = plt.gcf()
fig.set_size_inches(12, 8*.7)
#plt.xscale('log')
fig.savefig('/Users/nbaya/Desktop/upsampling_h2_observed_smoking_batches_all.png',dpi=300)

#####
# Combined h2 across batches
#####

combined_mean_h2 = np.mean(h2_20160.filter(regex=('observed')),axis=1)
combined_h2_se = np.mean(h2_20160.filter(regex=('_se'))**2,axis=1)**(1/2)
combined_h2_se_mean = (np.mean(h2_20160.filter(regex=('_se'))**2,axis=1)/h2_20160.filter(regex=('_se')).shape[1])**(1/2)


plt.axhline(y=h2_ref,color='k',alpha=1, linewidth=0.5) # h2part reference
#plt.axhline(y=0.073778,color='k',alpha=1, linewidth=0.5) # set 300 reference
plt.plot(h2_20160['n'],combined_mean_h2,'.-',color=[0, 0, 1])
plt.fill_between(h2_20160['n'],list(map(float,combined_mean_h2-2*combined_h2_se)),
                 list(map(float,combined_mean_h2+2*combined_h2_se)),alpha = 0.1,
                 color=[0, 0.5, 1])
plt.fill_between(h2_20160['n'],list(map(float,combined_mean_h2-2*combined_h2_se_mean)),
                 list(map(float,combined_mean_h2+2*combined_h2_se_mean)),alpha = 0.2,
                 color=[0, 0, 1])
plt.legend(['full data h2 ref ref\n(%f)' % h2_ref,'combined h2','prediction interval CI','CI for mean'],loc=1)
plt.axvline(x=25e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=50e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=100e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.ylabel('h2_observed')
plt.xlabel('n_non_missing')
plt.ylim([h2_ref-0.1,h2_ref+0.1])
plt.xlim([0,np.max(h2_20160['n'])*1.01])
plt.text(15e3,(h2_ref-0.1*0.8),'n = 25k')
plt.text(40e3,(h2_ref-0.1*0.8),'n = 50k')
plt.text(88.5e3,(h2_ref-0.1*0.8),'n = 100k')
plt.title('Combined h2_observed of 20160 ("Ever Smoked")\n('+str(n_batches)+' batches)')
fig = plt.gcf()
fig.set_size_inches(12, 8*.7)
fig.savefig('/Users/nbaya/Desktop/upsampling_'+phen+'_h2_observed_'+str(n_batches)+'_batches.png',dpi=300)




plt.axhline(y=h2_ref,color='k',alpha=1, linewidth=0.5) # h2part reference
#plt.axhline(y=0.474646,color='k',alpha=1, linewidth=0.5) # set 300 reference
plt.plot(h2_20160['n'],h2_20160.filter(regex='observed'))
for i in range(1,6):
    plt.fill_between(h2_20160['n'],h2_50['h2_observed_'+str(i)]-2*h2_20160['h2_se_'+str(i)],
                     h2_20160['h2_observed_'+str(i)]+2*h2_20160['h2_se_'+str(i)],alpha = 0.2)
for i in range(1,51):
    plt.fill_between(h2_20160['n'],h2_20160['h2_observed_1_'+str(i)]-2*h2_20160['h2_se_1_'+str(i)],
                     h2_20160['h2_observed_1_'+str(i)]+2*h2_20160['h2_se_1_'+str(i)],alpha = 0.2)
plt.axvline(x=25e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=50e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=100e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.ylim([h2_ref-0.1,h2_ref+0.1])
plt.xlim([0,np.max(h2_20160['n'])*1.01])
plt.text(15e3,(h2_ref-0.1*0.8),'n = 25k')
plt.text(40e3,(h2_ref-0.1*0.8),'n = 50k')
plt.text(88.5e3,(h2_ref-0.1*0.8),'n = 100k')
plt.title('h2_observed of height\n(%i batches)' % n_batches)
plt.legend(['full data h2 ref ref\n(0.485)'],loc=1)
fig = plt.gcf()
fig.set_size_inches(12, 8*.7)

#plt.errorbar(h2_20160['n'], combined_mean_h2, 2*combined_h2_se,
#             linestyle='None', marker='.',alpha=1,elinewidth=0.5)
#plt.ylabel('h2_observed')
#plt.xlabel('n_non_missing')
#plt.legend(['combined h2'],loc=1)
#plt.axvline(x=25e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
#plt.axvline(x=50e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
#plt.axvline(x=100e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
##plt.ylim([0,1])
#plt.ylim([0,0.2])
#plt.text(15e3,0.425,'n = 25k')
#plt.text(40e3,0.425,'n = 50k')
#plt.text(88.5e3,0.425,'n = 100k')
#plt.title('Combined h2_observed of smoking phenotype\n(4 batches)')
#fig = plt.gcf()
#fig.set_size_inches(12, 8*.7)

###############################################################################
"""
╔═════════════════════════════╗
║ h2part results for diabetes ║
╚═════════════════════════════╝
Batches: 1.1-1.50
phen code: 2443
"""
#Read h2part results file for diabetes phenotype (2443)

phen = '2443'
batch='1'

#for i in range(1,51):
#    os.system('gsutil cp gs://nbaya/split/meta_split/h2part/ukbb31063.h2part_results.'+
#              phen+'.downsample_n150.batch_1.'+str(i)+'.tsv.gz ~/Documents/lab/ukbb-sexdiff/h2part/')

#os.system('gsutil cp gs://nbaya/split/meta_split/h2part/ukbb31063.h2part_results.2443.downsample_n150.batch_1.1_set150.tsv.gz ~/Documents/lab/ukbb-sexdiff/h2part/')

h2_2443_batch1 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.'+phen+'.downsample_n150.batch_1.1.tsv.gz', 
                              sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20].copy()[['n','h2_observed','h2_observed_se']]
#h2_2443_batch1_set150 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.2443.downsample_n150.batch_1.1_set150.tsv.gz', 
#                              sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20].copy()[['n','h2_observed','h2_observed_se']]
#h2_2443_batch1 = h2_2443_batch1.append(h2_2443_batch1_set150.iloc[1,:].copy()[['n','h2_observed','h2_observed_se']],ignore_index=True) #only use when showing only batch 1 subbatch results

h2_2443 = h2_2443_batch1.rename(axis='columns',mapper={"h2_observed":"h2_observed_1_1","h2_observed_se":"h2_se_1_1"})
    

for i in range(2,51):
    other = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.'+phen+'.downsample_n150.batch_'+batch+'.'+str(i)+'.tsv.gz', sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20]
#    other = other.append(h2_2443_batch1_set150.iloc[1,:].copy()[['n','h2_observed','h2_observed_se']],ignore_index=True) #only use when showing only batch 1 subbatch results
    h2_2443[['h2_observed_1_'+str(i),'h2_se_1_'+str(i)]] = other[['h2_observed','h2_observed_se']]
    
#for i in range(300,301):
#    other = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.'+phen+'.downsample.batch_'+batch+'.'+str(i)+'.tsv.gz', sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20]
#    h2_20160[['h2_observed_'+str(i),'h2_se_'+str(i)]] = other[['h2_observed','h2_observed_se']]

h2 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/rg_sex/ukbb31063.both_sexes.h2part_results.phesant.tsv.gz', 
                              sep='\t',compression='gzip').rename(index=str,columns={'phenotype':'phen'}).iloc[:,0:20]
h2_ref = h2[h2['phen']==phen].h2_observed[0] # full data ref

n_batches = int((h2_2443.shape[1]-1)/2)

plt.axhline(y=h2_ref,color='k',alpha=1, linewidth=0.5) 
#plt.axhline(y=0.044821,color='k',alpha=1, linewidth=0.5) # full data ref
plt.plot(h2_2443['n'],h2_2443.filter(regex=('observed')),'.-')
plt.ylim([-0.4, 0.4])
for i in range(1,51):
    plt.fill_between(h2_2443['n'],h2_2443['h2_observed_1_'+str(i)]-2*h2_2443['h2_se_1_'+str(i)],
                 h2_2443['h2_observed_1_'+str(i)]+2*h2_2443['h2_se_1_'+str(i)],alpha = 0.2)
plt.ylabel('h2_observed')
plt.xlabel('n_non_missing')
plt.legend(['full data h2 ref ref\n(%f)' % h2_ref],loc=1)
plt.axvline(x=25e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=50e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=100e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.ylim([h2_ref-0.1,h2_ref+0.1])
plt.xlim([0,np.max(h2_2443['n'])*1.01])
plt.text(15e3,(h2_ref-0.1*0.8),'n = 25k')
plt.text(40e3,(h2_ref-0.1*0.8),'n = 50k')
plt.text(88.5e3,(h2_ref-0.1*0.8),'n = 100k')
plt.title('comparison of h2_observed between batches\nfor phen 2443 (diabetes)')
fig = plt.gcf()
fig.set_size_inches(12, 8*.7)
#plt.xscale('log')
fig.savefig('/Users/nbaya/Desktop/upsampling_'+phen+'_h2_observed_'+str(n_batches)+'_batches.png',dpi=300)

#####
#Combined h2 across batches
#####

combined_mean_h2 = np.mean(h2_2443.filter(regex=('observed')),axis=1)
combined_h2_se = np.mean(h2_2443.filter(regex=('_se'))**2,axis=1)**(1/2)
combined_h2_se_mean = (np.mean(h2_2443.filter(regex=('_se'))**2,axis=1)/h2_2443.filter(regex=('_se')).shape[1])**(1/2)

plt.axhline(y=h2_ref,color='k',alpha=1, linewidth=0.5) 
plt.plot(h2_2443['n'],combined_mean_h2,'.-',color=[0, 0, 1])
plt.fill_between(h2_2443['n'],list(map(float,combined_mean_h2-2*combined_h2_se)),
                 list(map(float,combined_mean_h2+2*combined_h2_se)),alpha = 0.1,
                 color=[0, 0.5, 1])
plt.fill_between(h2_2443['n'],list(map(float,combined_mean_h2-2*combined_h2_se_mean)),
                 list(map(float,combined_mean_h2+2*combined_h2_se_mean)),alpha = 0.2,
                 color=[0, 0, 1])
plt.legend(['full data h2 ref\n(%f)' % h2_ref,'combined h2','prediction interval CI','CI for mean'],loc=1)
plt.axvline(x=25e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=50e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=100e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.ylabel('h2_observed')
plt.xlabel('n_non_missing')
plt.ylim([h2_ref-0.1,h2_ref+0.1])
plt.xlim([0,np.max(h2_2443['n'])*1.01])
plt.text(15e3,(h2_ref-0.1*0.8),'n = 25k')
plt.text(40e3,(h2_ref-0.1*0.8),'n = 50k')
plt.text(88.5e3,(h2_ref-0.1*0.8),'n = 100k')
n_batches = int((h2_2443.shape[1]-1)/2)
plt.title('Combined h2_observed of 2443 ("Diabetes diagnosed by doctor")\n(%i batches)' % n_batches)
fig = plt.gcf()
fig.set_size_inches(12, 8*.7)
fig.savefig('/Users/nbaya/Desktop/upsampling_'+phen+'_h2_observed_'+str(n_batches)+'batches_combined.png',dpi=300)



###############################################################################
"""
╔══════════════════════════════════╗
║ h2part results for leg impedance ║
╚══════════════════════════════════╝
Batches: 1.1-1.50
phen code: 23107
"""

#Read h2part results file for leg impedance phenotype (23107)

phen = '23107'
batch='1'

#for i in range(1,51):
#    os.system('gsutil cp gs://nbaya/split/meta_split/h2part/ukbb31063.h2part_results.'+
#              phen+'.downsample_n300.batch_1.'+str(i)+'.tsv.gz ~/Documents/lab/ukbb-sexdiff/h2part/')

#os.system('gsutil cp gs://nbaya/split/meta_split/h2part/ukbb31063.h2part_results.23107.downsample_n300.batch_1.1_set300.tsv.gz ~/Documents/lab/ukbb-sexdiff/h2part/')


h2_23107_batch1 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.'+phen+'.downsample_n300.batch_1.1.tsv.gz', 
                              sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20].copy()[['n','h2_observed','h2_observed_se']]
h2_23107_batch1_set300 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.23107.downsample_n300.batch_1.1_set300.tsv.gz', 
                              sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20].copy()[['n','h2_observed','h2_observed_se']]
h2_23107_batch1 = h2_23107_batch1.append(h2_23107_batch1_set300.iloc[1,:].copy()[['n','h2_observed','h2_observed_se']],ignore_index=True) #only use when showing only batch 1 subbatch results

h2_23107 = h2_23107_batch1.rename(axis='columns',mapper={"h2_observed":"h2_observed_1_1","h2_observed_se":"h2_se_1_1"})
    

for i in range(2,51):
    other = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.'+phen+'.downsample_n300.batch_'+batch+'.'+str(i)+'.tsv.gz', sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20]
#    other = other.append(h2_23107_batch1_set300.iloc[1,:].copy()[['n','h2_observed','h2_observed_se']],ignore_index=True) #only use when showing only batch 1 subbatch results
    h2_23107[['h2_observed_1_'+str(i),'h2_se_1_'+str(i)]] = other[['h2_observed','h2_observed_se']]
    
#for i in range(300,301):
#    other = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.'+phen+'.downsample.batch_'+batch+'.'+str(i)+'.tsv.gz', sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20]
#    h2_20160[['h2_observed_'+str(i),'h2_se_'+str(i)]] = other[['h2_observed','h2_observed_se']]
h2 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/rg_sex/ukbb31063.both_sexes.h2part_results.phesant.tsv.gz', 
                              sep='\t',compression='gzip').rename(index=str,columns={'phenotype':'phen'}).iloc[:,0:20]
h2_ref = h2[h2['phen']==phen+'_irnt'].h2_observed[0] # full data ref

n_batches = int((h2_23107.shape[1]-1)/2)

#####
#Combined h2 across batches
#####

combined_mean_h2 = np.mean(h2_23107.filter(regex=('observed')),axis=1)
combined_h2_se = np.mean(h2_23107.filter(regex=('_se'))**2,axis=1)**(1/2)
combined_h2_se_mean = (np.mean(h2_23107.filter(regex=('_se'))**2,axis=1)/h2_23107.filter(regex=('_se')).shape[1])**(1/2)


plt.axhline(y=h2_ref,color='k',alpha=1, linewidth=0.5) 
plt.plot(h2_23107['n'],combined_mean_h2,'.-',color=[0, 0, 1])
plt.fill_between(h2_23107['n'],list(map(float,combined_mean_h2-2*combined_h2_se)),
                 list(map(float,combined_mean_h2+2*combined_h2_se)),alpha = 0.1,
                 color=[0, 0.5, 1])
plt.fill_between(h2_23107['n'],list(map(float,combined_mean_h2-2*combined_h2_se_mean)),
                 list(map(float,combined_mean_h2+2*combined_h2_se_mean)),alpha = 0.2,
                 color=[0, 0, 1])
plt.legend(['full data h2 ref\n(%f)' % h2_ref,'combined h2','prediction interval CI','CI for mean'],loc=1)
plt.axvline(x=25e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=50e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=100e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.ylabel('h2_observed')
plt.xlabel('n_non_missing')
plt.ylim([h2_ref-0.1,h2_ref+0.1])
plt.xlim([0,np.max(h2_23107['n'])*1.01])
plt.text(15e3,(h2_ref-0.1*0.8),'n = 25k')
plt.text(40e3,(h2_ref-0.1*0.8),'n = 50k')
plt.text(88.5e3,(h2_ref-0.1*0.8),'n = 100k')
n_batches = int((h2_23107.shape[1]-1)/2)
plt.title('Combined h2_observed of 23107 ("Impedance of leg (right)")\n(%i batches)' % n_batches)
fig = plt.gcf()
fig.set_size_inches(12, 8*.7)
fig.savefig('/Users/nbaya/Desktop/upsampling_'+phen+'_h2_observed_'+str(n_batches)+'batches_combined.png',dpi=300)



###############################################################################
"""
╔════════════════════════════════════════════╗
║ h2part results for Qual: None of the above ║
╚════════════════════════════════════════════╝
Batches: 1.1-1.50
"""

#Read h2part results file for Qualifications: None of the above (6138_100)

phen = '6138_100'
batch='1'

for i in range(1,51):
    if not os.path.isfile('/Users/nbaya/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.'+
              phen+'.downsample_n150.batch_1.'+str(i)+'.tsv.gz'):
        os.system('gsutil cp gs://nbaya/split/meta_split/h2part/ukbb31063.h2part_results.'+
              phen+'.downsample_n150.batch_1.'+str(i)+'.tsv.gz ~/Documents/lab/ukbb-sexdiff/h2part/')

h2_6138_100_batch1 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.'+phen+
                                 '.downsample_n150.batch_1.1.tsv.gz', sep='\t',compression='gzip'
                                 ).sort_values(by='n').reset_index().iloc[:,1:20].copy()[['n','h2_observed','h2_observed_se']]
#h2_6138_100_batch1_set150 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.6138_100.downsample_n150.batch_1.1_set150.tsv.gz', 
#                              sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20].copy()[['n','h2_observed','h2_observed_se']]
#h2_6138_100_batch1 = h2_6138_100_batch1.append(h2_6138_100_batch1_set150.iloc[1,:].copy()[['n','h2_observed','h2_observed_se']],ignore_index=True) #only use when showing only batch 1 subbatch results

h2_6138_100 = h2_6138_100_batch1.rename(axis='columns',mapper={"h2_observed":"h2_observed_1_1","h2_observed_se":"h2_se_1_1"})
    

for i in range(2,51):
    if os.path.isfile('/Users/nbaya/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.'+phen+'.downsample_n150.batch_'+batch+'.'+str(i)+'.tsv.gz'):
        other = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.'+phen+'.downsample_n150.batch_'+batch+'.'+str(i)+'.tsv.gz', sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20]
    #    other = other.append(h2_6138_100_batch1_set150.iloc[1,:].copy()[['n','h2_observed','h2_observed_se']],ignore_index=True) #only use when showing only batch 1 subbatch results
        h2_6138_100[['h2_observed_1_'+str(i),'h2_se_1_'+str(i)]] = other[['h2_observed','h2_observed_se']]
    
h2 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/rg_sex/ukbb31063.both_sexes.h2part_results.phesant.tsv.gz', 
                              sep='\t',compression='gzip').rename(index=str,columns={'phenotype':'phen'}).iloc[:,0:20]
h2_ref = h2[h2['phen']==phen].h2_observed[0] # full data ref
n_batches = int((h2_6138_100.shape[1]-1)/2)


#####
#Combined h2 across batches
#####

combined_mean_h2 = np.mean(h2_6138_100.filter(regex=('observed')),axis=1)
combined_h2_se = np.mean(h2_6138_100.filter(regex=('_se'))**2,axis=1)**(1/2)
combined_h2_se_mean = (np.mean(h2_6138_100.filter(regex=('_se'))**2,axis=1)/h2_6138_100.filter(regex=('_se')).shape[1])**(1/2)


plt.axhline(y=h2_ref,color='k',alpha=1, linewidth=0.5) 
plt.plot(h2_6138_100['n'],combined_mean_h2,'.-',color=[0, 0, 1])
plt.fill_between(h2_6138_100['n'],list(map(float,combined_mean_h2-2*combined_h2_se)),
                 list(map(float,combined_mean_h2+2*combined_h2_se)),alpha = 0.1,
                 color=[0, 0.5, 1])
plt.fill_between(h2_6138_100['n'],list(map(float,combined_mean_h2-2*combined_h2_se_mean)),
                 list(map(float,combined_mean_h2+2*combined_h2_se_mean)),alpha = 0.2,
                 color=[0, 0, 1])
plt.legend(['full data h2 ref\n(%f)' % h2_ref,'combined h2','prediction interval CI','CI for mean'],loc=1)
plt.axvline(x=25e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=50e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=100e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.ylabel('h2_observed')
plt.xlabel('n_non_missing')
plt.ylim([h2_ref-0.1,h2_ref+0.1])
plt.xlim([0,np.max(h2_6138_100['n'])*1.01])
plt.text(15e3,(h2_ref-0.1*0.8),'n = 25k')
plt.text(40e3,(h2_ref-0.1*0.8),'n = 50k')
plt.text(88.5e3,(h2_ref-0.1*0.8),'n = 100k')
plt.title('Combined h2_observed of 6138_100 ('+str(h2[h2['phen']==phen].description[0])+')\n(%i batches)' % n_batches)
fig = plt.gcf()
fig.set_size_inches(12, 8*.7)
fig.savefig('/Users/nbaya/Desktop/upsampling_'+phen+'_h2_observed_'+str(n_batches)+'_batches_combined.png',dpi=300)


###############################################################################
"""
╔═══════════════════════════╗
║ h2part results for 50_raw ║
╚═══════════════════════════╝
Batches: 1.1-1.50
"""

#Read h2part results file for Qualifications: None of the above (6138_100)

phen = '50_raw'
n_chunks = 300

cloud_wd = 'gs://nbaya/split/meta_split/h2part/'
wd = '/Users/nbaya/Documents/lab/ukbb-sexdiff/h2part/'
for i in range(1,51):
    filename = 'ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'.batch_1.'+str(i)+'.tsv.gz'
    if not os.path.isfile(wd+filename):
        os.system('gsutil cp '+cloud_wd+filename+' '+wd)

h2_50_raw = pd.read_csv('~/Documents/lab/ukbb-sexdiff/h2part/ukbb31063.h2part_results.'+phen+
                                 '.downsample_n'+str(n_chunks)+'.batch_1.1.tsv.gz', sep='\t',compression='gzip'
                                 ).sort_values(by='n').reset_index().iloc[:,1:20].copy(
                                         )[['n','h2_observed','h2_observed_se']].rename(
                                                 axis='columns',mapper={"h2_observed":"h2_observed_1_1","h2_observed_se":"h2_se_1_1"})

for i in range(2,51):
    filename = 'ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'.batch_1.'+str(i)+'.tsv.gz'
    if os.path.isfile(wd+filename):
        print(i)
        other = pd.read_csv(wd+filename, sep='\t',compression='gzip').sort_values(by='n').reset_index().iloc[:,1:20]
    #    other = other.append(h2_50_raw_batch1_set150.iloc[1,:].copy()[['n','h2_observed','h2_observed_se']],ignore_index=True) #only use when showing only batch 1 subbatch results
        h2_50_raw[['h2_observed_1_'+str(i),'h2_se_1_'+str(i)]] = other[['h2_observed','h2_observed_se']]
    
h2 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/rg_sex/ukbb31063.both_sexes.h2part_results.phesant.tsv.gz', 
                              sep='\t',compression='gzip').rename(index=str,columns={'phenotype':'phen'}).iloc[:,0:20]
h2_ref = h2[h2['phen']==phen].h2_observed[0] # full data ref
n_batches = int((h2_50_raw.shape[1]-1)/2)


#####
#Combined h2 across batches
#####

combined_mean_h2 = np.mean(h2_50_raw.filter(regex=('observed')),axis=1)
combined_h2_se = np.mean(h2_50_raw.filter(regex=('_se'))**2,axis=1)**(1/2)
combined_h2_se_mean = (np.mean(h2_50_raw.filter(regex=('_se'))**2,axis=1)/h2_50_raw.filter(regex=('_se')).shape[1])**(1/2)


plt.axhline(y=h2_ref,color='k',alpha=1, linewidth=0.5) 
plt.plot(h2_50_raw['n'],combined_mean_h2,'.-',color=[0, 0, 1])
plt.fill_between(h2_50_raw['n'],list(map(float,combined_mean_h2-2*combined_h2_se)),
                 list(map(float,combined_mean_h2+2*combined_h2_se)),alpha = 0.1,
                 color=[0, 0.5, 1])
plt.fill_between(h2_50_raw['n'],list(map(float,combined_mean_h2-2*combined_h2_se_mean)),
                 list(map(float,combined_mean_h2+2*combined_h2_se_mean)),alpha = 0.2,
                 color=[0, 0, 1])
plt.legend(['full data h2 ref\n(%f)' % h2_ref,'combined h2','prediction interval CI','CI for mean'],loc=1)
plt.axvline(x=25e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=50e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.axvline(x=100e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
plt.ylabel('h2_observed')
plt.xlabel('n_non_missing')
plt.ylim([h2_ref-0.1,h2_ref+0.1])
plt.xlim([0,np.max(h2_50_raw['n'])*1.01])
plt.text(15e3,(h2_ref-0.1*0.8),'n = 25k')
plt.text(40e3,(h2_ref-0.1*0.8),'n = 50k')
plt.text(88.5e3,(h2_ref-0.1*0.8),'n = 100k')
plt.title('Combined h2_observed of 50_raw ('+str(h2[h2['phen']==phen].description[0])+')\n(%i batches)' % n_batches)
fig = plt.gcf()
fig.set_size_inches(12, 8*.7)
fig.savefig('/Users/nbaya/Desktop/upsampling_'+phen+'_h2_observed_batches_combined.png',dpi=300)

###############################################################################
"""
╔═════════════════════════╗
║ Read genomicSEM results ║
╚═════════════════════════╝
Available phenotypes: 50, 30100
"""
gsem_rg_meta_50 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/rg_sex/genomicSEM_rg_meta_50.tsv', sep='\t') #n chunks = 300
gsem_rg_meta_30100 = pd.read_csv('~/Documents/lab/ukbb-sexdiff/rg_sex/genomicSEM_rg_meta_30100.tsv', sep='\t')
gsem_rg_nonmeta_50 = pd.DataFrame([0.9972281, 0.9982607, 1.0014674, 0.9901949, 0.9875092, 0.9795326, 0.9883829], columns=['x'])
gsem_rg_meta_n150_50 = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/genomicSEM_rg_meta_n150_50.tsv', sep='\t')

#only compare '50' and '30100' results from genomicSEM
sns.kdeplot(gsem_rg_meta_50.x, color='#ff7f0e')
sns.kdeplot(gsem_rg_meta_30100.x, color='#2ca02c')
plt.xlabel('rg')
plt.ylabel('density')
plt.legend(['50 rg','30100 rg'])
plt.title('Comparison of genomicSEM rg distr. \n 50 vs. 30100')
fig = plt.gcf()
fig.set_size_inches(6,4)
fig.savefig('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/plots/genomicSEM_50_30100_comparison.png',dpi=300)

#compare stock method to genomicSEM for '50'
sns.kdeplot(rg_df[rg_df.p1.str.contains('50_meta_A_batch')].rg, color='#ff7f0e')
sns.kdeplot(gsem_rg_meta_50.x, color='#ff7f0e', linestyle='--')
plt.legend(['stock ldsc','genomicSEM'])
plt.title('Comparison of stock ldsc to genomicSEM\nfor 50')
fig = plt.gcf()
fig.set_size_inches(6,4)
fig.savefig('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/plots/genomicSEM_50_comparison.png',dpi=300)

#compare stock method to genomicSEM for '30100'
sns.kdeplot(rg_df[rg_df.p1.str.contains('30100')].rg, color='#2ca02c')
sns.kdeplot(gsem_rg_meta_30100.x, color='#2ca02c', linestyle='--')
plt.legend(['stock ldsc','genomicSEM'])
plt.title('Comparison of stock ldsc to genomicSEM\nfor 30100')
fig = plt.gcf()
fig.set_size_inches(6,4)
fig.savefig('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/plots/genomicSEM_30100_comparison.png',dpi=300)

#compare '50' and '30100' for stock ldsc and genomicSEM
sns.kdeplot(rg_df[rg_df.p1.str.contains('50_meta_A_batch')].rg, color='#ff7f0e')
sns.kdeplot(gsem_rg_meta_50.x, color='#ff7f0e', linestyle='--')
sns.kdeplot(rg_df[rg_df.p1.str.contains('30100')].rg, color='#2ca02c')
sns.kdeplot(gsem_rg_meta_30100.x, color='#2ca02c', linestyle='--')  
plt.legend(['50 stock','50 genomicSEM','30100 stock', '30100 genomicSEM'])
plt.title('Comparison of stock ldsc to genomicSEM\nfor 50 and 30100')
fig = plt.gcf()
fig.set_size_inches(6,4)
fig.savefig('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/plots/stock_genomicSEM_50_30100_comparison.png',dpi=300)

#compare nonmeta to stock method to genomicSEM for '50'
sns.kdeplot(rg_split_50.rg, color='#ff7f0e')
sns.kdeplot(gsem_rg_nonmeta_50.x, color='#ff7f0e', linestyle='--')
plt.legend(['stock ldsc', 'genomicSEM', ])
plt.title('Comparison of stock ldsc to genomicSEM\nfor non-meta-analysis version of 50')
fig = plt.gcf()
fig.set_size_inches(6,4)
fig.savefig('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/plots/stock_genomicSEM_nonmeta_50_comparison.png',dpi=300)

#compare n300 genomicSEM to n150 genomicSEM for '50'
sns.kdeplot(gsem_rg_meta_50.x, color='#ff7f0e')
sns.kdeplot(gsem_rg_meta_n150_50.x, color='#ff7f0e', linestyle='--')
plt.legend(['n_chunks = 300', 'n_chunks = 150', ])
plt.title('Comparison of 300 to 150 chunks using genomicSEM\nfor meta-analysis version of 50')
fig = plt.gcf()
fig.set_size_inches(6,4)
fig.savefig('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/plots/genomicSEM_meta_50_n300_n150_comparison.png',dpi=300)

###############################################################################
"""
╔══════════════════════════════════╗
║ h2part results for any phenotype ║
╚══════════════════════════════════╝
Available phenotypes: 50, 20160, 2443, 23107, 6138_100, 50_raw, 30100, 50 n150
Batches: (1) and 1.1-1.50
"""
#Format:       {phenkey:[phen,  n_chunks, hasbatch1, haslastset, isirnt, isnested]} #hasbatch1 applies only to the first phenotypes we ran through (50, 20160). They used different random seeds than the batch 1 subbatches.
phendict = {       '50':['50',               300, True,  True,  True, float('NaN')], #there are also n150 results (fewer chunks, more samples in each chunk) for '50' (50_irnt) but only for the full set (i.e. set 150)
                '20160':['20160',            300, True,  False, False, float('NaN')],
                 '2443':['2443',             150, False, True,  False, float('NaN')],
                '23107':['23107',            300, False, True,  True, float('NaN')],
             '6138_100':['6138_100',         150, False, False, False, float('NaN')],
               '50_raw':['50_raw',           300, False, False, False, float('NaN')],
                '30100':['30100',            300, False, False, True, float('NaN')],
              '50_n150':['50',               150, False, True,  True, float('NaN')],
           '50_sim_inf':['50_sim_inf',       300, False, False, True, float('NaN')],
      '50_sim_inf_n150':['50_sim_inf',       150, False, False, True, float('NaN')],
    '50_sim_inf_h2_0.1':['50_sim_inf_h2_0.1',300, False, False, True, True],
       '50_nestedFalse':['50',               300, False, False, True, False],
           '50_raw_res':['50_raw_res',       300, False, False, False, True]
       }

def geth2batches(phenkey, phendict):
    cloud_wd = 'gs://nbaya/split/meta_split/h2part/'
    wd = '/Users/nbaya/Documents/lab/ukbb-sexdiff/h2part/downsample/'
    phen = phendict[phenkey][0]
    n_chunks = phendict[phenkey][1]
    hasbatch1 = phendict[phenkey][2]
    haslastset = phendict[phenkey][3]
    isnested = phendict[phenkey][5]
    
    #Download files if necessary
    if hasbatch1: #if the phenotype has a batch 1 result (which didn't use the subbatch format)
        filename='ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'.batch_1.tsv.gz'
        if not os.path.isfile(wd+filename):
            if os.system('gsutil cp '+cloud_wd+filename+' '+wd)==0:
                print('Downloaded '+filename+' \n        to '+wd) 
    for i in range(1,51):
        if isnested is True: 
            filename='ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'_nestedTrue.batch_1.'+str(i)+'.tsv.gz'
        elif isnested is False:
            filename='ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'_nestedFalse.batch_1.'+str(i)+'.tsv.gz'
        else: #include phens with filenames made before specifying nested=True/False. All phens are nested if not specified.
            filename='ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'.batch_1.'+str(i)+'.tsv.gz'
        if not os.path.isfile(wd+filename):
            if os.system('gsutil cp '+cloud_wd+filename+' '+wd)==0:
                print('Downloaded '+filename+' \n        to '+wd) 
    if haslastset: #if phen has last set (either set 300 or set 150 if n_chunks is 300 or 150 respectively)
        filename='ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'.batch_1.1_set'+str(n_chunks)+'.tsv.gz'
        if not os.path.isfile(wd+filename):
            if os.system('gsutil cp '+cloud_wd+filename+' '+wd)==0:
                print('Downloaded '+filename+' \n        to '+wd) 
    
    #Read into dataframes
    if isnested is True: 
        filename='ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'_nestedTrue.batch_1.1.tsv.gz'
    elif isnested is False:
        filename='ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'_nestedFalse.batch_1.1.tsv.gz'
    else: #include phens with filenames made before specifying nested=True/False. All phens are nested if not specified.
        filename='ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'.batch_1.1.tsv.gz'

    batch1_1 = pd.read_csv(wd+filename, sep='\t',compression='gzip').sort_values(by='n').reset_index()[['n','h2_observed','h2_observed_se','intercept']]
    
    if haslastset: #if phen has last set (either set 300 or set 150 if n_chunks is 300 or 150 respectively)
        filename='ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'.batch_1.1_set'+str(n_chunks)+'.tsv.gz'
        lastset = pd.read_csv(wd+filename, sep='\t',compression='gzip').sort_values(by='n').reset_index()[['n','h2_observed','h2_observed_se','intercept']]
        batch1_1 = batch1_1.append(lastset.iloc[1,:][['n','h2_observed','h2_observed_se','intercept']],ignore_index=True) #only use when showing only batch 1 subbatch results

    batch1_1 = batch1_1.rename(axis='columns',mapper={"h2_observed":"h2_observed_1_1","h2_observed_se":"h2_se_1_1"})

    if hasbatch1: #if the phenotype has a batch 1 result (only '50' and '20160')
        h2_batches = pd.read_csv(wd+'ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'.batch_1.tsv.gz',
                                   sep='\t',compression='gzip').sort_values(by='n').reset_index()[['n','h2_observed','h2_observed_se','intercept']]
        if haslastset: #if phen has last set (either set 300 or set 150 if n_chunks is 300 or 150 respectively)
            filename='ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'.batch_1.1_set'+str(n_chunks)+'.tsv.gz'
            lastset = pd.read_csv(wd+filename, sep='\t',compression='gzip').sort_values(by='n').reset_index()[['n','h2_observed','h2_observed_se','intercept']]
            h2_batches = h2_batches.append(lastset.iloc[1,:][['n','h2_observed','h2_observed_se','intercept']],ignore_index=True) #only use when showing only batch 1 subbatch results
        h2_batches = h2_batches.rename(axis='columns',mapper={"h2_observed":"h2_observed_1","h2_observed_se":"h2_se_1"})
        
        h2_batches[['h2_observed_1_1','h2_se_1_1']] = batch1_1[['h2_observed_1_1','h2_se_1_1']]
    else:
        h2_batches = batch1_1    
    
    for i in range(2,51):
        if isnested is True: 
            filename='ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'_nestedTrue.batch_1.'+str(i)+'.tsv.gz'
        elif isnested is False:
            filename='ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'_nestedFalse.batch_1.'+str(i)+'.tsv.gz'
        else: #include phens with filenames made before specifying nested=True/False. All phens are nested if not specified.
            filename='ukbb31063.h2part_results.'+phen+'.downsample_n'+str(n_chunks)+'.batch_1.'+str(i)+'.tsv.gz'
        if os.path.isfile(wd+filename):
            temp = pd.read_csv(wd+filename, sep='\t',compression='gzip').sort_values(by='n').reset_index()
            if haslastset:
                temp = temp.append(lastset.iloc[1,:].copy()[['n','h2_observed','h2_observed_se','intercept']],ignore_index=True) #only use when showing only batch 1 subbatch results
            h2_batches[['h2_observed_1_'+str(i),'h2_se_1_'+str(i)]] = temp[['h2_observed','h2_observed_se']]
            
        else:
            print('Missing: '+filename)

    return h2_batches
    
def ploth2upsampling(phenkey, phendict, plotlastset, uselastsetref, saveplot):
    wd = '/Users/nbaya/Documents/lab/ukbb-sexdiff/h2part/'
    h2_batches = geth2batches(phenkey, phendict)
    phen = phendict[phenkey][0]
    n_chunks = phendict[phenkey][1]
    haslastset = phendict[phenkey][3]
    isirnt = phendict[phenkey][4]
    isnested = phendict[phenkey][5]
    
    #Get h2 reference and n_batches
    h2 = pd.read_csv('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/ukbb31063.both_sexes.h2part_results.phesant.tsv.gz', 
                                      sep='\t',compression='gzip').rename(index=str,columns={'phenotype':'phen'}).iloc[:,0:20]
    if uselastsetref:
        if haslastset:
            h2_ref = h2_batches.filter(regex='observed').iloc[h2_batches.shape[0]-1,0] #grab last set h2 for the first h2_observed column
        else: 
            print('Phenotype '+phen+' does not have set '+str(n_chunks)+' (last set)')
            print('Plotting without set '+str(n_chunks)+' (last set)')
            uselastsetref = False
    if not uselastsetref: #use another if instead of else to include phenotypes that had uselastsetref set to false in the previous 'if' block
        if isirnt: #if the phenotype is irnt
            if phen == '50_sim_inf':
                h2_ref = h2[h2['phen']=='50_irnt'].h2_observed[0] # full data ref
            elif phen == '50_sim_inf_h2_0.1':
                h2_ref = 0.1
            else:
                h2_ref = h2[h2['phen']==phen+'_irnt'].h2_observed[0] # full data ref    
        else:
            if phen == '50_raw_res':
                h2_ref = h2[h2['phen']=='50_irnt'].h2_observed[0] # full data ref
            else:
                h2_ref = h2[h2['phen']==phen].h2_observed[0] # full data ref
    n_batches = int((h2_batches.shape[1]-1)/2)
    
    #Plot combined h2
    combined_mean_h2 = np.mean(h2_batches.filter(regex=('observed')),axis=1)
    combined_h2_se = np.mean(h2_batches.filter(regex=('_se'))**2,axis=1)**(1/2)
    combined_h2_se_mean = (np.mean(h2_batches.filter(regex=('_se'))**2,axis=1)/h2_batches.filter(regex=('_se')).shape[1])**(1/2)
    
    
    plt.axhline(y=h2_ref,color='k',alpha=1, linewidth=0.5) 
    plt.plot(h2_batches['n'],combined_mean_h2,'.-',color=[0, 0, 1])
    plt.fill_between(h2_batches['n'],list(map(float,combined_mean_h2-2*combined_h2_se)),
                     list(map(float,combined_mean_h2+2*combined_h2_se)),alpha = 0.1,
                     color=[0, 0.5, 1])
    plt.fill_between(h2_batches['n'],list(map(float,combined_mean_h2-2*combined_h2_se_mean)),
                     list(map(float,combined_mean_h2+2*combined_h2_se_mean)),alpha = 0.2,
                     color=[0, 0, 1])
    if uselastsetref:
        plt.legend(['set '+str(n_chunks)+' (full set) h2 ref\n(%f)' % h2_ref,'combined h2','prediction interval CI','CI for mean'],loc=1)
    else:
        plt.legend(['full data h2 ref\n(%f)' % h2_ref,'combined h2','prediction interval CI','CI for mean'],loc=1)
    plt.axvline(x=25e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
    plt.axvline(x=50e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
    plt.axvline(x=100e3,color='r',alpha=1, linewidth=0.5,linestyle='--')
    plt.ylabel('h2_observed')
    plt.xlabel('n_non_missing')
    ylim_param = 0.2 #default: 0.1
    plt.ylim([h2_ref-ylim_param,h2_ref+ylim_param])
    if plotlastset or not haslastset: 
        plt.xlim([0,np.max(h2_batches['n'])*1.01])
    else:
        plt.xlim([0,np.max(h2_batches[h2_batches['n']<np.max(h2_batches['n'])]['n'])*1.02])
    plt.text(15e3,(h2_ref-ylim_param*0.8),'n = 25k')
    plt.text(40e3,(h2_ref-ylim_param*0.8),'n = 50k')
    plt.text(88.5e3,(h2_ref-ylim_param*0.8),'n = 100k')
    if isirnt: #if the phenotype is irnt
        if 'sim' in phen: 
            plt.title('Combined h2_observed of '+phen+' (sim h2='+str(h2_ref)+')\n('+str(n_batches)+' batches, n chunks = '+str(n_chunks)+', nested ='+str(isnested)+')')
        else:
            plt.title('Combined h2_observed of '+phen+' ('+str(h2[h2['phen']==phen+'_irnt'].description[0])+')\n('+str(n_batches)+' batches, n chunks = '+str(n_chunks)+', nested ='+str(isnested)+')')
    else:
        if 'res' in phen:
            plt.title('Combined h2_observed of '+phen+' ('+str(h2[h2['phen']=='50_raw'].description[0])+')\n('+str(n_batches)+' batches, n chunks = '+str(n_chunks)+', nested ='+str(isnested)+')')
        else:
            plt.title('Combined h2_observed of '+phen+' ('+str(h2[h2['phen']==phen].description[0])+')\n('+str(n_batches)+' batches, n chunks = '+str(n_chunks)+', nested ='+str(isnested)+')')
    fig = plt.gcf()
    fig.set_size_inches(12, 8*.7)
    if saveplot:
        if plotlastset:
            if uselastsetref:
                fig.savefig(wd+'plots/'+'upsampling_'+phen+'_n'+str(n_chunks)+'_nested'+str(isnested)+'_h2_observed_'+str(n_batches)+'batches_combined_includelastset_lastsetref.png',dpi=300)
            else:
                fig.savefig(wd+'plots/'+'upsampling_'+phen+'_n'+str(n_chunks)+'_nested'+str(isnested)+'_h2_observed_'+str(n_batches)+'batches_combined_includelastset.png',dpi=300)
        else:
            fig.savefig(wd+'plots/'+'upsampling_'+phen+'_n'+str(n_chunks)+'_nested'+str(isnested)+'_h2_observed_'+str(n_batches)+'batches_combined.png',dpi=300)

h2_phen = geth2batches('50_raw_res',phendict)
ploth2upsampling('50_raw_res',phendict,plotlastset=False,uselastsetref=False,saveplot=False)
