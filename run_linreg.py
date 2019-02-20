#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 14:47:48 2018

Used to generate residuals for a given phenotype in the UKB.
Linreg does not use dosage.

@author: nbaya
"""

import hail as hl
import numpy as np

phen='50_raw'
batch='1'
variant_set = 'hm3'

mt0 = hl.read_matrix_table('gs://nbaya/split/ukb31063.'+variant_set+'_variants.gwas_samples_repart.mt')

if phen == '50_raw':
    phen_tb_all = hl.import_table('gs://nbaya/ukb31063.50_raw.tsv.bgz',missing='',impute=True,types={'s': hl.tstr}).rename({phen: 'phen'})
else:
    phen_tb_all = hl.import_table('gs://phenotype_31063/ukb31063.phesant_phenotypes.both_sexes.tsv.bgz',
                                  missing='',impute=True,types={'"userId"': hl.tstr}).rename({ '"userId"': 's', '"'+phen+'"': 'phen'})

#Remove withdrawn samples
withdrawn = hl.import_table('gs://nbaya/w31063_20181016.csv',missing='',no_header=True)
withdrawn_set = set(withdrawn.f0.take(withdrawn.count()))
phen_tb_all.describe()
phen_tb_all = phen_tb_all.filter(hl.literal(withdrawn_set).contains(phen_tb_all['s']),keep=False) 

phen_tb = phen_tb_all.select(phen_tb_all['s'], phen_tb_all['phen'])
phen_tb = phen_tb.key_by('s')
phen_tb.count() #expect to be 361194 - 50 = 361144 for every phenotype before filtering, minus withdrawn samples 

mt1 = mt0.annotate_cols(phen_str = hl.str(phen_tb[mt0.s]['phen']).replace('\"',''))
mt1 = mt1.filter_cols(mt1.phen_str == '',keep=False)

if phen_tb.phen.dtype == hl.dtype('bool'):
    mt1 = mt1.annotate_cols(phen = hl.bool(mt1.phen_str)).drop('phen_str')
else:
    mt1 = mt1.annotate_cols(phen = hl.float64(mt1.phen_str)).drop('phen_str')


n_samples = mt1.count_cols()
print('\n>>> N samples = '+str(n_samples)+' <<<') #expect n samples to match n_non_missing from phenotypes.both_sexes.tsv
mt2 = mt1.add_col_index()

mt = mt1.rename({'dosage': 'x', 'phen': 'y'})
    
mt = mt1.rename({'dosage': 'x', 'phen': 'y'})
    
ht = mt.cols()
    
cov_list = [ ht['isFemale'], ht['age'], ht['age_squared'], ht['age_isFemale'],
            ht['age_squared_isFemale'] ]+ [ht['PC{:}'.format(i)] for i in range(1, 21)] 


linreg_results = ht.aggregate(hl.agg.linreg(y=ht.y, x = [1] + cov_list))
ht = ht.annotate(linreg = linreg_results)

ht = ht.annotate(y_pred =   linreg_results.beta[0]                +linreg_results.beta[1]*ht.isFemale    +linreg_results.beta[2]*ht.age+
                            linreg_results.beta[3] *ht.age_squared+linreg_results.beta[4]*ht.age_isFemale+linreg_results.beta[5]*ht.age_squared_isFemale+
                            linreg_results.beta[6] *ht.PC1        +linreg_results.beta[7]*ht.PC2         +linreg_results.beta[8]*ht.PC3+
                            linreg_results.beta[9] *ht.PC4        +linreg_results.beta[10]*ht.PC5        +linreg_results.beta[11]*ht.PC6+
                            linreg_results.beta[12]*ht.PC7        +linreg_results.beta[13]*ht.PC8        +linreg_results.beta[14]*ht.PC9+
                            linreg_results.beta[15]*ht.PC10       +linreg_results.beta[16]*ht.PC11       +linreg_results.beta[17]*ht.PC12+
                            linreg_results.beta[18]*ht.PC13       +linreg_results.beta[19]*ht.PC14       +linreg_results.beta[20]*ht.PC15+
                            linreg_results.beta[21]*ht.PC16       +linreg_results.beta[22]*ht.PC17       +linreg_results.beta[23]*ht.PC18+
                            linreg_results.beta[24]*ht.PC19       +linreg_results.beta[25]*ht.PC20)

ht = ht.annotate(res = ht.y-ht.y_pred)

ht.write('gs://nbaya/split/'+phen+'_linreg.ht')