#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 16 16:19:00 2019

Runs a GWAS on a matrix table, given the dependent/independent variables and
list of covariates. 
Options:    list of covariates for linreg, list of variables to pass through,
            path to save


@author: nbaya
"""

import hail as hl

def normalize_genotypes(mt, genotype):
    '''Normalizes genotypes'''
    print('\rNormalizing genotypes...')
    mt1 = mt.annotate_entries(__gt = genotype)
    mt2 = mt1.annotate_rows(__gt_stats = hl.agg.stats(mt1.__gt))
    return mt2.annotate_entries(__norm_gt = (mt2.__gt-mt2.__gt_stats.mean)/mt2.__gt_stats.stdev)   

def gwas(mt, x, y, cov_list=[], with_intercept=True, normalize_x=True, pass_through=[],
         path_to_save=None, is_std_cov_list=True):
    '''Runs GWAS'''
    header = '\n****************************************\n'
    header += 'Running GWAS\n'
    header += 'Covariates list: {}\n'.format(cov_list if not is_std_cov_list else 'Standard 26 covariates')
    header += 'With intercept: {}\n'.format(with_intercept)
    header += 'Normalizing X: {}\n'.format(normalize_x)
    header += '' if pass_through is [] else 'Variables to pass through: {}\n'.format(pass_through)
    header += '' if path_to_save is None else 'Saving to: {}\n'.format(path_to_save)
    header += '****************************************'
    print(header)
    
    
    mt = mt._annotate_all(col_exprs={'y':y},
                           entry_exprs={'x':x})
    if normalize_x:
        mt = normalize_genotypes(mt, mt.x) 
        mt = mt.annotate_entries(x = mt.__norm_gt).drop('__norm_gt')
    
    if is_std_cov_list:
        cov_list = ['isFemale','age','age_squared','age_isFemale',
                    'age_squared_isFemale']+['PC{:}'.format(i) for i in range(1, 21)]
        
    if str in list(map(lambda x: type(x),cov_list)):
        cov_list = list(map(lambda x: mt[x] if type(x) is str else x,cov_list))
        
    cov_list = ([1] if with_intercept else [])+cov_list

    gwas_ht = hl.linear_regression_rows(y=mt.y,
                                        x=mt.x,
                                        covariates=cov_list,
                                        pass_through = pass_through)
    gwas_ht = gwas_ht.annotate_globals(with_intercept = with_intercept)
    
    if path_to_save is not None:
        gwas_ht.write(path_to_save)
        
    return gwas_ht
    
if __name__ == '__main__':
    gs_bucket = 'gs://nbaya/ldscsim/'
    mt0 = hl.read_matrix_table(gs_bucket+'hm3.50_sim_h2_0.08.mt/')
    ht0 = mt0.cols()
    mt  = hl.read_matrix_table(gs_bucket+'ukbb31063.GT.chr22.MAF_0.05.white_british.mt/')
    mt = mt.annotate_cols(PC1 = ht0[mt.s]['PC1'])
    
    gwas_ht = gwas(mt, mt.GT.n_alt_alleles(),)
    
    
