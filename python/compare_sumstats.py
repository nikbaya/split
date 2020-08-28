#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 13:15:22 2018

Used to generate QQ plots of rg split sumstats.

@author: nbaya
"""

import hail as hl
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

splitA = hl.import_table('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/sumstats/50_sim_inf_meta_A_n300_batch_1_s0.tsv.bgz', impute=True)
splitB = hl.import_table('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/sumstats/50_sim_inf_meta_B_n300_batch_1_s0.tsv.bgz', impute=True)

splitA = hl.import_table('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/sumstats/vds1_5.tsv.gz', impute=True)
splitB = hl.import_table('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/sumstats/vds2_5.tsv.gz', impute=True)

splitA = hl.import_table('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/sumstats/50_sim_inf_meta_A_n300_batch_1_s0_old.tsv.bgz', impute=True)
splitB = hl.import_table('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/sumstats/50_sim_inf_meta_B_n300_batch_1_s0_old.tsv.bgz', impute=True)

splitA = hl.import_table('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/sumstats/20160_meta_A_batch_1_s0.tsv.bgz', impute=True)
splitB = hl.import_table('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/sumstats/20160_meta_B_batch_1_s0.tsv.bgz', impute=True)

splitA = hl.import_table('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/sumstats/50_sim_inf_sample_A_batch_1.1_set84.tsv.bgz', impute=True)
splitB = hl.import_table('/Users/nbaya/Documents/lab/ukbb-sexdiff/rg_sex/sumstats/50_sim_inf_sample_A_batch_1.2_set84.tsv.bgz', impute=True)

df_A = splitA.to_pandas()
df_B = splitB.to_pandas()

mean_chi2_A = np.mean(df_A.Z*df_A.Z)
mean_chi2_B = np.mean(df_B.Z*df_B.Z)

df_A['P'] = stats.norm.sf(abs(df_A.Z))*2
df_B['P'] = stats.norm.sf(abs(df_B.Z))*2

plt.subplot(2,1,1)
plt.plot(-np.log10(np.linspace(1,1/df_A.shape[0],df_A.shape[0])), -np.log10(df_A.sort_values(by='P',ascending=False).P),'o', alpha=0.5)
plt.plot([0,10],[0,10],'k--')
plt.xlim([0,np.max(-np.log10(np.linspace(1,1/df_A.shape[0],df_A.shape[0])))*1.05])
plt.ylim([0,np.max(-np.log10(df_A.sort_values(by='P',ascending=False).P))*1.05])

plt.subplot(2,1,2)
plt.plot(-np.log10(np.linspace(1,1/df_B.shape[0],df_B.shape[0])), -np.log10(df_B.sort_values(by='P',ascending=False).P),'o', alpha=0.5)
plt.plot([0,10],[0,10],'k--')
plt.xlim([0,np.max(-np.log10(np.linspace(1,1/df_B.shape[0],df_B.shape[0])))*1.05])
plt.ylim([0,np.max(-np.log10(df_B.sort_values(by='P',ascending=False).P))*1.05])
fig = plt.gcf()
fig.set_size_inches(8, 12)
