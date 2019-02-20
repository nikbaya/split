#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 08:37:36 2018

@author: nbaya
"""

import hail as hl
import numpy as np
import datetime
#from multiprocessing import Pool


phen = '50'
desc = 'height'
batch = '1'

#Globals
n = 300 #number of subgroups

"""
###################
Part 1: Generate subset
More precisely, this determines the order in which group IDs are added to the subset
###################
"""
groups = list(range(n))
randstate = np.random.RandomState(int(batch))
randstate.shuffle(groups)

"""
###################
Part 2: Calculate meta summary statistics (via inverse-variance weighting meta-analysis)
###################
"""
print('Starting Part 5: Splitting into populations A and B, calculatng summary statistics')
print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))

gmt = hl.read_matrix_table('gs://nbaya/split/meta_split/ukb31063.hm3_'+phen+'_gmt_batch_'+batch+'.mt')
gmt = gmt.add_col_index()
gmt = gmt.key_cols_by('group_id')
gmt = gmt.rename({'rsid': 'SNP'})

for subset_i in range(10):  #restricted range (< n) to start with
    subset_i += 1
    pi = np.array(['B']*n)
    pi[groups[0:subset_i]] = ['A']*(subset_i) #label the groups that belong to the subset
    pi = pi.tolist()
    gmt_lab = gmt.annotate_cols(label = hl.literal(pi)[hl.int32(gmt.col_idx)])
    mt = gmt_lab.group_cols_by(gmt_lab.label).aggregate(
            unnorm_meta_beta=hl.agg.sum(gmt_lab.beta / gmt_lab.standard_error ** 2),
            inv_se2 = hl.agg.sum(1 / gmt_lab.standard_error ** 2))
    