#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 09:46:37 2018

@author: nbaya
"""

import hail as hl

phen='50'
variant_set='hm3'
n_chunks=300
batch='1'

mt = hl.read_matrix_table('gs://nbaya/split/ukb31063.'+variant_set+'_variants.gwas_samples_'+phen+'_grouped'+str(n_chunks)+'_batch_'+batch+'.mt') 

mt1 = mt.annotate_entries(gt = hl.int(hl.int(mt.dosage*3/2)*2/3))

mt1 = mt1.annotate_entries(GT = hl.call(mt1.gt))

hl.identity_by_descent(mt1).write('gs://nbaya/split/ibd.'+variant_set+'_variants.'+phen+'_grouped'+str(n_chunks)+'_batch_'+batch+'.ht')

