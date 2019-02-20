#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 15:04:54 2018

@author: nbaya
"""

import json
import gzip
import pandas as pd
import subprocess



phen = '20160'
desc = 'smoking'


i_start = 1     # iteration index to start at (change to prevent overwriting previous samples)
i_n = 4         # number of iterations

for sample in range(i_start,i_start+i_n):
    for i in ["A","B"]:
        filename = phen+'_vds'+ i + '_s' + str(sample)
        dlpath_cloud = 'gs://nbaya/split/'
        path_local = '/Users/nbaya/Documents/lab/ukbb-sexdiff/split/'
    
#        print('Starting conversion for sample '+ str(sample))
#        subprocess.call(['gsutil','cp',dlpath_cloud + filename + '.tsv.bgz', path_local])
#        
#        with gzip.open(path_local + filename + '.tsv.bgz', 'r') as f:
#            f.readline()
#            rows1 = []
#            rows2 = []
#            for line in f:
#                v, va = line.strip().split(b'\t')
#                json_line2 = json.loads(va)
#                rows1.append(v.decode('utf-8'))
#                rows2.append(json_line2)
#    
#        rsid = list(map(lambda x: x['rsid'], rows2))
#        A1 = list(map(lambda x: rows1[x].split(':')[3], range(len(rows1))))
#        A2 = list(map(lambda x: rows1[x].split(':')[2], range(len(rows1))))
#    
#        df_linreg = pd.DataFrame(list(map(lambda x: (x['linreg']), rows2)))
#        N = df_linreg['nCompleteSamples']
#        Z = list(map(lambda x: x[0], df_linreg['tstat']))
#    
#        df = pd.DataFrame({'SNP': rsid, 'A1': A1, 'A2': A2, 'N':N, 'Z': Z})
#    
#        file_path_local = path_local + filename + '.tsv.gz'
#        df.to_csv(file_path_local, sep='\t', compression='gzip', index=False)
#        path_cloud = 'gs://nbaya/split/test/'+filename+'.tsv.gz'
#        subprocess.call(['gsutil','cp',file_path_local, path_cloud])
#        print('Finished conversion for sample '+ str(sample))
        
    h2part_df = pd.read_csv(path_local+desc+'.h2part.tsv',delimiter='\t')
    h2part_df.at[0,'female_file'] = phen+'_vdsA_s'+str(sample)+'.tsv.gz'
    h2part_df.at[0,'male_file'] = phen+'_vdsB_s'+str(sample)+'.tsv.gz'
    h2part_df.to_csv(path_local+desc+'.h2part.s'+str(sample)+'.tsv',sep='\t')
    subprocess.call(['gsutil','cp',path_local+desc+'.h2part.s'+str(sample)+'.tsv','gs://nbaya/rg_sex/'])


