#!/usr/bin/env python2

###
# load packages
print "Loading packages..."
###
import datetime
import numpy as np
import random
from hail import *
hc = HailContext()

phen = '20160'
desc = 'smoking' #one word description of the phenotype


print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))

vds = hc.read('gs://nbaya/split/ukb31063.hm3_variants.gwas_samples.vds')

phen_table = hc.read_table("gs://nbaya/split/phen_tables/"+phen+'_'+desc+".kt")

cov = hc.read_table("gs://phenotype_31063/hail/0.1/ukb31063.gwas_covariates.both_sexes.kt")

print "Annotating vds with covariates and phenotype..."
vds = vds.annotate_samples_table(cov, root = 'sa.covariates')
vds = vds.annotate_samples_table(phen_table, root = 'sa.'+desc)
print "Finished annotating"

cov_list = ['sa.covariates.isFemale',
           'sa.covariates.age',
           'sa.covariates.age_squared',
           'sa.covariates.age_isFemale',
           'sa.covariates.age_squared_isFemale'] + \
['sa.covariates.PC{:}'.format(i) for i in xrange(1,21)]

all_ID = vds.sample_ids    # list of all sample IDs
all_i = np.arange(len(all_ID))  # array of all indices

i_start = 5     # iteration index to start at (change to prevent overwriting previous samples)
                # check gs://nbaya/split/ for previous index numbers
i_n = 1         # number of iterations

for sample in np.arange(i_start,i_start+i_n):
    print "\nStarting sample "+str(sample)
    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))

    random.shuffle(all_i)  # randomly shuffle all indices
    
    vdsA_i = all_i[:len(all_i)/2]   # take first half of randomly shuffled indices
    vdsB_i = all_i[len(all_i)/2:]   # take second half of randomly shuffled indices
    
    vdsA_ID = [all_ID[i] for i in vdsA_i]  # list of sample IDs for first half
    vdsB_ID = [all_ID[i] for i in vdsB_i]  # list of sample IDs for first half
    
    vdsA = vds.filter_samples_list(vdsA_ID)
    vdsB = vds.filter_samples_list(vdsB_ID)
    
    vdsA_result = vdsA.linreg3(ys = ['sa.'+desc], covariates=cov_list, root = 'va.linreg', use_dosages = True)
    vdsB_result = vdsB.linreg3(ys = ['sa.'+desc], covariates=cov_list, root = 'va.linreg', use_dosages = True)
    
    vdsA_result.variants_table().export('gs://nbaya/split/'+phen+'_vdsA_s'+str(sample)+'.tsv.bgz')
    vdsB_result.variants_table().export('gs://nbaya/split/'+phen+'_vdsB_s'+str(sample)+'.tsv.bgz')
    
print "\nFinished"
print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
