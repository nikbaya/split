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

print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))

vds = hc.read('gs://nbaya/split/ukb31063.hm3_variants.gwas_samples.vds')

phen = '20160'
desc = 'smoking'


phen_str = '"'+phen+'"'
phen_tb = hc.import_table('gs://phenotype_31063/ukb31063.phesant_phenotypes.both_sexes.tsv.bgz', types = {phen_str: TString()}, missing = "")
phen_tb = phen_tb.select(['"userId"',phen_str])
phen_tb = phen_tb.rename({ '"userId"': 's', phen_str: desc})
phen_tb = phen_tb.key_by('s')

phen_tb = phen_tb.annotate("smoking = smoking.replace('\"','').toFloat()")        #check phenotype
phen_tb.show()

cov = hc.read_table("gs://phenotype_31063/hail/0.1/ukb31063.gwas_covariates.both_sexes.kt")

vds = vds.annotate_samples_table(cov, root = 'sa.covariates')
vds = vds.annotate_samples_table(phen_tb, root = 'sa.'+desc)

print "Exporting to tsv..."
vds.samples_table().flatten().export('gs://nbaya/split/phen_tables/ukb31063.hm3_variants.gwas_samples_'+phen+'_'+desc+'.tsv')

print "\nFinished"
print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))



