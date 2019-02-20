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

vds = hc.read('gs://nbaya/ukb31063.hm3_variants.gwas_samples.vds')

kt_height = hc.read_table("gs://nbaya/split/height.kt")

cov = hc.read_table("gs://phenotype_31063/hail/0.1/ukb31063.gwas_covariates.both_sexes.kt")

print "Annotating vds with covariates and height phenotype..."
vds = vds.annotate_samples_table(cov, root = 'sa.covariates')
vds = vds.annotate_samples_table(kt_height, root = 'sa.height')
print "Finished annotating"

print "Exporting to vcf..."
vds.export_vcf('gs://nbaya/split/ukb31063.hm3_variants.gwas_samples_height.vcf.bgz')

print "\nFinished"
print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))


