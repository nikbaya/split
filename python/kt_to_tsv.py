#!/usr/bin/env python2

from hail import *
hc = HailContext()

kt = hc.read_table('gs://ukb31063-mega-gwas/ldsc/ld_ref_panel/hm3.r3.b37.auto_bi_af.ukbb_gwas_qcpos.no_mhc.kt')

print kt.count()

kt = kt.select('v')

print kt.count()

kt.export('gs://nbaya/split/hapmap3_variants.tsv')

