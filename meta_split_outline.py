#!/usr/bin/env python3

ht = mt.cols()
rows keyed by sample with y: Float, cov: Array[Float], group_id: Int

mt:
columns keyed by sample with y: Float, cov: Array[Float], group_id: Int
rows keyed by locus, alleles
entry field x: Float

meta:
gmt = (mt.group_cols_by(mt.group_id)
        .aggregate(linreg = agg.linreg(y=mt.y, x=hl.array([mt.x, 1]).extend(mt.cov))))
gmt.select_entries(beta=mt.linreg.beta[0], standard_error=mt.linreg.standard_error[0]).write(PATH)

To do meta analysis:

you need a map from group_id to A or B: create an array pi of length 350 of 0s and 1s.

loop over pi:

pi = hl.literal(pi)
gmt = hl.read_matrix_table(PATH)
ht = gmt.group_cols_by(pi[gmt.group_id]).aggregate(unnorm_meta_beta=agg.sum(gmt.beta / gmt.se ** 2), inv_se2 = agg.sum(1 / gmt.se ** 2)).make_table()

consider replacing variant with row index at the beginning.
you could export to TSV

a = np.array(hl.array([ht.A.unnorm_meta_beta / ht.A.inv_se2, hl.sqrt(1 / ht.A.inv_se2), ht.B.unnorm_meta_beta / ht.B.inv_se2, hl.sqrt(1 / ht.B.inv_se2)]).collect())
a.save(...)