import pandas as pd
import hail as hl

variant_set = 'hm3'

mt0 = hl.read_matrix_table('gs://nbaya/split/ukb31063.'+variant_set+'_variants.gwas_samples_repart.mt')

df = mt0.rows().to_pandas()

df['chr'] = df['locus.contig'].astype(int)

for chr in range(1,23):
    print(f'chr: {chr}')
    ld = pd.read_csv(f'/home/nbaya/ldscores/{chr}.l2.ldscore.gz',compression='gzip',sep='\t')
    # ld_test = pd.read_csv(f'/home/nbaya/ldsc/test_{chr}.l2.ldscore.gz',compression='gzip',sep='\t')
    tmp = df.loc[df.chr==chr]

    tmp.loc[:,'SNP'] = tmp.apply(lambda row: f'{row["chr"]}:{row["locus.position"]}:{row["alleles"][0]}:{row["alleles"][1]}',axis=1)    
    tmp1 = tmp.drop(['locus.contig','locus.position','alleles','varid','chr'],axis=1)

    tmp1 = ld.merge(tmp1,left_on='SNP',right_on='SNP')
    tmp1['SNP'] = tmp1['rsid']
    tmp1 = tmp1.drop('rsid',axis=1)
    tmp1.to_csv(f'/home/nbaya/ldsc/{chr}.l2.ldscore.gz',compression='gzip',sep='\t',index=False)
