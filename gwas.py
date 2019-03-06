def gwas(mt, x, y, cov_list=[], with_intercept=True, pass_through=[], path_to_save=None, 
         normalize_x=True, is_std_cov_list=False):
    '''Runs GWAS'''
    
    mt = mt._annotate_all(col_exprs={'y':y},
                           entry_exprs={'x':x})
    if normalize_x:
        mt = normalize_genotypes(mt, mt.x) 
        mt = mt.annotate_entries(x = mt.__norm_gt).drop('__norm_gt')
    
    if is_std_cov_list:
        cov_list = ['isFemale','age','age_squared','age_isFemale',
                    'age_squared_isFemale']+['PC{:}'.format(i) for i in range(1, 21)]
        
    if str in list(map(lambda x: type(x),cov_list)):
        cov_list = list(map(lambda x: mt[x] if type(x) is str else x,cov_list))
        
    cov_list = ([1] if with_intercept else [])+cov_list
    
    print(pass_through)

    gwas_ht = hl.linear_regression_rows(y=mt.y,
                                        x=mt.x,
                                        covariates=cov_list,
                                        pass_through = ['rsid']+pass_through)
    
    gwas_ht = gwas_ht.annotate_globals(with_intercept = with_intercept)
        
    gwas_ht = gwas_ht.rename({'rsid':'SNP'}).key_by('SNP')
        
    gwas_ht = gwas_ht.select(Z = gwas_ht.t_stat,
                             N = gwas_ht.n)
    
    ss_template = hl.read_table('gs://nbaya/rg_sex/hm3.sumstats_template.ht') # sumstats template as a hail table
    ss_template = ss_template.key_by('SNP')
    ss_template = ss_template.annotate(N = n_samples)
        
    ss = ss_template.annotate(Z = gwas_ht[ss_template.SNP].Z,
                              N = gwas_ht[ss_template.SNP].N)
    
    if path_to_save is not None:
        ss.export(path_to_save)
        
    return ss
