def gwas(mt, x, y, cov_list=[], with_intercept=True, pass_through=[], path_to_save=None, 
         normalize_x=False, is_std_cov_list=False):
    '''Runs GWAS'''
    
    mt = mt._annotate_all(col_exprs={'__y':y},
                           entry_exprs={'__x':x})
    if normalize_x:
        mt = mt.annotate_rows(__gt_stats = hl.agg.stats(mt.__x))
        mt = mt.annotate_entries(__x= (mt.__x-mt.__gt_stats.mean)/mt.__gt_stats.stdev) 
        mt = mt.drop('__gt_stats')
    
    if is_std_cov_list:
        cov_list = ['isFemale','age','age_squared','age_isFemale',
                    'age_squared_isFemale']+['PC{:}'.format(i) for i in range(1, 21)]
        
    if str in list(map(lambda x: type(x),cov_list)):
        cov_list = list(map(lambda x: mt[x] if type(x) is str else x,cov_list))
        
    cov_list = ([1] if with_intercept else [])+cov_list
    
    print(pass_through)

    gwas_ht = hl.linear_regression_rows(y=mt.__y,
                                        x=mt.__x,
                                        covariates=cov_list,
                                        pass_through = ['rsid']+pass_through)
    
    gwas_ht = gwas_ht.annotate_globals(with_intercept = with_intercept)
        
    gwas_ht = gwas_ht.rename({'rsid':'SNP'}).key_by('SNP')
        
    gwas_ht = gwas_ht.select(Z = gwas_ht.t_stat,
                             N = gwas_ht.n)
    
    ss_template = hl.read_table('gs://nbaya/rg_sex/hm3.sumstats_template.ht') # sumstats template as a hail table
    ss_template = ss_template.key_by('SNP')
        
    ss = ss_template.annotate(Z = gwas_ht[ss_template.SNP].Z,
                              N = gwas_ht[ss_template.SNP].N)
    
    if path_to_save is not None:
        ss.export(path_to_save)
        
    return ss
