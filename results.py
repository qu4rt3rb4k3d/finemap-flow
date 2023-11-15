import numpy
import pandas

def get_results_df(snp_dfs, posts, ss_dfs, af_field):
    fm_dfs = [snp_df[['rsid', 'mean', 'prob']].copy() for snp_df in snp_dfs]
    results_dfs = []
    for fm_df, post, ss_df in zip(fm_dfs, posts, ss_dfs):
        fm_df['mean'] = numpy.array(fm_df['mean']) * (1 - post[0])
        fm_df['prob'] = numpy.array(fm_df['prob']) * (1 - post[0])
        ss_df = ss_df.sort_values(by='variant')
        fm_df = fm_df.sort_values(by='rsid')
        fm_df['af'] = ss_df[af_field].to_list()
        results_dfs.append(fm_df)
    return pandas.concat(results_dfs)
