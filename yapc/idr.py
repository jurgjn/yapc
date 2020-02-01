
import math

import numpy as np
import pandas as pd

import idr
import idr.idr

def call_idr(scores_rep1, scores_rep2):
    """Runs an IDR analysis by directly calling IDR model fitting routines.

    The IDR module does not seem to have a public, documented API. This wrapper
    is based on reading the IDR source code. Explicit references to the
    originating lines are given in the comments.
    """

    # https://github.com/nboley/idr/blob/74665e73bffb689a440948640c386b1188eea1e3/idr/idr.py#L295-L296
    rank1 = np.lexsort((np.random.random(len(scores_rep1)), scores_rep1)).argsort()
    rank2 = np.lexsort((np.random.random(len(scores_rep2)), scores_rep2)).argsort()

    # https://github.com/nboley/idr/blob/74665e73bffb689a440948640c386b1188eea1e3/idr/idr.py#L298-L299
    # https://github.com/nboley/idr/blob/74665e73bffb689a440948640c386b1188eea1e3/idr/idr.py#L408
    # https://github.com/nboley/idr/blob/74665e73bffb689a440948640c386b1188eea1e3/idr/idr.py#L858-L865
    localIDRs = idr.idr.fit_model_and_calc_local_idr(r1=np.array(rank1, dtype=np.int), r2=np.array(rank2, dtype=np.int))

    # https://github.com/nboley/idr/blob/74665e73bffb689a440948640c386b1188eea1e3/idr/idr.py#L883
    IDRs = idr.idr.calc_global_IDR(localIDRs)

    # https://github.com/nboley/idr/blob/74665e73bffb689a440948640c386b1188eea1e3/idr/idr.py#L334
    def scaledIDR(IDR): return -math.log10(max(1e-5, IDR))
    scaledIDRs = pd.Series([* map(scaledIDR, IDRs) ])
    return scaledIDRs

def call_idr_system(chroms, starts, ends, scores_rep1, scores_rep2, prefix_out, condition):
    """Runs an IDR analysis using the command line script.

    This writes "mock" input files (with identical peak coordinates), and calls
    the IDR script as an external process. Not actively used, but kept around
    for sanity checking.
    """

    fp_rep1 = prefix_out + '_peaksidr_%(condition)s_rep1.bed' % locals()
    fp_rep2 = prefix_out + '_peaksidr_%(condition)s_rep2.bed' % locals()
    fp_iout = prefix_out + '_peaksidr_%(condition)s.bed' % locals()

    df_ = pd.DataFrame()
    df_['chrom'] = chroms
    df_['start'] = starts
    df_['end'] = ends
    df_['name'] = '.'

    # https://github.com/nboley/idr/blob/74665e73bffb689a440948640c386b1188eea1e3/idr/idr.py#L287-L299
    # https://github.com/nboley/idr/blob/74665e73bffb689a440948640c386b1188eea1e3/idr/idr.py#L845
    # https://github.com/nboley/idr/blob/74665e73bffb689a440948640c386b1188eea1e3/idr/idr.py#L858-L865
    df_['score_rep1'] = pd.Series(scores_rep1).rank()
    df_['score_rep2'] = pd.Series(scores_rep2).rank()
    df_['strand'] = '.'

    # https://github.com/nboley/idr/blob/74665e73bffb689a440948640c386b1188eea1e3/idr/idr.py#L64
    # values don't seem to be used, but seem to be required to be present in the input...
    df_['col6'] = 0 # signalValue?
    df_['col7'] = 0 # pValue?
    df_['col8'] = 0 # qValue?

    df_[['chrom', 'start', 'end', 'name', 'score_rep1', 'strand', 'col6', 'col7', 'col8']].to_csv(fp_rep1, header=False, index=False, sep='\t')
    df_[['chrom', 'start', 'end', 'name', 'score_rep2', 'strand', 'col6', 'col7', 'col8']].to_csv(fp_rep2, header=False, index=False, sep='\t')

    idr_args = '--input-file-type bed --rank score --peak-merge-method max --plot'
    idr_cmd = 'idr %(idr_args)s --samples %(fp_rep1)s %(fp_rep2)s --output-file %(fp_iout)s' % locals()
    print('Running idr on condition %(condition)s:' % locals())
    print(idr_cmd)
    os.system(idr_cmd)

    col_prefix = condition
    names_ = [
        'chrom', 'concave_start', 'concave_end', 'name',
        '%(col_prefix)s_scaledIDR' % locals(), 'strand',
        '%(col_prefix)s_localIDR' % locals(), '%(col_prefix)s_globalIDR' % locals(),
        '%(col_prefix)s_rep1_start' % locals(), '%(col_prefix)s_rep1_end' % locals(), '%(col_prefix)s_rep1_score' % locals(),
        '%(col_prefix)s_rep2_start' % locals(), '%(col_prefix)s_rep2_end' % locals(), '%(col_prefix)s_rep2_score' % locals(),
    ]
    df_ = pd.read_csv(fp_iout, names=names_, sep='\t').sort_values(['chrom', 'concave_start', 'concave_end']).reset_index(drop=True)
    return df_['%(col_prefix)s_globalIDR' % locals()]
