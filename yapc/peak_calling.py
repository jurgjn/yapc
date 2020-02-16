import argparse
import collections
import csv
import itertools
import logging
import os
import math
import sys

from pprint import pprint

import numpy as np
import pandas as pd
import pyBigWig

import yapc
import yapc.idr
import yapc.utils

def prepare_second_derivative_kernel(width, times):
    """
    Return composite kernel that calculates the smoothed second derivative via a single convolution with the input signal.

    Weights of the second derivative kernel: https://doi.org/10.1090/S0025-5718-1988-0935077-0
    Rationale for smoothing 3 times (by default): https://terpconnect.umd.edu/~toh/spectrum/Differentiation.html#Smoothing
    """
    kernel = [-1.0/560, 8.0/315, -1.0/5, 8.0/5, -205.0/72, 8.0/5, -1.0/5, 8.0/315, -1.0/560]
    rolling_mean_kernel = np.ones(width)/float(width)
    for i in range(times):
        kernel = np.convolve(kernel, rolling_mean_kernel)
    ADHOC_LARGE_NUMBER = 1e6 # scale to have better visualisation of d2smooth in IGV
    return ADHOC_LARGE_NUMBER*kernel[::-1]

def coverage(l_fp_inp, fp_out):
    print('calculating mean of all input signal:')
    l_fh_inp = [pyBigWig.open(fp_inp) for fp_inp in l_fp_inp]
    assert pyBigWig.numpy == 1
    assert all([fh_inp.isBigWig() for fh_inp in l_fh_inp])

    chroms_sizes_ = dict() # Iterate over the *union* of all chromosomes, as empty chromosomes might be missing in (some) files (e.g. Ensembl CB4 short sequences)
    for fh_inp in l_fh_inp:
        for chrom, size in fh_inp.chroms().items():
            chroms_sizes_[chrom] = max(size, chroms_sizes_.get(chrom, 0))
    chroms_sizes = collections.OrderedDict(sorted(chroms_sizes_.items(), key=lambda chrom_size: bytearray(chrom_size[0], encoding='ascii')))

    fh_out = pyBigWig.open(fp_out, 'w')
    fh_out.addHeader([(chrom, size) for (chrom, size) in chroms_sizes.items()])

    for chrom, size in itertools.islice(chroms_sizes.items(), None):
        l_val = []
        print(chrom, size)
        for fh_inp in l_fh_inp:
            #val = fill_gap(fh_inp.values(chrom, 0, size, numpy=True)) # This is doing something else.
            try:
                val = fh_inp.values(chrom, 0, size, numpy=True)
            except RuntimeError:
                print('RuntimeError suggests empty chromosome; returning flat signal')
                val = np.zeros(size)
            #print(val.shape)
            l_val.append(val)
        fh_out.addEntries(chrom, 0, values=np.mean(l_val, axis=0), span=1, step=1)

    # Close files
    for fh_inp in l_fh_inp:
        fh_inp.close()
    fh_out.close()

def d2smooth(fp_inp, fp_out, kernel):
    fh_inp = pyBigWig.open(fp_inp)
    fh_out = pyBigWig.open(fp_out, 'w')
    fh_out.addHeader([(chrom, size) for (chrom, size) in fh_inp.chroms().items()])

    for chrom, size in itertools.islice(fh_inp.chroms().items(), None):
        print(chrom, size)
        val = fh_inp.values(chrom, 0, size, numpy=True)
        val_d2 = np.convolve(val, kernel, mode='same')
        fh_out.addEntries(chrom, 0, values=val_d2, span=1, step=1)

    fh_inp.close()
    fh_out.close()

# https://doi.org/10.1101/gr.153668.112
# For analyses where a single TSS position was required, we considered the distribution of cap 5' ends
# within the TIC, and selected the position with the most tags (the mode). In the case of a tie (two or more
# positions with the same number of tags), we selected the mode closest to the median of the TIC.
def nanargmax_median(a):
    a_max = np.nanmax(a)
    a_max_indices = np.flatnonzero(a == a_max)
    if len(a_max_indices) == 0:
        return 0
    return a_max_indices[len(a_max_indices) // 2]

assert nanargmax_median([1,2,3,3,3,2,1]) == 3
assert nanargmax_median([1,2,2,1]) == 2
assert nanargmax_median([1]) == 0

def find_concave_regions_chrom(d2y, chrom='.', tol=1e-10):
    s = np.where(np.diff((d2y < tol).astype(int))==1)[0] + 1
    e = np.where(np.diff((d2y < tol).astype(int))==-1)[0] + 1
    if len(s) == 0:
        print('No raw concave regions on chrom %s' % (chrom,))
        return pd.DataFrame(columns=['chrom', 'concave_start', 'concave_end', 'mode'])

    if s[0] > e[0]:
        s = np.insert(s, 0, 0)
    if s[-1] > e[-1]:
        e = np.insert(e, len(e), len(d2y))
    v = [-np.mean(d2y[i:j]) for i,j in zip(s,e)]
    m = [i + nanargmax_median(-d2y[i:j]) for i,j in zip(s,e)]

    df_regions = pd.DataFrame(collections.OrderedDict([
        ('chrom', [chrom]*len(s)),
        ('concave_start', s),
        ('concave_end', e),
        #('val', v),
        ('mode', m),
    ]))

    print('%d raw concave regions on chrom %s' % (len(df_regions), chrom))
    return df_regions

def find_concave_regions(fp_inp, fp_out_tsv, fp_out_bed, chroms=None):
    fh_inp = pyBigWig.open(fp_inp)
    if chroms is None:
        chroms_sizes = fh_inp.chroms().items()
    else:
        chroms_sizes = [(chrom, size) for chrom, size in fh_inp.chroms().items() if chrom in chroms]

    df_regions = pd.concat([
        find_concave_regions_chrom(d2y=fh_inp.values(chrom, 0, size, numpy=True), chrom=chrom)
        for chrom, size in itertools.islice(chroms_sizes, None)], axis=0, ignore_index=True)
    fh_inp.close()
    print('%d peaks total' % (len(df_regions),))
    return df_regions

def scores(fp_inp, chroms, starts, ends, kernel):
    """
    Calculate scores for a list of regions, and return an iterator of the scores.

    chroms, starts, ends -- lists of genomic coordinates specifying the regions
    fp_inp -- name of the BigWig file containing the raw coverage signal
    kernel -- specifies the kernel to convolve the raw signal with (this should match the kernel used to define the concave regions)

    Scores are currently defined as the reverse of the mean smoothed second derivative within the region.
    """
    fh_inp = pyBigWig.open(fp_inp)
    d_chrom_d2smooth = {}
    for chrom, size in itertools.islice(fh_inp.chroms().items(), None):
        if not(chrom in chroms):
            # DEBUG: skip convolution of a particular chromosome if it won't be used in scoring
            continue
        val = fh_inp.values(chrom, 0, size, numpy=True)
        val_d2 = np.convolve(val, kernel, mode='same')
        d_chrom_d2smooth[chrom] = val_d2

    for chrom, start, end in zip(chroms, starts, ends):
        try:
            yield -np.mean(d_chrom_d2smooth[chrom][start:end])
        except KeyError: # empty chromosome
            yield -np.inf

def check_recycle(label, fp, recycle):
    # File does not exist
    if not os.path.isfile(fp):
        logging.info('[compute]\t%s: %s' % (label, fp))
        return True

    # File exists and recycle enabled
    elif recycle:
        logging.info('[recycle]\t%s: %s' % (label, fp))
        return False

    # File exists and recycle not enabled
    else:
        os.remove(fp) # delete existing file
        logging.info('[recompute]\t%s: %s' % (label, fp))
        return True

def main(df_samples, prefix_out, kernel, min_concave_region_width, truncate_idr_input, fixed_peak_halfwidth, pseudoreplicates, recycle, chroms_subset):

    # fp_coverage
    fp_coverage = prefix_out + '_coverage.bw'
    if check_recycle('pooled coverage track', fp_coverage, recycle):
        coverage(df_samples.query('~pseudoreplicate')['fp'].tolist(), fp_coverage)

    # fp_d2smooth
    fp_d2smooth = prefix_out + '_d2smooth.bw'
    if check_recycle('d2smooth track', fp_d2smooth, recycle):
        d2smooth(fp_coverage, fp_d2smooth, kernel)

    # fp_peaksall -- all raw peaks, scored by D2
    fp_peaksall_tsv = prefix_out + '_peaksall.tsv'
    fp_peaksall_bed = prefix_out + '_peaksall.bed'

    # Use df_samples index values to identify throughout
    df_samples.set_index((f"{condition}_{sample}" for condition, sample in zip(df_samples['condition'], df_samples['sample'])), inplace=True)

    if check_recycle('raw concave regions', fp_peaksall_tsv, recycle):
        # call concave regions
        df_regions = find_concave_regions(fp_d2smooth, fp_peaksall_tsv, fp_peaksall_bed, chroms=chroms_subset)

        # score concave regions
        for i, r in itertools.islice(df_samples.iterrows(), None):
            fp_ = r['fp']
            col_ = f'{i}_score'
            print('scoring peaks for: %s' % (i,))
            df_regions[col_] = [*scores(fp_, df_regions.chrom.tolist(), df_regions.concave_start.tolist(), df_regions.concave_end.tolist(), kernel)]

        # write all concave regions to output
        df_regions.to_csv(fp_peaksall_tsv, header=True, index=False, sep='\t')
        yapc.utils.write_gffbed(fp_peaksall_bed,
            chrom = df_regions['chrom'],
            start = df_regions['concave_start'],
            end = df_regions['concave_end'],
            attr = df_regions[[f'{i}_score' for i in df_samples.index]],
            thickStart = df_regions['mode'],
            thickEnd = df_regions['mode'] + 1,
        )

    df_regions = pd.read_csv(fp_peaksall_tsv, sep='\t')
    logging.info('%d raw concave regions found' % (len(df_regions),))

    # discard narrow peaks
    logging.info('%d raw peaks' % (len(df_regions),))
    df_regions = df_regions.query('(concave_end - concave_start) >= @min_concave_region_width').reset_index(drop=True)
    logging.info('%d peaks after discarding narrow concave regions (min_concave_region_width=%d)' % (len(df_regions), min_concave_region_width))

    # rank scores & truncate peak list based on "best" rank => 100,000 peaks
    score_cols = [f'{i}_score' for i in df_samples.index]
    df_regions['min_rank'] = df_regions[score_cols].rank(axis=0, ascending=False).min(axis=1)
    df_regions_idr = df_regions.sort_values('min_rank').head(truncate_idr_input).sort_values(['chrom', 'concave_start', 'concave_end']).reset_index(drop=True)
    logging.info('%d peaks after filtering based on best rank across all conditions (truncate_idr_input=%d)' % (len(df_regions_idr), truncate_idr_input))

    str_ = df_regions[score_cols].head().to_string(line_width=1000, index=False)
    logging.info('curvature index-scores:\n%(str_)s' % locals())

    # Final peak boundaries -- either defined by the concave region or fixed at flank_len bases either side of the mode
    if fixed_peak_halfwidth is None:
        df_regions_idr['start'] = df_regions_idr['concave_start']
        df_regions_idr['end'] = df_regions_idr['concave_end']
    else:
        df_regions_idr['start'] = df_regions_idr['mode'] - fixed_peak_halfwidth
        df_regions_idr['end'] = df_regions_idr['mode'] + fixed_peak_halfwidth + 1

    # Final output of all raw peaks subject to IDR analyses
    fp_regions_idr = prefix_out + '.tsv'
    if check_recycle('IDR analyses', fp_regions_idr, recycle):
        # run IDR using scores from each condition
        for condition in df_samples['condition'].unique():
            col_prefix = condition
            (rep1, rep2) = df_samples.query("condition == @condition").index.values
            col_rep1_score = f"{rep1}_score"
            col_rep2_score = f"{rep2}_score"
            print(col_rep1_score)
            if not pseudoreplicates:
                #IDRs_ = yapc.idr.call_idr_system(df_regions_idr['chrom'], df_regions_idr['concave_start'], df_regions_idr['concave_end'], df_regions_idr['%(condition)s_rep1_score' % locals()], df_regions_idr['%(condition)s_rep2_score' % locals()], prefix_out, condition)
                IDRs_ = yapc.idr.call_idr(df_regions_idr[col_rep1_score], df_regions_idr[col_rep2_score])
                df_regions_idr['%(col_prefix)s_globalIDR' % locals()] = IDRs_
            else:
                # TODO pseudoreplicate code is not updated to sample sheet refactoring
                th_select = 0.001
                th_globalIDR = -math.log(th_select, 10)
                IDRs_rep = run_idr_call(df_regions_idr['%(condition)s_rep1_score' % locals()], df_regions_idr['%(condition)s_rep2_score' % locals()])
                IDRs_prp = run_idr_call(df_regions_idr['%(condition)s_prp1_score' % locals()], df_regions_idr['%(condition)s_prp2_score' % locals()])

                n_rep = sum(IDRs_rep >= th_globalIDR)
                n_prp = sum(IDRs_prp >= th_globalIDR)

                if n_rep >= n_prp:
                    logging.info('condition %s: n_rep=%d, n_prp=%d (at IDR=%.3f) => using replicates' % (condition, n_rep, n_prp, th_select))
                    df_regions_idr['%(col_prefix)s_globalIDR' % locals()] = IDRs_rep
                else:
                    logging.info('condition %s: n_rep=%d, n_prp=%d (at IDR=%.3f) => using pseudoreplicates' % (condition, n_rep, n_prp, th_select))
                    df_regions_idr['%(col_prefix)s_globalIDR' % locals()] = IDRs_prp

        df_regions_idr.to_csv(fp_regions_idr, header=True, index=False, sep='\t')
        logging.info('Wrote all %d IDR-scored peaks to %s' % (len(df_regions_idr), fp_regions_idr))

    df_regions_idr = pd.read_csv(fp_regions_idr, sep='\t')
    logging.info('%d IDR-scored peaks' % (len(df_regions),))

    logging.info('globalIDR:')
    col_ = ['%(col_prefix)s_globalIDR' % locals() for col_prefix in df_samples['condition'].unique()]
    logging.info(df_regions_idr[col_].head().to_string(line_width=1000, index=False))

    # Write final .bed-files at a range of thresholds
    col_ = ['%(col_prefix)s_globalIDR' % locals() for col_prefix in df_samples['condition'].unique()]
    df_regions_idr['max_globalIDR'] = df_regions_idr[col_].max(axis=1)
    for th in [0.001, 0.005, 0.01, 0.05, 0.1, 0.2]:
        th_globalIDR = -math.log(th, 10)
        fp_regions_idr_th = prefix_out + '_%(th)s.bed' % locals()
        df_regions_idr_th = df_regions_idr.query('max_globalIDR >= @th_globalIDR').reset_index(drop=True)
        yapc.utils.write_gffbed(fp_regions_idr_th,
            chrom = df_regions_idr_th['chrom'],
            start = df_regions_idr_th['start'],
            end = df_regions_idr_th['end'],
            attr = df_regions_idr_th[col_],
        )
        logging.info('Wrote %d peaks at IDR=%.3f to %s' % (len(df_regions_idr_th), th, fp_regions_idr_th))
