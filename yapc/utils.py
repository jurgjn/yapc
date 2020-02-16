
import csv
import itertools
import os
import os.path

import pandas as pd

import pyBigWig

def makedirsp(fp):
    try:
        os.makedirs(fp)
    except:
        if not(os.path.isdir(fp)):
            raise

def is_bw(fp):
    if not os.path.isfile(fp):
        return False
    with pyBigWig.open(fp) as fh:
        f = fh.isBigWig()
    return f

def write_gffbed(fp,
                 chrom, start, end, name='', attr=None, score=0, strand='.',
                 thickStart=None, thickEnd=None, itemRgb='#0072b2',
                 trackline='#track gffTags=on', v=False):
    df = pd.DataFrame()
    def force_iter(l_):
        try:
            l = list(l_)
            if (len(l) > 1) and not(type(l_) is str):
                return l
            else:
                return list(itertools.repeat(l_, len(df)))
        except TypeError:
            return list(itertools.repeat(l_, len(df)))
    #return list(l) if hasattr(l, '__iter__') else list(itertools.repeat(l, len(df)))
    df['chrom'] = list(chrom)
    df['start'] = force_iter(start)
    df['end'] = force_iter(end)
    def pack_row(ir): return (";".join([("%s=%s" % (k, v)).replace(" ", "%20") for k, v in zip(ir[1].index, ir[1])]))
    attr_ = pd.concat([pd.DataFrame({'Name': force_iter(name)}), attr], axis=1)
    df['name'] = list(map(pack_row, attr_.iterrows()))
    df['score'] = force_iter(score)
    df['strand'] = force_iter(strand)

    if not(thickStart is None):
        df['thickStart'] = force_iter(thickStart)
    else:
        df['thickStart'] = df['start'].copy().tolist()

    if not(thickEnd is None):
        df['thickEnd'] = force_iter(thickEnd)
    else:
        df['thickEnd'] = df['end'].copy().tolist()

    df['itemRgb'] = force_iter(itemRgb)
    with open(fp, 'w') as fh:
        print(trackline, file=fh)
        df.sort_values(['chrom', 'start', 'end']).to_csv(fh, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)
    if v: return df
