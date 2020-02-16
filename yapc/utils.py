
import os
import os.path

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
