#!/usr/bin/env python

"""Produce a /<dataname>T dataframe (shape->[<nr of anchors>,<2>])
<dataname>T["label_code"] codes of labels
<dataname>T["label_name"] names of labels

binary lables (thresholded at 0) supported so far
"""


import logging
#import pdb
import argparse

import pandas as pd

from dimer import archive

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Taylor Lab (odenas@emory.edu)")
    parser.add_argument("input", type=archive.dset_path, help=archive.DSPEC_MSG)
    parser.add_argument("--labels", nargs='+', type=str, help="label names in order")

    opt = parser.parse_args()

    with pd.get_store(archive.archname(opt.input)) as store:
        key = archive.basename(opt.input)
        y = store["%s/Y" % key]

        codel = y.map(int)
        coden = codel.map(lambda v: opt.labels[v])

        t = pd.DataFrame({"label_code": codel, 'label_name': coden})
        store["%s/T" % key] = t
