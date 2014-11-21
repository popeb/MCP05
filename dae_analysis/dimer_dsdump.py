#!/usr/bin/env python

"""dump a dataset into a set of .npy files"""

import os
import sys
import logging
import argparse
from operator import itemgetter

import numpy as np

from dimer import archive, data, ops, argutils

logging.basicConfig(level=logging.INFO)
lg = logging.getLogger()

def base_name(dspath):
    """remove the extension from the archivename
    and combine with dataset name

    e.g., 'ds.h5:ciao' results in 'ds_ciao'

    :param str dspath: path to dataset (see dimer.archive.DSPEC_MSG
    :rtype: string
    """

    arch_base = os.path.splitext(archive.archname(dspath))[0]
    ds_base = archive.split(dspath)[1]

    return "_".join((arch_base, ds_base))

def out_name(dspath, dscomp):
    """combine the base name (specific to the archive name) with
    with the dataset component

    e.g., the 'X' component of 'ds.h5:ciao' results in 'ds_ciao_X'
    e.g., 'gene_names' of 'ds.h5:ciao' results in 'ds_ciao_gene_names'
    """

    return "_".join((base_name(dspath), dscomp))

def dump_index(ofd, idx):
    for i in idx:
        ofd.write("%s\n" % str(i))

if __name__ != "__main__":
    lg.error("this is a script. Do not import")
    sys.exit(1)

parser = argutils.parser_inst(__doc__)
parser.add_argument("dataset", type=argutils.ArgTypes.is_dsarch,
                    help="Dataset " + archive.DSPEC_MSG)
parser.add_argument("raw", action='store_true', default=False,
                    help="Select the raw dataset, if true")

opt = parser.parse_args()
ds = data.AnchorDataset._from_archive(opt.dataset, opt.raw)

## dump ds/X ##
lg.info("dumping X to %s_XXX ...", base_name(opt.dataset))

# data
np.save(out_name(opt.dataset, "X.npy"), ds.X)

# attributes
with open(out_name(opt.dataset, "items.txt"), 'w') as ofd:
    dump_index(ofd, ds.pX.items)

with open(out_name(opt.dataset, "major_axis.txt"), 'w') as ofd:
    dump_index(ofd, ds.pX.major_axis)

with open(out_name(opt.dataset, "minor_axis.txt"), 'w') as ofd:
    dump_index(ofd, ds.pX.minor_axis)

if ds.Y is None:
    lg.info("no Y and no T component. DONE")
    sys.exit(0)

## dump ds/Y ##
lg.info("dumping Y to %s_XXX ...", base_name(opt.dataset))

# data
np.save(out_name(opt.dataset, "Y.npy"), ds.Y)

if ds.T is None:
    lg.info("no T component. DONE")
    sys.exit(0)

## dump ds/T ##
lg.info("dumping T to %s_XXX ...", base_name(opt.dataset))

# data
np.save(out_name(opt.dataset, "T_class.npy"), ds.T)

with open(out_name(opt.dataset, "T_names.txt"), 'w') as ofd:
    dump_index(ofd, ds.dfT['label_name'])
