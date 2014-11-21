#!/usr/bin/env python

"infer segments from model and data"

import sys
import logging
import os
from itertools import izip

import numpy as np
from dimer.archive import basename
from dimer.data import AnnotatedDataset
from dimer.argutils import ArgTypes, parser_inst

from sklearn.decomposition import PCA
from sklearn.preprocessing import scale

rng = np.random.RandomState()

logging.basicConfig(level=logging.INFO)
lg = logging.getLogger()

if __name__ != '__main__':
    lg.error("this is a script do not import")
    sys.exit(1)


parser = parser_inst(__doc__)
parser.add_argument("input", type=ArgTypes.is_dsarch, help="Input dataset.")
parser.add_argument("components", nargs="+", type=ArgTypes.is_uint,
                    help="Size of the output")
opt = parser.parse_args()

ds = AnnotatedDataset._from_archive(opt.input, False)
ds.flatten()
## TODO: subtract training_mean, divide by training_sd

parse_reduced_row = lambda r: " ".join(map(lambda v: "%.1f" % v, r))
for pca_comp in opt.components:
    model = PCA(n_components=pca_comp)
    out_dt = model.fit_transform(scale(ds.X))
    out_fname = "%s_%d_%d.txt" % (os.path.splitext(basename(opt.input))[0],
                                  ds.X.shape[1], pca_comp)
    lg.info("%d -> %d ... saving to: %s ", ds.X.shape[1], pca_comp,
            out_fname)
    with open(out_fname, 'w') as ofd:
        for i in xrange(ds.X.shape[0]):
            print >>ofd, "%s\t%s" % (ds.label_names[int(ds.Y[i])],
                                     parse_reduced_row(out_dt[i]))
