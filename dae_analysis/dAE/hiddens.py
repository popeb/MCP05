#!/usr/bin/env python
"""uniform sized layers for the given dataset"""


import sys, os, re
import argparse
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

def get_nN(arch_fn, n):
    P = re.compile(".*W(\d+)_B(\d+).*.h5")
    #ch, natsg, w_, b_ = os.path.splitext(arch_fn)[0].split("_")
    w, b = map(int, P.match(arch_fn).groups()) #int(w_[1:]), int(b_[1:])
    N = n * w / b
    return n, N

parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    epilog="Taylor Lab (odenas@emory.edu)")
parser.add_argument("input", type=str, help="HDF5 archive of dataset")
parser.add_argument("nlev", type=int, help="number of layers")
parser.add_argument("--nin", type=int, default=31, help="number of input tracks")
parser.add_argument("--nout", type=int, default=31, help="number of output tracks")
opt = parser.parse_args()

log.warning("assuming %d in-tracks and %d out-tracks", opt.nin, opt.nout)
s = sum(get_nN(opt.input, opt.nin))
print " ".join(map(str, reversed(range(opt.nout, s, s/(1 + opt.nlev))))[1:])
