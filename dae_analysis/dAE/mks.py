#!/usr/bin/env python

"infer segments from model and data"

import sys
import logging
import os
from itertools import izip

import numpy as np
from bx.bbi.bigwig_file import BigWigFile
from theano.tensor.shared_randomstreams import RandomStreams
from dimer import archive, human_readable
from dimer.data import AnnotatedDataset
from dimer.nnet import autoencoder
from dimer.argutils import ArgTypes, parser_inst
from dimer.nnet import autoencoder
from theano.tensor.shared_randomstreams import RandomStreams


rng = np.random.RandomState()
thrng = RandomStreams(rng.randint(100000))

logging.basicConfig(level=logging.INFO)
lg = logging.getLogger()
logging.getLogger("dimer.nnet").setLevel(logging.INFO)
#logging.getLogger("dimer.nnet.autoencoder").setLevel(logging.DEBUG)
#logging.getLogger("dimer.nnet.config_spec").setLevel(logging.DEBUG)
#logging.getLogger("dimer.archive").setLevel(logging.INFO)


def load_model(inp, trid, input_data_type):
    from dimer.nnet import autoencoder
    from theano.tensor.shared_randomstreams import RandomStreams

    rng = np.random.RandomState()
    thrng = RandomStreams(rng.randint(100000))

    arch, dsp = archive.split(inp)
    model = autoencoder.AEStack._from_archive(arch, os.path.join(dsp, trid),
                                              rng, thrng, input_data_type)
    return model


def repr_batched_iter(model, X, layer, batch_size, dry_run):
    W = model[0].get_weights()[0]
    Wo = model[layer].get_weights()[0]
#    assert W.shape[0] % X.shape[1] == 0, ("%d, %d" % (W.shape[0], X.shape[0]))

    for i in range(0, X.shape[0], batch_size):
        bstart, bend = i, min(i + batch_size, X.shape[0])
        yield (np.zeros((bend - bstart, Wo.shape[1])) if dry_run
               else model.compute_state(layer, X[bstart:bend]))

def named_repr_iter(model, dataset, batch_size, layer, dry_run):
    biter = repr_batched_iter(model, dataset.X, layer, batch_size, dry_run)
    for bidx, x in enumerate(biter):
        names = dataset.pX.items[bidx * batch_size:(bidx+1) * batch_size]
        assert names.shape[0] == x.shape[0] == batch_size
        for i in range(batch_size):
            yield "\t".join((names[i],
                             " ".join(map(lambda v: "%.1f" % v, x[i]))))

if __name__ != '__main__':
    lg.error("this is a script do not import")
    sys.exit(1)


parser = parser_inst(__doc__)
parser.add_argument("model_arch", type=ArgTypes.is_dsarch,
                    help="Archive model was trained on." + archive.DSPEC_MSG)
parser.add_argument("trid", type=str, help="Train id")
parser.add_argument("input", type=ArgTypes.is_dsarch,
                    help="Data to compute." + archive.DSPEC_MSG)
parser.add_argument("--layer", type=ArgTypes.is_uint, default=4,
                    help="metalabels are repr. of this layer")
parser.add_argument("--batch_size", type=ArgTypes.is_pint, default=100,
                    help="batch size")
parser.add_argument("--dry_run", action='store_true', default=False,
                    help="Dry run")
opt = parser.parse_args()


#training_dataset = AnnotatedDataset._from_archive(opt.input, True)
dataset = AnnotatedDataset._from_archive(opt.input, False)
dataset.flatten()
## TODO: subtract training_mean, divide by training_sd

model = load_model(opt.model_arch, opt.trid, dataset.X.dtype)
#metalab_iter = repr_iter(model, dataset.X, opt.layer, opt.batch_size, opt.dry_run)
for line in named_repr_iter(model, dataset, opt.batch_size, opt.layer,
        opt.dry_run):
    print line
