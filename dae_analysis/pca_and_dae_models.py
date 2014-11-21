import sys
import random
from dimer import data
from dimer.nnet import config_spec
from utils import *
from matplotlib import pyplot as plt
import cPickle


## DATASET ##
full_ds_spec = config_spec.ExpSpec("W500000_B200.h5", "full", "")
ds = data.AnchorDataset._from_archive(full_ds_spec.ds_path, True)
ds.flatten()
SMPL_IDX = random.sample(xrange(ds.X.shape[0]), 22)
SMPL_IDX = range(ds.X.shape[0])

L = [0, 1, 2, 3]
P = [2, 5, 10, 100, 1000]

kmtr = train(ds, SMPL_IDX, pca_sizes=P, layers=L,
            cluster_methods=[km], keep_data=True)
score_data = map(lambda (n, m, x, y): (n, -m.score(x)), kmtr)

with open('km_scores.pkl', 'w') as ofd:
    cPickle.dump(pd.Series(data=map(itemgetter(1), score_data),
                           index=map(itemgetter(0), score_data)),
                 ofd)

tr = (train(ds, SMPL_IDX, pca_sizes=P, layers=L, cluster_methods=[ward])
      + map(lambda (name, model, x, y): (name, model, None, y), kmtr))

with open('tr.pkl', 'w') as ofd:
    cPickle.dump(tr, ofd)
