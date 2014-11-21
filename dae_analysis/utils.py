# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os
import logging
import inspect
from operator import itemgetter, attrgetter
from collections import namedtuple, OrderedDict
from itertools import permutations, product

import numpy as np
import pandas as pd

from sklearn import metrics
from sklearn.cluster import KMeans, Ward

logging.basicConfig(level=logging.INFO)
lg = logging.getLogger()


class CR(namedtuple("CR_tup", "pca precision recall f1")):
    @classmethod
    def _from_pred(cls, y_pred, y_true, pca=0, red=False):
        if red:
            y_pred = cls.reduce_labels(y_pred)
            y_true = cls.reduce_labels(y_true)
        metrics_f = metrics.precision_recall_fscore_support
        rep = metrics_f(y_true, cls.relabel(y_true, y_pred),
                        average='weighted', pos_label=None)
        return cls._make((pca, rep[0], rep[1], rep[2]))

    @classmethod
    def _toDF(cls, cr_lst):
        idx, cols = CR._fields[0], CR._fields[1:]
        return pd.DataFrame(np.array(map(attrgetter(*cols), cr_lst)),
                            index=map(attrgetter(idx), cr_lst),
                            columns=cols)

    @classmethod
    def relabel(cls, y_true, y_pred):
        """relabel `y_pred` so that it has the highest agreement with `y_true`

        e.g., for y_pred = [0, 0, 0, 1] and y_true = [1, 1, 1, 0]
        the classifier is pretty good, just the cluster code is off"""

        target_classes = np.unique(np.concatenate((y_true, y_pred)))
        target_classes = sorted(target_classes.tolist())
        metrics_f = metrics.precision_recall_fscore_support

        instances = []
        for p in permutations(target_classes, len(target_classes)):
            mp = dict(zip(p, target_classes))
            new_y_pred = map(lambda v: mp[v], y_pred)

            rep = metrics_f(y_true, new_y_pred,
                            average='weighted', pos_label=None)
            instances.append((mp, cls._make((0, rep[0], rep[1], rep[2]))))

        mp, inst = sorted(instances, key=lambda t: t[1].f1)[-1]
        return np.array(map(lambda v: mp[v], y_pred))

    @classmethod
    def confmat(cls, ds, y_pred, y_true, red=False):
        if red:
            y_pred = cls.reduce_labels(y_pred)
            y_true = cls.reduce_labels(y_true)
            label_names = ["E", "TTR/L"]
        else:
            label_names = ds.label_names
        M = metrics.confusion_matrix(y_true, cls.relabel(y_true, y_pred))
        col_names = map(lambda s: "%s (pred)" % s, label_names)
        row_names = map(lambda s: "%s (true)" % s, label_names)
        return pd.DataFrame(data=OrderedDict(map(lambda (i, c): (c, M[:, i]),
                                             enumerate(col_names))),
                            index=row_names)

    @staticmethod
    def reduce_labels(y):
        "merge the classes 1 and 2 together"

        yj = y.copy()
        yj[y == 2] = 1
        return yj


def _repr_paths(d, cond_f, cond_f_i):
    to_path = lambda s: os.path.join(d, s)
    return filter(cond_f_i, map(to_path, filter(cond_f, os.listdir(d))))


def dAE_repr_paths(layer=None):
    """get the paths of dAE representations dumped in 'dAE'"""

    cond_f = lambda s: s.startswith("tr__") and s.endswith(".repr")
    cond_f_i = lambda s: (True if layer is None
                          else s.endswith(".%d.repr" % layer))
    return _repr_paths("./dAE", cond_f, cond_f_i)


def PCA_repr_paths(sz=None, dsname="full"):
    cond_f = lambda s: s.startswith(dsname) and s.endswith(".txt")
    cond_f_i = lambda s: True if sz is None else s.endswith("_%d.txt" % sz)
    return _repr_paths("./PCA", cond_f, cond_f_i)


def load_repr(ds, fpath):
    """load pre-computed representations.
    each line is: <label_name><tab><representation>"""

    cache_fpath = "%s.npz" % fpath

    if os.path.isfile(cache_fpath):
        with open(cache_fpath) as ofd:
            npzfile = np.load(ofd)
            lg.debug("loaded representations %s from %s",
                     str(npzfile['reduced_data'].shape), cache_fpath)
            return npzfile['y_true'], npzfile['reduced_data']

    lab2code = lambda l: ds.label_names.index(l.split("\t")[0].split("_")[0])
    rep2code = lambda l: map(float, l.split("\t")[1].split())

    with open(fpath) as fd:
        A = map(lambda l: (lab2code(l), rep2code(l)), fd)
    y_true = np.array(map(itemgetter(0), A))
    reduced_data = np.array(map(itemgetter(1), A))
    lg.debug("loaded representations %s from %s",
             str(reduced_data.shape), fpath)
    ## save cache before returning
    with open(cache_fpath, 'w') as ofd:
        np.savez(ofd, y_true=y_true, reduced_data=reduced_data)
    return y_true, reduced_data


def ward(X, n_clust):
    "H"

    ward = Ward(n_clusters=n_clust)
    ward.fit(X)
    return ward


def km(X, n_clust=3):
    "KM"

    kme = KMeans(init="k-means++", n_clusters=n_clust,
                 n_init=5, n_jobs=-1, verbose=0)
    kme.fit(X)
    return kme


def train(ds, sample_idx, pca_sizes=[2, 10], layers=[3],
          cluster_methods=[km, ward], n_clust=3, keep_data=False):
    def pca_loader(sz, sample_idx=sample_idx):
        y, x = load_repr(ds, PCA_repr_paths(sz=sz)[0])
        return y[sample_idx], x[sample_idx]

    def dAE_loader(layer, sample_idx=sample_idx):
        layer_deep_data = map(lambda p: load_repr(ds, p),
                              dAE_repr_paths(layer=layer))

        yt = np.array(zip(*layer_deep_data)[0]).mean(axis=0)
        X = np.array(zip(*layer_deep_data)[1]).mean(axis=0)

        return yt[sample_idx], X[sample_idx]

    ## PCA
    pca_models = []
    for (pca_size, clst_f) in product(pca_sizes, cluster_methods):
        yt, X = pca_loader(pca_size)
        model_name = "PCA%d_%s" % (pca_size, inspect.getdoc(clst_f))
        lg.info("training %s on X%s", model_name, str(X.shape))
        pca_models.append((model_name, clst_f(X, n_clust=n_clust),
                           X if keep_data else None, yt))

    ## RAW
    raw_models = []
    for clst_f in cluster_methods:
        X, yt = ds.X[sample_idx], ds.Y[sample_idx]
        model_name = "RAW_%s" % inspect.getdoc(clst_f)
        lg.info("training %s on X%s", model_name, str(X.shape))
        raw_models.append((model_name, clst_f(X, n_clust=n_clust),
                           X if keep_data else None, yt))

    ## dAE
    dae_models = []
    for (layer, clst_f) in product(layers, cluster_methods):
        yt, X = dAE_loader(layer)
        model_name = "dAE%d_%s" % (layer, inspect.getdoc(clst_f))
        lg.info("training %s on X%s", model_name, str(X.shape))
        dae_models.append((model_name, clst_f(X, n_clust=n_clust),
                           X if keep_data else None, yt))


    return pca_models + dae_models + raw_models


def predict(model_results):
    P = []
    for model_name, model, X, yt in model_results:
        P.append((model_name,
                  pd.DataFrame({"y_true": yt,
                                "y_pred": model.labels_})))
    return pd.Panel(OrderedDict(P))


def analyze(P):
    analyze_f = lambda cn: CR._from_pred(P[cn]["y_pred"].values,
                                         P[cn]["y_true"].values, cn)
    return CR._toDF(map(analyze_f, P.items))


def scatter_by_label(plt, X, yt, label_names):
    """E.g.,
    yt, X = load_repr(ds, PCA_repr_paths(2)[0])
    scatter_by_label(ds, X, yt)
    """

    f, ax = plt.subplots()
    for i, c in enumerate('rgb'):
        ax.scatter(X[yt == i, 0], X[yt == i, 1], color=c, alpha=0.5,
                   label=label_names[i], s=2)
    plt.legend()


if __name__ == "__main__":
    from sklearn.preprocessing import scale
    from dimer.nnet import config_spec
    from dimer import data

    full_ds_spec = config_spec.ExpSpec("W500000_B200.h5", "full", "")
    ds = data.AnchorDataset._from_archive(full_ds_spec.ds_path, True)
    ds.flatten()
    y_true, X = load_repr(ds, dAE_repr_paths(layer=3)[0])

    km_sc = km(scale(X))
    km_r = km(X)
