"""Comparison of SplAdder results (objects, pickles, HDF5 files, text files, result directories)."""

import gzip
import os
import pickle

import h5py
import numpy as np
import scipy.sparse


def deep_equal(a, b, path='root'):
    """Structural comparison including dtypes; returns None if equal, else a description of the first difference."""

    if scipy.sparse.issparse(a) or scipy.sparse.issparse(b):
        return deep_equal(a.toarray(), b.toarray(), path)
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        if not (isinstance(a, np.ndarray) and isinstance(b, np.ndarray)):
            return '%s: %s vs %s' % (path, type(a), type(b))
        if a.dtype != b.dtype or a.shape != b.shape:
            return '%s: dtype/shape %s %s vs %s %s' % (path, a.dtype, a.shape, b.dtype, b.shape)
        if a.dtype == object:
            for i, (x, y) in enumerate(zip(a.ravel(), b.ravel())):
                d = deep_equal(x, y, '%s[%d]' % (path, i))
                if d:
                    return d
            return None
        return None if np.array_equal(a, b, equal_nan=(a.dtype.kind == 'f')) else '%s: values differ' % path
    if isinstance(a, dict):
        if not isinstance(b, dict) or set(a) != set(b):
            return '%s: keys differ' % path
        for k in a:
            d = deep_equal(a[k], b[k], '%s[%r]' % (path, k))
            if d:
                return d
        return None
    if isinstance(a, (list, tuple)):
        if type(a) != type(b) or len(a) != len(b):
            return '%s: length/type differ' % path
        for i, (x, y) in enumerate(zip(a, b)):
            d = deep_equal(x, y, '%s[%d]' % (path, i))
            if d:
                return d
        return None
    if hasattr(a, '__dict__') and not isinstance(a, type):
        if type(a) != type(b):
            return '%s: type differs' % path
        return deep_equal(vars(a), vars(b), path)
    return None if a == b else '%s: %r vs %r' % (path, a, b)


def compare_dirs(a, b):
    """Compares all result files of two output directories; returns a list of differences."""

    diffs = []
    files_a = sorted(os.path.relpath(os.path.join(r, f), a) for r, _, fs in os.walk(a) for f in fs if not r.startswith(os.path.join(a, 'tmp')))
    files_b = sorted(os.path.relpath(os.path.join(r, f), b) for r, _, fs in os.walk(b) for f in fs if not r.startswith(os.path.join(b, 'tmp')))
    if files_a != files_b:
        diffs.append('file lists differ: %s' % sorted(set(files_a) ^ set(files_b)))
    for f in sorted(set(files_a) & set(files_b)):
        pa, pb = os.path.join(a, f), os.path.join(b, f)
        if f.endswith('.pickle') or '.quick_' in f:
            d = deep_equal(pickle.load(open(pa, 'rb')), pickle.load(open(pb, 'rb')), f)
        elif f.endswith('.hdf5'):
            with h5py.File(pa, 'r') as x, h5py.File(pb, 'r') as y:
                d = None if set(x.keys()) == set(y.keys()) else '%s: datasets differ' % f
                for k in (x.keys() if d is None else []):
                    d = d or deep_equal(x[k][()], y[k][()], '%s:%s' % (f, k))
        else:
            op = gzip.open if f.endswith('.gz') else open
            d = None if op(pa, 'rb').read() == op(pb, 'rb').read() else '%s: content differs' % f
        if d:
            diffs.append(d)
    return diffs
