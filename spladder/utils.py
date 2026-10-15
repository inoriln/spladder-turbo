import contextlib
import gc
import os
import pickle
import uuid

import numpy as np

def jit_cached(fn):
    """Compiles the function with numba and caches the machine code on disk, so that it is not compiled again in every process."""

    from numba import jit
    try:
        return jit(nopython=True, cache=True)(fn)
    except RuntimeError:
        ### no writable location for the cache
        return jit(nopython=True)(fn)


def load_pickle(fname, **kwargs):
    """pickle.load of a file with very many small objects (the splice graphs of all genes) with the garbage collector paused.

    Creating that many objects triggers collections of an ever growing heap, which about doubles the loading time.
    There is no garbage in the loaded data, so the objects are exactly the same."""

    enabled = gc.isenabled()
    gc.disable()
    try:
        with open(fname, 'rb') as fd:
            return pickle.load(fd, **kwargs)
    finally:
        if enabled:
            gc.enable()


def dump_pickle(obj, fname):
    """pickle.dump to a file with the garbage collector paused (saves about a third of the time).

    The file is written in one piece: a process reading the file while it is written gets the old or the complete new one."""

    target = os.path.realpath(fname)
    ### (processes on different hosts writing to the same directory can have the same process id)
    fname_tmp = '%s.tmp%s' % (target, uuid.uuid4().hex[:12])
    enabled = gc.isenabled()
    gc.disable()
    try:
        try:
            fd = open(fname_tmp, 'wb')
        except OSError:
            ### no place for a temporary file: write the file itself
            with open(target, 'wb') as fd:
                pickle.dump(obj, fd, -1)
            return
        try:
            with fd:
                pickle.dump(obj, fd, -1)
            os.replace(fname_tmp, target)
        except BaseException:
            if os.path.exists(fname_tmp):
                os.remove(fname_tmp)
            raise
    finally:
        if enabled:
            gc.enable()


### the file loaded last by load_pickle_shared: (identity of the file, content)
_SHARED = [None, None]

def load_pickle_shared(fname, **kwargs):
    """Like load_pickle, but a file that was loaded by this process before and has not changed since is not read again.

    All callers get the same object: it must only be read (or copied by forked processes)."""

    st = os.stat(fname)
    key = (os.path.abspath(fname), st.st_mtime_ns, st.st_size, st.st_ino, tuple(sorted(kwargs.items())))
    if _SHARED[0] != key:
        ### only one file at a time: the graphs are large
        _SHARED[0], _SHARED[1] = None, None
        _SHARED[1] = load_pickle(fname, **kwargs)
        _SHARED[0] = key
    return _SHARED[1]


@contextlib.contextmanager
def frozen_heap():
    """Processes forked inside of this context do not collect (and thereby copy) the objects the parent process holds."""

    gc.freeze()
    try:
        yield
    finally:
        gc.unfreeze()


def isequal(A, B):
    
    if A is None or B is None:
        return False

    if not np.all(A.shape == B.shape):
        return False
    else:
        return np.all(A == B)

def issubset(A, B):
    # function found = issubset(A,B)
    #
    # returns true if A is a subset of B, where 
    # both A and B are vectors with 1 where elements exists and 0 otherwise.
    return np.all(np.where(A)[0] == np.isin(np.where(A)[0], np.where(B)[0]))


#def intersect_rows(array1, array2, index = None):
#    """Return intersection of rows"""
#
#    if (array1.shape[0] == 0):
#        if index == True:
#            return (array1, np.zeros((0,)), np.zeros((0,)))
#        else:
#            return array1
#    if (array2.shape[0] == 0):
#        if index == True:
#            return (array2, np.zeros((0,)), np.zeros((0,)))
#        else:
#            return array2
#
#    array1_v = array1.view([('', array1.dtype)] * array1.shape[1])
#    array2_v = array2.view([('', array2.dtype)] * array2.shape[1])
#    array_i = np.intersect1d(array1_v, array2_v)
#
#    if index == True:
#        a1_i = np.where(np.isin(array1_v, array_i))[0]
#        a2_i = np.where(np.isin(array2_v, array_i))[0]
#        return (array_i.view(array1.dtype).reshape(array_i.shape[0], array1.shape[1]), a1_i, a2_i)
#    else:
#        return array_i.view(array1.dtype).reshape(array_i.shape[0], array1.shape[1])


def unique_rows(array, index = None):
    """Make array unique by rows"""

    if array.shape[0] == 0:
        if index == True:
            return (array, [])
        else:
            return (array)

    if len(array.shape) == 1:
        if index == True:
            return (array, [0])
        else:
            return array

    (array_s, s_idx) = sort_rows(array, True)
    ### keep the first of every group of identical (neighbouring) rows
    keep = np.ones((array.shape[0],), dtype='bool')
    keep[1:] = ~np.all(array_s[1:, :] == array_s[:-1, :], axis=1)
    k_idx = np.where(keep)[0]

    if index == True:
        return (array[s_idx[k_idx], :], s_idx[k_idx])
    else:
        return array[s_idx[k_idx], :]


def intersect_rows(array1, array2, index=False):
    """Return array with rows that intersect between array1 and array2"""

    tmp1 = np.array(['-'.join(array1[i, :].astype('str')) for i in range(array1.shape[0])])
    tmp2 = np.array(['-'.join(array2[i, :].astype('str')) for i in range(array2.shape[0])])
    
    idx = np.where(np.isin(tmp1, tmp2))[0]
    if index:
        idx2 = np.where(np.isin(tmp2, tmp1))[0]

    if index:
        return (array1[idx, :], idx, idx2)
    else:
        return (array1[idx, :], None, None)

def sort_rows(array, index = None):
    """Sort array by rows"""

    ### empty array
    if array.shape[0] == 0:
        if index == True:
            return (array, [])
        else:
            return (array)

    ### only one row
    if len(array.shape) == 1:
        if index == True:
            return (array, [0])
        else:
            return (array)

    ### more than one row
    s_idx = np.lexsort([array[:, -i] for i in range(1, array.shape[1] + 1)])

    if index == True:
        return (array[s_idx, :], s_idx)
    else:
        return array[s_idx, :]


def ismember(element, array, rows=False):
    """Check if element is member of array"""

    if rows:
        return np.any(np.all(array == element, axis=1))
    else:
        return np.all([element[i] in array for i in element.shape[0]])

        ### TODO I think this is not quite right yet

def replace_sub_matrix(mat_in, idx, mat_put):
    """Replaces the values in mat_in in rows and cols idx with values of mat_put"""
    
    assert((idx.shape[0] * idx.shape[0]) == mat_put.ravel().shape[0])

    np.put(mat_in, np.ravel_multi_index([[x for x in idx for _ in idx], [x for _ in idx for x in idx]], (mat_in.shape[0], mat_in.shape[1])), mat_put.ravel())

    return mat_in

def flatten(l):
    
    r = []
    for i in l:
        if isinstance(i, list):
            r.extend(flatten(i))
        else:
            r.append(i)
    return r
