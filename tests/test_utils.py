"""Optimized helper functions must reproduce the original implementations exactly."""

import os
import random

import numpy as np
import pytest

from spladder import utils
from spladder.utils import sort_rows


def reference_unique_rows(array, index=None):
    """Verbatim original implementation of utils.unique_rows."""

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
    tmp = [False]
    tmp.extend([np.all(array_s[i-1, :] == array_s[i, :]) for i in range(1, array.shape[0])])
    k_idx = np.where(~np.array(tmp, dtype='bool'))[0]

    if index == True:
        return (array[s_idx[k_idx], :], s_idx[k_idx])
    else:
        return array[s_idx[k_idx], :]


def reference_ismember(element, array, rows=False):
    """Verbatim original implementation of utils.ismember (rows=True)."""

    if rows:
        return np.any([np.all(array[x, :] == element) for x in range(array.shape[0])])


def _same(a, b):
    if isinstance(a, tuple):
        return isinstance(b, tuple) and len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    if isinstance(a, list):
        return isinstance(b, list) and a == b
    return isinstance(b, np.ndarray) and a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b)


@pytest.mark.parametrize('dtype', ['int64', 'uint32', 'int32', 'float64'])
@pytest.mark.parametrize('index', [None, True])
def test_unique_rows_matches_original(dtype, index):
    rng = random.Random(3)
    for n in list(range(0, 25)) * 20:
        ncols = rng.choice([1, 2, 3, 4])
        a = np.array([[rng.randint(0, 3) for _ in range(ncols)] for _ in range(n)], dtype=dtype).reshape(n, ncols)
        assert _same(reference_unique_rows(a, index), utils.unique_rows(a, index)), a
    ### one dimensional input is returned as it is
    a = np.arange(5, dtype=dtype)
    assert _same(reference_unique_rows(a, index), utils.unique_rows(a, index))


@pytest.mark.parametrize('dtype', ['int64', 'uint32'])
def test_ismember_rows_matches_original(dtype):
    rng = random.Random(4)
    for n in list(range(0, 15)) * 30:
        a = np.array([[rng.randint(0, 4), rng.randint(0, 4)] for _ in range(n)], dtype=dtype).reshape(n, 2)
        element = [np.int64(rng.randint(0, 4)), np.uint32(rng.randint(0, 4))]
        exp, got = reference_ismember(element, a, rows=True), utils.ismember(element, a, rows=True)
        assert exp == got and type(exp) == type(got), (a, element)


def test_load_pickle_gives_the_object_and_restores_the_collector(tmp_path):
    import gc
    import pickle
    obj = {'a': [np.arange(5), 'x'], 'b': (1, 2.5)}
    fname = str(tmp_path / 'o.pickle')
    pickle.dump(obj, open(fname, 'wb'), -1)
    for enabled in (True, False):
        (gc.enable if enabled else gc.disable)()
        try:
            loaded = utils.load_pickle(fname)
            assert gc.isenabled() == enabled
        finally:
            gc.enable()
        assert loaded['b'] == obj['b'] and np.array_equal(loaded['a'][0], obj['a'][0])
    ### the collector is on again, even if loading fails
    open(fname, 'wb').write(b'not a pickle')
    with pytest.raises(Exception):
        utils.load_pickle(fname)
    assert gc.isenabled()


def test_shared_pickle_is_loaded_once_per_file_and_version(tmp_path):
    import pickle
    f1, f2 = str(tmp_path / 'one.pickle'), str(tmp_path / 'two.pickle')
    pickle.dump([1, 2, 3], open(f1, 'wb'))
    pickle.dump([4], open(f2, 'wb'))
    first = utils.load_pickle_shared(f1)
    assert utils.load_pickle_shared(f1) is first                 # not read again
    assert utils.load_pickle_shared(f2) == [4]                   # another file
    assert utils.load_pickle_shared(f1) is not first             # only one file is kept
    pickle.dump([1, 2, 3, 4, 5, 6], open(f1, 'wb'))              # a changed file is read again
    assert utils.load_pickle_shared(f1) == [1, 2, 3, 4, 5, 6]
    ### options of the unpickler are part of the identity
    assert utils.load_pickle_shared(f1, encoding='latin1') == [1, 2, 3, 4, 5, 6]


def test_frozen_heap_unfreezes():
    import gc
    before = gc.get_freeze_count()
    with utils.frozen_heap():
        assert gc.get_freeze_count() > before
    assert gc.get_freeze_count() == before


def test_dump_pickle_writes_the_file_in_one_piece(tmp_path):
    import gc
    import pickle
    obj = {'a': np.arange(10), 'b': ['x', (1, 2)]}
    fname = str(tmp_path / 'o.pickle')
    utils.dump_pickle(obj, fname)
    loaded = pickle.load(open(fname, 'rb'))
    assert np.array_equal(loaded['a'], obj['a']) and loaded['b'] == obj['b']
    assert os.listdir(str(tmp_path)) == ['o.pickle']              # no temporary file left
    assert gc.isenabled()
    ### an existing file is replaced, also through a link
    target = tmp_path / 'target.pickle'
    pickle.dump('old', open(str(target), 'wb'))
    link = tmp_path / 'link.pickle'
    os.symlink(str(target), str(link))
    utils.dump_pickle([1, 2, 3], str(link))
    assert os.path.islink(str(link)) and pickle.load(open(str(target), 'rb')) == [1, 2, 3]
    ### failing leaves the old file and no temporary file, and the collector on
    class Unpicklable:
        def __reduce__(self):
            raise ValueError('no')
    with pytest.raises(ValueError):
        utils.dump_pickle(Unpicklable(), str(target))
    assert pickle.load(open(str(target), 'rb')) == [1, 2, 3]
    assert sorted(os.listdir(str(tmp_path))) == ['link.pickle', 'o.pickle', 'target.pickle']
    assert gc.isenabled()


def test_dump_pickle_of_concurrent_writers_is_never_seen_partially(tmp_path):
    """Processes writing one file while another reads it, even if they have the same process id (on different hosts, in containers)."""
    import pickle
    import traceback
    fname = str(tmp_path / 'shared.pickle')
    size, rounds, writers = 300000, 25, 3
    utils.dump_pickle(np.zeros(size, dtype=np.int64), fname)
    pids = []
    for w in range(writers):
        pid = os.fork()
        if pid == 0:
            status = 1
            try:
                os.getpid = lambda: 4242
                for r in range(rounds):
                    utils.dump_pickle(np.full(size, w * 1000 + r, dtype=np.int64), fname)
                status = 0
            except BaseException:
                traceback.print_exc()
            finally:
                os._exit(status)
        pids.append(pid)
    running, seen = set(pids), 0
    try:
        while running:
            data = pickle.load(open(fname, 'rb'))                 # a partial file raises here
            assert data.shape == (size,) and (data == data[0]).all()
            seen += 1
            for pid in list(running):
                done, status = os.waitpid(pid, os.WNOHANG)
                if done:
                    running.discard(pid)
                    assert status == 0, 'a writer failed'
    finally:
        for pid in running:
            os.kill(pid, 9)
            os.waitpid(pid, 0)
    assert seen > 0
    assert os.listdir(str(tmp_path)) == ['shared.pickle']          # no temporary file left
