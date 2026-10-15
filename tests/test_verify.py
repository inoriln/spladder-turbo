"""Optimized helpers of the event verification must reproduce the original behaviour exactly.

(The verification of whole events is tested in test_parallel.py::test_event_verification_is_unchanged.)
"""

import random

import numpy as np
import pytest


@pytest.mark.parametrize('shape', ['flat', 'column'])
@pytest.mark.parametrize('ordered', [True, False])
def test_gene_rows_match_a_scan_of_the_ids(shape, ordered):
    from spladder.alt_splice import verify
    rng = random.Random(7)
    for n in list(range(0, 40)) * 5:
        ids = np.array(sorted(rng.randint(0, 9) for _ in range(n)), dtype='int')
        if not ordered:
            ids = np.array(rng.sample(list(ids), len(ids)), dtype='int')
        ids = ids.reshape(n, 1) if shape == 'column' else ids
        is_sorted = verify._is_sorted_ids(ids)
        assert is_sorted == bool(n < 2 or np.all(ids.reshape(-1)[1:] >= ids.reshape(-1)[:-1]))
        for g in range(-1, 11):
            exp = np.where(ids == g)[0]
            rows, index = verify._gene_rows(ids, g, is_sorted)
            assert index.dtype == exp.dtype and np.array_equal(index, exp)
            data = np.arange(n * 3).reshape(n, 3)
            assert np.array_equal(data[rows, :], data[exp, :])


def test_count_index_follows_the_count_file(tmp_path):
    """The quick index files are written once, and not used any more when the counts have been replaced."""
    import os
    import h5py
    from spladder.alt_splice import verify

    def write_counts(fname, offset):
        with h5py.File(fname, 'w') as out:
            out.create_dataset('gene_ids_segs', data=np.arange(5).reshape(5, 1) + offset)
            out.create_dataset('gene_ids_edges', data=np.arange(3).reshape(3, 1) + offset)
            out.create_dataset('edge_idx', data=np.arange(3) + offset)

    fname = str(tmp_path / 'g.count.hdf5')
    write_counts(fname, 0)
    first = verify.load_count_index(fname)
    assert [x.shape for x in first] == [(5, 1), (3, 1), (3,)]
    quick = sorted(f for f in os.listdir(str(tmp_path)) if '.quick_' in f)
    assert quick == ['g.count.hdf5.quick_edge_idx', 'g.count.hdf5.quick_ids_edges', 'g.count.hdf5.quick_ids_segs']

    ### used as long as the counts are the same: the count file is not opened (a corrupt one would fail)
    stamp = os.stat(fname).st_mtime_ns
    open(fname + '.quick_edge_idx', 'ab').close()
    os.utime(fname + '.quick_edge_idx', ns=(stamp + 10 ** 9, stamp + 10 ** 9))
    again = verify.load_count_index(fname)
    assert all(np.array_equal(a, b) for a, b in zip(first, again))

    ### replaced counts: the old index must not be used
    write_counts(fname, 100)
    os.utime(fname, ns=(stamp + 5 * 10 ** 9, stamp + 5 * 10 ** 9))
    for f in quick:
        os.utime(str(tmp_path / f), ns=(stamp + 10 ** 9, stamp + 10 ** 9))
    fresh = verify.load_count_index(fname)
    assert [int(x.min()) for x in fresh] == [100, 100, 100]
    ### (and the index written then is used again)
    assert [int(x.min()) for x in verify.load_count_index(fname)] == [100, 100, 100]

    ### without the count file the error is the one of opening it
    os.remove(fname)
    with pytest.raises(OSError):
        verify.load_count_index(fname)
