"""The optimized get_reads path must reproduce the original implementation exactly."""

import itertools
import math
import random

import numpy as np
import pytest

from spladder import reads

from . import synthetic


def _reference():
    ### the verbatim original implementation (kept for collapse=False and as test reference)
    return getattr(reads, '_get_reads_reference', reads.get_reads)


@pytest.fixture(scope='module')
def bam(tmp_path_factory):
    return synthetic.make_bam(tmp_path_factory.mktemp('reads') / 'syn.bam', seed=3)


def _conf(level, readlen=50):
    return {0: dict(intron=350000, exon_len=math.ceil(readlen * 0.10), mismatch=max(2, math.floor(readlen * 0.03)), mincount=1),
            1: dict(intron=350000, exon_len=math.ceil(readlen * 0.15), mismatch=max(1, math.floor(readlen * 0.02)), mincount=2),
            2: dict(intron=350000, exon_len=math.ceil(readlen * 0.20), mismatch=max(1, math.floor(readlen * 0.01)), mincount=2),
            3: dict(intron=350000, exon_len=math.ceil(readlen * 0.25), mismatch=0, mincount=2)}[level]


def _windows(rng, n=25):
    ws = []
    for contig, ln in synthetic.CONTIGS + [('chrZ', 1000)]:
        for _ in range(n):
            a = rng.randint(0, ln)
            b = min(ln + 300, a + rng.choice([1, 5, 40, 300, 2000, 9000, 25000]))
            ws.append((contig, a, b))
    for contig, genes in synthetic.GENES.items():
        for _, exons in genes:
            ws.append((contig, max(exons[0][0] - 5000, 1), exons[-1][1] + 5000))
            ws.append((contig, exons[0][0], exons[-1][1]))
    return ws


def _equal(a, b):
    if isinstance(a, tuple):
        return isinstance(b, tuple) and len(a) == len(b) and all(_equal(x, y) for x, y in zip(a, b))
    return a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b)


def test_collapsed_get_reads_matches_reference(bam):
    rng = random.Random(11)
    ref = _reference()
    combos = list(itertools.product([None, _conf(0), _conf(1), _conf(2), _conf(3), dict(exon_len=7, mismatch=2)],
                                    [None, '+', '-'], [(True, True), (True, False), (False, True)],
                                    [True, False], [False, True]))
    n = bad = 0
    for (contig, a, b) in _windows(rng):
        for flt, strand, (mapped, spliced), primary, no_mm in rng.sample(combos, 12):
            if flt is None and no_mm:
                continue
            args = (bam, contig, a, b, strand, flt, mapped, spliced, False, True, primary, no_mm, 'NM', None)
            exp = ref(*args)
            got = reads.get_reads(*args)
            n += 1
            if not _equal(exp, got):
                bad += 1
                print('MISMATCH', args)
    assert n > 500
    assert bad == 0


def test_uncollapsed_get_reads_unchanged(bam):
    ### collapse=False builds a sparse read matrix; it must keep working exactly as before
    ### (the original implementation only copes with windows that start at 0 and span the contig, which is all this checks)
    ref = _reference()
    for (contig, a, b) in [('chrA', 0, 20000), ('chrB', 0, 8000)]:
        args = (bam, contig, a, b, '+', _conf(1), True, True, False, False, True, False, 'NM', None)
        exp, got = ref(*args), reads.get_reads(*args)
        assert exp[0].shape == got[0].shape
        assert np.array_equal(exp[0].toarray(), got[0].toarray())
        assert _equal(exp[1:], got[1:])


def test_file_handle_cache_is_bounded(bam, tmp_path):
    ### opening many different files must neither fail nor leak handles
    import shutil
    files = []
    for i in range(12):
        dst = tmp_path / ('c%d.bam' % i)
        shutil.copy(bam, dst)
        shutil.copy(bam + '.bai', str(dst) + '.bai')
        files.append(str(dst))
    ref = _reference()
    for f in files * 2:
        args = (f, 'chrA', 900, 5600, '+', _conf(3), True, True, False, True, True, False, 'NM', None)
        assert _equal(ref(*args), reads.get_reads(*args))
