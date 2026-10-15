"""The single-pass read evidence must reproduce the per-window implementation exactly.

Reference = the original code path: the verbatim original get_reads, no evidence scope.
Candidate = reads.add_reads_from_bam inside `evidence.evidence_scope()`.
Every array is compared including dtype and shape.
"""

import itertools
import math
import random
import shutil
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pytest

from spladder import evidence, reads

from . import synthetic


def _conf(level, readlen=50):
    return {0: dict(intron=350000, exon_len=math.ceil(readlen * 0.10), mismatch=max(2, math.floor(readlen * 0.03)), mincount=1),
            1: dict(intron=350000, exon_len=math.ceil(readlen * 0.15), mismatch=max(1, math.floor(readlen * 0.02)), mincount=2),
            2: dict(intron=350000, exon_len=math.ceil(readlen * 0.20), mismatch=max(1, math.floor(readlen * 0.01)), mincount=2),
            3: dict(intron=350000, exon_len=math.ceil(readlen * 0.25), mismatch=0, mincount=2)}[level]


FILTERS = [None, _conf(0), _conf(1), _conf(2), _conf(3), dict(exon_len=7, mismatch=2), dict(exon_len=30, mismatch=3, mincount=5)]
TYPES = [['exon_track'], ['intron_list'], ['exon_track', 'intron_list'], ['intron_list', 'exon_track']]


@pytest.fixture(scope='module')
def bam(tmp_path_factory):
    return synthetic.make_bam(tmp_path_factory.mktemp('evidence') / 'syn.bam', seed=5)


def _block(contig, start, stop, strand, numpy_ints=False):
    if numpy_ints:
        start, stop = np.int64(start), np.int64(stop)
    return np.array([SimpleNamespace(chr=contig, start=start, stop=stop, strand=strand)], dtype='object')


def _legacy(*args, **kwargs):
    """The original code path."""
    with mock.patch.object(reads, 'get_reads', getattr(reads, '_get_reads_reference', reads.get_reads)):
        return reads.add_reads_from_bam(*args, use_cache=False, **kwargs)


def _equal(a, b):
    if isinstance(a, tuple):
        return isinstance(b, tuple) and len(a) == len(b) and all(_equal(x, y) for x, y in zip(a, b))
    if isinstance(a, list):
        return isinstance(b, list) and len(a) == len(b) and all(_equal(x, y) for x, y in zip(a, b))
    return isinstance(b, np.ndarray) and a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b)


def _windows(rng, n=14):
    ws = []
    for contig, ln in synthetic.CONTIGS + [('chrZ', 1000)]:
        for _ in range(n):
            a = rng.randint(1, ln)
            b = min(ln + 300, a + rng.choice([1, 5, 40, 300, 2000, 9000, 25000]))
            ws.append((contig, a, b))
    for contig, genes in synthetic.GENES.items():
        for _, exons in genes:
            ws.append((contig, max(exons[0][0] - 5000, 1), exons[-1][1] + 5000))
            ws.append((contig, exons[0][0], exons[-1][1]))
            for x, y in zip(exons[:-1], exons[1:]):
                ws.append((contig, x[1], y[0]))
    return ws


def _random_calls(rng, windows, per_window):
    for (contig, a, b) in windows:
        for _ in range(per_window):
            yield dict(window=(contig, a, b),
                       strand=rng.choice([None, '+', '-']),
                       types=rng.choice(TYPES),
                       filter=rng.choice(FILTERS),
                       unstranded=rng.choice([True, False]),
                       primary_only=rng.choice([True, False]),
                       no_mm=rng.choice([True, False]),
                       numpy_ints=rng.choice([True, False]))


def _call(c, bam, **extra):
    contig, a, b = c['window']
    return reads.add_reads_from_bam(_block(contig, a, b, c['strand'], c['numpy_ints']), bam, c['types'], c['filter'], False,
                                    c['primary_only'], c['no_mm'], unstranded=c['unstranded'], mm_tag='NM', cram_ref=None, **extra)


def _call_legacy(c, bam):
    contig, a, b = c['window']
    return _legacy(_block(contig, a, b, c['strand'], c['numpy_ints']), bam, c['types'], c['filter'], False,
                   c['primary_only'], c['no_mm'], unstranded=c['unstranded'], mm_tag='NM', cram_ref=None)


def _run_equivalence(bam, calls):
    bad = []
    with evidence.evidence_scope() as scope:
        for c in calls:
            if c['filter'] is None and c['no_mm']:
                continue
            if not _equal(_call_legacy(c, bam), _call(c, bam)):
                bad.append(c)
    assert not bad, 'first mismatches: %s' % bad[:3]
    return scope


def test_equivalence_random_calls(bam):
    rng = random.Random(1)
    scope = _run_equivalence(bam, list(_random_calls(rng, _windows(rng), 14)))
    assert scope.stats['hits'] > 500
    assert scope.stats['fallbacks'] == 0


def test_equivalence_exhaustive_on_gene_windows(bam):
    ### every parameter combination on a few windows with plenty of spliced reads
    calls = []
    for window in [('chrA', 500, 6500), ('chrB', 100, 7000), ('chrA', 11000, 18000)]:
        for strand, types, flt, unstranded, primary_only, no_mm in itertools.product([None, '+', '-'], TYPES, FILTERS, [True, False], [True, False], [True, False]):
            calls.append(dict(window=window, strand=strand, types=types, filter=flt, unstranded=unstranded,
                              primary_only=primary_only, no_mm=no_mm, numpy_ints=False))
    _run_equivalence(bam, calls)


def _all_introns(path):
    """All introns present in the file, taken from a window spanning the complete contigs (original code path)."""
    found = set()
    for contig, ln in synthetic.CONTIGS:
        out = _legacy(_block(contig, 1, ln + 1, None), path, ['intron_list'], None, False, True, True, unstranded=True)
        found.update((contig, int(a), int(b)) for a, b, _ in out[0])
    return sorted(found)


def test_window_boundaries_at_intron_ends(bam):
    ### windows starting / ending exactly at (or one base next to) intron ends: strict vs. non strict comparisons matter here
    introns = _all_introns(bam)
    assert len(introns) > 50
    rng = random.Random(8)
    calls = []
    for (contig, a, b) in rng.sample(introns, 40):
        for d1, d2 in itertools.product([-1, 0, 1], [-1, 0, 1]):
            for (ws, we) in [(a + d1, b + d2), (a + d1, b + 60 + d2), (a - 60 + d1, b + d2)]:
                if ws >= 1 and we > ws:
                    calls.append(dict(window=(contig, ws, we), strand=rng.choice([None, '+', '-']), types=['intron_list'],
                                      filter=rng.choice([None, _conf(1), _conf(3)]), unstranded=rng.choice([True, False]),
                                      primary_only=True, no_mm=True, numpy_ints=rng.choice([True, False])))
    _run_equivalence(bam, calls)


def test_one_build_per_contig_and_filter(bam):
    builds = []
    original = evidence.build_contig_evidence

    def spy(*args, **kwargs):
        builds.append(args[:2])
        return original(*args, **kwargs)

    rng = random.Random(2)
    flt = _conf(3)
    with mock.patch.object(evidence, 'build_contig_evidence', spy):
        with evidence.evidence_scope() as scope:
            for _ in range(60):
                contig = rng.choice(['chrA', 'chrB'])
                a = rng.randint(1, 5000)
                reads.add_reads_from_bam(_block(contig, a, a + rng.randint(10, 3000), rng.choice(['+', '-'])), bam, ['exon_track', 'intron_list'], flt,
                                         False, True, False, unstranded=False, mm_tag='NM', cram_ref=None)
    assert sorted(builds) == [(bam, 'chrA'), (bam, 'chrB')]
    assert scope.stats['builds'] == 2
    assert scope.stats['hits'] == 60


def test_outside_a_scope_nothing_is_cached(bam):
    assert evidence.current_scope() is None
    with mock.patch.object(evidence, 'build_contig_evidence', side_effect=AssertionError('must not be called')):
        out = reads.add_reads_from_bam(_block('chrA', 900, 5600, '+'), bam, ['exon_track', 'intron_list'], _conf(3), False, True, False,
                                       unstranded=False, mm_tag='NM', cram_ref=None)
    assert _equal(out, _legacy(_block('chrA', 900, 5600, '+'), bam, ['exon_track', 'intron_list'], _conf(3), False, True, False,
                               unstranded=False, mm_tag='NM', cram_ref=None))


def test_use_cache_false_bypasses_the_scope(bam):
    with evidence.evidence_scope() as scope:
        reads.add_reads_from_bam(_block('chrA', 900, 5600, '+'), bam, ['exon_track'], _conf(3), False, True, False, use_cache=False)
    assert scope.stats['hits'] == 0 and scope.stats['builds'] == 0


def test_scope_is_released(bam):
    with evidence.evidence_scope() as scope:
        reads.add_reads_from_bam(_block('chrA', 900, 5600, '+'), bam, ['exon_track'], _conf(3), False, True, False)
        assert evidence.current_scope() is scope
        assert scope.nbytes > 0
    assert evidence.current_scope() is None
    assert scope.nbytes == 0


def test_environment_variable_disables_evidence(bam, monkeypatch):
    monkeypatch.setenv('SPLADDER_NO_EVIDENCE', '1')
    with evidence.evidence_scope() as scope:
        reads.add_reads_from_bam(_block('chrA', 900, 5600, '+'), bam, ['exon_track'], _conf(3), False, True, False)
    assert scope.stats['hits'] == 0 and scope.stats['builds'] == 0


def test_memory_cap_evicts_but_stays_exact(bam):
    rng = random.Random(4)
    calls = [c for c in _random_calls(rng, _windows(rng, 6), 6) if not (c['filter'] is None and c['no_mm'])]
    bad = []
    with evidence.evidence_scope(max_bytes=1) as scope:
        for c in calls:
            if not _equal(_call_legacy(c, bam), _call(c, bam)):
                bad.append(c)
    assert not bad
    assert scope.stats['evictions'] > 0


@pytest.mark.parametrize('kwargs', [dict(strand='.'),                       # unsupported strand value
                                    dict(var_aware=True),                   # variation aware filtering
                                    dict(types=['exon_track', 'exon_track']),
                                    dict(types=['intron_track'])])
def test_unsupported_requests_fall_back(bam, kwargs):
    strand = kwargs.get('strand', '+')
    types = kwargs.get('types', ['exon_track', 'intron_list'])
    var_aware = kwargs.get('var_aware', False)
    flt = dict(exon_len=7, mismatch=2) if var_aware else _conf(3)
    if var_aware:
        ### the synthetic reads carry no XM/XG tags: only the fall-back decision is of interest here
        with evidence.evidence_scope() as scope:
            with pytest.raises(KeyError):
                reads.add_reads_from_bam(_block('chrA', 900, 5600, strand), bam, types, flt, var_aware, True, False)
        assert scope.stats['hits'] == 0
        return
    args = (_block('chrA', 900, 5600, strand), bam, types, flt, var_aware, True, False)
    with evidence.evidence_scope() as scope:
        got = reads.add_reads_from_bam(*args, unstranded=False)
    assert scope.stats['hits'] == 0
    assert _equal(_legacy(*args, unstranded=False), got)


def test_several_alignment_files_fall_back(bam, tmp_path):
    other = str(tmp_path / 'other.bam')
    shutil.copy(bam, other)
    shutil.copy(bam + '.bai', other + '.bai')
    args = (_block('chrA', 900, 5600, '+'), [bam, other], ['exon_track', 'intron_list'], _conf(3), False, True, False)
    with evidence.evidence_scope() as scope:
        got = reads.add_reads_from_bam(*args, unstranded=False)
    assert scope.stats['hits'] == 0
    assert _equal(_legacy(*args, unstranded=False), got)


def test_mito_contig_and_unknown_contig_stay_empty(bam):
    with evidence.evidence_scope():
        for contig in ['MT', 'chrE', 'chrZ']:
            got = reads.add_reads_from_bam(_block(contig, 100, 2500, '+'), bam, ['exon_track', 'intron_list'], _conf(3), False, True, False, unstranded=False)
            exp = _legacy(_block(contig, 100, 2500, '+'), bam, ['exon_track', 'intron_list'], _conf(3), False, True, False, unstranded=False)
            assert _equal(exp, got)
            ### the original implementation blanks the MT contig whatever the alignments say
            if contig == 'MT':
                assert not got[0].any() and got[1][0].shape == (0, 3)


def test_missing_nm_tag_keeps_original_behaviour(tmp_path):
    ### reads without NM in one region: windows elsewhere work, windows touching the region abort exactly like before
    path = synthetic.make_bam(tmp_path / 'nonm.bam', seed=9, nm_missing_region=('chrA', 7000, 9000))
    flt = _conf(3)
    with evidence.evidence_scope() as scope:
        ok = reads.add_reads_from_bam(_block('chrA', 12000, 17500, '+'), path, ['exon_track', 'intron_list'], flt, False, True, False, unstranded=False)
        assert _equal(ok, _legacy(_block('chrA', 12000, 17500, '+'), path, ['exon_track', 'intron_list'], flt, False, True, False, unstranded=False))
        with pytest.raises(SystemExit):
            reads.add_reads_from_bam(_block('chrA', 6500, 9500, '-'), path, ['exon_track', 'intron_list'], flt, False, True, False, unstranded=False)
        ### without the mismatch filter the tag is not needed at all
        ok2 = reads.add_reads_from_bam(_block('chrA', 6500, 9500, '-'), path, ['exon_track', 'intron_list'], flt, False, True, True, unstranded=False)
        assert _equal(ok2, _legacy(_block('chrA', 6500, 9500, '-'), path, ['exon_track', 'intron_list'], flt, False, True, True, unstranded=False))
    assert scope.stats['fallbacks'] > 0
    with pytest.raises(SystemExit):
        _legacy(_block('chrA', 6500, 9500, '-'), path, ['exon_track', 'intron_list'], flt, False, True, False, unstranded=False)


def test_unindexed_file_falls_back(bam, tmp_path):
    ### pysam would return the reads of all contigs for a contig wide fetch of a file without index
    unindexed = str(tmp_path / 'unindexed.bam')
    shutil.copy(bam, unindexed)
    args = (_block('chrA', 900, 5600, '+'), unindexed, ['exon_track', 'intron_list'], _conf(3), False, True, False)
    with evidence.evidence_scope() as scope:
        with pytest.raises(Exception) as got:
            reads.add_reads_from_bam(*args, unstranded=False)
    with pytest.raises(Exception) as exp:
        _legacy(*args, unstranded=False)
    assert type(exp.value) == type(got.value)
    assert scope.stats['hits'] == 0


def _evidence_equal(a, b):
    if a is None or b is None:
        return a is b
    for (pa, va), (pb, vb) in zip(a.cov, b.cov):
        if not (_equal(pa, pb) and _equal(va, vb)):
            return False
    for strand in ('+', '-', None):
        for (ia, sa), (ib, sb) in zip(a.introns[strand], b.introns[strand]):
            if not (_equal(ia, ib) and _equal(sa, sb)):
                return False
    return True


CONTIGS = ['chrA', 'chrB', 'MT', 'chrE', 'chrZ']


def test_prefetch_gives_the_same_evidence_as_building_on_demand(bam):
    spec = evidence.make_spec(bam, _conf(3), True, False, 'NM', None)
    with evidence.evidence_scope() as lazy:
        expected = {c: lazy.contig_evidence(spec, c) for c in CONTIGS}
    with evidence.evidence_scope() as scope:
        scope.prefetch(spec, CONTIGS, 3)
        assert scope.stats['builds'] == len(CONTIGS)
        got = {c: scope.contig_evidence(spec, c) for c in CONTIGS}
        assert scope.stats['builds'] == len(CONTIGS)          # everything was prefetched
    assert all(_evidence_equal(expected[c], got[c]) for c in CONTIGS)


def test_prefetched_evidence_answers_requests_like_the_original(bam):
    rng = random.Random(6)
    flt = _conf(3)
    spec = evidence.make_spec(bam, flt, True, False, 'NM', None)
    with evidence.evidence_scope() as scope:
        scope.prefetch(spec, CONTIGS, 2)
        for _ in range(60):
            contig = rng.choice(CONTIGS)
            a = rng.randint(1, 6000)
            args = (_block(contig, a, a + rng.randint(5, 4000), rng.choice(['+', '-'])), bam, ['exon_track', 'intron_list'], flt, False, True, False)
            assert _equal(_legacy(*args, unstranded=False), reads.add_reads_from_bam(*args, unstranded=False))
    assert scope.stats['builds'] == len(CONTIGS) and scope.stats['hits'] == 60


def test_prefetch_remembers_contigs_that_cannot_be_served(tmp_path):
    path = synthetic.make_bam(tmp_path / 'nonm.bam', seed=9, nm_missing_region=('chrA', 7000, 9000))
    spec = evidence.make_spec(path, _conf(3), True, False, 'NM', None)
    with evidence.evidence_scope() as scope:
        scope.prefetch(spec, ['chrA', 'chrB'], 2)
        assert scope.stats['builds'] == 1                       # chrB; chrA has reads without NM tag
        got = reads.add_reads_from_bam(_block('chrA', 12000, 17500, '+'), path, ['exon_track'], _conf(3), False, True, False)
        assert _equal(got, _legacy(_block('chrA', 12000, 17500, '+'), path, ['exon_track'], _conf(3), False, True, False))
    assert scope.stats['fallbacks'] == 1 and scope.stats['builds'] == 1


def test_prefetch_does_nothing_for_a_single_contig_or_worker(bam):
    spec = evidence.make_spec(bam, _conf(3), True, False, 'NM', None)
    with evidence.evidence_scope() as scope:
        scope.prefetch(spec, ['chrA'], 4)
        scope.prefetch(spec, CONTIGS, 1)
    assert scope.stats['builds'] == 0
