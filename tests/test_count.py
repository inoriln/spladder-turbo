"""The per gene counting of segments and edges must reproduce the original loops exactly."""

import random

import numpy as np
import pytest

from spladder import count
from spladder.classes.segmentgraph import Segmentgraph

from . import reference
from .test_gene import _random_gene


def _case(rng, n_exons, coord_range, density, intron_style):
    gene = _random_gene(rng, n_exons, coord_range, density, rng.choice([0.0, 0.4]), True)
    gene.segmentgraph = Segmentgraph(gene)
    segs = gene.segmentgraph.segments
    if segs.shape[1] == 0:
        return None
    gene.start = segs.ravel().min()
    gene.stop = segs.ravel().max()
    length = int(gene.stop - gene.start)
    ### integer valued coverage as float, like the read evidence
    tracks = np.array([[rng.choice([0, 0, 1, 2, 5, 17, 300]) for _ in range(length)]], dtype='float64')
    ### introns: some match edges of the segment graph exactly, some do not, some occur twice (like strand specific duplicates)
    k, l = np.where(gene.segmentgraph.seg_edges == 1)
    rows = []
    for m in range(k.shape[0]):
        if rng.random() < 0.6:
            rows.append([segs[1, k[m]], segs[0, l[m]], rng.randint(1, 50)])
            if rng.random() < 0.3:
                rows.append([segs[1, k[m]], segs[0, l[m]], rng.randint(1, 50)])
    for _ in range(rng.randint(0, 5)):
        rows.append([rng.randint(0, coord_range), rng.randint(0, coord_range), rng.randint(1, 9)])
    if intron_style == 'empty':
        rows = []
    introns = np.array(rows, dtype=intron_style if intron_style in ('uint32', 'int64') else 'uint32').reshape(len(rows), 3)
    if introns.shape[0]:
        introns = introns[np.lexsort((introns[:, 1], introns[:, 0]))]
    return gene, tracks, introns


@pytest.mark.parametrize('intron_style', ['uint32', 'int64', 'empty'])
@pytest.mark.parametrize('density', [0.0, 0.3, 0.7])
def test_count_gene_matches_original(density, intron_style):
    rng = random.Random(int(density * 10) + len(intron_style))
    done = matched_float = 0
    for n in list(range(1, 14)) * 40:
        case = _case(rng, n, rng.choice([20, 60, 300]), density, intron_style)
        if case is None:
            continue
        gene, tracks, introns = case
        exp = reference.count_gene_reference(gene, tracks, introns)
        got = count._count_gene(gene, tracks, introns, 0)
        for attr in ('segments', 'seg_pos', 'edges'):
            a, b = getattr(exp, attr), getattr(got, attr)
            assert a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b), (attr, a, b)
        done += 1
        matched_float += exp.edges.dtype == np.float64 and exp.edges.shape[0] > 0
    assert done > 100
    if intron_style == 'uint32' and density > 0:
        assert matched_float > 10      # matched introns turn the edge counts into floats, which has to be reproduced


def _random_counts_table(rng, n_samples, n_genes, edge_dtypes, with_edges):
    """counts_tmp[sample, gene]: Counts objects; all samples count the same segments and edges of a gene, the dtypes of the edge counts differ."""
    from spladder.classes.counts import Counts
    table = np.zeros((n_samples, n_genes), dtype='object')
    for c in range(n_genes):
        n_seg = rng.choice([0, 1, 2, 3, 5, 8])
        m = rng.choice([0, 0, 1, 2, 4]) if with_edges else 0
        flat = np.sort(np.array([rng.randint(0, 50) for _ in range(m)]))
        for s in range(n_samples):
            counts = Counts(n_seg)
            counts.segments[:] = [rng.choice([0, 1.5, 7, 100]) for _ in range(n_seg)]
            counts.seg_pos[:] = [rng.randint(0, 40) for _ in range(n_seg)]
            if m:
                dtype = rng.choice(edge_dtypes)
                counts.edges = np.c_[flat, np.array([rng.randint(0, 99) for _ in range(m)])].astype(dtype)
            table[s, c] = counts
    return table


def _stacked(table, function):
    """What count_graph_coverage_wrapper makes of the pieces gathered by `function`."""
    counts = dict(segments=[], seg_pos=[], gene_ids_segs=[], edges=[], gene_ids_edges=[])
    function(table, counts)
    for key in counts:
        counts[key] = np.vstack(counts[key]) if len(counts[key]) > 0 else counts[key]
    return counts


@pytest.mark.parametrize('edge_dtypes', [['int64'], ['float64'], ['int64', 'float64']])
@pytest.mark.parametrize('n_samples', [1, 2, 3])
@pytest.mark.parametrize('with_edges', [True, False])
def test_counts_of_all_genes_are_put_together_like_gene_by_gene(n_samples, edge_dtypes, with_edges):
    rng = random.Random(n_samples + len(edge_dtypes))
    for n_genes in (1, 2, 7, 40):
        table = _random_counts_table(rng, n_samples, n_genes, edge_dtypes, with_edges)
        exp = _stacked(table, count._collect_counts_by_gene)
        got = _stacked(table, count._collect_counts)
        assert sorted(exp) == sorted(got)
        for key in exp:
            a, b = exp[key], got[key]
            if isinstance(a, list):
                assert isinstance(b, list) and a == b == []
            else:
                assert a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b), key


def test_unusual_counts_are_put_together_gene_by_gene():
    """Samples that did not count the same segments of a gene are not mixed up by the matrix version."""
    rng = random.Random(1)
    table = _random_counts_table(rng, 2, 6, ['int64'], True)
    table[1, 3].segments = np.zeros(table[0, 3].segments.shape[0] + 1)
    table[1, 3].seg_pos = np.zeros(table[0, 3].segments.shape[0] + 1)
    calls = []
    original = count._collect_counts_by_gene
    count._collect_counts_by_gene = lambda *args: calls.append(1) or original(*args)
    try:
        counts = dict(segments=[], seg_pos=[], gene_ids_segs=[], edges=[], gene_ids_edges=[])
        with pytest.raises(ValueError):
            count._collect_counts(table, counts)
    finally:
        count._collect_counts_by_gene = original
    assert calls == [1]
