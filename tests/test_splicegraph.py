"""Splicegraph operations must reproduce the original implementations exactly."""

import copy
import random

import numpy as np
import pytest

from spladder.classes.splicegraph import Splicegraph
from spladder.utils import sort_rows


def reference_uniquify(self):
    """Verbatim original implementation of Splicegraph.uniquify."""

    self.sort()
    (s_tmp, s_idx) = sort_rows(self.vertices.T, index=True)
    self.vertices = s_tmp.T
    self.edges = self.edges[s_idx, :][:, s_idx]
    self.terminals = self.terminals[:, s_idx]

    rm_idx = []
    for j in range(1, self.vertices.shape[1]):
        if np.all(self.vertices[:, j-1] == self.vertices[:, j]):
            self.edges[:, j] = self.edges[:, j-1] | self.edges[:, j]
            self.edges[j, :] = self.edges[j-1, :] | self.edges[j, :]
            rm_idx.append(j - 1)

    keep_idx = np.where(~np.isin(np.array(np.arange(self.vertices.shape[1])), rm_idx))[0]
    self.vertices = self.vertices[:, keep_idx]
    self.edges = self.edges[keep_idx, :][:, keep_idx]
    self.terminals = self.terminals[:, keep_idx]


def random_splicegraph(rng, n, coord_range, symmetric):
    sg = Splicegraph()
    starts = np.array([rng.randint(0, coord_range) for _ in range(n)], dtype='int')
    stops = starts + np.array([rng.randint(1, coord_range) for _ in range(n)], dtype='int')
    sg.vertices = np.vstack([starts, stops]).reshape(2, n)
    edges = np.array([[rng.randint(0, 1) if rng.random() < 0.3 else 0 for _ in range(n)] for _ in range(n)], dtype='int').reshape(n, n)
    if symmetric:
        edges = np.triu(edges, 1)
        edges = edges + edges.T
    sg.edges = edges
    sg.terminals = np.array([[rng.randint(0, 1) for _ in range(n)] for _ in range(2)], dtype='int').reshape(2, n)
    return sg


def _same(a, b):
    return (a.vertices.dtype == b.vertices.dtype and a.vertices.shape == b.vertices.shape and np.array_equal(a.vertices, b.vertices) and
            a.edges.dtype == b.edges.dtype and a.edges.shape == b.edges.shape and np.array_equal(a.edges, b.edges) and
            a.terminals.dtype == b.terminals.dtype and a.terminals.shape == b.terminals.shape and np.array_equal(a.terminals, b.terminals))


@pytest.mark.parametrize('symmetric', [True, False])
@pytest.mark.parametrize('coord_range', [3, 6, 40])
def test_uniquify_matches_original(symmetric, coord_range):
    rng = random.Random(coord_range * 10 + int(symmetric))
    for n in list(range(0, 13)) * 60:
        sg = random_splicegraph(rng, n, coord_range, symmetric)
        exp, got = copy.deepcopy(sg), copy.deepcopy(sg)
        reference_uniquify(exp)
        got.uniquify()
        assert _same(exp, got), (n, sg.vertices, sg.edges)


def test_uniquify_returns_fresh_arrays():
    ### callers keep references to the vertices and modify the graph afterwards: no aliasing with the arrays before the call
    rng = random.Random(5)
    sg = random_splicegraph(rng, 8, 40, True)
    before = (sg.vertices, sg.edges, sg.terminals)
    sg.uniquify()
    for old, new in zip(before, (sg.vertices, sg.edges, sg.terminals)):
        assert not np.shares_memory(old, new)


def _random_transcripts(rng):
    """Transcripts of a gene as arrays of exons: few distinct exons, so that exons occur in several transcripts; single exon transcripts, repeated exons, identical transcripts."""
    pool = [(a, a + rng.randint(1, 30)) for a in rng.sample(range(0, 400, 7), rng.randint(2, 9))]
    transcripts = []
    for _ in range(rng.randint(0, 6)):
        if rng.random() < 0.2:
            exons = [rng.choice(pool)]
        else:
            n = rng.randint(2, 6)
            exons = sorted(rng.choice(pool) for _ in range(n)) if rng.random() < 0.7 else [rng.choice(pool) for _ in range(n)]
        transcripts.append(np.array(exons, dtype='int').reshape(len(exons), 2))
    return transcripts


def test_splicegraph_from_gene_matches_original():
    import random
    from spladder.classes.gene import Gene
    from spladder.classes.splicegraph import Splicegraph
    from . import reference
    rng = random.Random(7)
    sizes = set()
    for _ in range(3000):
        gene = Gene(name='G', start=0, stop=10, chr='chr1', strand='+')
        gene.transcripts = ['t%d' % i for i in range(0)]
        for i, exons in enumerate(_random_transcripts(rng)):
            gene.transcripts.append('t%d' % i)
            gene.exons.append(exons)
        expected = Splicegraph()
        reference.splicegraph_from_gene_reference(expected, gene)
        got = Splicegraph(gene)
        for attr in ('vertices', 'edges', 'terminals'):
            a, b = getattr(expected, attr), getattr(got, attr)
            assert a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b), (attr, a, b, gene.exons)
        sizes.add(got.vertices.shape[1])
    assert max(sizes) >= 5 and 0 in sizes


def test_add_exon_matches_original():
    import random
    from spladder.classes.gene import Gene
    from . import reference
    rng = random.Random(3)
    for _ in range(300):
        a, b = Gene(name='G', start=0, stop=10, chr='chr1', strand='+'), Gene(name='G', start=0, stop=10, chr='chr1', strand='+')
        for _ in range(rng.randint(0, 12)):
            idx = min(len(a.exons) + rng.choice([0, 0, 0, 1]), len(a.exons)) if rng.random() < 0.3 else rng.randint(0, max(len(a.exons) - 1, 0))
            idx = idx if idx <= len(a.exons) else len(a.exons)
            exon = np.array([rng.randint(0, 100), rng.randint(100, 200)], dtype='int')
            a.add_exon(exon, idx)
            reference.gene_add_exon_reference(b, exon, idx)
        assert len(a.exons) == len(b.exons)
        for x, y in zip(a.exons, b.exons):
            assert x.dtype == y.dtype and x.shape == y.shape and np.array_equal(x, y)
