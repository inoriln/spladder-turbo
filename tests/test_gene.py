"""Gene.label_alt must label exactly like the original implementation."""

import random

import numpy as np
import pytest

from spladder.classes.gene import Gene

from . import reference


def _random_gene(rng, n, coord_range, density, shared, symmetric=True):
    """Random splice graph with exons sorted by (start, stop); `shared` is the chance to reuse the start or end of an earlier exon."""

    starts, stops = [], []
    for _ in range(n):
        if starts and rng.random() < shared:
            j = rng.randrange(len(starts))
            s, e = starts[j], stops[j]
            if rng.random() < 0.5:
                e = s + rng.randint(1, coord_range // 3 + 1)
            else:
                s = max(e - rng.randint(1, coord_range // 3 + 1), 0)
        else:
            s = rng.randint(0, coord_range)
            e = s + rng.randint(1, coord_range // 3 + 1)
        starts.append(s)
        stops.append(e)
    order = sorted(range(n), key=lambda i: (starts[i], stops[i]))
    gene = Gene(name='G', start=0, stop=coord_range * 2, chr='chr1', strand='+')
    gene.splicegraph.vertices = np.array([[starts[i] for i in order], [stops[i] for i in order]], dtype='int').reshape(2, n)
    edges = np.triu(np.array([[1 if rng.random() < density else 0 for _ in range(n)] for _ in range(n)], dtype='int').reshape(n, n), 1)
    ### the original code only reads the upper triangle of the adjacency matrix
    gene.splicegraph.edges = edges + edges.T if symmetric else edges + np.tril(np.array([[1 if rng.random() < density else 0 for _ in range(n)] for _ in range(n)], dtype='int').reshape(n, n), -1)
    return gene


@pytest.mark.parametrize('symmetric', [True, False])
@pytest.mark.parametrize('density', [0.0, 0.15, 0.4, 0.8])
@pytest.mark.parametrize('shared', [0.0, 0.4, 0.8])
def test_label_alt_matches_original(density, shared, symmetric):
    rng = random.Random(int(density * 100) + int(shared * 10))
    seen = set()
    for n in list(range(0, 13)) * 50:
        coord_range = rng.choice([20, 60, 300])
        gene = _random_gene(rng, n, coord_range, density, shared, symmetric)
        ref = _random_gene(random.Random(0), 0, 1, 0, 0)
        ref.splicegraph.vertices = gene.splicegraph.vertices.copy()
        ref.splicegraph.edges = gene.splicegraph.edges.copy()
        reference.label_alt_reference(ref)
        gene.label_alt()
        assert (ref.is_alt, ref.is_alt_spliced) == (gene.is_alt, gene.is_alt_spliced), (gene.splicegraph.vertices, gene.splicegraph.edges)
        assert type(ref.is_alt) == type(gene.is_alt) and type(ref.is_alt_spliced) == type(gene.is_alt_spliced)
        seen.add((gene.is_alt, gene.is_alt_spliced))
    ### the random graphs cover the different outcomes
    assert len(seen) >= 2


def _outcome(function):
    try:
        return ('ok', function())
    except Exception as e:
        return ('error', type(e))


@pytest.mark.parametrize('symmetric', [True, False])
@pytest.mark.parametrize('density', [0.0, 0.2, 0.6])
def test_segmentgraph_matches_original(density, symmetric):
    from spladder.classes.segmentgraph import Segmentgraph
    rng = random.Random(int(density * 10) + int(symmetric))
    ok = errors = 0
    for n in list(range(0, 15)) * 40:
        coord_range = rng.choice([20, 60, 300])
        gene = _random_gene(rng, n, coord_range, density, rng.choice([0.0, 0.4]), symmetric)
        if rng.random() < 0.1 and n > 0:
            ### a zero length exon matches no segment
            k = rng.randrange(n)
            gene.splicegraph.vertices[1, k] = gene.splicegraph.vertices[0, k]
        ref = Segmentgraph()
        exp = _outcome(lambda: reference.segmentgraph_from_gene_reference(ref, gene))
        got = _outcome(lambda: Segmentgraph(gene))
        assert exp[0] == got[0], (exp, got, gene.splicegraph.vertices, gene.splicegraph.edges)
        if exp[0] == 'error':
            assert exp[1] == got[1]
            errors += 1
            continue
        ok += 1
        new = got[1]
        for attr in ('segments', 'seg_match', 'seg_edges'):
            a, b = getattr(ref, attr), getattr(new, attr)
            assert a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b), (attr, gene.splicegraph.vertices, gene.splicegraph.edges)
    assert ok > 100


def _assert_same_state(a, b):
    """Two genes in the same sparse or dense state hold the very same arrays (dtype, shape, values)."""
    assert sorted(a.__dict__) == sorted(b.__dict__)
    for key in a.__dict__:
        x, y = getattr(a, key), getattr(b, key)
        if isinstance(x, np.ndarray):
            assert isinstance(y, np.ndarray) and x.dtype == y.dtype and x.shape == y.shape and np.array_equal(x, y), key
        elif key == 'splicegraph':
            assert sorted(x.__dict__) == sorted(y.__dict__)
            for k2 in x.__dict__:
                v, w = getattr(x, k2), getattr(y, k2)
                if isinstance(v, np.ndarray):
                    assert isinstance(w, np.ndarray) and v.dtype == w.dtype and v.shape == w.shape and np.array_equal(v, w), k2
                else:
                    assert v is None and w is None or v == w
        elif key not in ('segmentgraph', 'exons', 'transcripts', 'introns_anno'):
            assert x == y, key


def _random_matrix(rng, n, dtype, density):
    m = np.array([[rng.randint(0, 9) if rng.random() < density else 0 for _ in range(n)] for _ in range(n)]).reshape(n, n)
    if dtype == 'bool':
        return m.astype('bool')
    if dtype == 'float64':
        return m.astype('float64') / 2
    return m.astype(dtype)


@pytest.mark.parametrize('dtype', ['int64', 'int32', 'uint8', 'bool', 'float64'])
@pytest.mark.parametrize('density', [0.0, 0.2, 1.0])
def test_sparse_conversion_matches_scipy(dtype, density):
    import copy
    rng = random.Random(len(dtype) + int(density * 10))
    for n in list(range(0, 14)) * 10:
        with_count = rng.random() < 0.6
        gene = Gene(name='G', start=0, stop=10, chr='chr1', strand='+')
        gene.splicegraph.vertices = np.zeros((2, n), dtype='int')
        gene.splicegraph.edges = _random_matrix(rng, n, dtype, density)
        if with_count:
            gene.edge_count = _random_matrix(rng, n, 'int64', density)
        exp, got = copy.deepcopy(gene), copy.deepcopy(gene)
        reference.to_sparse_reference(exp)
        got.to_sparse()
        _assert_same_state(exp, got)
        ### the restored dense matrices are the original ones, whichever implementation wrote the sparse form
        for sparse, restore in ((exp, 'new'), (got, 'reference')):
            back = copy.deepcopy(sparse)
            if restore == 'new':
                back.from_sparse()
            else:
                reference.from_sparse_reference(back)
            assert back.splicegraph.edges.dtype == gene.splicegraph.edges.dtype and np.array_equal(back.splicegraph.edges, gene.splicegraph.edges)
            if with_count:
                assert back.edge_count.dtype == gene.edge_count.dtype and np.array_equal(back.edge_count, gene.edge_count)
            _assert_same_state(back, gene)
        ### both implementations take no sparse matrices in the dense state again
        with pytest.raises(AttributeError):
            got.to_sparse()


def test_sparse_conversion_of_unusual_matrices_matches_scipy():
    """Matrices that are not plain two-dimensional arrays go through scipy and give what they gave before."""
    import copy
    import scipy.sparse as spsp
    for edges in (spsp.csr_matrix(np.eye(3, dtype='int')), spsp.coo_matrix(np.eye(3, dtype='int'))):
        gene = Gene(name='G', start=0, stop=10, chr='chr1', strand='+')
        gene.splicegraph.edges = edges
        exp, got = copy.deepcopy(gene), copy.deepcopy(gene)
        reference.to_sparse_reference(exp)
        got.to_sparse()
        _assert_same_state(exp, got)


def _adversarial_gene(rng):
    """Random graph with many coinciding coordinates, zero length and sometimes negative length exons, any pattern of edges, several dtypes."""
    n = rng.choice([0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 14])
    coord = rng.choice([4, 8, 15, 40, 200])
    starts = [rng.randint(0, coord) for _ in range(n)]
    low = -1 if rng.random() < 0.15 else 0
    stops = [s + rng.randint(low, coord // 3 + 1) for s in starts]
    if rng.random() < 0.5:
        order = sorted(range(n), key=lambda i: (starts[i], stops[i]))
        starts, stops = [starts[i] for i in order], [stops[i] for i in order]
    density = rng.choice([0.0, 0.05, 0.15, 0.3, 0.6, 1.0])
    values = rng.choice([[0, 1], [0, 1, 2]])
    edges = np.array([[rng.choice(values) if rng.random() < density else 0 for _ in range(n)] for _ in range(n)]).reshape(n, n)
    if rng.random() < 0.5:
        edges = np.triu(edges, 1)
        edges = edges + edges.T
    elif rng.random() < 0.3:
        edges = np.triu(edges, 1)
    if max(values) == 1 and rng.random() < 0.2:
        edges = edges.astype(bool)
    gene = Gene(name='G', start=0, stop=10, chr='chr1', strand='+')
    gene.splicegraph.vertices = np.array([starts, stops], dtype='int64').reshape(2, n)
    gene.splicegraph.edges = edges
    return gene


def _label_with(method, gene):
    ref = Gene(name='G', start=0, stop=10, chr='chr1', strand='+')
    ref.splicegraph.vertices = gene.splicegraph.vertices.copy()
    ref.splicegraph.edges = gene.splicegraph.edges.copy()
    return _outcome(lambda: (method(ref), ref.is_alt, ref.is_alt_spliced, type(ref.is_alt), type(ref.is_alt_spliced))[1:])


def test_label_alt_matches_original_on_adversarial_graphs(monkeypatch):
    """The compiled version, the general version and the version for graphs the compiled one cannot handle all label like the original."""
    rng = random.Random(1234)
    kinds = set()
    for _ in range(6000):
        gene = _adversarial_gene(rng)
        expected = _label_with(lambda g: reference.label_alt_reference(g), gene)
        kinds.add(expected)
        assert _label_with(Gene.label_alt, gene) == expected, (gene.splicegraph.vertices, gene.splicegraph.edges)
        assert _label_with(Gene._label_alt_general, gene) == expected, (gene.splicegraph.vertices, gene.splicegraph.edges)
        assert _label_with(Gene._label_alt_original, gene) == expected
    assert len(kinds) >= 4          # all outcomes occur, including failures
    ### the version for large graphs
    monkeypatch.setattr('spladder.classes.gene._MAX_COMPILED_EXONS', 2)
    rng = random.Random(99)
    for _ in range(1500):
        gene = _adversarial_gene(rng)
        assert _label_with(Gene.label_alt, gene) == _label_with(lambda g: reference.label_alt_reference(g), gene), (gene.splicegraph.vertices, gene.splicegraph.edges)
