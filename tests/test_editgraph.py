"""Vectorized helpers of the splice graph editing must reproduce the original loops exactly."""

import random

import numpy as np
import pytest

from spladder import editgraph
from spladder.utils import unique_rows


def reference_all_introns(introns, vertices, edges):
    """The original loop of insert_cassette_exons that adds the introns implied by the splice graph."""

    all_introns = introns[:, :2]
    for k in range(edges.shape[0]):
        for l in range(k+1, edges.shape[0]):
            if edges[k, l] == 1:
                all_introns = np.r_[all_introns, np.array([[vertices[1, k], vertices[0, l]]])]
    return unique_rows(all_introns)


def reference_candidate_pairs(all_introns, vertices):
    """The conditions of the original pair loop of insert_cassette_exons that precede the coverage test."""

    pairs = []
    for k in range(all_introns.shape[0]):
        for l in range(k + 1, all_introns.shape[0]):
            if all_introns[k, 1] >= all_introns[l, 0]:
                continue
            ### only take intron pair, if outer ends are supported by current exons
            if (not all_introns[k, 0] in vertices[1, :]) or (not all_introns[l, 1] in vertices[0, :]):
                continue
            curr_exon = [all_introns[k, 1], all_introns[l, 0]]
            ### do not allow curr_exon to overlap existing exon
            if np.sum((vertices[0, :] < curr_exon[1]) & (vertices[1, :] > curr_exon[0])) > 0:
                continue
            pairs.append((k, l))
    return pairs


def _random_case(rng, coord_range, dtype):
    nv = rng.randint(0, 9)
    starts = [rng.randint(0, coord_range) for _ in range(nv)]
    ### short exons leave room for new exons in between
    vertices = np.array([starts, [s + rng.randint(1, max(coord_range // 8, 2)) for s in starts]], dtype='int64').reshape(2, nv)
    edges = np.triu(np.array([[1 if rng.random() < 0.3 else 0 for _ in range(nv)] for _ in range(nv)], dtype='int').reshape(nv, nv), 1)
    edges = edges + edges.T
    ni = rng.randint(0, 16)
    base = np.array([[rng.randint(0, coord_range), 0] for _ in range(ni)], dtype='int64').reshape(ni, 2)
    base[:, 1] = base[:, 0] + np.array([rng.randint(1, coord_range) for _ in range(ni)], dtype='int64')
    ### one end of an intron often coincides with an exon boundary, the other one is a novel position
    for i in range(ni):
        if nv:
            if rng.random() < 0.5:
                base[i, 0] = vertices[1, rng.randrange(nv)]
            elif rng.random() < 0.5:
                base[i, 1] = vertices[0, rng.randrange(nv)]
    return base.astype(dtype), vertices, edges


@pytest.mark.parametrize('dtype', ['uint32', 'int64'])
@pytest.mark.parametrize('coord_range', [10, 30, 120])
def test_splicegraph_introns_match_original_loop(dtype, coord_range):
    rng = random.Random(coord_range)
    for _ in range(600):
        introns, vertices, edges = _random_case(rng, coord_range, dtype)
        exp = reference_all_introns(introns, vertices, edges)
        got = editgraph._all_introns(introns, vertices, edges)
        assert exp.dtype == got.dtype and exp.shape == got.shape and np.array_equal(exp, got), (introns, vertices, edges)


@pytest.mark.parametrize('dtype', ['uint32', 'int64'])
@pytest.mark.parametrize('coord_range', [10, 30, 120])
def test_cassette_candidate_pairs_match_original_loop(dtype, coord_range):
    rng = random.Random(coord_range + 1)
    n_pairs = 0
    for _ in range(600):
        introns, vertices, edges = _random_case(rng, coord_range, dtype)
        all_introns = reference_all_introns(introns, vertices, edges)
        exp = reference_candidate_pairs(all_introns, vertices)
        got = [(int(k), int(l)) for k, l in editgraph._cassette_candidate_pairs(all_introns, vertices)]
        assert exp == got, (all_introns, vertices)
        n_pairs += len(exp)
    assert n_pairs > 100    # the random cases do contain plenty of candidate pairs

