"""Merging the splice graphs of several samples: the union of the genes, exons and edges, with the number of samples that hold an edge."""

import os
import pickle
import random
from types import SimpleNamespace

import numpy as np
import pytest

from spladder.classes.gene import Gene
from spladder.merge import merge_genes_by_splicegraph


def make_gene(name, exons, edge_pairs, position):
    """A gene whose splice graph has the given exons ((start, stop) pairs) and edges (pairs of positions in `exons`)."""
    gene = Gene(name=name, start=position, stop=position + 10000, chr='chr1', strand='+')
    n = len(exons)
    gene.splicegraph.vertices = np.array([[e[0] for e in exons], [e[1] for e in exons]], dtype='int')
    gene.splicegraph.edges = np.zeros((n, n), dtype='int')
    for a, b in edge_pairs:
        gene.splicegraph.edges[a, b] = gene.splicegraph.edges[b, a] = 1
    gene.splicegraph.terminals = np.zeros((2, n), dtype='int')
    gene.exons = [np.array(exons)]
    gene.transcripts = ['t']
    return gene


def coordinates_of_edges(gene):
    """{(exon, exon): number of samples} of the merged gene, by exon coordinates."""
    gene.from_sparse()
    v = gene.splicegraph.vertices
    counts = {}
    for a in range(v.shape[1]):
        for b in range(a + 1, v.shape[1]):
            if gene.splicegraph.edges[a, b]:
                counts[(tuple(v[:, a]), tuple(v[:, b]))] = int(gene.edge_count[a, b])
    return counts


def random_samples(rng, n_samples):
    """Per sample a few of the genes G0...G9, every gene with some of the exons of a common pool."""
    names = ['G%d' % i for i in range(10)]
    pool = {name: [(1000 * i + 100 * k, 1000 * i + 100 * k + 60) for k in range(6)] for i, name in enumerate(names)}
    samples = []
    for _ in range(n_samples):
        genes = []
        for i, name in enumerate(names):
            if rng.random() < 0.6:
                idx = sorted(rng.sample(range(6), rng.randint(2, 5)))
                exons = [pool[name][k] for k in idx]
                pairs = [(a, b) for a in range(len(exons)) for b in range(a + 1, len(exons)) if rng.random() < 0.5 or b == a + 1]
                genes.append((name, exons, pairs, 1000 * i))
        samples.append(genes)
    return samples


@pytest.mark.parametrize('n_samples', [2, 3, 4])
@pytest.mark.parametrize('seed', range(8))
def test_merge_of_samples_with_different_genes(tmp_path, seed, n_samples):
    rng = random.Random(seed * 10 + n_samples)
    samples = random_samples(rng, n_samples)
    files = []
    for si, genes in enumerate(samples):
        objs = np.array([make_gene(*g) for g in genes], dtype='object')
        for g in objs:
            g.to_sparse()
        fname = str(tmp_path / ('s%d.pickle' % si))
        pickle.dump((objs, {}), open(fname, 'wb'))
        files.append(fname)
    out = str(tmp_path / 'merged.pickle')
    merge_genes_by_splicegraph(SimpleNamespace(do_merge_all=False, outdir=str(tmp_path), confidence=3), merge_list=files, fn_out=out)
    merged, _ = pickle.load(open(out, 'rb'))

    ### what has to come out: every gene of any sample once, with the edges of all samples and the number of samples holding them
    expected = {}
    for genes in samples:
        for name, exons, pairs, position in genes:
            edges = expected.setdefault(name, {})
            for a, b in pairs:
                key = (exons[a], exons[b])
                edges[key] = edges.get(key, 0) + 1
    assert sorted(g.name for g in merged) == sorted(expected)
    for gene in merged:
        got = coordinates_of_edges(gene)
        assert got == expected[gene.name], gene.name
