"""The compiled test for intron retentions must decide exactly like the original loops, and the splice graphs must come out the same."""

import copy
import random

import numpy as np
import pysam
import pytest

from spladder import editgraph, evidence, settings
from spladder import spladder as spladder_cli
from spladder.classes.gene import Gene
from spladder.editgraph import insert_intron_retentions

from . import reference
from .compare import deep_equal

CONTIG_LEN = 400000


def random_genes_and_reads(rng, n_genes):
    """Genes with random exon structures and reads with a random depth and extent in each exon and intron."""
    genes = []
    reads = []
    origin = 1000
    for gi in range(n_genes):
        n_exons = rng.randint(2, 6)
        exons = []
        pos = origin
        for _ in range(n_exons):
            pos += rng.choice([rng.randint(20, 200), rng.randint(200, 2000)])
            ln = rng.choice([rng.randint(30, 120), rng.randint(120, 600)])
            exons.append((pos, pos + ln))
            pos += ln
        end = pos
        gene = Gene(name='G%d' % gi, start=exons[0][0], stop=exons[-1][1], chr='c1', strand=rng.choice(['+', '-']))
        gene.splicegraph.vertices = np.array([[e[0] for e in exons], [e[1] for e in exons]], dtype='int64')
        edges = np.zeros((n_exons, n_exons), dtype='int64')
        for k in range(n_exons - 1):
            edges[k, k + 1] = edges[k + 1, k] = 1
        for _ in range(rng.randint(0, 2)):
            k = rng.randint(0, n_exons - 2)
            l = rng.randint(k + 2, n_exons) if k + 2 < n_exons else None
            if l is not None and l < n_exons:
                edges[k, l] = edges[l, k] = 1
        gene.splicegraph.edges = edges
        gene.splicegraph.terminals = np.zeros((2, n_exons), dtype='int64')
        genes.append(gene)
        ### reads in the exons and the introns: the depth and the extent differ from region to region
        n_regions = n_exons
        base = rng.choice([8, 30, 60, 120])
        retaining = rng.random() < 0.6
        for r, (a, b) in enumerate([(e[0], e[1]) for e in exons] + [(exons[k][1], exons[k + 1][0]) for k in range(n_exons - 1)]):
            if r < n_regions:
                depth = int(base * rng.uniform(0.7, 1.3))
                extent = 1.0
            elif retaining:
                depth = int(base * rng.choice([0, 0.15, 0.3, 0.6, 0.9, 1.1, 1.6, 2.5]))
                extent = rng.choice([0.5, 0.85, 0.95, 1.0, 1.0, 1.0])
            else:
                depth = int(base * rng.choice([0, 0, 0.1]))
                extent = rng.choice([0.3, 1.0])
            lo = a + int((b - a) * (1 - extent) * rng.random())
            hi = min(b, lo + max(int((b - a) * extent), 1))
            for _ in range(depth * max(1, (hi - lo) // 50)):
                ln = rng.randint(20, 60)
                if hi - lo <= ln:
                    continue
                reads.append((rng.randint(lo, hi - ln), ln))
        origin = end + 3000
    return genes, reads


def write_bam(path, reads, rng):
    header = {'HD': {'VN': '1.4', 'SO': 'coordinate'}, 'SQ': [{'SN': 'c1', 'LN': CONTIG_LEN}]}
    reads = sorted(reads)
    with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
        for i, (pos, ln) in enumerate(reads):
            a = pysam.AlignedSegment(out.header)
            a.query_name = 'r%d' % i
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = 60
            a.flag = 0
            a.cigartuples = [(0, ln)]
            a.query_sequence = 'A' * ln
            a.set_tag('NM', rng.choice([0, 0, 0, 1]))
            out.write(a)
    pysam.index(str(path))
    return str(path)


def make_options(tmp_path, bam, confidence):
    gtf = tmp_path / 'a.gtf'
    gtf.write_text('')
    (mode, options) = spladder_cli.parse_options(['spladder', 'build', '-a', str(gtf), '-o', str(tmp_path / 'out'), '-b', bam, '--readlen', '50', '-c', str(confidence)])
    options = settings.parse_args(options)
    return settings.set_confidence_level(options)


@pytest.mark.parametrize('confidence', [0, 3])
@pytest.mark.parametrize('use_evidence', [False, True])
def test_retentions_match_original(tmp_path, monkeypatch, confidence, use_evidence):
    rng = random.Random(10 + confidence)
    genes, reads = random_genes_and_reads(rng, 200)
    bam = write_bam(tmp_path / 'r.bam', reads, rng)
    options = make_options(tmp_path, bam, confidence)
    ### the genes are sorted by contig, strand does not matter
    expected_genes = np.array(copy.deepcopy(genes), dtype='object')
    got_genes = np.array(copy.deepcopy(genes), dtype='object')

    def run(function, gs):
        opts = copy.deepcopy(options)
        if use_evidence:
            with evidence.evidence_scope() as scope:
                evidence.prefetch_contigs(scope, bam, ['c1'], opts.read_filter, opts)
                return function(gs, bam, opts)
        return function(gs, bam, opts)

    exp, exp_inserted = run(reference.insert_intron_retentions_reference, expected_genes)
    calls = []
    kernel = editgraph._retentions_of_gene
    monkeypatch.setattr(editgraph, '_retentions_of_gene', lambda *args: calls.append(1) or kernel(*args))
    got, got_inserted = run(insert_intron_retentions, got_genes)
    assert len(calls) >= 150            # the compiled test decided for (almost) all genes
    ### the test only means something, if retentions are inserted and some are not
    assert exp_inserted > 10
    assert any(e.splicegraph.vertices.shape[1] != g.splicegraph.vertices.shape[1] for e, g in zip(genes, exp))
    assert exp_inserted == got_inserted
    assert deep_equal(list(exp), list(got)) is None
