"""The optimized graph editing steps must give exactly the splice graphs of the original implementations.

A synthetic annotation lacks some of the exons the (synthetic) reads originate from, so that the
steps that insert exons / introns / retentions into the graph really have something to do.
"""

import pickle

import pytest

from spladder import spladder as spladder_cli
from spladder.core import gen_graphs

from . import reference, synthetic
from .compare import deep_equal

### transcripts of the annotation: some of the exons that are expressed in the synthetic reads are missing
ANNOTATION = [
    ('G1', 'chrA', '+', [[(1000, 1300), (3200, 3500), (5000, 5400)]]),                         # exon 2000-2200 is not annotated
    ('G2', 'chrA', '-', [[(7000, 7300), (8500, 8800), (9500, 9900)], [(7000, 7300), (9500, 9900)]]),
    ('G3', 'chrA', '+', [[(12000, 12400), (13500, 13700), (17000, 17300)]]),                   # exon 15000-15500 is not annotated
    ('G4', 'chrA', '+', [[(18000, 18100), (18130, 18230), (18270, 18400)]]),
    ('G5', 'chrB', '-', [[(500, 800), (3000, 3300), (4500, 4900), (6000, 6400)]]),             # exon 1500-1700 is not annotated
]


def write_gtf(path):
    with open(path, 'w') as fd:
        for gene, contig, strand, transcripts in ANNOTATION:
            start = min(e[0] for t in transcripts for e in t) + 1
            stop = max(e[1] for t in transcripts for e in t)
            fd.write('%s\tsyn\tgene\t%d\t%d\t.\t%s\t.\tgene_id "%s"; gene_type "protein_coding"; gene_name "%s";\n' % (contig, start, stop, strand, gene, gene))
            for ti, exons in enumerate(transcripts):
                for (a, b) in exons:
                    fd.write('%s\tsyn\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s.t%d"; gene_type "protein_coding"; gene_name "%s";\n' % (contig, a + 1, b, strand, gene, gene, ti, gene))


@pytest.fixture(scope='module')
def inputs(tmp_path_factory):
    d = tmp_path_factory.mktemp('stages')
    bam = synthetic.make_bam(d / 'syn.bam', seed=21, n_reads=9000)
    gtf = str(d / 'annotation.gtf')
    write_gtf(gtf)
    return bam, gtf


def build_graph(bam, gtf, outdir, **extra):
    args = ['spladder', 'build', '-a', gtf, '-o', str(outdir), '-b', bam, '--merge-strat', 'single', '--no-extract-ase', '--no-quantify-graph',
            '--readlen', '50', '-c', str(extra.get('confidence', 3)), '-v']
    spladder_cli.main(args)
    return pickle.load(open(str(outdir / 'spladder' / ('genes_graph_conf%d.syn.pickle' % extra.get('confidence', 3))), 'rb'))


@pytest.mark.parametrize('confidence', [3, 1])
def test_cassette_exon_insertion_is_unchanged(inputs, tmp_path, monkeypatch, capsys, confidence):
    bam, gtf = inputs
    monkeypatch.setattr(gen_graphs, 'insert_cassette_exons', reference.insert_cassette_exons_reference)
    expected = build_graph(bam, gtf, tmp_path / 'reference', confidence=confidence)
    monkeypatch.undo()
    got = build_graph(bam, gtf, tmp_path / 'optimized', confidence=confidence)
    capsys.readouterr()
    ### the scenario is only meaningful, if exons are really inserted
    assert expected[1]['cassette_exon'] > 0
    assert deep_equal(expected, got) is None


def build_multi(outdir, gtf, bams, parallel_procs):
    spladder_cli.main(['spladder', 'build', '-a', gtf, '-o', str(outdir), '-b', ','.join(bams), '--merge-strat', 'merge_graphs', '--extract-ase',
                       '--readlen', '50', '--parallel', str(parallel_procs), '-v'])


@pytest.mark.parametrize('n_samples,procs', [(1, 3), (2, 2), (2, 5)])
def test_parallel_prefetch_and_nested_parallelism_are_exact(tmp_path, capsys, n_samples, procs):
    ### several contigs: the evidence of the contigs is prefetched in parallel (one sample: all processes, several: the share of each sample)
    from .compare import compare_dirs
    bams = [synthetic.make_bam(tmp_path / ('s%d.bam' % i), seed=40 + i, n_reads=5000) for i in range(n_samples)]
    gtf = str(tmp_path / 'annotation.gtf')
    write_gtf(gtf)
    build_multi(tmp_path / 'serial', gtf, bams, 1)
    build_multi(tmp_path / 'parallel', gtf, bams, procs)
    capsys.readouterr()
    assert compare_dirs(str(tmp_path / 'serial'), str(tmp_path / 'parallel')) == []


def test_prefetch_covers_all_requests_of_the_pipeline(tmp_path, capsys):
    ### with enough processes the contigs are read in parallel in advance: no contig is read on demand in the main process any more
    from unittest import mock
    from spladder import evidence
    bam = synthetic.make_bam(tmp_path / 's.bam', seed=50, n_reads=5000)
    gtf = str(tmp_path / 'annotation.gtf')
    write_gtf(gtf)
    on_demand = []
    original = evidence.build_contig_evidence

    def spy(*args, **kwargs):
        on_demand.append(args[:2])
        return original(*args, **kwargs)

    with mock.patch.object(evidence, 'build_contig_evidence', spy):
        build_multi(tmp_path / 'out', gtf, [bam], 3)
    capsys.readouterr()
    assert on_demand == []
