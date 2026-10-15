"""Editing the splice graphs of different contigs in different processes must give exactly the graphs of the serial run
and leave the settings the way the serial run leaves them."""

import pickle

import pytest

from spladder import spladder as spladder_cli
from spladder import settings
from spladder.core import gen_graphs
from spladder.core.spladdercore import spladder_core
from spladder.spladder_prep import prep_annotation

from . import synthetic
from .compare import compare_dirs, deep_equal

N_EXTRA = 8
### exons 1500-1700 of this gene are missing in the annotation: reads from them make the editing insert something
EXTRA_GENE = [(500, 800), (3000, 3300), (4500, 4900), (6000, 6400)]


@pytest.fixture
def many_contigs(tmp_path, monkeypatch):
    extra = ['chrX%d' % i for i in range(N_EXTRA)]
    monkeypatch.setattr(synthetic, 'CONTIGS', [('chrA', 20000), ('chrB', 8000), ('MT', 3000)] + [(c, 8000) for c in extra])
    genes = dict(synthetic.GENES)
    for c in extra:
        genes[c] = synthetic.GENES['chrB']
    monkeypatch.setattr(synthetic, 'GENES', genes)
    bam = synthetic.make_bam(tmp_path / 's.bam', seed=77, n_reads=14000, contig_choices=['chrA'] * 3 + ['chrB'] * 2 + ['MT'] + extra * 2)
    gtf = str(tmp_path / 'annotation.gtf')
    with open(gtf, 'w') as fd:
        for gene, contig, strand, exons in [('GA', 'chrA', '+', [(1000, 1300), (3200, 3500), (5000, 5400)]), ('GB', 'chrB', '-', EXTRA_GENE)] + \
                                           [('GX%d' % i, c, '-', EXTRA_GENE) for i, c in enumerate(extra)]:
            fd.write('%s\tsyn\tgene\t%d\t%d\t.\t%s\t.\tgene_id "%s"; gene_type "protein_coding"; gene_name "%s";\n' % (contig, exons[0][0] + 1, exons[-1][1], strand, gene, gene))
            for (a, b) in exons:
                fd.write('%s\tsyn\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s.t0"; gene_type "protein_coding"; gene_name "%s";\n' % (contig, a + 1, b, strand, gene, gene, gene))
    return bam, gtf


def build(gtf, bam, outdir, procs):
    spladder_cli.main(['spladder', 'build', '-a', gtf, '-o', str(outdir), '-b', bam, '--merge-strat', 'merge_graphs', '--extract-ase', '--readlen', '50',
                       '--parallel', str(procs), '-v'])


def graph_of_one_sample(gtf, bam, outdir, procs):
    """The graphs as the first step of the pipeline computes them, and the settings that step leaves behind."""
    (mode, options) = spladder_cli.parse_options(['spladder', 'build', '-a', gtf, '-o', str(outdir), '-b', bam, '--merge-strat', 'single', '--readlen', '50',
                                                  '--parallel', str(procs), '-v'])
    options = settings.parse_args(options)
    options = settings.set_confidence_level(options)
    options = prep_annotation(options)
    spladder_core(bam, str(outdir / 'graph.pickle'), options)
    return pickle.load(open(str(outdir / 'graph.pickle'), 'rb')), options


@pytest.mark.parametrize('procs', [2, 3])
def test_graphs_of_many_contigs_are_exact(many_contigs, tmp_path, capsys, monkeypatch, procs):
    ### more contigs than shards: several contigs are edited by the same process
    bam, gtf = many_contigs
    calls = []
    original = gen_graphs._augment_genes_in_shards
    monkeypatch.setattr(gen_graphs, '_augment_genes_in_shards', lambda *args: calls.append(1) or original(*args))
    build(gtf, bam, tmp_path / 'serial', 1)
    assert calls == []
    build(gtf, bam, tmp_path / 'sharded', procs)
    capsys.readouterr()
    assert calls == [1]
    assert compare_dirs(str(tmp_path / 'serial'), str(tmp_path / 'sharded')) == []


def test_settings_after_editing_in_shards_are_the_serial_ones(many_contigs, tmp_path, capsys, monkeypatch):
    bam, gtf = many_contigs
    (graphs_serial, inserted_serial), opts_serial = graph_of_one_sample(gtf, bam, tmp_path / 'serial', 1)
    calls = []
    original = gen_graphs._augment_genes_in_shards
    monkeypatch.setattr(gen_graphs, '_augment_genes_in_shards', lambda *args: calls.append(1) or original(*args))
    (graphs_sharded, inserted_sharded), opts_sharded = graph_of_one_sample(gtf, bam, tmp_path / 'sharded', 3)
    capsys.readouterr()
    assert calls == [1]
    ### the scenario is only meaningful, if something was inserted
    assert inserted_serial['intron_retention'] + inserted_serial['cassette_exon'] + inserted_serial['exon_skip'] + inserted_serial['alt_53_prime'] > 0
    assert deep_equal(inserted_serial, inserted_sharded) is None
    assert deep_equal(graphs_serial, graphs_sharded) is None
    ### what the next sample starts with
    assert opts_serial.read_filter == opts_sharded.read_filter
    assert (opts_serial.read_filter is opts_serial.intron_retention['read_filter']) == (opts_sharded.read_filter is opts_sharded.intron_retention['read_filter'])
    assert opts_serial.intron_retention['read_filter'] == opts_sharded.intron_retention['read_filter']
    for attr in ('insert_intron_iterations', 'debug'):
        assert getattr(opts_serial, attr, None) == getattr(opts_sharded, attr, None)
    assert opts_serial.chrm_lookup == opts_sharded.chrm_lookup
