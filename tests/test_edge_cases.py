"""Unusual but valid inputs and situations: the pipeline must finish (or fail the way it should), not crash on a coding error."""

import os
import random
import shutil
import signal
import subprocess
import sys
import time

import pysam
import pytest

from spladder import spladder as spladder_cli

CONTIGS = [('c1', 30000), ('c2', 30000)]
G1 = ('G1', 'c1', '+', [[(1000, 1300), (3000, 3300), (5000, 5400)], [(1000, 1300), (5000, 5400)]])
G2 = ('G2', 'c1', '-', [[(8000, 8300), (9500, 9800), (11000, 11400)]])


def write_gtf(path, genes):
    with open(str(path), 'w') as fd:
        for name, contig, strand, transcripts in genes:
            start = min(e[0] for t in transcripts for e in t) + 1
            stop = max(e[1] for t in transcripts for e in t)
            fd.write('%s\tsyn\tgene\t%d\t%d\t.\t%s\t.\tgene_id "%s"; gene_type "protein_coding"; gene_name "%s";\n' % (contig, start, stop, strand, name, name))
            for ti, exons in enumerate(transcripts):
                for (a, b) in exons:
                    fd.write('%s\tsyn\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s.t%d"; gene_type "protein_coding"; gene_name "%s";\n' % (contig, a + 1, b, strand, name, name, ti, name))


def write_bam(path, reads, contigs=CONTIGS):
    header = {'HD': {'VN': '1.4', 'SO': 'coordinate'}, 'SQ': [{'SN': n, 'LN': l} for n, l in contigs]}
    tid = {n: i for i, (n, _) in enumerate(contigs)}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
        for i, (contig, pos, cigar, strand) in enumerate(sorted(reads, key=lambda r: (tid[r[0]], r[1]))):
            a = pysam.AlignedSegment(out.header)
            a.query_name = 'r%d' % i
            a.reference_id = tid[contig]
            a.reference_start = pos
            a.mapping_quality = 60
            a.flag = 0
            a.cigartuples = cigar
            a.query_sequence = 'A' * sum(l for o, l in cigar if o in (0, 1))
            a.set_tag('NM', 0)
            if any(o == 3 for o, _ in cigar):
                a.set_tag('XS', strand, value_type='A')
            out.write(a)
    pysam.index(str(path))
    return str(path)


def spliced_reads(contig, exons, strand, n, rng):
    reads = []
    for _ in range(n):
        j = rng.randrange(len(exons) - 1)
        left, right = rng.randint(10, 40), rng.randint(10, 40)
        reads.append((contig, exons[j][1] - left, [(0, left), (3, exons[j + 1][0] - exons[j][1]), (0, right)], strand))
    return reads


def build(tmp_path, genes, reads, parallel_procs=1, n_samples=2, extra=()):
    gtf = tmp_path / 'a.gtf'
    write_gtf(gtf, genes)
    bams = [write_bam(tmp_path / ('s%d.bam' % i), reads) for i in range(n_samples)]
    spladder_cli.main(['spladder', 'build', '-a', str(gtf), '-o', str(tmp_path / 'out'), '-b', ','.join(bams), '--merge-strat', 'merge_graphs', '--extract-ase',
                       '--readlen', '50', '--parallel', str(parallel_procs), '-v'] + list(extra))
    return tmp_path / 'out'


def test_annotation_without_any_edge(tmp_path, capsys):
    ### only single exon genes: no segment graph has an edge (writing the counts used to fail)
    genes = [('S1', 'c1', '+', [[(1000, 1300)]]), ('S2', 'c2', '-', [[(2000, 2500)]])]
    reads = [('c1', 1010, [(0, 50)], '+')] * 30 + [('c2', 2010, [(0, 50)], '-')] * 30
    out = build(tmp_path, genes, reads)
    capsys.readouterr()
    assert os.path.exists(str(out / 'spladder' / 'genes_graph_conf3.merge_graphs.count.hdf5'))


def test_gene_on_an_unknown_strand_is_read_as_plus_strand(tmp_path, capsys):
    ### (the warning about it raised an error)
    genes = [G1, ('GD', 'c2', '.', [[(14000, 14200), (15000, 15200)]])]
    rng = random.Random(1)
    reads = spliced_reads('c1', [(1000, 1300), (3000, 3300), (5000, 5400)], '+', 100, rng)
    with pytest.warns(UserWarning, match='strand of gene was provided as \\. - automatically set to \\+'):
        build(tmp_path, genes, reads)
    capsys.readouterr()


def test_verbose_run_with_processes_and_nothing_to_detect(tmp_path, capsys):
    ### no gene is alternatively spliced: the progress of the event detection had nothing to show (division by zero)
    genes = [('GI', 'c1', '+', [[(1000, 1300), (3000, 3300)], [(1000, 1300), (3000, 3300)]])]
    reads = spliced_reads('c1', [(1000, 1300), (3000, 3300)], '+', 100, random.Random(2))
    out = build(tmp_path, genes, reads, parallel_procs=2)
    capsys.readouterr()
    assert os.path.exists(str(out / 'merge_graphs_exon_skip_C3.pickle'))


def test_progress_of_nothing_is_not_shown(capsys):
    from spladder.helpers import log_progress
    log_progress(0, 0)
    log_progress(4, 0)
    assert capsys.readouterr().out == ''
    log_progress(3, 10)           # the first call does not have to be the one for the first item
    assert '4 / 10' in capsys.readouterr().out


def test_empty_alignment_file(tmp_path, capsys):
    out = build(tmp_path, [G1, G2], [], parallel_procs=2)
    capsys.readouterr()
    assert os.path.exists(str(out / 'spladder' / 'genes_graph_conf3.merge_graphs.pickle'))


def test_graph_plot_with_a_highlighted_range():
    ### (the highlight used a module that was not imported)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from spladder.viz.graph import plot_graph
    vertices = np.array([[100, 300, 500], [150, 350, 560]])
    edges = np.array([[0, 1, 1], [0, 0, 1], [0, 0, 0]])
    fig, ax = plt.subplots()
    try:
        plot_graph(vertices, edges, ax, highlight=(120, 400))
        assert len(ax.patches) == 3 + 1                    # the exons and the highlighted range
    finally:
        plt.close(fig)


def test_event_ids_of_an_unknown_event_type_are_refused_with_an_error():
    ### (the error class used for that did not exist)
    import numpy as np
    from types import SimpleNamespace
    from spladder.alt_splice import quantify
    IN = {'gene_idx': np.array([0]), 'gene_chr': np.array(['chr1']), 'event_pos': np.arange(8).reshape(1, 8)}
    options = SimpleNamespace(verbose=False)
    ids = quantify.get_event_ids(IN, 'exon_skip', np.array([0]), options)
    assert len(ids) == 2 and ids[0].shape == (1,)
    with pytest.raises(ValueError, match='not known'):
        quantify.get_event_ids(IN, 'no_such_event', np.array([0]), options)


def test_size_factors_by_upper_quartile():
    ### (the percentile came from scipy, which is not imported any more)
    import numpy as np
    from types import SimpleNamespace
    from spladder.count import get_size_factors
    counts = np.array([[1., 10.], [2., 20.], [3., 30.], [4., 40.], [5., 50.]])
    factors = get_size_factors(counts, SimpleNamespace(verbose=False), kind='uq')
    assert np.allclose(factors, [4 + 1, 40 + 1])


def test_gene_information_of_the_events_from_a_given_graph_file(tmp_path):
    ### (the option spladderfile: the file was opened in text mode, the mode went to pickle.load)
    import pickle
    from types import SimpleNamespace
    import h5py
    import numpy as np
    from spladder.alt_splice import analyze
    genes = np.array([SimpleNamespace(name='G1', chr='c1', strand='+', start=10, stop=90), SimpleNamespace(name='G2', chr='c2', strand='-', start=5, stop=50)], dtype='object')
    fname = str(tmp_path / 'graph.pickle')
    with open(fname, 'wb') as fd:
        pickle.dump((genes, {}), fd)
    options = SimpleNamespace(spladderfile=fname, samples=np.array(['s1', 's2']))
    with h5py.File(str(tmp_path / 'out.hdf5'), 'w') as out:
        analyze._prepare_count_hdf5(options, out, ['e1_cov', 'e2_cov'])
        assert [x.decode() if isinstance(x, bytes) else x for x in out['gene_names'][:]] == ['G1', 'G2']
        assert out['gene_pos'][:].tolist() == [[10, 90], [5, 50]]


def test_preparing_sparse_alignment_files(tmp_path):
    ### (the preparation used pysam and two helpers it did not import)
    from types import SimpleNamespace
    from spladder import spladder_prep
    from . import synthetic
    bam = synthetic.make_bam(tmp_path / 'syn.bam', seed=5, n_reads=500)
    options = SimpleNamespace(bam_fnames=[bam], parallel=1, verbose=False, primary_only=False, var_aware=False, mm_tag='NM', confidence=3,
                              read_filter={'intron': 20000, 'exon_len': 8, 'mismatch': 1, 'mincount': 1}, tmpdir=str(tmp_path), cram_ref=None)
    spladder_prep.prep_sparse_bam_full(options)
    assert os.path.exists(str(tmp_path / 'syn.hdf5'))
    del options.chrm_lookup
    spladder_prep.prep_sparse_bam_filtered(options)
    assert os.path.exists(str(tmp_path / 'syn.conf_3.filt.hdf5'))
