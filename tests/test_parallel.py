"""Processing samples in parallel must give exactly the results of the serial run."""

import os
import pickle
import shutil
import signal
import sys
import time
import types

import pytest

from spladder import helpers, parallel, settings, spladder as spladder_cli

from .compare import compare_dirs, deep_equal

DATA = os.path.join(os.path.dirname(__file__), 'testcase_events', 'data')


def build_args(outdir, gtf, n_samples, parallel_procs, extra=()):
    bams = ','.join(os.path.join(DATA, 'align', 'testcase_events_1_sample%d.bam' % (i + 1)) for i in range(n_samples))
    return ['spladder', 'build', '-a', gtf, '-o', str(outdir), '-b', bams, '--event-types', 'exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip',
            '--merge-strat', 'merge_graphs', '--extract-ase', '--readlen', '50', '--output-txt', '--output-gff3', '--output-conf-icgc', '--output-conf-tcga',
            '--parallel', str(parallel_procs), '-v'] + list(extra)


@pytest.fixture()
def gtf(tmp_path):
    ### a private copy: the parsed annotation is cached next to the GTF file
    dst = tmp_path / 'annotation.gtf'
    shutil.copy(os.path.join(DATA, 'testcase_events_spladder.gtf'), dst)
    return str(dst)


def run_build(outdir, gtf, n_samples, parallel_procs, extra=()):
    spladder_cli.main(build_args(outdir, gtf, n_samples, parallel_procs, extra))


def test_parallel_matches_serial(tmp_path, gtf, capsys):
    run_build(tmp_path / 'serial', gtf, 8, 1)
    run_build(tmp_path / 'parallel', gtf, 8, 4)
    assert compare_dirs(str(tmp_path / 'serial'), str(tmp_path / 'parallel')) == []


def test_parallel_matches_serial_when_the_read_filter_changes(tmp_path, gtf, monkeypatch, capsys):
    ### a low feasibility limit makes the intron filter tighten for some samples: in the serial run every sample
    ### starts with the filter state left by its predecessors, which the parallel run has to reproduce exactly
    monkeypatch.setattr(helpers, 'MAX_FEASIBLE_INTRONS', 5)
    run_build(tmp_path / 'serial', gtf, 8, 1)
    serial_log = capsys.readouterr().out
    assert serial_log.count('[feasibility]') > 0
    parallel.stats['reruns'] = 0
    run_build(tmp_path / 'parallel', gtf, 8, 4)
    assert compare_dirs(str(tmp_path / 'serial'), str(tmp_path / 'parallel')) == []
    ### the speculatively started samples had to be repeated at least once
    assert parallel.stats['reruns'] > 0


def test_parallel_resume_skips_existing_samples(tmp_path, gtf, capsys):
    run_build(tmp_path / 'serial', gtf, 6, 1)
    ### pre-compute two of the sample graphs, then let the parallel run do the rest
    run_build(tmp_path / 'resumed', gtf, 2, 1, extra=['--no-extract-ase', '--merge-strat', 'single', '--no-quantify-graph'])
    for f in os.listdir(str(tmp_path / 'resumed' / 'spladder')):
        if f.endswith('.count.hdf5'):
            os.remove(str(tmp_path / 'resumed' / 'spladder' / f))
    run_build(tmp_path / 'resumed', gtf, 6, 3)
    ### labels are file names, so the pre-computed graphs are the ones of the first two samples
    a = pickle.load(open(str(tmp_path / 'serial' / 'spladder' / 'genes_graph_conf3.testcase_events_1_sample3.pickle'), 'rb'))
    b = pickle.load(open(str(tmp_path / 'resumed' / 'spladder' / 'genes_graph_conf3.testcase_events_1_sample3.pickle'), 'rb'))
    assert deep_equal(a, b) is None


def test_failing_sample_aborts_like_the_serial_run(tmp_path, gtf, monkeypatch, capsys):
    from spladder.core import spladdercore
    real = spladdercore.spladder_core

    def failing(bam_fnames, out_fname, options):
        if 'sample3.' in str(bam_fnames):
            print('aborting sample 3')
            raise SystemExit(7)
        return real(bam_fnames, out_fname, options)

    import spladder.spladder_build as build
    monkeypatch.setattr(build, 'spladder_core', failing)
    for procs, name in [(1, 'serial'), (4, 'parallel')]:
        with pytest.raises(SystemExit) as e:
            run_build(tmp_path / name, gtf, 6, procs)
        assert e.value.code == 7
    done = lambda name: sorted(f for f in os.listdir(str(tmp_path / name / 'spladder')) if f.startswith('genes_graph'))
    assert done('serial') == done('parallel') == ['genes_graph_conf3.testcase_events_1_sample%d.pickle' % i for i in (1, 2)]
    assert not [f for f in os.listdir(str(tmp_path / 'parallel' / 'spladder')) if '.tmp' in f]


def test_read_filter_state_round_trip():
    options = types.SimpleNamespace(readlen=50, confidence=3, intron_retention=dict())
    options = settings.set_confidence_level(options)
    state = parallel.snapshot_read_filter_state(options)
    assert state[2] is True                                   # aliased right after set_confidence_level
    options.read_filter['exon_len'] += 4                      # in-place change as done by make_introns_feasible
    assert options.intron_retention['read_filter']['exon_len'] == options.read_filter['exon_len']
    parallel.apply_read_filter_state(options, state)
    assert options.read_filter['exon_len'] == state[0]['exon_len']
    assert options.intron_retention['read_filter'] is options.read_filter

    unaliased = ({'exon_len': 13, 'mismatch': 0, 'mincount': 2}, {'exon_len': 20, 'mismatch': 0, 'mincount': 2}, False)
    parallel.apply_read_filter_state(options, unaliased)
    assert options.read_filter is not options.intron_retention['read_filter']
    assert parallel.snapshot_read_filter_state(options) == unaliased


def test_event_verification_is_unchanged(tmp_path, gtf, monkeypatch, capsys):
    ### verify_wrapper reads the counts of the genes from the HDF5 file and collects the results of the samples
    from spladder.alt_splice import verify
    from . import reference
    monkeypatch.setattr(verify, 'verify_wrapper', reference.verify_wrapper_reference)
    run_build(tmp_path / 'reference', gtf, 6, 1)
    monkeypatch.undo()
    run_build(tmp_path / 'optimized', gtf, 6, 1)
    assert compare_dirs(str(tmp_path / 'reference'), str(tmp_path / 'optimized')) == []


def test_events_are_collected_while_the_graph_is_quantified(tmp_path, gtf, monkeypatch, capsys):
    ### the collection of the events overlaps with the quantification in a parallel run and gives the very same files
    run_build(tmp_path / 'serial', gtf, 4, 1)
    serial_log = capsys.readouterr().out
    jobs = []
    original = parallel.Background

    class Recording(original):
        def __init__(self, function, *args):
            jobs.append(function.__name__)
            super().__init__(function, *args)

    monkeypatch.setattr(parallel, 'Background', Recording)
    run_build(tmp_path / 'parallel', gtf, 4, 3)
    parallel_log = capsys.readouterr().out
    ### one job for each group of event types (exon skips, multiple exon skips, the rest)
    assert jobs == ['collect_events'] * 3
    assert compare_dirs(str(tmp_path / 'serial'), str(tmp_path / 'parallel')) == []
    ### what the job prints is shown after the output of the quantification, as in the serial run
    assert serial_log.index('Quantifying gene expression') < serial_log.index('saving exon skips')
    assert parallel_log.index('Quantifying gene expression') < parallel_log.index('saving exon skips')
    assert not parallel._BACKGROUND


def test_background_job_reports_output_and_errors(capsys):
    def noisy(a, b):
        print('hello', a + b)
        import sys
        sys.stderr.write('warning\n')

    job = parallel.Background(noisy, 1, 2)
    assert capsys.readouterr().out == ''            # nothing before join
    job.join()
    captured = capsys.readouterr()
    assert captured.out == 'hello 3\n' and captured.err == 'warning\n'

    def failing():
        raise ValueError('broken')

    job = parallel.Background(failing)
    with pytest.raises(RuntimeError, match='broken'):
        job.join()
    assert not parallel._BACKGROUND


def test_background_job_is_ended_when_the_main_job_fails(capsys):
    import time

    def slow():
        time.sleep(30)

    job = parallel.Background(slow)
    process = job._process
    parallel.abort_background()
    assert not process.is_alive()
    assert not parallel._BACKGROUND


def test_event_types_are_verified_at_the_same_time(tmp_path, gtf, monkeypatch, capsys):
    from spladder import spladder_build
    run_build(tmp_path / 'serial', gtf, 3, 1)
    calls = []
    original = parallel.run_shards

    def recording(function, shared_data, n_shards, workers, daemon=True):
        calls.append((function.__name__, n_shards, workers, daemon))
        return original(function, shared_data, n_shards, workers, daemon)

    monkeypatch.setattr(parallel, 'run_shards', recording)
    run_build(tmp_path / 'parallel', gtf, 3, 2)
    capsys.readouterr()
    ### the six event types, each of them with processes of its own (they start worker processes)
    assert ('_analyze_task', 6, 6, False) in calls
    assert compare_dirs(str(tmp_path / 'serial'), str(tmp_path / 'parallel')) == []
    ### nothing is left behind next to the results
    assert not [f for f in os.listdir(str(tmp_path / 'parallel') + '/spladder') if '.tmp' in f]


### -------------------------------------------------------------------------------------------------------------
### workers must not outlive the main process

_LIFECYCLE_SCRIPT = r'''
import os, sys, time
sys.path.insert(0, {repo!r})
from spladder import parallel

marker, flag = sys.argv[1], sys.argv[2]

def inner(x):
    time.sleep(120)

def task(x):
    ### a worker that has workers of its own (like a sample worker with the processes of its contigs)
    for _ in parallel.ordered_map(inner, [1, 2], 2, daemon=True):
        pass

with parallel.sigterm_unwinds():
    try:
        for _ in parallel.ordered_map(task, [1, 2], 2, daemon=False):
            pass
    finally:
        open(flag, 'w').write('cleaned up')
'''


def _processes_with(marker):
    import subprocess
    out = subprocess.run(['pgrep', '-f', marker], capture_output=True, text=True).stdout.split()
    return [int(p) for p in out]


@pytest.mark.skipif(not sys.platform.startswith('linux') or shutil.which('pgrep') is None, reason='needs Linux and pgrep')
@pytest.mark.parametrize('sig,cleans_up', [(signal.SIGTERM, True), (signal.SIGINT, True), (signal.SIGKILL, False)])
def test_workers_end_with_the_main_process(tmp_path, sig, cleans_up):
    import subprocess
    marker = 'lifecycle-%s' % tmp_path.name
    flag = str(tmp_path / 'flag')
    script = _LIFECYCLE_SCRIPT.format(repo=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ### (a new session, like a job of a scheduler, and with the default action for SIGINT, which a background job of a shell would not have)
    main = subprocess.Popen([sys.executable, '-c', script, marker, flag], start_new_session=True)
    try:
        deadline = time.time() + 30
        while len(_processes_with(marker)) < 1 + 2 + 4 and time.time() < deadline:
            time.sleep(0.2)
        assert len(_processes_with(marker)) == 1 + 2 + 4       # the main process, its two workers, their four workers
        main.send_signal(sig)
        deadline = time.time() + 15
        while _processes_with(marker) and time.time() < deadline:
            time.sleep(0.2)
        assert _processes_with(marker) == []
        main.wait(timeout=10)
        assert os.path.exists(flag) == cleans_up
    finally:
        for pid in _processes_with(marker):
            os.kill(pid, signal.SIGKILL)
        main.kill()
        main.wait()


def test_temporary_files_of_workers_that_were_ended_are_removed(tmp_path):
    ### a worker that is ended while it writes leaves the file dump_pickle was writing: <graph>.part<pid of the main process>.tmp<random>
    from types import SimpleNamespace
    from spladder.utils import dump_pickle
    options = SimpleNamespace(parallel=2, read_filter={'mismatch': 2}, intron_retention={'read_filter': None}, insert_ir=False, logfile='-',
                              samples=['a', 'b'], bam_fnames=['a.bam', 'b.bam'])
    graph_fname = lambda idx: str(tmp_path / ('g_%s.pickle' % options.samples[idx]))
    leftover = graph_fname(1) + '.part%i.tmp0123456789ab' % os.getpid()
    of_another_run = graph_fname(1) + '.part1.tmp0123456789ab'
    for f in (leftover, of_another_run):
        open(f, 'wb').write(b'x')
    parallel.build_sample_graphs(options, [0, 1], graph_fname, lambda bam, out, opts: dump_pickle(bam, out))
    assert sorted(os.listdir(str(tmp_path))) == ['g_a.pickle', 'g_b.pickle', 'g_b.pickle.part1.tmp0123456789ab']
