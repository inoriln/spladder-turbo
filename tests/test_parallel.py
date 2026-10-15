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
