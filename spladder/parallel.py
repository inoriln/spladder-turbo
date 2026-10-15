"""Processing of samples in parallel processes.

The graphs of the individual samples (and the counting of the graph on the individual samples) do not depend
on each other, with one exception: the serial implementation shares the read filter settings between the
samples.  If the intron filter has to be tightened for genes with too many introns (helpers.make_introns_feasible),
the tightened filter is also used for the samples that follow, and, because of the way the settings are
handed over between the processing steps, the settings of the first sample are not exactly those of the others.

To give exactly the results of the serial run in any case, the graphs are computed speculatively: every
sample starts from the settings it has in the serial run if no sample changes them (which is almost always the
case).  The results are then accepted in sample order and only used, if the settings a sample really started
with in the serial run are the ones it was started with.  In the rare case that a sample changes the settings,
the affected samples are computed again, in order, in the main process.
"""

import contextlib
import gc
import io
import multiprocessing as mp
import multiprocessing.connection
import os
import re
import signal
import sys
import threading
import traceback

import numpy as np

from .utils import frozen_heap

### number of speculatively computed samples that had to be computed again (for tests and diagnostics)
stats = {'reruns': 0}

### the work description inherited by the forked worker processes
_TASK = None


def snapshot_read_filter_state(options):
    """The state of the read filter settings: (filter values, values of the intron retention filter, whether both are the same object)."""

    rf = options.read_filter
    ir = options.intron_retention.get('read_filter', None)
    return (dict(rf), None if ir is None else dict(ir), ir is rf)


def apply_read_filter_state(options, state):
    rf_values, ir_values, aliased = state
    options.read_filter = dict(rf_values)
    if ir_values is not None:
        options.intron_retention['read_filter'] = options.read_filter if aliased else dict(ir_values)


def _next_state_without_changes(options, state):
    """State after a sample that does not change the filter values: only the intron retention step decouples the two filters."""

    rf_values, ir_values, aliased = state
    if ir_values is None:
        return state
    return (rf_values, rf_values if aliased else ir_values, aliased and not options.insert_ir)


def _fork_context():
    """The workers inherit the (large) data structures of the main process by fork; where this is not safe, everything runs serially."""

    if not sys.platform.startswith('linux'):
        return None
    try:
        return mp.get_context('fork')
    except ValueError:
        return None


def can_fork():
    return _fork_context() is not None


def can_run_parallel(options, n_tasks):
    return options.parallel > 1 and n_tasks > 1 and can_fork()


def _ignore_sigint():
    """An interrupt is handled by the main process of a pool, which terminates its workers."""

    signal.signal(signal.SIGINT, signal.SIG_IGN)


def _die_with_parent(parent_pid):
    """Linux: asks the kernel to send SIGTERM to this process when the process that started it ends, however it ends.

    Without this, the workers of a main process that is killed (or terminated by a job scheduler) keep running to the end
    of their tasks, and so do the processes they have started themselves: every worker does this, so it is a chain."""

    try:
        import ctypes
        libc = ctypes.CDLL(None, use_errno=True)
        libc.prctl(1, int(signal.SIGTERM), 0, 0, 0)         # PR_SET_PDEATHSIG
    except Exception:
        return
    if os.getppid() != parent_pid:
        ### the parent has ended before this was set up
        os._exit(143)


def _worker_start(parent_pid):
    """The first thing every forked worker does."""

    ### (the main process may handle SIGTERM itself: a worker stops at once)
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    _die_with_parent(parent_pid)


def _pool_worker_start(parent_pid):
    _ignore_sigint()
    _worker_start(parent_pid)


def _child(function, task, conn, parent_pid):
    _worker_start(parent_pid)
    try:
        conn.send(('ok', function(task)))
    except BaseException:
        conn.send(('error', traceback.format_exc()))
    finally:
        conn.close()


def ordered_map(function, tasks, workers, daemon=True):
    """Yields function(task) for all tasks in the order of the tasks, computed in up to `workers` forked processes.

    In contrast to multiprocessing.Pool, a worker that dies (e.g., because it was killed for running out of
    memory) raises an error instead of blocking forever.  Workers that start processes themselves must not be daemons."""

    ctx = _fork_context()
    pending = list(enumerate(tasks))
    running = {}
    done = {}
    next_idx = 0
    ### the garbage collector of a forked process would touch (and thereby copy) all the objects of the parent
    gc.freeze()
    try:
        while next_idx < len(tasks):
            while pending and len(running) < workers:
                i, task = pending.pop(0)
                receiver, sender = ctx.Pipe(duplex=False)
                process = ctx.Process(target=_child, args=(function, task, sender, os.getpid()), daemon=daemon)
                process.start()
                sender.close()
                running[i] = (process, receiver)
            if next_idx in done:
                yield done.pop(next_idx)
                next_idx += 1
                continue
            multiprocessing.connection.wait([r for _, r in running.values()] + [p.sentinel for p, _ in running.values()])
            for i, (process, receiver) in list(running.items()):
                if receiver.poll():
                    status, value = receiver.recv()
                    process.join()
                    del running[i]
                    if status != 'ok':
                        raise RuntimeError('error in worker process:\n%s' % value)
                    done[i] = value
                elif not process.is_alive():
                    raise RuntimeError('worker process died unexpectedly (exit code %s), probably because it ran out of memory' % process.exitcode)
    finally:
        gc.unfreeze()
        for process, receiver in running.values():
            if process.is_alive():
                process.terminate()
            process.join()
            receiver.close()


### the data inherited by the worker processes of forked_pool
_SHARED = None


def shared():
    """The `shared` data handed to forked_pool (to be used in the functions that run in the worker processes)."""

    return _SHARED


def forked_pool(workers, shared_data):
    """multiprocessing pool whose worker processes see `shared_data` (any Python objects) as it is while they are started.

    In contrast to the arguments of the tasks, nothing has to be pickled and sent to the workers: for large data (the
    splice graphs of all genes) this is what makes the difference.  Returns None, where processes cannot be forked, so
    that the caller has to send the data with the tasks.  Release the data with release_forked_pool afterwards."""

    global _SHARED

    ctx = _fork_context()
    if ctx is None:
        return None
    _SHARED = shared_data
    ### the garbage collector of a forked process would touch (and thereby copy) all the objects of the parent
    with frozen_heap():
        return ctx.Pool(processes=workers, initializer=_pool_worker_start, initargs=(os.getpid(),))


def release_forked_pool(pool):
    global _SHARED

    pool.terminate()
    pool.join()
    _SHARED = None


def balanced_shards(contigs, workers, per_worker=3):
    """Puts the contigs ((name, array of positions) pairs) together to shards of about the same size, the biggest first.

    Returns the sorted positions of each shard: about `per_worker` shards for every worker balance the load without making the shards too small."""

    n_shards = max(min(len(contigs), per_worker * workers), 1)
    shards = [[] for _ in range(n_shards)]
    load = [0] * n_shards
    for name, idx in sorted(contigs, key=lambda c: -c[1].shape[0]):
        k = load.index(min(load))
        shards[k].append(idx)
        load[k] += idx.shape[0]
    order = sorted((k for k in range(n_shards) if shards[k]), key=lambda k: -load[k])
    return [np.sort(np.concatenate(shards[k])) for k in order]


def _run_shard(task):
    function, number = task
    result = []
    outcome, out, err = _run_captured(lambda: result.append(function(number)))
    return (outcome, result[0] if result else None, out, err)


def run_shards(function, shared_data, n_shards, workers, daemon=True):
    """Yields function(number) for number in range(n_shards) in forked processes (up to `workers` at a time), in order of the numbers.

    The function finds `shared_data` with shared().  What the shards print is shown by the main process, errors are raised there.
    Shards that start processes of their own must not be daemons."""

    global _SHARED

    _SHARED = shared_data
    ### receiving the results creates very many objects: no need to look for garbage in between
    gc_enabled = gc.isenabled()
    gc.disable()
    try:
        for (outcome, value, out, err) in ordered_map(_run_shard, [(function, number) for number in range(n_shards)], workers, daemon=daemon):
            _replay(outcome, out, err)
            yield value
    finally:
        _SHARED = None
        if gc_enabled:
            gc.enable()


def _run_captured(function, *args):
    """Runs function(*args) with the console output captured; returns (outcome, stdout, stderr)."""

    out, err = io.StringIO(), io.StringIO()
    outcome = None
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            function(*args)
        except SystemExit as e:
            outcome = ('exit', e.code)
        except Exception:
            outcome = ('error', traceback.format_exc())
    return outcome, out.getvalue(), err.getvalue()


def _replay(outcome, out, err):
    """Writes the captured output and raises what was raised in the worker."""

    sys.stdout.write(out)
    sys.stderr.write(err)
    sys.stdout.flush()
    if outcome is not None:
        if outcome[0] == 'exit':
            sys.exit(outcome[1])
        raise RuntimeError('error in worker process:\n%s' % outcome[1])


def _graph_worker(task):
    idx, state_in = task
    options, spladder_core, graph_fname, tag, inner = _TASK
    ### the processes that are not needed for the samples read the contigs of a sample in parallel
    options.parallel = inner
    if options.logfile != '-':
        options.logfile = '%s.%s' % (options.logfile, options.samples[idx])
    apply_read_filter_state(options, state_in)
    tmp_fname = '%s%s' % (graph_fname(idx), tag)
    outcome, out, err = _run_captured(spladder_core, options.bam_fnames[idx], tmp_fname, options)
    return (idx, state_in, snapshot_read_filter_state(options), tmp_fname, outcome, out, err)


def build_sample_graphs(options, todo, graph_fname, spladder_core):
    """Computes the graphs of the samples in `todo` in parallel and writes them to graph_fname(idx), exactly as the serial loop would."""

    global _TASK

    state0 = snapshot_read_filter_state(options)
    guess = _next_state_without_changes(options, state0)
    tasks = [(idx, state0 if i == 0 else guess) for i, idx in enumerate(todo)]

    tag = '.part%i' % os.getpid()
    workers = min(options.parallel, len(todo))
    inner = max(options.parallel // workers, 1)
    _TASK = (options, spladder_core, graph_fname, tag, inner)
    expected = state0
    try:
        for (idx, state_in, state_out, tmp_fname, outcome, out, err) in ordered_map(_graph_worker, tasks, workers, daemon=(inner == 1)):
            final_fname = graph_fname(idx)
            if state_in != expected:
                ### the speculation was wrong for this sample: compute it as the serial run does
                stats['reruns'] += 1
                if os.path.exists(tmp_fname):
                    os.remove(tmp_fname)
                apply_read_filter_state(options, expected)
                spladder_core(options.bam_fnames[idx], final_fname, options)
                expected = snapshot_read_filter_state(options)
                continue
            if outcome is None:
                os.replace(tmp_fname, final_fname)
                out = out.replace(tmp_fname, final_fname)
            elif os.path.exists(tmp_fname):
                os.remove(tmp_fname)
            _replay(outcome, out, err)
            expected = state_out
    finally:
        _TASK = None
        ### leave no temporary files of samples that were not used (nor the file a worker was writing when it was ended: dump_pickle's .tmp<random>)
        directory = os.path.dirname(graph_fname(todo[0]))
        temporary = re.compile(re.escape(tag) + r'(\.tmp[0-9a-f]{12})?$')
        for f in os.listdir(directory):
            if temporary.search(f):
                os.remove(os.path.join(directory, f))
    apply_read_filter_state(options, expected)


### the jobs started by Background that have not been joined yet
_BACKGROUND = []


def _background_main(function, args, conn, parent_pid):
    _worker_start(parent_pid)
    outcome, out, err = _run_captured(function, *args)
    try:
        conn.send((outcome, out, err))
    finally:
        conn.close()


class Background:
    """Runs function(*args) in a forked process while the caller does something else; join() waits for it.

    What the function prints is shown by join(), so that the output of the caller and of the job do not mix; an
    error is raised there.  The job may use processes of its own."""

    def __init__(self, function, *args):
        ctx = _fork_context()
        self._receiver, sender = ctx.Pipe(duplex=False)
        ### no daemon: the job starts processes itself
        self._process = ctx.Process(target=_background_main, args=(function, args, sender, os.getpid()), daemon=False)
        self._process.start()
        sender.close()
        _BACKGROUND.append(self)

    def join(self):
        try:
            outcome, out, err = self._receiver.recv()
        except EOFError:
            self.abort()
            raise RuntimeError('background process died unexpectedly (exit code %s), probably because it ran out of memory' % self._process.exitcode)
        self._process.join()
        self._close()
        _replay(outcome, out, err)

    def abort(self):
        if self._process.is_alive():
            self._process.terminate()
        self._process.join()
        self._close()

    def _close(self):
        self._receiver.close()
        if self in _BACKGROUND:
            _BACKGROUND.remove(self)


def _raise_on_sigterm(signum, frame):
    raise SystemExit(128 + signum)


@contextlib.contextmanager
def sigterm_unwinds():
    """While active, SIGTERM (e.g., from a job scheduler) raises SystemExit in the main thread, so that the workers are ended and the
    temporary files are removed; with the default action the process would be gone with all of that left behind.  A handler of the
    caller is left alone."""

    previous = None
    try:
        if threading.current_thread() is threading.main_thread() and signal.getsignal(signal.SIGTERM) == signal.SIG_DFL:
            previous = signal.signal(signal.SIGTERM, _raise_on_sigterm)
    except ValueError:
        previous = None
    try:
        yield
    finally:
        if previous is not None:
            signal.signal(signal.SIGTERM, previous)


def abort_background():
    """Ends the jobs that are still running (when the caller does not need their results any more, e.g., after an error)."""

    for job in list(_BACKGROUND):
        job.abort()


def _count_worker(s_idx):
    genes, bam_fnames, options, count_one, inner = _TASK
    options.parallel = inner
    result = []
    outcome, out, err = _run_captured(lambda: result.append(count_one(genes, bam_fnames[s_idx], options)))
    return (result[0] if result else None, outcome, out, err)


def count_samples(genes, bam_fnames, options, count_one):
    """count_one(genes, bam_fname, options) for all samples in parallel; the results are returned in sample order."""

    global _TASK

    workers = min(options.parallel, len(bam_fnames))
    inner = max(options.parallel // workers, 1)
    _TASK = (genes, bam_fnames, options, count_one, inner)
    results = []
    ### receiving the counts creates very many objects: no need to look for garbage in between
    gc_enabled = gc.isenabled()
    gc.disable()
    try:
        for (counts, outcome, out, err) in ordered_map(_count_worker, range(len(bam_fnames)), workers, daemon=(inner == 1)):
            _replay(outcome, out, err)
            results.append(counts)
    finally:
        _TASK = None
        if gc_enabled:
            gc.enable()
    return results
