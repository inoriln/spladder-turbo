"""Helpers shared by the validation tools.

The tools run two programs on the same input and compare what they write: SplAdder 3.1.1 ("upstream", a source tree of the original
package) and spladder-turbo (this repository). Both are run with the same Python interpreter, which needs the dependencies of
SplAdder (numpy, numba, pysam, h5py, ...). The program is chosen with PYTHONPATH, so it does not matter whether one of them, or
both, are installed in that environment.
"""

import os
import shutil
import subprocess
import sys
import time
from collections import namedtuple

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from tests.compare import compare_dirs  # noqa: E402

### the linear algebra libraries only multiply small matrices: threads burn CPU and distort timings
THREAD_LIMITS = {'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'NUMBA_NUM_THREADS': '1'}

Run = namedtuple('Run', ['returncode', 'seconds', 'tail'])


def check_source(path, what):
    """Stops with a message, if `path` is not a directory that contains the package spladder."""

    if not os.path.isfile(os.path.join(path, 'spladder', '__init__.py')):
        sys.exit('%s: "%s" is not a directory that contains the package spladder' % (what, path))


class Program:
    """SplAdder or spladder-turbo: a source tree (the directory that contains the package `spladder`) and an interpreter."""

    def __init__(self, name, src, python=None, key=None):
        self.name = name
        self.key = key or name.lower().replace(' ', '_')
        self.src = os.path.abspath(src)
        self.python = python or sys.executable

    def env(self, extra=None):
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONUNBUFFERED='1')
        env['PYTHONPATH'] = self.src + (os.pathsep + os.environ['PYTHONPATH'] if os.environ.get('PYTHONPATH') else '')
        env.update(THREAD_LIMITS)
        env.update(extra or {})
        return env

    def build_command(self, annotation, outdir, bams, extra=()):
        return [self.python, '-m', 'spladder.spladder', 'build', '-a', annotation, '-o', outdir, '-b', ','.join(bams)] + list(extra)

    def run(self, argv, log, cwd=None, timeout=None, env_extra=None):
        """Runs argv with this program on the path; the output (also the error output) goes to the file `log`."""

        start = time.time()
        with open(log, 'w') as fd:
            try:
                returncode = subprocess.run(argv, cwd=cwd, env=self.env(env_extra), stdout=fd, stderr=subprocess.STDOUT, timeout=timeout).returncode
            except subprocess.TimeoutExpired:
                returncode = -1
        return Run(returncode, time.time() - start, tail_of(log))


def tail_of(path, n=8):
    with open(path, errors='replace') as fd:
        lines = fd.read().replace('\r', '\n').splitlines()
    return '\n'.join(lines[-n:])


def link_annotation(annotation, directory):
    """The annotation as a file of its own in `directory` (a link, not a copy): SplAdder writes the converted annotation next to it."""

    os.makedirs(directory, exist_ok=True)
    target = os.path.join(directory, os.path.basename(annotation))
    for stale in (target, target + '.pickle'):
        if os.path.lexists(stale):
            os.remove(stale)
    try:
        os.symlink(os.path.abspath(annotation), target)
    except OSError:
        shutil.copy(annotation, target)
    return target


def write_gtf(path, genes):
    """Annotation for small synthetic data sets.

    genes: [(name, contig, strand, [transcript, ...])], a transcript is a list of (start, end) exons, 0-based and half open."""

    with open(path, 'w') as fd:
        for name, contig, strand, transcripts in genes:
            start = min(e[0] for t in transcripts for e in t) + 1
            stop = max(e[1] for t in transcripts for e in t)
            fd.write('%s\tsyn\tgene\t%d\t%d\t.\t%s\t.\tgene_id "%s"; gene_type "protein_coding"; gene_name "%s";\n' % (contig, start, stop, strand, name, name))
            for ti, exons in enumerate(transcripts):
                for (a, b) in exons:
                    fd.write('%s\tsyn\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s"; transcript_id "%s.t%d"; gene_type "protein_coding"; gene_name "%s";\n'
                             % (contig, a + 1, b, strand, name, name, ti, name))


def same_array(a, b):
    return a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b)


def same_evidence(a, b):
    """True if two results of build_contig_evidence are the same: values, dtypes and shapes of coverage and introns.

    (Other results, e.g. None or the text that stands for 'cannot be built', are the same if they are equal.)"""

    if a is None or b is None or isinstance(a, str) or isinstance(b, str):
        return a == b
    if len(a.cov) != len(b.cov) or set(a.introns) != set(b.introns):
        return False
    if not all(same_array(x, y) for pa, pb in zip(a.cov, b.cov) for x, y in zip(pa, pb)):
        return False
    for key in a.introns:
        if len(a.introns[key]) != len(b.introns[key]):
            return False
        if not all(same_array(x, y) for pa, pb in zip(a.introns[key], b.introns[key]) for x, y in zip(pa, pb)):
            return False
    return True


def compare_outputs(dir_a, dir_b):
    """(list of differences, number of files) of two output directories of `spladder build`."""

    count = sum(len(files) for root, _, files in os.walk(dir_a) if os.path.relpath(root, dir_a).split(os.sep)[0] != 'tmp')
    return compare_dirs(dir_a, dir_b), count


def verdict(upstream, turbo, diffs, n_files):
    """(status, text) for the two runs and the differences between their outputs.

    status: 'identical', 'both failed', 'exit codes differ' or 'different'."""

    if upstream.returncode != turbo.returncode:
        return 'exit codes differ', 'exit codes differ: SplAdder %d, spladder-turbo %d' % (upstream.returncode, turbo.returncode)
    if upstream.returncode != 0:
        return 'both failed', 'both programs failed with exit code %d (this does not show that their results are the same)' % upstream.returncode
    if diffs:
        return 'different', 'DIFFERENT: %s' % '; '.join(str(d)[:200] for d in diffs[:3])
    return 'identical', 'IDENTICAL (%d files)' % n_files
