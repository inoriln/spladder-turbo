#!/usr/bin/env python
"""Run SplAdder (upstream) and spladder-turbo on the same data and compare everything they write.

Example:

    python validation/verify_against_upstream.py --upstream-src ~/src/spladder-3.1.1 \\
        --annotation genes.gtf --bams a.bam,b.bam -- --parallel 4

Everything after `--` is given to `spladder build` of both programs, the same options for both. The tool needs

* a source tree of SplAdder 3.1.1 (the directory that contains the package `spladder`, e.g. a clone of
  https://github.com/ratschlab/spladder at the tag v3.1.1),
* a Python interpreter with the dependencies of SplAdder (default: the one that runs this script). It is not important which of
  the two programs is installed in it, if any: the program is chosen with PYTHONPATH.

Each program gets its own output directory and its own link to the annotation (SplAdder converts the annotation to a .pickle file
next to it). The outputs are compared with tests/compare.py: text files by content, HDF5 files by dataset, pickle files by object.
Use a small part of your data (for example one chromosome of a few samples): SplAdder is slow.

Exit code: 0 if the outputs are identical, 1 if they differ (or the programs end differently), 2 if both programs fail (which does not
say anything about the results, e.g. when an input is missing).
"""

import argparse
import os
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402


def split_arguments(argv):
    """(arguments of this tool, options for spladder build): they are separated by `--`."""

    if '--' in argv:
        i = argv.index('--')
        return argv[:i], argv[i + 1:]
    return argv, []


def parse_arguments(argv):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], epilog='Options for `spladder build` follow after `--`.')
    parser.add_argument('--upstream-src', required=True, metavar='DIR', help='source tree of SplAdder 3.1.1 (contains the package spladder)')
    parser.add_argument('--turbo-src', default=common.REPO, metavar='DIR', help='source tree of spladder-turbo [this repository]')
    parser.add_argument('--annotation', '-a', required=True, metavar='FILE', help='annotation (GTF or GFF3)')
    parser.add_argument('--bams', '-b', required=True, metavar='A.bam,B.bam', help='sorted, indexed alignment files, comma separated')
    parser.add_argument('--python', default=sys.executable, metavar='PYTHON', help='interpreter that has the dependencies [the current one]')
    parser.add_argument('--workdir', metavar='DIR', help='where the outputs and logs are written [a new temporary directory]')
    parser.add_argument('--concurrent', action='store_true', help='run the two programs at the same time (the timings are then less reliable)')
    parser.add_argument('--timeout', type=float, metavar='SECONDS', help='stop a program that runs longer')
    return parser.parse_args(argv)


def run_program(program, args, workdir, extra):
    directory = os.path.join(workdir, program.key)
    annotation = common.link_annotation(args.annotation, os.path.join(directory, 'annotation'))
    outdir = os.path.join(directory, 'out')
    argv = program.build_command(annotation, outdir, [os.path.abspath(b) for b in args.bams.split(',')], extra)
    print('%s: %s' % (program.name, ' '.join(argv)), flush=True)
    return program.run(argv, os.path.join(directory, 'log.txt'), cwd=directory, timeout=args.timeout)


def main(argv=None):
    argv, extra = split_arguments(sys.argv[1:] if argv is None else argv)
    args = parse_arguments(argv)
    common.check_source(args.upstream_src, 'upstream')
    common.check_source(args.turbo_src, 'spladder-turbo')
    for bam in args.bams.split(','):
        if not os.path.exists(bam):
            sys.exit('alignment file %s does not exist' % bam)
    workdir = os.path.abspath(args.workdir) if args.workdir else tempfile.mkdtemp(prefix='verify_')
    os.makedirs(workdir, exist_ok=True)

    upstream = common.Program('SplAdder', args.upstream_src, args.python, key='upstream')
    turbo = common.Program('spladder-turbo', args.turbo_src, args.python, key='turbo')
    if args.concurrent:
        with ThreadPoolExecutor(2) as pool:
            run_upstream, run_turbo = pool.map(lambda p: run_program(p, args, workdir, extra), (upstream, turbo))
    else:
        run_upstream = run_program(upstream, args, workdir, extra)
        run_turbo = run_program(turbo, args, workdir, extra)

    print()
    for program, run in ((upstream, run_upstream), (turbo, run_turbo)):
        print('%-16s exit code %d, %.1f s' % (program.name, run.returncode, run.seconds))
        if run.returncode != 0:
            print('  last lines of the log (%s):\n    %s' % (os.path.join(workdir, program.key, 'log.txt'), run.tail.replace('\n', '\n    ')))
    if run_upstream.returncode == 0 and run_turbo.returncode == 0 and run_turbo.seconds > 0:
        print('speed-up: %.1f x' % (run_upstream.seconds / run_turbo.seconds))

    diffs, n_files = [], 0
    if run_upstream.returncode == 0 and run_turbo.returncode == 0:
        diffs, n_files = common.compare_outputs(os.path.join(workdir, 'upstream', 'out'), os.path.join(workdir, 'turbo', 'out'))
    status, text = common.verdict(run_upstream, run_turbo, diffs, n_files)
    print('\nResult: %s\nOutputs and logs: %s' % (text, workdir))
    return {'identical': 0, 'both failed': 2}.get(status, 1)


if __name__ == '__main__':
    sys.exit(main())
