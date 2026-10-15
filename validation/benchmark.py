#!/usr/bin/env python
"""Time SplAdder (upstream) and spladder-turbo on the same data: wall-clock time, CPU time and the memory of the largest process.

Example:

    python validation/benchmark.py --upstream-src ~/src/spladder-3.1.1 --annotation genes.gtf --bams a.bam,b.bam \\
        --repeat 3 --discard-first -- --parallel 8

Everything after `--` is given to `spladder build` of both programs. The tool needs the same things as verify_against_upstream.py.

* Every run starts from nothing: a new output directory and (unless --reuse-annotation is given) the annotation is parsed again.
* CPU time is the user and system time of the program and of all the processes that it started. The memory is the resident
  set size of the largest of these processes (not the sum).
* The first run of spladder-turbo compiles numba functions, which takes about half a minute; --discard-first makes one extra
  run per program that is not reported.
* The linear algebra libraries are limited to one thread for both programs, which is the setting that is most favourable to SplAdder.
* The programs run one after the other. Do not run other work on the machine while measuring, and say on what machine
  you measured: the report starts with a description of the machine and the software.
"""

import argparse
import json
import os
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402


def machine_description(python):
    cpu = ''
    try:
        with open('/proc/cpuinfo') as fd:
            cpu = next((line.split(':', 1)[1].strip() for line in fd if line.startswith('model name')), '')
    except OSError:
        pass
    versions = subprocess.run([python, '-c', 'import sys, numpy, numba, pysam; print("Python %s, numpy %s, numba %s, pysam %s" % (sys.version.split()[0], numpy.__version__, numba.__version__, pysam.__version__))'],
                              capture_output=True, text=True).stdout.strip()
    return ['%s, %s' % (platform.system(), platform.release()), '%s, %d logical CPUs' % (cpu or platform.processor() or 'unknown processor', os.cpu_count()),
            versions or 'versions of the packages: unknown', time.strftime('measured on %Y-%m-%d')]


def timed_run(program, argv, log, cwd):
    """Runs argv; returns (exit code, wall-clock seconds, CPU seconds, largest resident set size in GB)."""

    start = time.time()
    with open(log, 'w') as fd:
        process = subprocess.Popen(argv, cwd=cwd, env=program.env(), stdout=fd, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(process.pid, 0)
    wall = time.time() - start
    return os.waitstatus_to_exitcode(status), wall, usage.ru_utime + usage.ru_stime, usage.ru_maxrss / 1024.0 ** 2


def parse_arguments(argv):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], epilog='Options for `spladder build` follow after `--`.')
    parser.add_argument('--upstream-src', metavar='DIR', help='source tree of SplAdder 3.1.1 (contains the package spladder); without it only spladder-turbo is timed')
    parser.add_argument('--turbo-src', default=common.REPO, metavar='DIR', help='source tree of spladder-turbo [this repository]')
    parser.add_argument('--annotation', '-a', required=True, metavar='FILE')
    parser.add_argument('--bams', '-b', required=True, metavar='A.bam,B.bam')
    parser.add_argument('--python', default=sys.executable, metavar='PYTHON', help='interpreter that has the dependencies [the current one]')
    parser.add_argument('--workdir', metavar='DIR', help='outputs, logs and results.json [a new temporary directory]')
    parser.add_argument('--repeat', type=int, default=1, metavar='N', help='runs per program [1]')
    parser.add_argument('--discard-first', action='store_true', help='one more run per program that is not reported (compilation, file cache)')
    parser.add_argument('--reuse-annotation', action='store_true', help='keep the converted annotation between the runs of a program')
    return parser.parse_args(argv)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    extra = []
    if '--' in argv:
        i = argv.index('--')
        argv, extra = argv[:i], argv[i + 1:]
    args = parse_arguments(argv)
    common.check_source(args.turbo_src, 'spladder-turbo')
    programs = []
    if args.upstream_src:
        common.check_source(args.upstream_src, 'upstream')
        programs.append(common.Program('SplAdder', args.upstream_src, args.python, key='upstream'))
    programs.append(common.Program('spladder-turbo', args.turbo_src, args.python, key='turbo'))
    bams = [os.path.abspath(b) for b in args.bams.split(',')]
    workdir = os.path.abspath(args.workdir) if args.workdir else tempfile.mkdtemp(prefix='benchmark_')
    os.makedirs(workdir, exist_ok=True)

    print('\n'.join(machine_description(args.python)))
    print('annotation: %s, alignment files: %d, options: %s\n' % (os.path.basename(args.annotation), len(bams), ' '.join(extra) or '(none)'))
    results = {}
    for program in programs:
        directory = os.path.join(workdir, program.key)
        annotation = common.link_annotation(args.annotation, os.path.join(directory, 'annotation'))
        runs = []
        for number in range(-1 if args.discard_first else 0, args.repeat):
            outdir = os.path.join(directory, 'out')
            shutil.rmtree(outdir, ignore_errors=True)
            if not args.reuse_annotation and os.path.exists(annotation + '.pickle'):
                os.remove(annotation + '.pickle')
            command = program.build_command(annotation, outdir, bams, extra)
            if number == (-1 if args.discard_first else 0):
                print('%s: %s' % (program.name, ' '.join(command)), flush=True)
            returncode, wall, cpu, memory = timed_run(program, command, os.path.join(directory, 'log.txt'), directory)
            label = 'discarded' if number < 0 else 'run %d' % (number + 1)
            print('  %-10s exit %d   wall-clock %9.1f s   CPU %9.1f s   largest process %.2f GB' % (label, returncode, wall, cpu, memory), flush=True)
            if returncode != 0:
                sys.exit('%s failed; see %s' % (program.name, os.path.join(directory, 'log.txt')))
            if number >= 0:
                runs.append({'wall': wall, 'cpu': cpu, 'memory_gb': memory})
        results[program.key] = runs

    print('\n| | wall-clock (median) | CPU time (median) | largest process |')
    print('|---|---|---|---|')
    summary = {}
    for program in programs:
        runs = results[program.key]
        summary[program.key] = {k: statistics.median(r[k] for r in runs) for k in ('wall', 'cpu')}
        summary[program.key]['memory_gb'] = max(r['memory_gb'] for r in runs)
        print('| %s | %.1f s (%.1f-%.1f) | %.1f s | %.2f GB |' % (program.name, summary[program.key]['wall'], min(r['wall'] for r in runs), max(r['wall'] for r in runs),
                                                                     summary[program.key]['cpu'], summary[program.key]['memory_gb']))
    if 'upstream' in summary:
        print('\nspeed-up: %.1f x wall-clock, %.1f x CPU time' % (summary['upstream']['wall'] / summary['turbo']['wall'], summary['upstream']['cpu'] / summary['turbo']['cpu']))
    with open(os.path.join(workdir, 'results.json'), 'w') as fd:
        json.dump({'machine': machine_description(args.python), 'options': extra, 'runs': results, 'summary': summary}, fd, indent=1)
    print('results: %s' % os.path.join(workdir, 'results.json'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
