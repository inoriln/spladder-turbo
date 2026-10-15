#!/usr/bin/env python
"""Many combinations of options (and multi-step workflows): SplAdder and spladder-turbo must write identical outputs for each of them.

Every configuration is run with both programs on the same data and all outputs are compared. The 42 configurations are

* 20 sets of options for `spladder build`: confidence levels, merge strategies, sparse alignment files, filters, event types, the
  output formats, ...;
* the same sets again with `--parallel 4` (16 of them);
* 6 workflows of several commands: a cohort (one graph per sample, merging, quantification per sample and collected, event
  detection), the same with chunked merging and with `--parallel 3`, one sample, two samples, two samples with sparse files.

Example:

    python validation/config_matrix.py --upstream-src ~/src/spladder-3.1.1 --annotation genes.gtf --bams a.bam,b.bam,c.bam --jobs 4

Use a small part of your data (for example one chromosome of three samples): SplAdder is slow, and a configuration takes about as
long as one run of SplAdder. At least three alignment files are needed (the workflows use the first three). The times in the table
are only a rough guide when several configurations run at once (--jobs).

Exit code: 0 if all configurations are identical, 1 if one is not, 2 if there is no difference, but a configuration failed in both programs.
"""

import argparse
import os
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402

OUTPUT_FORMATS = ['--output-txt', '--output-gff3', '--output-struc', '--output-struc-conf', '--output-bed', '--output-conf-bed', '--output-conf-tcga', '--output-conf-icgc']

CONFIGS = {
    'default':      [],
    'conf0':        ['-c', '0'],
    'conf1':        ['-c', '1'],
    'conf2':        ['-c', '2'],
    'single':       ['--merge-strat', 'single'],
    'merge_bams':   ['--merge-strat', 'merge_bams'],
    'merge_all':    ['--merge-strat', 'merge_all'],
    'sparse':       ['--sparse-bam'],
    'parallel3':    ['--parallel', '3'],
    'ignore_mm':    ['--ignore-mismatches'],
    'all_align':    ['--no-primary-only'],
    'validate_sg':  ['--validate-sg'],
    'event_subset': ['--event-types', 'exon_skip,intron_retention'],
    'all_outputs':  OUTPUT_FORMATS,
    'no_insert_ir': ['--no-insert-ir'],
    'no_insert_es': ['--no-insert-es'],
    'no_insert_ni': ['--no-insert-ni'],
    'iter2_se':     ['-I', '2', '--remove-se'],
    'readlen30':    ['-n', '30'],
    'qmode_single': ['--qmode', 'single'],
}

### the same settings with several processes (contigs are edited and counted in parallel, events are collected while quantifying)
for _name in ['default', 'conf0', 'conf1', 'merge_bams', 'merge_all', 'sparse', 'ignore_mm', 'all_align', 'validate_sg', 'event_subset', 'all_outputs',
              'no_insert_ir', 'no_insert_es', 'no_insert_ni', 'iter2_se', 'readlen30']:
    CONFIGS[_name + '_p4'] = CONFIGS[_name] + ['--parallel', '4']

NO_QUANT = ['--no-extract-ase', '--no-quantify-graph']
SINGLE_GRAPHS = [([i], ['--merge-strat', 'single'] + NO_QUANT) for i in range(3)]
MERGE = [('all', ['--merge-strat', 'merge_graphs'] + NO_QUANT)]
QUANTIFY = [([i], ['--merge-strat', 'merge_graphs', '--no-extract-ase', '--quantify-graph', '--qmode', 'single']) for i in range(3)] + \
           [('all', ['--merge-strat', 'merge_graphs', '--no-extract-ase', '--quantify-graph', '--qmode', 'collect'])]
EVENTS = [('all', ['--merge-strat', 'merge_graphs', '--event-types', 'exon_skip,intron_retention'])]
CHUNKED_MERGE = [('all', ['--merge-strat', 'merge_graphs'] + NO_QUANT + ['--chunksize', '2', '--chunked-merge'] + levels)
                 for levels in (['1', '2', '0', '2'], ['1', '2', '2', '4'], ['2', '2', '0', '2'])]

### a workflow is a list of steps (samples, options): the samples are 'all' or a list of indices
WORKFLOWS = {
    'wf_cohort':         SINGLE_GRAPHS + MERGE + QUANTIFY + EVENTS,
    'wf_cohort_chunked': SINGLE_GRAPHS + CHUNKED_MERGE + QUANTIFY + EVENTS,
    'wf_cohort_par3':    [(samples, options + ['--parallel', '3']) for samples, options in SINGLE_GRAPHS + MERGE + QUANTIFY + EVENTS],
    'wf_one_sample':     [([0], ['--merge-strat', 'single'])],
    'wf_two_samples':    [([0, 1], ['--merge-strat', 'merge_graphs'])],
    'wf_two_sparse':     [([0, 1], ['--merge-strat', 'merge_graphs', '--sparse-bam', '--parallel', '2'])],
}


def steps_of(name):
    if name in WORKFLOWS:
        return WORKFLOWS[name]
    options = CONFIGS[name]
    return [('all', ([] if '--merge-strat' in options else ['--merge-strat', 'merge_graphs']) + options)]


def link_alignments(bams, directory):
    """Links to the alignment files and their indices in `directory`: SplAdder writes the sparse representation next to them."""

    os.makedirs(directory, exist_ok=True)
    links = []
    for bam in bams:
        base = os.path.splitext(bam)[0]
        for source in (bam, bam + '.bai', bam + '.csi', base + '.bai'):
            if os.path.exists(source):
                target = os.path.join(directory, os.path.basename(source))
                if not os.path.lexists(target):
                    os.symlink(source, target)
        links.append(os.path.join(directory, os.path.basename(bam)))
    return links


def run_workflow(program, name, args, workdir, bams):
    """Runs the steps of a configuration; returns (Run of the last step or the failed one, total seconds, output directory)."""

    work = os.path.join(workdir, name, program.key)
    outdir = os.path.join(work, 'out')
    annotation = common.link_annotation(args.annotation, os.path.join(work, 'annotation'))
    links = link_alignments(bams, os.path.join(work, 'alignments'))
    start = time.time()
    run = None
    for number, (samples, options) in enumerate(steps_of(name), 1):
        selected = links if samples == 'all' else [links[i] for i in samples]
        command = program.build_command(annotation, outdir, selected, ['-n', str(args.readlen), '--extract-ase', '-v'] + options)
        run = program.run(command, os.path.join(work, 'step%d.log' % number), cwd=work, timeout=args.timeout)
        if run.returncode != 0:
            break
    return run, time.time() - start, outdir


def run_configuration(task):
    name, args, workdir, bams, upstream, turbo = task
    (run_up, seconds_up, out_up), (run_tb, seconds_tb, out_tb) = (run_workflow(p, name, args, workdir, bams) for p in (upstream, turbo))
    diffs, n_files = [], 0
    if run_up.returncode == 0 and run_tb.returncode == 0:
        diffs, n_files = common.compare_outputs(out_up, out_tb)
    status, text = common.verdict(run_up, run_tb, diffs, n_files)
    return name, status, text, seconds_up, seconds_tb


def parse_arguments(argv):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--upstream-src', metavar='DIR', help='source tree of SplAdder 3.1.1 (contains the package spladder); required, except for --list')
    parser.add_argument('--turbo-src', default=common.REPO, metavar='DIR', help='source tree of spladder-turbo [this repository]')
    parser.add_argument('--annotation', '-a', metavar='FILE', help='annotation (GTF or GFF3); required, except for --list')
    parser.add_argument('--bams', '-b', metavar='A.bam,B.bam,C.bam', help='at least three sorted, indexed alignment files, comma separated; required, except for --list')
    parser.add_argument('--readlen', type=int, default=50, metavar='N', help='read length that is given to every run (-n) [50]')
    parser.add_argument('--python', default=sys.executable, metavar='PYTHON', help='interpreter that has the dependencies [the current one]')
    parser.add_argument('--configs', metavar='NAME,NAME', help='only these configurations (--list shows them) [all]')
    parser.add_argument('--list', action='store_true', help='show the configurations and exit')
    parser.add_argument('--jobs', type=int, default=2, metavar='N', help='configurations at the same time [2]')
    parser.add_argument('--timeout', type=float, metavar='SECONDS', help='stop a command that runs longer')
    parser.add_argument('--workdir', metavar='DIR', help='outputs and logs [a new temporary directory]')
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_arguments(sys.argv[1:] if argv is None else argv)
    available = list(CONFIGS) + list(WORKFLOWS)
    if args.list:
        print('\n'.join(available))
        return 0
    if not (args.upstream_src and args.annotation and args.bams):
        sys.exit('--upstream-src, --annotation and --bams are required')
    selected = args.configs.split(',') if args.configs else available
    unknown = [name for name in selected if name not in available]
    if unknown:
        sys.exit('unknown configuration(s): %s' % ', '.join(unknown))
    bams = [os.path.abspath(b) for b in args.bams.split(',')]
    if len(bams) < 3:
        sys.exit('at least three alignment files are needed')
    if len({os.path.basename(b) for b in bams}) < len(bams):
        sys.exit('the alignment files must have different names (they are linked into one directory)')
    missing = [b for b in bams if not os.path.exists(b)]
    if missing or not os.path.exists(args.annotation):
        sys.exit('does not exist: %s' % ', '.join(missing + ([] if os.path.exists(args.annotation) else [args.annotation])))
    common.check_source(args.upstream_src, 'upstream')
    common.check_source(args.turbo_src, 'spladder-turbo')
    workdir = os.path.abspath(args.workdir) if args.workdir else tempfile.mkdtemp(prefix='config_matrix_')
    upstream = common.Program('SplAdder', args.upstream_src, args.python, key='upstream')
    turbo = common.Program('spladder-turbo', args.turbo_src, args.python, key='turbo')

    tasks = [(name, args, workdir, bams, upstream, turbo) for name in selected]
    print('%-18s %-60s %10s %10s' % ('configuration', 'result', 'SplAdder', 'turbo'), flush=True)
    counts = {}
    with ThreadPoolExecutor(args.jobs) as pool:
        for name, status, text, seconds_up, seconds_tb in pool.map(run_configuration, tasks):
            counts[status] = counts.get(status, 0) + 1
            print('%-18s %-60s %9.0fs %9.0fs' % (name, text[:60], seconds_up, seconds_tb), flush=True)
    print('\n%d configurations: %s   (work directory %s)' % (len(tasks), ', '.join('%d %s' % (n, s) for s, n in sorted(counts.items())), workdir))
    if counts.get('different', 0) or counts.get('exit codes differ', 0):
        return 1
    return 2 if counts.get('both failed', 0) else 0


if __name__ == '__main__':
    sys.exit(main())
