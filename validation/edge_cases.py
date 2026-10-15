#!/usr/bin/env python
"""Unusual inputs: SplAdder and spladder-turbo must end in the same way.

14 small data sets (an empty alignment file, a gene on an unknown strand, an annotation without any splice junction, an alignment
file without index, ...) are run through both programs in three ways: serial, `--parallel 3` and `--parallel 3 -c 1`. The result
of each run is one of

* identical       both programs finish and write identical outputs;
* same failure    both programs stop with the same error;
* expected        SplAdder stops with a coding error that spladder-turbo has fixed (see FIXED_DEFECTS), spladder-turbo finishes;
* PROBLEM         anything else: different outputs, different exit codes, a different error.

Example:

    python validation/edge_cases.py --upstream-src ~/src/spladder-3.1.1 --jobs 8

Exit code: 0 if there is no PROBLEM, 1 otherwise.
"""

import argparse
import os
import random
import sys
import tempfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402

BUILD_OPTIONS = ['--merge-strat', 'merge_graphs', '--extract-ase', '--readlen', '50', '-v']
MODES = [('serial', []), ('parallel', ['--parallel', '3']), ('parallel -c 1', ['--parallel', '3', '-c', '1'])]

### SplAdder 3.1.1 stops on these inputs because of coding errors, which spladder-turbo has fixed (docs/PERFORMANCE_CHANGES.md, section 9)
### case, runs, defect, exception of SplAdder
FIXED_DEFECTS = [('unstranded_dot_gene', 'all', 'B7', 'TypeError'),            # a gene on strand '.'
                 ('single_exon_genes', 'all', 'B8', 'AttributeError'),         # an annotation without any edge
                 ('identical_transcripts', 'parallel', 'B10', 'ZeroDivisionError'),    # nothing to verify, with -v and --parallel
                 ('many_small_contigs', 'parallel', 'B10', 'ZeroDivisionError')]


def write_bam(path, contigs, reads, with_nm=True, index=True, xs=True):
    """reads: (contig, position, CIGAR, strand)"""

    header = {'HD': {'VN': '1.4', 'SO': 'coordinate'}, 'SQ': [{'SN': n, 'LN': l} for n, l in contigs]}
    tid = {n: i for i, (n, _) in enumerate(contigs)}
    with pysam.AlignmentFile(path, 'wb', header=header) as out:
        for i, (contig, pos, cigar, strand) in enumerate(sorted(reads, key=lambda r: (tid[r[0]], r[1]))):
            a = pysam.AlignedSegment(out.header)
            a.query_name = 'r%d' % i
            a.reference_id = tid[contig]
            a.reference_start = pos
            a.mapping_quality = 60
            a.flag = 0
            a.cigartuples = cigar
            a.query_sequence = 'A' * sum(length for op, length in cigar if op in (0, 1))
            if with_nm:
                a.set_tag('NM', 0)
            if xs and any(op == 3 for op, _ in cigar):
                a.set_tag('XS', strand, value_type='A')
            out.write(a)
    if index:
        pysam.index(path)


def spliced_reads(contig, exons, strand, n, rng):
    reads = []
    for _ in range(n):
        j = rng.randrange(len(exons) - 1)
        left, right = rng.randint(10, 40), rng.randint(10, 40)
        reads.append((contig, exons[j][1] - left, [(0, left), (3, exons[j + 1][0] - exons[j][1]), (0, right)], strand))
    for e in exons:
        for _ in range(n // 2):
            length = rng.randint(30, 60)
            if e[1] - e[0] > length + 2:
                reads.append((contig, rng.randint(e[0], e[1] - length), [(0, length)], strand))
    return reads


def build_cases():
    """{name: (genes, contigs, reads, options of write_bam)}"""

    rng = random.Random(1)
    contigs = [('c1', 20000), ('c2', 20000)]
    exons1 = [(1000, 1300), (3000, 3300), (5000, 5400)]
    exons2 = [(8000, 8300), (9500, 9800), (11000, 11400)]
    g1 = ('G1', 'c1', '+', [exons1, [(1000, 1300), (5000, 5400)]])
    g2 = ('G2', 'c1', '-', [exons2])
    reads = spliced_reads('c1', exons1, '+', 200, rng) + spliced_reads('c1', exons2, '-', 200, rng)

    cases = {}

    def case(name, genes, case_contigs, case_reads, **options):
        cases[name] = (genes, case_contigs, case_reads, options)

    case('normal', [g1, g2], contigs, reads)
    case('empty_bam', [g1, g2], contigs, [])
    case('contig_missing_in_bam', [g1, g2, ('G3', 'cX', '+', [[(100, 300), (700, 900)]])], contigs, reads)
    case('unstranded_dot_gene', [g1, ('GD', 'c1', '.', [[(14000, 14200), (15000, 15200)]])], contigs, reads)
    case('single_exon_genes', [('S1', 'c1', '+', [[(1000, 1300)]]), ('S2', 'c2', '-', [[(2000, 2500)]])], contigs,
         [('c1', 1010, [(0, 50)], '+')] * 30 + [('c2', 2010, [(0, 50)], '-')] * 30)
    case('only_unspliced_reads', [g1, g2], contigs, [r for r in reads if len(r[2]) == 1])
    case('no_nm_tag', [g1, g2], contigs, reads, with_nm=False)
    case('no_xs_tag', [g1, g2], contigs, reads, xs=False)
    case('one_gene_one_contig', [g1], [('c1', 20000)], [r for r in reads if r[1] < 6000])
    case('mt_contig', [g1, ('GM', 'MT', '+', [[(100, 300), (500, 700)]])], contigs + [('MT', 5000)], reads + [('MT', 110, [(0, 50)], '+')] * 20)
    many = [('G%d' % i, 'k%d' % i, '+', [[(500, 800), (1500, 1800)]]) for i in range(12)]
    case('many_small_contigs', many, [('k%d' % i, 5000) for i in range(12)],
         [r for i in range(12) for r in spliced_reads('k%d' % i, [(500, 800), (1500, 1800)], '+', 40, rng)])
    case('identical_transcripts', [('GI', 'c1', '+', [[(1000, 1300), (3000, 3300)], [(1000, 1300), (3000, 3300)]])], contigs,
         spliced_reads('c1', [(1000, 1300), (3000, 3300)], '+', 100, rng))
    case('overlapping_genes', [g1, ('GO', 'c1', '+', [[(1100, 1250), (3050, 3250), (5100, 5300)]])], contigs, reads)
    case('no_index', [g1, g2], contigs, reads, index=False)
    return cases


def fixed_defect(case, mode):
    """(defect, exception) if SplAdder is known to stop on this case in this mode, else None."""

    for name, runs, defect, exception in FIXED_DEFECTS:
        if name == case and (runs == 'all' or (runs == 'parallel' and mode != 'serial')):
            return defect, exception
    return None


def last_line(run, *paths):
    """The last line of the log, with the paths of the run replaced, so that the lines of two programs can be compared."""

    lines = [line for line in run.tail.splitlines() if line.strip()]
    text = lines[-1] if lines else ''
    for i, path in enumerate(paths):
        text = text.replace(path, '<path%d>' % i)
    return text


def judge(case, mode, upstream, turbo, diffs, n_files, replaced):
    """(status, text) of one case in one mode; status is 'identical', 'same failure', 'expected' or 'PROBLEM'."""

    known = fixed_defect(case, mode)
    if known:
        defect, exception = known
        if upstream.returncode != 0 and turbo.returncode == 0 and exception in upstream.tail:
            return 'expected', 'SplAdder stops with %s, spladder-turbo finishes (%s)' % (exception, defect)
        return 'PROBLEM', 'a known defect of SplAdder (%s, %s) was expected, but SplAdder ended with %d and spladder-turbo with %d' % (
            defect, exception, upstream.returncode, turbo.returncode)
    status, text = common.verdict(upstream, turbo, diffs, n_files)
    if status == 'identical':
        return 'identical', text
    if status == 'both failed':
        lines = last_line(upstream, *replaced[0]), last_line(turbo, *replaced[1])
        if lines[0] == lines[1]:
            return 'same failure', 'both stop with exit code %d: %s' % (upstream.returncode, lines[0][:120])
        return 'PROBLEM', 'both stop with exit code %d, with different errors: "%s" / "%s"' % (upstream.returncode, lines[0][:100], lines[1][:100])
    return 'PROBLEM', text


def run_case(task):
    directory, case, mode, options, gtf, bams, upstream, turbo, timeout = task
    runs, outdirs = [], []
    for program in (upstream, turbo):
        work = os.path.join(directory, mode.replace(' ', '_'), program.key)
        annotation = common.link_annotation(gtf, os.path.join(work, 'annotation'))
        outdir = os.path.join(work, 'out')
        runs.append(program.run(program.build_command(annotation, outdir, bams, BUILD_OPTIONS + options), os.path.join(work, 'log.txt'), cwd=work, timeout=timeout))
        outdirs.append(outdir)
    diffs, n_files = [], 0
    if runs[0].returncode == 0 and runs[1].returncode == 0:
        diffs, n_files = common.compare_outputs(*outdirs)
    replaced = [(outdir, program.src) for outdir, program in zip(outdirs, (upstream, turbo))]
    return judge(case, mode, runs[0], runs[1], diffs, n_files, replaced)


def parse_arguments(argv):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--upstream-src', metavar='DIR', help='source tree of SplAdder 3.1.1 (contains the package spladder); required, except for --list')
    parser.add_argument('--turbo-src', default=common.REPO, metavar='DIR', help='source tree of spladder-turbo [this repository]')
    parser.add_argument('--python', default=sys.executable, metavar='PYTHON', help='interpreter that has the dependencies [the current one]')
    parser.add_argument('--cases', metavar='NAME,NAME', help='only these cases (--list shows them) [all]')
    parser.add_argument('--list', action='store_true', help='show the cases and exit')
    parser.add_argument('--jobs', type=int, default=4, metavar='N', help='runs at the same time [4]')
    parser.add_argument('--timeout', type=float, default=900, metavar='SECONDS', help='stop a run that takes longer [900]')
    parser.add_argument('--workdir', metavar='DIR', help='data, outputs and logs [a new temporary directory]')
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_arguments(sys.argv[1:] if argv is None else argv)
    cases = build_cases()
    if args.list:
        print('\n'.join(cases))
        return 0
    selected = args.cases.split(',') if args.cases else list(cases)
    unknown = [name for name in selected if name not in cases]
    if unknown:
        sys.exit('unknown case(s): %s' % ', '.join(unknown))
    if not args.upstream_src:
        sys.exit('--upstream-src is required')
    common.check_source(args.upstream_src, 'upstream')
    common.check_source(args.turbo_src, 'spladder-turbo')
    workdir = os.path.abspath(args.workdir) if args.workdir else tempfile.mkdtemp(prefix='edge_cases_')
    upstream = common.Program('SplAdder', args.upstream_src, args.python, key='upstream')
    turbo = common.Program('spladder-turbo', args.turbo_src, args.python, key='turbo')

    tasks = []
    for name in selected:
        genes, contigs, reads, options = cases[name]
        directory = os.path.join(workdir, name)
        os.makedirs(directory, exist_ok=True)
        gtf = os.path.join(directory, 'a.gtf')
        common.write_gtf(gtf, genes)
        bams = [os.path.join(directory, 's%d.bam' % i) for i in (1, 2)]
        for bam in bams:
            write_bam(bam, contigs, reads, **options)
        tasks.extend((directory, name, mode, extra, gtf, bams, upstream, turbo, args.timeout) for mode, extra in MODES)

    with ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(run_case, tasks))
    counts = Counter()
    for task, (status, text) in zip(tasks, results):
        counts[status] += 1
        print('%-24s %-14s %-13s %s' % (task[1], task[2], status, text))
    print('\n%d runs: %s   (work directory %s)' % (len(tasks), ', '.join('%d %s' % (counts[s], s) for s in ('identical', 'same failure', 'expected', 'PROBLEM')), workdir))
    return 1 if counts['PROBLEM'] else 0


if __name__ == '__main__':
    sys.exit(main())
