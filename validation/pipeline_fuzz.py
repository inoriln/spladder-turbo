#!/usr/bin/env python
"""Differential fuzz of the whole pipeline: random annotations and reads are run through SplAdder and spladder-turbo.

For every seed a small data set is generated: 2 to 6 contigs with 2 to 6 genes each, 2 or 3 samples, and reads from the annotated
transcripts and from variants that are not in the annotation: novel cassette exons, retained introns, alternative 5' and 3' sites,
novel junctions that skip two exons, mutually exclusive exons, single exon transcripts. SplAdder runs once, spladder-turbo runs
serially and with --parallel 3 (see --modes); every output must be identical.

Example:

    python validation/pipeline_fuzz.py --upstream-src ~/src/spladder-3.1.1 --first-seed 0 --count 10

Exit code: 0 if no problem was found, 1 otherwise.
"""

import argparse
import os
import random
import shutil
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402

BUILD_OPTIONS = ['--merge-strat', 'merge_graphs', '--extract-ase', '--readlen', '60', '-v']


def build_dataset(seed, directory):
    """Writes the annotation a.gtf and the samples s0.bam, s1.bam, ... to `directory`; returns their paths."""

    rng = random.Random(seed)
    n_contigs = rng.randint(2, 6)
    contigs = [('ct%d' % i, 400000) for i in range(n_contigs)]
    genes = []        # (name, contig, strand, annotated transcripts, transcripts of the reads, introns with reads inside, exons)
    for ci, (cname, clen) in enumerate(contigs):
        pos = 5000
        for gi in range(rng.randint(2, 6)):
            n_ex = rng.randint(2, 6)
            exons = []
            p = pos
            for _ in range(n_ex):
                p += rng.choice([rng.randint(200, 800), rng.randint(800, 4000)])
                ln = rng.choice([rng.randint(60, 150), rng.randint(150, 500)])
                exons.append((p, p + ln))
                p += ln
            strand = rng.choice(['+', '-'])
            annotated = [exons]
            if n_ex > 3 and rng.random() < 0.6:
                skip = rng.randrange(1, n_ex - 1)
                annotated.append([e for k, e in enumerate(exons) if k != skip])
            if rng.random() < 0.2:
                annotated.append(exons[:1])        # a single exon transcript
            if n_ex >= 4 and rng.random() < 0.4:   # two mutually exclusive exons
                i = rng.randrange(1, n_ex - 2)
                # (without the transcript that has both: the two exons are only used alternatively)
                annotated[:] = [[e for k, e in enumerate(exons) if k != i + 1], [e for k, e in enumerate(exons) if k != i]] + annotated[1:]
            # transcripts of the reads: the annotated ones plus novel variants
            true = [list(t) for t in annotated if len(t) > 1]
            if rng.random() < 0.5:        # novel cassette exon inside the largest intron
                k = max(range(n_ex - 1), key=lambda i: exons[i + 1][0] - exons[i][1])
                a, b = exons[k][1], exons[k + 1][0]
                if b - a > 400:
                    start = a + rng.randint(100, (b - a) // 2)
                    mid = (start, start + rng.randint(60, 150))
                    t = list(exons)
                    t.insert(k + 1, mid)
                    true.append(t)
            if rng.random() < 0.5:        # alternative 5' / 3' site
                k = rng.randrange(n_ex - 1)
                t = list(exons)
                t[k] = (t[k][0], max(t[k][0] + 40, t[k][1] - rng.randint(15, 40)))
                true.append(t)
            if rng.random() < 0.4:        # novel junction skipping two exons
                if n_ex >= 4:
                    k = rng.randrange(n_ex - 3)
                    true.append(exons[:k + 1] + exons[k + 3:])
            retained = [k for k in range(n_ex - 1) if rng.random() < 0.3]
            genes.append((('G%d_%d' % (ci, gi)), cname, strand, annotated, true, retained, exons))
            pos = p + rng.randint(2000, 6000)

    gtf = os.path.join(directory, 'a.gtf')
    common.write_gtf(gtf, [(name, cname, strand, annotated) for name, cname, strand, annotated, _, _, _ in genes])

    bams = []
    for si in range(rng.randint(2, 3)):
        header = {'HD': {'VN': '1.4', 'SO': 'coordinate'}, 'SQ': [{'SN': n, 'LN': l} for n, l in contigs]}
        tid = {n: i for i, (n, _) in enumerate(contigs)}
        reads = []
        for name, cname, strand, annotated, true, retained, exons in genes:
            depth = rng.choice([10, 40, 120, 300])
            for iso in true:
                w = rng.choice([0.3, 1.0, 2.0])
                for _ in range(int(depth * w)):
                    k = rng.randrange(len(iso) - 1)
                    left = min(rng.choice([rng.randint(5, 14), rng.randint(15, 60)]), iso[k][1] - iso[k][0])
                    right = min(rng.choice([rng.randint(5, 14), rng.randint(15, 60)]), iso[k + 1][1] - iso[k + 1][0])
                    cigar = [(0, left), (3, iso[k + 1][0] - iso[k][1]), (0, right)]
                    reads.append((cname, iso[k][1] - left, cigar, strand))
                for e in iso:
                    for _ in range(int(depth * w) // 2):
                        ln = rng.randint(30, 70)
                        if e[1] - e[0] > ln + 2:
                            reads.append((cname, rng.randint(e[0], e[1] - ln), [(0, ln)], strand))
            for k in retained:
                a, b = exons[k][1], exons[k + 1][0]
                for _ in range(int(depth * rng.choice([0.2, 0.6, 1.0]) * (b - a) / 60)):
                    ln = rng.randint(30, 70)
                    if b - a > ln + 2:
                        reads.append((cname, rng.randint(a, b - ln), [(0, ln)], strand))
        reads.sort(key=lambda r: (tid[r[0]], r[1]))
        path = os.path.join(directory, 's%d.bam' % si)
        with pysam.AlignmentFile(path, 'wb', header=header) as out:
            for i, (c, pos, cigar, strand) in enumerate(reads):
                a = pysam.AlignedSegment(out.header)
                a.query_name = 'r%d' % i
                a.reference_id = tid[c]
                a.reference_start = pos
                a.mapping_quality = 60
                a.flag = rng.choice([0, 0, 16, 256 if rng.random() < 0.05 else 0])
                a.cigartuples = cigar
                a.query_sequence = 'A' * sum(ln for op, ln in cigar if op in (0, 1))
                a.set_tag('NM', rng.choice([0, 0, 0, 1, 2, 4]))
                if any(op == 3 for op, _ in cigar) and rng.random() < 0.9:
                    a.set_tag('XS', strand if rng.random() < 0.9 else ('+' if strand == '-' else '-'), value_type='A')
                out.write(a)
        pysam.index(path)
        bams.append(path)
    return gtf, bams


def run_build(program, gtf, bams, directory, parallel):
    annotation = common.link_annotation(gtf, os.path.join(directory, 'annotation'))
    outdir = os.path.join(directory, 'out')
    extra = BUILD_OPTIONS + (['--parallel', str(parallel)] if parallel > 1 else [])
    return program.run(program.build_command(annotation, outdir, bams, extra), os.path.join(directory, 'log.txt'), cwd=directory, timeout=1800), outdir


def parse_arguments(argv):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--upstream-src', required=True, metavar='DIR', help='source tree of SplAdder 3.1.1 (contains the package spladder)')
    parser.add_argument('--turbo-src', default=common.REPO, metavar='DIR', help='source tree of spladder-turbo [this repository]')
    parser.add_argument('--python', default=sys.executable, metavar='PYTHON', help='interpreter that has the dependencies [the current one]')
    parser.add_argument('--first-seed', type=int, default=0, metavar='N', help='first seed [0]')
    parser.add_argument('--count', type=int, default=5, metavar='N', help='number of seeds [5]')
    parser.add_argument('--modes', default='1,3', metavar='N,N', help='values of --parallel for spladder-turbo, 1 is serial [1,3]')
    parser.add_argument('--workdir', metavar='DIR', help='data sets, outputs and logs [a new temporary directory]')
    parser.add_argument('--cleanup', action='store_true', help='delete the files of a seed when it is identical')
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_arguments(sys.argv[1:] if argv is None else argv)
    common.check_source(args.upstream_src, 'upstream')
    common.check_source(args.turbo_src, 'spladder-turbo')
    modes = [int(m) for m in args.modes.split(',')]
    workdir = os.path.abspath(args.workdir) if args.workdir else tempfile.mkdtemp(prefix='pipeline_fuzz_')
    os.makedirs(workdir, exist_ok=True)
    upstream = common.Program('SplAdder', args.upstream_src, args.python)
    turbo = common.Program('spladder-turbo', args.turbo_src, args.python)

    problems = 0
    for seed in range(args.first_seed, args.first_seed + args.count):
        directory = os.path.join(workdir, 'seed%d' % seed)
        os.makedirs(directory, exist_ok=True)
        gtf, bams = build_dataset(seed, directory)
        run_up, out_up = run_build(upstream, gtf, bams, os.path.join(directory, 'upstream'), 1)
        parts, seed_ok = [], True
        for parallel in modes:
            run_tb, out_tb = run_build(turbo, gtf, bams, os.path.join(directory, 'turbo_p%d' % parallel), parallel)
            diffs, n_files = ([], 0)
            if run_up.returncode == 0 and run_tb.returncode == 0:
                diffs, n_files = common.compare_outputs(out_up, out_tb)
            status, text = common.verdict(run_up, run_tb, diffs, n_files)
            parts.append('turbo --parallel %d: %s' % (parallel, text))
            seed_ok = seed_ok and status == 'identical'
        problems += not seed_ok
        print('seed %-5d %s' % (seed, ' | '.join(parts)), flush=True)
        if seed_ok and args.cleanup:
            shutil.rmtree(directory, ignore_errors=True)
    print('problems: %d   (work directory %s)' % (problems, workdir))
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
