#!/usr/bin/env python
"""Fuzz test of the compiled BAM reader (spladder/bamscan.py) against the pysam loop that it replaces.

Random BAM files are generated with every CIGAR operation, many flags, several types of the NM tag, XS tags and other tags. For each
contig and for six combinations of the read filters, the evidence that is extracted by the pysam loop (SPLADDER_NO_BAMSCAN=1) is
compared with that of the default path, which uses the compiled reader whenever it can: coverage, introns, dtypes and shapes.

* mode `clean`: only records that the compiled reader is meant to decode itself; every comparison must be identical and the reader
  must not give up.
* mode `dirty`: also records with unusual tags (a float NM, a string XS, no NM); the reader hands those contigs over to pysam, and
  the result must still be identical.

Example:

    python validation/bam_fuzz.py --files 100 --procs 8 --mode clean

Exit code: 0 if no difference was found, 1 otherwise. Needs only this repository and its dependencies.
"""

import argparse
import multiprocessing as mp
import os
import random
import shutil
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402
from spladder import bamscan, evidence  # noqa: E402

SPECS = [dict(primary_only=True, mismatch=2, exon_len=9, mm_tag='NM'), dict(primary_only=False, mismatch=0, exon_len=30, mm_tag='NM'),
         dict(primary_only=True, mismatch=None, exon_len=None, mm_tag=None), dict(primary_only=False, mismatch=1, exon_len=None, mm_tag='NM'),
         dict(primary_only=True, mismatch=None, exon_len=12, mm_tag=None), dict(primary_only=True, mismatch=300, exon_len=1, mm_tag='NM')]


def make_bam(path, seed, clean):
    rng = random.Random(seed)
    contigs = [('c%d' % i, rng.choice([2000, 60000, 300000])) for i in range(rng.randint(1, 4))]
    header = {'HD': {'VN': '1.4', 'SO': 'coordinate'}, 'SQ': [{'SN': n, 'LN': l} for n, l in contigs]}
    reads = []
    for i in range(rng.randint(0, 3000)):
        tid = rng.randrange(len(contigs))
        n_ops = rng.choice([1, 1, 1, 2, 3, 3, 5, 8])
        ops = []
        for _ in range(n_ops):
            op = rng.choice([0, 0, 0, 0, 3, 3, 1, 2, 4, 5, 6, 7, 8])
            length = rng.choice([0, 1, 5, 20, 50, 200, 1500]) if rng.random() < 0.9 else rng.choice([0, 1, 2])
            ops.append((op, length))
        ref_len = sum(length for op, length in ops if op in (0, 2, 3, 7, 8))
        pos = rng.randint(0, max(contigs[tid][1] - ref_len - 1, 0))
        if pos + ref_len > contigs[tid][1]:
            continue
        flag = rng.choice([0, 0, 0, 16, 256, 272, 2048, 1024, 4, 99, 147, 83, 163, 512, 0])
        reads.append((tid, pos, ops, flag, i))
    reads.sort(key=lambda r: (r[0], r[1]))
    with pysam.AlignmentFile(path, 'wb', header=header) as out:
        for tid, pos, ops, flag, i in reads:
            a = pysam.AlignedSegment(out.header)
            a.query_name = 'r%d' % i
            a.reference_id = tid
            a.reference_start = pos
            a.mapping_quality = rng.choice([0, 1, 60, 255])
            a.flag = flag
            a.cigartuples = ops
            qlen = sum(length for op, length in ops if op in (0, 1, 4, 7, 8))
            a.query_sequence = 'ACGT' * (qlen // 4) + 'A' * (qlen % 4)
            u = rng.random()
            if clean or u < 0.93:
                value = rng.choice([0, 0, 0, 1, 2, 3, 5, 9, 40, 300, 70000])
                kind = rng.choice([None, None, 'c', 'C', 's', 'S', 'i', 'I']) if value < 120 else rng.choice(['S', 'i', 'I', None])
                if kind == 'c' and value > 127:
                    kind = 'C'
                if kind == 's' and value > 32767:
                    kind = 'S'
                if kind == 'S' and value > 65535:
                    kind = 'i'
                if kind is None:
                    a.set_tag('NM', value)
                else:
                    a.set_tag('NM', value, value_type=kind)
            elif u < 0.95:
                a.set_tag('NM', 1.5, value_type='f')           # unusual type: handed over
            # (otherwise there is no NM at all)
            for tag in ('NH', 'HI', 'AS'):
                if rng.random() < 0.5:
                    a.set_tag(tag, rng.randint(0, 100))
            if rng.random() < 0.3:
                a.set_tag('MD', '10A5', value_type='Z')
            w = rng.random()
            if w < 0.3:
                a.set_tag('XS', rng.choice(['+', '-', '.']), value_type='A')
            elif w < 0.33 and not clean:
                a.set_tag('XS', '+', value_type='Z')           # unusual type: handed over
            if rng.random() < 0.2:
                a.set_tag('ZZ', 'text', value_type='Z')
            out.write(a)
    pysam.index(path)


def run_one(arguments):
    seed, clean = arguments
    directory = tempfile.mkdtemp(prefix='bamfuzz')
    try:
        path = os.path.join(directory, 'f.bam')
        make_bam(path, seed, clean)
        comparisons = differences = handed_over = 0
        names = [s['SN'] for s in pysam.AlignmentFile(path).header.to_dict()['SQ']]
        for contig in names:
            for spec_values in SPECS:
                spec = evidence.Spec(path, None, spec_values['primary_only'], spec_values['mismatch'], spec_values['exon_len'], spec_values['mm_tag'])
                results = []
                for disable in ('1', ''):
                    os.environ['SPLADDER_NO_BAMSCAN'] = disable
                    os.environ['SPLADDER_DEBUG_BAMSCAN'] = '1'
                    try:
                        results.append(evidence.build_contig_evidence(path, contig, spec))
                    except evidence.Unsupported:
                        results.append('unsupported')
                comparisons += 1
                if not common.same_evidence(results[0], results[1]):
                    differences += 1
                    print('DIFFERENT: seed %d contig %s filters %s' % (seed, contig, spec_values), flush=True)
                try:
                    tid = pysam.AlignmentFile(path).gettid(contig)
                    bamscan.scan_contig(path, tid, len(names), spec_values['primary_only'], spec_values['mm_tag'], spec_values['mismatch'], spec_values['exon_len'])
                except bamscan.NotSupported:
                    handed_over += 1
        return comparisons, differences, handed_over
    finally:
        shutil.rmtree(directory, ignore_errors=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--files', type=int, default=50, metavar='N', help='number of random BAM files [50]')
    parser.add_argument('--procs', type=int, default=4, metavar='N', help='processes [4]')
    parser.add_argument('--mode', choices=['clean', 'dirty'], default='clean', help='see above [clean]')
    parser.add_argument('--first-seed', type=int, default=5000, metavar='N', help='seed of the first file [5000]')
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)
    with mp.Pool(args.procs) as pool:
        results = pool.map(run_one, [(seed, args.mode == 'clean') for seed in range(args.first_seed, args.first_seed + args.files)])
    comparisons, differences, handed_over = (sum(r[i] for r in results) for i in range(3))
    print('%s BAM files: %d   comparisons: %d   DIFFERENT: %d   handed over by the compiled reader: %d' % (args.mode, args.files, comparisons, differences, handed_over))
    problems = differences + (handed_over if args.mode == 'clean' else 0)
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
