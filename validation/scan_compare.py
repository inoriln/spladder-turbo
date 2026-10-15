#!/usr/bin/env python
"""The compiled BAM reader against the pysam loop on real alignment files.

For every contig that has alignments, and for three settings of the read filters (those used to build graphs, those used to count,
and one that takes all alignments), the evidence that spladder-turbo extracts from the reads (coverage and introns) is built
twice: with the pysam loop that SplAdder uses (SPLADDER_NO_BAMSCAN=1) and with the compiled reader. Values, dtypes and shapes
must be identical. The tool also reports the CPU time of both and the number of cases in which the compiled reader handed the
contig over to pysam (an alignment file with something it does not decode itself; the result is then the pysam result).

Example:

    python validation/scan_compare.py --bams a.bam,b.bam --procs 16

Exit code: 0 if everything is identical, 1 otherwise.
"""

import argparse
import multiprocessing as mp
import os
import sys
import time

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402
from spladder import bamscan, evidence  # noqa: E402

SETTINGS = {'graph': dict(primary_only=True, mismatch=2, exon_len=9, mm_tag='NM'),
            'count': dict(primary_only=True, mismatch=None, exon_len=None, mm_tag=None),
            'all': dict(primary_only=False, mismatch=1, exon_len=13, mm_tag='NM')}


def compare_one(task):
    """[(file, contig, setting, identical, handed over, seconds of pysam, seconds of the compiled reader)] for the settings of one contig."""

    fname, contig = task
    rows = []
    for name, values in SETTINGS.items():
        spec = evidence.Spec(fname, None, values['primary_only'], values['mismatch'], values['exon_len'], values['mm_tag'])
        os.environ['SPLADDER_NO_BAMSCAN'] = '1'
        start = time.time()
        reference = evidence.build_contig_evidence(fname, contig, spec)
        seconds_pysam = time.time() - start
        os.environ['SPLADDER_NO_BAMSCAN'] = ''
        os.environ['SPLADDER_DEBUG_BAMSCAN'] = '1'                  # a failure of the compiled reader is an error, not a fall-back
        start = time.time()
        result = evidence.build_contig_evidence(fname, contig, spec)
        seconds_compiled = time.time() - start
        with pysam.AlignmentFile(fname) as infile:
            tid, n_ref = infile.gettid(contig), infile.nreferences
        try:
            bamscan.scan_contig(fname, tid, n_ref, values['primary_only'], values['mm_tag'], values['mismatch'], values['exon_len'])
            handed_over = False
        except bamscan.NotSupported:
            handed_over = True
        rows.append((os.path.basename(fname), contig, name, common.same_evidence(reference, result), handed_over, seconds_pysam, seconds_compiled))
    return rows


def contigs_with_reads(fname):
    with pysam.AlignmentFile(fname) as infile:
        return [s.contig for s in sorted(infile.get_index_statistics(), key=lambda s: -s.mapped) if s.mapped > 0]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--bams', '-b', required=True, metavar='A.bam,B.bam', help='sorted, indexed alignment files, comma separated')
    parser.add_argument('--contigs', metavar='NAME,NAME', help='only these contigs [all that have alignments]')
    parser.add_argument('--procs', type=int, default=4, metavar='N', help='processes [4]')
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)

    tasks = []
    for fname in args.bams.split(','):
        if not os.path.exists(fname):
            sys.exit('alignment file %s does not exist' % fname)
        contigs = args.contigs.split(',') if args.contigs else contigs_with_reads(fname)
        tasks.extend((os.path.abspath(fname), contig) for contig in contigs)
    ### (the contigs of a file come in the order of their size, the largest first: the long ones start early)
    with mp.Pool(args.procs) as pool:
        rows = [row for part in pool.imap_unordered(compare_one, tasks) for row in part]

    different = [r for r in rows if not r[3]]
    handed_over = sum(r[4] for r in rows)
    print('compared %d (file, contig, setting) combinations: %d DIFFERENT, %d handed over to pysam by the compiled reader' % (len(rows), len(different), handed_over))
    for r in different[:10]:
        print('  DIFFERENT: %s %s %s' % r[:3])
    total_pysam, total_compiled = sum(r[5] for r in rows), sum(r[6] for r in rows)
    if total_compiled > 0:
        print('CPU time, pysam loop: %.0f s, compiled reader: %.0f s (%.1f x)' % (total_pysam, total_compiled, total_pysam / total_compiled))
        for name in SETTINGS:
            a, b = sum(r[5] for r in rows if r[2] == name), sum(r[6] for r in rows if r[2] == name)
            print('  %-6s pysam %7.0f s   compiled %7.0f s   (%.1f x)' % (name, a, b, a / b if b else 0))
    return 1 if different else 0


if __name__ == '__main__':
    sys.exit(main())
