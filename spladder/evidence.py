# Copyright (c) 2026, Yihao Zhang. Part of spladder-turbo, a fork of SplAdder (see ATTRIBUTION.md).
# SPDX-License-Identifier: BSD-3-Clause
"""Single-pass read evidence for coverage and intron queries.

SplAdder asks for the read coverage and the spliced alignments of one gene after the other and does so in
several processing steps (intron list, cassette exons, intron retentions, new intron edges, counting).
Answering each request from the alignment file decodes every read once per overlapping gene and step.

Inside an `evidence_scope` the coverage and the introns of a complete contig are instead extracted in ONE
sequential pass, kept in a compact form, and every window request of `reads.add_reads_from_bam` is answered
from there.  The answers are identical to the ones obtained from the alignment file, including dtype and
shape (see tests/test_evidence.py), because

* the per-read filters are the ones of `reads.get_reads` and do not depend on the requested window, except for
  the strand check, which is resolved by keeping the reads apart by their XS tag,
* the coverage at a position only depends on the reads covering that position, which all overlap every
  window containing the position,
* an intron lying completely inside a window belongs to reads overlapping the window, so its read count is the
  same, whether it is counted over the window or over the contig.

Anything that cannot be answered this way (several alignment files, variation aware filters, unusual strand
values, alignments the original code would abort on, ...) is passed back to the per-window implementation,
so the behaviour in those cases is the original one.

The pass over a contig is made by bamscan (a compiled reader for indexed BAM files) if the file allows for it, and
with pysam otherwise; both give the same evidence.

Setting the environment variable SPLADDER_NO_EVIDENCE=1 disables the mechanism altogether, SPLADDER_NO_BAMSCAN=1
the compiled reader.
"""

import os

from array import array
from collections import OrderedDict, defaultdict, namedtuple
from contextlib import contextmanager

import numpy as np
import pysam

from . import bamscan
from . import parallel
from .utils import sort_rows

ENV_DISABLE = 'SPLADDER_NO_EVIDENCE'
ENV_MAX_GB = 'SPLADDER_EVIDENCE_MAX_GB'
ENV_NO_SCAN = 'SPLADDER_NO_BAMSCAN'
ENV_DEBUG_SCAN = 'SPLADDER_DEBUG_BAMSCAN'
DEFAULT_MAX_BYTES = 6 * 1024 ** 3

### read classes by the XS tag: no tag, '+', '-', any other value
_NONE, _PLUS, _MINUS, _OTHER = 0, 1, 2, 3

_SUPPORTED_TYPES = {('exon_track',), ('intron_list',), ('exon_track', 'intron_list'), ('intron_list', 'exon_track')}

_scopes = []


class Unsupported(Exception):
    """The evidence of this contig cannot be built; the per-window code path has to be used."""


### everything that determines which reads of a contig are taken into account
Spec = namedtuple('Spec', ['fname', 'cram_ref', 'primary_only', 'mismatch', 'exon_len', 'mm_tag'])


class ContigEvidence:
    """Coverage and introns of one contig for one Spec.

    cov     -- per read class a step function (positions, coverage from that position on)
    introns -- for the strands '+', '-' and None: ((introns_p, starts_p), (introns_m, starts_m)), the
               introns being sorted (n, 3) arrays [start, stop, count] with the start column duplicated
               as contiguous array for the binary search
    """

    __slots__ = ('cov', 'introns', 'nbytes')

    def __init__(self, cov, introns):
        self.cov = cov
        self.introns = introns
        seen = {}
        for pos, val in cov:
            seen[id(pos)] = pos
            seen[id(val)] = val
        for views in introns.values():
            for arr, starts in views:
                seen[id(arr)] = arr
                seen[id(starts)] = starts
        self.nbytes = sum(a.nbytes for a in seen.values())


def _empty_introns():
    return np.zeros((0, 3), dtype='uint32')


def _empty_evidence():
    empty = (np.zeros(0, dtype='int64'), np.zeros(0, dtype='int32'))
    e = _empty_introns()
    es = np.zeros(0, dtype='int64')
    return ContigEvidence((empty,) * 4, {s: ((e, es), (e, es)) for s in ('+', '-', None)})


def _step_function(starts, stops):
    """Breakpoints and coverage values from the starts and (exclusive) stops of aligned blocks."""

    if len(starts) == 0:
        return (np.zeros(0, dtype='int64'), np.zeros(0, dtype='int32'))
    s = np.frombuffer(starts, dtype='int64')
    e = np.frombuffer(stops, dtype='int64')
    pos = np.concatenate((s, e))
    delta = np.concatenate((np.ones(s.shape[0], dtype='int32'), -np.ones(e.shape[0], dtype='int32')))
    order = np.argsort(pos, kind='stable')
    pos = pos[order]
    delta = delta[order]
    first = np.flatnonzero(np.concatenate(([True], pos[1:] != pos[:-1])))
    change = np.add.reduceat(delta, first)
    value = np.cumsum(change, dtype='int32')
    keep = change != 0
    return (pos[first][keep], value[keep])


def _intron_array(counts):
    if len(counts) == 0:
        return _empty_introns()
    return sort_rows(np.array([[k[0], k[1], v] for k, v in counts.items()], dtype='uint32'))


def _merge_intron_arrays(arrays):
    """Adds up the counts of identical introns of several sorted intron arrays."""

    arrays = [a for a in arrays if a.shape[0] > 0]
    if len(arrays) == 0:
        return _empty_introns()
    if len(arrays) == 1:
        return arrays[0]
    a = np.concatenate(arrays).astype('int64')
    a = a[np.lexsort((a[:, 1], a[:, 0]))]
    first = np.flatnonzero(np.concatenate(([True], (a[1:, 0] != a[:-1, 0]) | (a[1:, 1] != a[:-1, 1]))))
    out = np.empty((first.shape[0], 3), dtype='uint32')
    out[:, 0] = a[first, 0]
    out[:, 1] = a[first, 1]
    out[:, 2] = np.add.reduceat(a[:, 2], first)
    return out


def _with_starts(arr):
    return (arr, np.ascontiguousarray(arr[:, 0]).astype('int64'))


def build_contig_evidence(fname, contig, spec):
    """One sequential pass over `contig` with the read filters of reads.get_reads.

    Raises Unsupported, if the original code path has to be used for this contig."""

    ### TODO THIS IS A HACK (the original implementation never looks at the MT contig)
    if contig == 'MT':
        return _empty_evidence()

    try:
        if fname.lower().endswith('.bam'):
            infile = pysam.AlignmentFile(fname, 'rb')
        elif fname.lower().endswith('.cram'):
            infile = pysam.AlignmentFile(fname, 'rc', reference_filename=spec.cram_ref, ignore_truncation=True)
        else:
            raise Unsupported('unknown alignment format')
    except Unsupported:
        raise
    except Exception as e:
        raise Unsupported(str(e))

    try:
        ### without an index fetch would return the reads of all contigs
        if not infile.has_index():
            raise Unsupported('alignment file is not indexed')
        tid = infile.gettid(contig)
        if tid < 0:
            return _empty_evidence()

        ### the compiled reader gives the very same answers, where it can be used
        scanned = None
        if os.environ.get(ENV_NO_SCAN, '') in ('', '0'):
            try:
                scanned = bamscan.scan_contig(fname, tid, infile.nreferences, spec.primary_only, spec.mm_tag, spec.mismatch, spec.exon_len)
            except bamscan.NotSupported:
                scanned = None
            except Exception:
                ### whatever goes wrong there, the general reader is still there
                if os.environ.get(ENV_DEBUG_SCAN, '') not in ('', '0'):
                    raise
                scanned = None

        if scanned is None:
            starts = [array('q') for _ in range(4)]
            stops = [array('q') for _ in range(4)]
            introns = [defaultdict(int) for _ in range(4)]
            primary_only = spec.primary_only
            check_nm = spec.mismatch is not None
            mismatch = spec.mismatch
            mm_tag = spec.mm_tag
            check_len = spec.exon_len is not None
            exon_len = spec.exon_len

            try:
                for read in infile.fetch(contig, until_eof=True):

                    flag = read.flag
                    if flag & 4:                          # unmapped
                        continue
                    if primary_only and (flag & 256):     # secondary
                        continue

                    cigar = read.cigartuples
                    is_spliced = False
                    for op, ln in cigar:
                        if op == 3:
                            is_spliced = True
                            break

                    if check_nm and mismatch < read.get_tag(mm_tag):
                        continue

                    if check_len and is_spliced:
                        ### min segment length: aligned length (M, D, =, X) of the shortest segment between introns
                        seg = 0
                        min_seg = -1
                        for op, ln in cigar:
                            if op == 3:
                                if min_seg < 0 or seg < min_seg:
                                    min_seg = seg
                                seg = 0
                            elif op == 0 or op == 2 or op == 7 or op == 8:
                                seg += ln
                        if min_seg < 0 or seg < min_seg:
                            min_seg = seg
                        if min_seg <= exon_len:
                            continue

                    cls = _NONE
                    if read.has_tag('XS'):
                        xs = read.get_tag('XS')
                        cls = _PLUS if xs == '+' else (_MINUS if xs == '-' else _OTHER)

                    st = starts[cls]
                    sp = stops[cls]
                    p = read.reference_start
                    if len(cigar) == 1 and cigar[0][0] == 0:
                        ### the most frequent case: a single aligned block
                        st.append(p)
                        sp.append(p + cigar[0][1])
                        continue
                    for op, ln in cigar:
                        if op == 3:
                            introns[cls][(p, p + ln)] += 1
                            p += ln
                        elif op == 0 or op == 2 or op == 7 or op == 8:
                            st.append(p)
                            sp.append(p + ln)
                            p += ln
            except (KeyError, TypeError, ValueError, OverflowError) as e:
                ### e.g. a missing NM tag: the original code decides whether this matters for the requested windows
                raise Unsupported(repr(e))
    finally:
        infile.close()

    if scanned is None:
        cov = tuple(_step_function(starts[c], stops[c]) for c in range(4))
        arr = [_intron_array(introns[c]) for c in range(4)]
    else:
        cov = tuple(_step_function(scanned[0][c][0], scanned[0][c][1]) for c in range(4))
        arr = scanned[1]
    ### the strand a window is requested for decides which reads were taken into account:
    ### reads without XS go to the plus-type introns, XS '-' reads to the minus-type introns
    p_plus = _with_starts(_merge_intron_arrays([arr[_NONE], arr[_PLUS]]))
    p_minus = _with_starts(arr[_NONE])
    p_all = _with_starts(_merge_intron_arrays([arr[_NONE], arr[_PLUS], arr[_OTHER]]))
    m_minus = _with_starts(arr[_MINUS])
    m_none = _with_starts(_empty_introns())
    return ContigEvidence(cov, {'+': (p_plus, m_none), '-': (p_minus, m_minus), None: (p_all, m_minus)})


def _coverage(ev, start, stop, strand):
    """Integer coverage of [start, stop) of the reads taken into account for `strand`."""

    out = np.zeros(stop - start, dtype='int64')
    if strand == '+':
        classes = (_NONE, _PLUS)
    elif strand == '-':
        classes = (_NONE, _MINUS)
    else:
        classes = (_NONE, _PLUS, _MINUS, _OTHER)
    for c in classes:
        pos, val = ev.cov[c]
        if pos.shape[0] == 0:
            continue
        i0 = int(np.searchsorted(pos, start, side='right'))
        i1 = int(np.searchsorted(pos, stop, side='left'))
        base = int(val[i0 - 1]) if i0 > 0 else 0
        if i1 <= i0:
            if base:
                out += base
            continue
        bounds = np.empty(i1 - i0 + 2, dtype='int64')
        bounds[0] = start
        bounds[1:-1] = pos[i0:i1]
        bounds[-1] = stop
        values = np.empty(i1 - i0 + 1, dtype='int64')
        values[0] = base
        values[1:] = val[i0:i1]
        out += np.repeat(values, np.diff(bounds))
    return out


def _introns_inside(view, start, stop, win_start, win_stop):
    """Introns of a (sorted array, starts) view completely inside the window; `win_*` are the original objects."""

    arr, starts = view
    if arr.shape[0] == 0:
        return _empty_introns()
    lo = int(np.searchsorted(starts, start, side='right'))
    hi = int(np.searchsorted(starts, stop, side='left'))
    sub = arr[lo:hi]
    return sub[(sub[:, 0] > win_start) & (sub[:, 1] < win_stop)]


class EvidenceScope:
    """Holds the contig evidence of the alignment files used while the scope is open."""

    def __init__(self, max_bytes=None, enabled=True):
        if max_bytes is None:
            max_bytes = int(float(os.environ.get(ENV_MAX_GB, DEFAULT_MAX_BYTES / 1024 ** 3)) * 1024 ** 3)
        self.max_bytes = max_bytes
        self.enabled = enabled
        self.nbytes = 0
        self.stats = {'builds': 0, 'hits': 0, 'fallbacks': 0, 'evictions': 0}
        self._store = OrderedDict()

    def contig_evidence(self, spec, contig):
        """Returns the ContigEvidence or None, if the per-window implementation has to be used."""

        key = (spec, contig)
        if key in self._store:
            self._store.move_to_end(key)
            return self._store[key]
        try:
            ev = build_contig_evidence(spec.fname, contig, spec)
            self.stats['builds'] += 1
        except Unsupported:
            ev = None
        self.put(key, ev)
        return ev

    def prefetch(self, spec, contigs, workers):
        """Builds the evidence of several contigs in parallel (the contigs with most reads first), if that pays off."""

        todo = [c for c in contigs if (spec, c) not in self._store]
        if workers < 2 or len(todo) < 2 or not self.enabled or not parallel.can_fork():
            return
        try:
            with pysam.AlignmentFile(spec.fname, 'rc' if spec.fname.lower().endswith('.cram') else 'rb', reference_filename=spec.cram_ref) as infile:
                reads = {s.contig: s.mapped for s in infile.get_index_statistics()}
            todo.sort(key=lambda c: -reads.get(c, 0))
        except Exception:
            pass
        for contig, ev in zip(todo, parallel.ordered_map(_build_task, [(spec.fname, c, spec) for c in todo], min(workers, len(todo)))):
            if ev is not None:
                self.stats['builds'] += 1
            self.put((spec, contig), ev)

    def put(self, key, ev):
        self._store[key] = ev
        if ev is None:
            return
        self.nbytes += ev.nbytes
        ### evict the least recently used evidence, but never the one that was just built
        while self.nbytes > self.max_bytes:
            victim = next((k for k, v in self._store.items() if v is not None and k != key), None)
            if victim is None:
                break
            self.nbytes -= self._store.pop(victim).nbytes
            self.stats['evictions'] += 1

    def release(self):
        self._store.clear()
        self.nbytes = 0


def _build_task(task):
    fname, contig, spec = task
    try:
        return build_contig_evidence(fname, contig, spec)
    except Unsupported:
        return None


def current_scope():
    return _scopes[-1] if _scopes else None


@contextmanager
def evidence_scope(max_bytes=None):
    """Answers the read requests of `reads.add_reads_from_bam` from single-pass evidence while the scope is open."""

    scope = EvidenceScope(max_bytes=max_bytes, enabled=os.environ.get(ENV_DISABLE, '') in ('', '0'))
    _scopes.append(scope)
    try:
        yield scope
    finally:
        _scopes.pop()
        scope.release()


def prefetch_contigs(scope, bam_fnames, contigs, read_filter, options, no_mm=False):
    """Reads the contigs of the alignment file of a sample in parallel, if the options allow for several processes."""

    if options.parallel < 2 or options.sparse_bam or options.var_aware:
        return
    fname = bam_fnames if isinstance(bam_fnames, str) else (bam_fnames[0] if len(bam_fnames) == 1 else None)
    if not isinstance(fname, str):
        return
    spec = make_spec(fname, read_filter, options.primary_only, no_mm, options.mm_tag, options.ref_genome)
    if spec is not None:
        scope.prefetch(spec, sorted(set(contigs)), options.parallel)


def make_spec(fname, filter, primary_only, no_mm, mm_tag, cram_ref):
    """Spec of the reads taken into account; None, if the requested filter cannot be answered from evidence."""

    if filter is None:
        return Spec(fname, cram_ref, bool(primary_only), None, None, None)
    if not isinstance(filter, dict) or 'maps' in filter or 'mismatch' not in filter or 'exon_len' not in filter:
        return None
    if filter['exon_len'] is None or (not no_mm and filter['mismatch'] is None):
        return None
    return Spec(fname, cram_ref, bool(primary_only), None if no_mm else filter['mismatch'], filter['exon_len'], None if no_mm else mm_tag)


def add_reads_from_evidence(blocks, filenames, types, filter, var_aware, primary_only, no_mm, unstranded, mm_tag, cram_ref):
    """Same result as the per-window code of reads.add_reads_from_bam, or NotImplemented."""

    scope = current_scope()
    if scope is None or not scope.enabled:
        return NotImplemented

    spec = None
    if tuple(types) in _SUPPORTED_TYPES and not var_aware and blocks.shape[0] == 1:
        if isinstance(filenames, str):
            fname = filenames
        elif len(filenames) == 1 and isinstance(filenames[0], str):
            fname = filenames[0]
        else:
            fname = None
        gg = blocks[0]
        if (fname is not None and gg.strand in ('+', '-', None) and isinstance(gg.chr, str)
                and isinstance(gg.start, (int, np.integer)) and isinstance(gg.stop, (int, np.integer))
                and gg.stop > gg.start and os.path.exists(fname)):
            spec = make_spec(fname, filter, primary_only, no_mm, mm_tag, cram_ref)
    if spec is None:
        scope.stats['fallbacks'] += 1
        return NotImplemented

    ev = scope.contig_evidence(spec, gg.chr)
    if ev is None:
        scope.stats['fallbacks'] += 1
        return NotImplemented
    scope.stats['hits'] += 1

    start = int(gg.start)
    stop = int(gg.stop)
    block_len = int(gg.stop - gg.start)
    tracks = np.zeros((0, block_len))
    intron_list = []
    for ttype in types:
        if ttype == 'exon_track':
            coverage = np.zeros((1, block_len))
            coverage += _coverage(ev, start, stop, gg.strand)
            tracks = np.r_[tracks, coverage]
        else:
            (introns_p, starts_p), (introns_m, starts_m) = ev.introns[gg.strand]
            introns_p = _introns_inside((introns_p, starts_p), start, stop, gg.start, gg.stop)
            introns_m = _introns_inside((introns_m, starts_m), start, stop, gg.start, gg.stop)
            if introns_p.shape[0] > 0 or introns_m.shape[0] > 0:
                if unstranded:
                    introns = sort_rows(np.r_[introns_p, introns_m])
                else:
                    introns = introns_m if gg.strand == '-' else introns_p
                if filter is not None and 'mincount' in filter:
                    take_idx = np.where(introns[:, 2] >= filter['mincount'])[0]
                    if take_idx.shape[0] > 0:
                        intron_list.append(introns[take_idx, :])
                    else:
                        intron_list.append(_empty_introns())
                else:
                    intron_list.append(introns)
            else:
                intron_list.append(_empty_introns())

    if len(types) == 1 and types[0] == 'intron_list':
        return intron_list
    elif 'intron_list' in types:
        return (tracks, intron_list)
    else:
        return tracks
