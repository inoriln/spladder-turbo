"""Reading the alignments of one contig of an indexed BAM file without creating an object per alignment.

evidence.build_contig_evidence needs a few facts about every alignment of a contig: where its aligned blocks and
its introns are, and whether it passes the read filters.  Asking pysam for these costs a Python object and
several attribute lookups per alignment, which is about four times the time it takes to decompress the file.  Here
the contig is located with the index (BAI), its BGZF blocks are decompressed, and the alignment records are decoded
by a compiled function, which applies exactly the filters of evidence.build_contig_evidence (same order, same
conditions).

This is an accelerator only: whenever something is not the usual (no BAI, a CRAM file, a tag of an unusual type,
a missing NM tag, a record that does not look right, ...), NotSupported is raised and the caller uses the general
reader, which gives the answers (or raises the errors) of the original implementation.
"""

import os
import struct
import zlib

import numpy as np

from .utils import jit_cached


class NotSupported(Exception):
    """The file or an alignment in it has to be read by the general reader."""


### results of the compiled function
_NEED_MORE, _DONE, _UNSUPPORTED, _CAPACITY = 0, 1, 2, 3

### the bin of the BAM index that holds the offsets of the first and the last alignment of a reference
_PSEUDO_BIN = 37450

### decompressed bytes handed to the compiled function at once
_WINDOW = 24 * 2 ** 20

### aligned blocks and introns the output arrays can hold at first
_INITIAL_CAPACITY = 2 ** 20


@jit_cached
def _u16(buf, p):
    return np.int64(buf[p]) | (np.int64(buf[p + 1]) << 8)


@jit_cached
def _u32(buf, p):
    return np.int64(buf[p]) | (np.int64(buf[p + 1]) << 8) | (np.int64(buf[p + 2]) << 16) | (np.int64(buf[p + 3]) << 24)


@jit_cached
def _i32(buf, p):
    v = _u32(buf, p)
    if v >= 2147483648:
        v -= 4294967296
    return v


@jit_cached
def _u64(buf, p):
    return _u32(buf, p) | (_u32(buf, p + 4) << 32)


@jit_cached
def _bai_ranges(buf, n_ref):
    """(ok, first, last) virtual offsets of the alignments of every reference: from the pseudo bin if the index has one, else from all chunks."""

    first = np.full(n_ref, -1, dtype=np.int64)
    last = np.full(n_ref, -1, dtype=np.int64)
    n = buf.shape[0]
    pos = 8
    for r in range(n_ref):
        if pos + 4 > n:
            return (False, first, last)
        n_bin = _i32(buf, pos)
        pos += 4
        lo = np.int64(-1)
        hi = np.int64(-1)
        pseudo_lo = np.int64(-1)
        pseudo_hi = np.int64(-1)
        for b in range(n_bin):
            if pos + 8 > n:
                return (False, first, last)
            bin_id = _u32(buf, pos)
            n_chunk = _i32(buf, pos + 4)
            pos += 8
            if n_chunk < 0 or pos + 16 * n_chunk > n:
                return (False, first, last)
            if bin_id == 37450:
                if n_chunk >= 1:
                    pseudo_lo = _u64(buf, pos)
                    pseudo_hi = _u64(buf, pos + 8)
            else:
                for c in range(n_chunk):
                    chunk_lo = _u64(buf, pos + 16 * c)
                    chunk_hi = _u64(buf, pos + 16 * c + 8)
                    if lo < 0 or chunk_lo < lo:
                        lo = chunk_lo
                    if chunk_hi > hi:
                        hi = chunk_hi
            pos += 16 * n_chunk
        if pos + 4 > n:
            return (False, first, last)
        n_intv = _i32(buf, pos)
        pos += 4
        if n_intv < 0 or pos + 8 * n_intv > n:
            return (False, first, last)
        pos += 8 * n_intv
        if pseudo_lo >= 0 and pseudo_hi >= 0:
            first[r] = pseudo_lo
            last[r] = pseudo_hi
        elif lo >= 0:
            first[r] = lo
            last[r] = hi
    return (True, first, last)


_INDEX_RANGES = {}


def _index_ranges(fname, n_ref):
    """Virtual offsets (first, end) of the alignments of every reference of the BAM file, from its BAI index."""

    for candidate in (fname + '.bai', os.path.splitext(fname)[0] + '.bai'):
        if os.path.isfile(candidate):
            break
    else:
        raise NotSupported('no BAI index')
    st = os.stat(candidate)
    key = (candidate, st.st_mtime_ns, st.st_size)
    if key not in _INDEX_RANGES:
        with open(candidate, 'rb') as fd:
            data = np.frombuffer(fd.read(), dtype=np.uint8)
        if data.shape[0] < 8 or bytes(data[:4]) != b'BAI\x01' or int(data[4:8].view('<i4')[0]) != n_ref:
            raise NotSupported('unusual BAI index')
        ok, first, last = _bai_ranges(data, n_ref)
        if not ok:
            raise NotSupported('unusual BAI index')
        _INDEX_RANGES.clear()
        _INDEX_RANGES[key] = (first, last)
    return _INDEX_RANGES[key]


@jit_cached
def _find_tag(buf, pos, end, t0, t1):
    """Index of the type byte of the first auxiliary field with the tag t0 t1; -1 if there is none, -2 if the fields are malformed."""

    while pos + 3 <= end:
        typ = buf[pos + 2]
        if buf[pos] == t0 and buf[pos + 1] == t1:
            ### the value of a number or a character has to be there completely
            need = 0
            if typ == 65 or typ == 99 or typ == 67:
                need = 1
            elif typ == 115 or typ == 83:
                need = 2
            elif typ == 105 or typ == 73 or typ == 102:
                need = 4
            if pos + 3 + need > end:
                return -2
            return pos + 2
        pos += 3
        if typ == 65 or typ == 99 or typ == 67:                     # A c C
            pos += 1
        elif typ == 115 or typ == 83:                               # s S
            pos += 2
        elif typ == 105 or typ == 73 or typ == 102:                 # i I f
            pos += 4
        elif typ == 90 or typ == 72:                                # Z H
            while pos < end and buf[pos] != 0:
                pos += 1
            pos += 1
        elif typ == 66:                                             # B
            if pos + 5 > end:
                return -2
            sub = buf[pos]
            count = _u32(buf, pos + 1)
            if sub == 99 or sub == 67:
                size = 1
            elif sub == 115 or sub == 83:
                size = 2
            elif sub == 105 or sub == 73 or sub == 102:
                size = 4
            else:
                return -2
            pos += 5 + count * size
        else:
            return -2
    if pos == end:
        return -1
    return -2


@jit_cached
def _scan(buf, tid, started, primary_only, check_nm, mismatch, check_len, exon_len, t0, t1,
          cls_out, start_out, stop_out, key_out, n_blk, n_key):
    """Decodes the complete alignment records at the start of buf.

    Mirrors the loop of evidence.build_contig_evidence: the aligned blocks (class by XS tag, start, stop) and the
    introns (class, start, stop packed into one number) of the alignments that pass the filters are appended to the
    output arrays.  Returns (status, bytes used, blocks, introns, started)."""

    n = buf.shape[0]
    pos = 0
    cap_blk = cls_out.shape[0]
    cap_key = key_out.shape[0]
    while True:
        if pos + 4 > n:
            return (_NEED_MORE, pos, n_blk, n_key, started)
        block_size = _i32(buf, pos)
        ### (no alignment record is anywhere near that large: a damaged file must not make us read the rest of the file)
        if block_size < 32 or block_size > 67108864:
            return (_UNSUPPORTED, pos, n_blk, n_key, started)
        end = pos + 4 + block_size
        if end > n:
            return (_NEED_MORE, pos, n_blk, n_key, started)

        rec = pos + 4
        refid = _i32(buf, rec)
        if refid != tid:
            ### records are sorted by reference: before the contig, skip; after it, done
            if started or refid < 0 or refid > tid:
                return (_DONE, pos, n_blk, n_key, started)
            pos = end
            continue
        started = True

        flag = _u16(buf, rec + 14)
        if flag & 4:
            pos = end
            continue
        if primary_only and (flag & 256):
            pos = end
            continue

        l_name = np.int64(buf[rec + 8])
        n_cigar = _u16(buf, rec + 12)
        l_seq = _i32(buf, rec + 16)
        cigar = rec + 32 + l_name
        aux = cigar + 4 * n_cigar + (l_seq + 1) // 2 + l_seq
        if l_seq < 0 or aux > end:
            return (_UNSUPPORTED, pos, n_blk, n_key, started)
        if n_cigar == 0:
            return (_UNSUPPORTED, pos, n_blk, n_key, started)
        ### an alignment with more than 65535 operations has `<length of read>S<length on reference>N` here and the operations
        ### in the CG tag, where the general reader (htslib) finds them
        if n_cigar == 2 and (_u32(buf, cigar) & 15) == 4 and (_u32(buf, cigar + 4) & 15) == 3:
            if _find_tag(buf, aux, end, 67, 71) != -1:                  # CG
                return (_UNSUPPORTED, pos, n_blk, n_key, started)

        is_spliced = False
        query_len = np.int64(0)
        for k in range(n_cigar):
            v = _u32(buf, cigar + 4 * k)
            op = v & 15
            if op == 3:
                is_spliced = True
            elif op == 0 or op == 1 or op == 4 or op == 7 or op == 8:
                query_len += v >> 4
        ### htslib refuses an alignment whose CIGAR and sequence have different lengths (a stored sequence is optional)
        if l_seq > 0 and query_len != l_seq:
            return (_UNSUPPORTED, pos, n_blk, n_key, started)

        if check_nm:
            at = _find_tag(buf, aux, end, t0, t1)
            if at < 0:
                return (_UNSUPPORTED, pos, n_blk, n_key, started)
            typ = buf[at]
            if typ == 99:                       # c
                nm = np.int64(buf[at + 1])
                if nm >= 128:
                    nm -= 256
            elif typ == 67:                     # C
                nm = np.int64(buf[at + 1])
            elif typ == 115:                    # s
                nm = _u16(buf, at + 1)
                if nm >= 32768:
                    nm -= 65536
            elif typ == 83:                     # S
                nm = _u16(buf, at + 1)
            elif typ == 105:                    # i
                nm = _i32(buf, at + 1)
            elif typ == 73:                     # I
                nm = _u32(buf, at + 1)
            else:
                return (_UNSUPPORTED, pos, n_blk, n_key, started)
            if mismatch < nm:
                pos = end
                continue

        if check_len and is_spliced:
            ### min segment length: aligned length (M, D, =, X) of the shortest segment between introns
            seg = 0
            min_seg = -1
            for k in range(n_cigar):
                v = _u32(buf, cigar + 4 * k)
                op = v & 15
                if op == 3:
                    if min_seg < 0 or seg < min_seg:
                        min_seg = seg
                    seg = 0
                elif op == 0 or op == 2 or op == 7 or op == 8:
                    seg += v >> 4
            if min_seg < 0 or seg < min_seg:
                min_seg = seg
            if min_seg <= exon_len:
                pos = end
                continue

        cls = 0
        at = _find_tag(buf, aux, end, 88, 83)           # XS
        if at == -2:
            return (_UNSUPPORTED, pos, n_blk, n_key, started)
        if at >= 0:
            if buf[at] != 65:                           # A
                return (_UNSUPPORTED, pos, n_blk, n_key, started)
            if buf[at + 1] == 43:                       # +
                cls = 1
            elif buf[at + 1] == 45:                     # -
                cls = 2
            else:
                cls = 3

        if n_blk + n_cigar > cap_blk or n_key + n_cigar > cap_key:
            return (_CAPACITY, pos, n_blk, n_key, started)

        p = _i32(buf, rec + 4)
        for k in range(n_cigar):
            v = _u32(buf, cigar + 4 * k)
            op = v & 15
            ln = v >> 4
            if op == 3:
                if p < 0 or p + ln >= 2147483648:
                    return (_UNSUPPORTED, pos, n_blk, n_key, started)
                key_out[n_key] = (np.uint64(cls) << np.uint64(62)) | (np.uint64(p) << np.uint64(31)) | np.uint64(p + ln)
                n_key += 1
                p += ln
            elif op == 0 or op == 2 or op == 7 or op == 8:
                cls_out[n_blk] = cls
                start_out[n_blk] = p
                stop_out[n_blk] = p + ln
                n_blk += 1
                p += ln
        pos = end


def _bgzf_blocks(fd, offset, last_offset, file_size):
    """Decompressed data of the BGZF blocks from the block at `offset` to the block at `last_offset` (inclusive)."""

    data = b''
    pos = 0
    next_read = offset
    while True:
        if len(data) - pos < 18:
            if next_read >= file_size:
                if len(data) - pos > 0:
                    raise NotSupported('truncated BGZF block')
                return
            data = data[pos:] + os.pread(fd, 4 * 2 ** 20, next_read)
            next_read += 4 * 2 ** 20
            pos = 0
            if len(data) < 18:
                raise NotSupported('truncated BGZF block')
        if data[pos] != 31 or data[pos + 1] != 139 or data[pos + 3] & 4 == 0:
            raise NotSupported('not a BGZF block')
        xlen = data[pos + 10] | (data[pos + 11] << 8)
        ### the size of the block is in the 'BC' subfield
        bsize = -1
        sub = pos + 12
        while sub < pos + 12 + xlen:
            if sub + 4 > len(data):
                break
            slen = data[sub + 2] | (data[sub + 3] << 8)
            if data[sub] == 66 and data[sub + 1] == 67 and slen == 2:
                bsize = (data[sub + 4] | (data[sub + 5] << 8)) + 1
            sub += 4 + slen
        if bsize < 0:
            if pos + 12 + xlen > len(data):
                data = data[pos:] + os.pread(fd, 4 * 2 ** 20, next_read)
                next_read += 4 * 2 ** 20
                pos = 0
                continue
            raise NotSupported('not a BGZF block')
        if pos + bsize > len(data):
            if next_read >= file_size:
                raise NotSupported('truncated BGZF block')
            data = data[pos:] + os.pread(fd, 4 * 2 ** 20, next_read)
            next_read += 4 * 2 ** 20
            pos = 0
            continue
        crc, isize = struct.unpack_from('<II', data, pos + bsize - 8)
        if isize > 65536:
            raise NotSupported('unusual BGZF block')
        raw = zlib.decompress(data[pos + 12 + xlen:pos + bsize - 8], -15, isize) if isize else b''
        ### (the general reader does not accept a damaged block either)
        if len(raw) != isize or zlib.crc32(raw) != crc:
            raise NotSupported('damaged BGZF block')
        block_offset = offset
        yield block_offset, raw
        offset += bsize
        pos += bsize
        if block_offset >= last_offset:
            return


def scan_contig(fname, tid, n_ref, primary_only, mm_tag, mismatch, exon_len):
    """Aligned blocks and introns of the alignments of reference `tid` that pass the filters.

    Returns (blocks, introns): blocks[c] = (starts, stops) of the aligned blocks of the alignments of class c (no XS tag, XS +, XS -,
    any other XS), introns[c] = sorted array of the introns (start, stop, number of alignments) of class c.  Raises NotSupported."""

    if not fname.lower().endswith('.bam'):
        raise NotSupported('not a BAM file')
    check_nm = mismatch is not None
    check_len = exon_len is not None
    for value in (mismatch, exon_len):
        if value is not None and (isinstance(value, bool) or not isinstance(value, (int, np.integer))):
            raise NotSupported('unusual filter values')
    if check_nm and (not isinstance(mm_tag, str) or len(mm_tag) != 2 or not mm_tag.isascii()):
        raise NotSupported('unusual mismatch tag')
    t0, t1 = (ord(mm_tag[0]), ord(mm_tag[1])) if check_nm else (0, 0)

    first, last = _index_ranges(fname, n_ref)
    if tid < 0 or tid >= n_ref:
        raise NotSupported('unknown reference')
    if first[tid] < 0:
        ### no alignment on this reference
        return ([(np.zeros(0, dtype='int64'),) * 2 for _ in range(4)], [np.zeros((0, 3), dtype='uint32') for _ in range(4)])
    begin_offset, end_offset = int(first[tid]), int(last[tid])

    cls_out = np.empty(_INITIAL_CAPACITY, dtype=np.int8)
    start_out = np.empty(_INITIAL_CAPACITY, dtype=np.int64)
    stop_out = np.empty(_INITIAL_CAPACITY, dtype=np.int64)
    key_out = np.empty(_INITIAL_CAPACITY, dtype=np.uint64)
    n_blk = 0
    n_key = 0
    started = False
    done = False
    buf = bytearray()

    fd = os.open(fname, os.O_RDONLY)
    try:
        file_size = os.fstat(fd).st_size
        skip = begin_offset & 0xffff
        pending_skip = skip
        reached_end = False
        for block_offset, raw in _bgzf_blocks(fd, begin_offset >> 16, end_offset >> 16, file_size):
            reached_end = block_offset >= (end_offset >> 16)
            if pending_skip:
                raw = raw[pending_skip:]
                pending_skip = 0
            buf += raw
            if len(buf) < _WINDOW and block_offset < (end_offset >> 16):
                continue
            while True:
                view = np.frombuffer(buf, dtype=np.uint8)
                (status, used, n_blk, n_key, started) = _scan(view, tid, started, bool(primary_only), check_nm, int(mismatch) if check_nm else 0, check_len,
                                                              int(exon_len) if check_len else 0, t0, t1, cls_out, start_out, stop_out, key_out, n_blk, n_key)
                del view
                del buf[:used]
                if status == _CAPACITY:
                    new_capacity = 2 * cls_out.shape[0]
                    cls_out = np.resize(cls_out, new_capacity)
                    start_out = np.resize(start_out, new_capacity)
                    stop_out = np.resize(stop_out, new_capacity)
                    key_out = np.resize(key_out, new_capacity)
                    continue
                if status == _UNSUPPORTED:
                    raise NotSupported('unusual alignment record')
                done = status == _DONE
                break
            if done:
                break
        else:
            ### all blocks used: the file has to reach the end of the contig, and nothing may be left of an incomplete record
            if not reached_end:
                raise NotSupported('the file ends before the alignments of the contig')
            if len(buf) > 0 and not done:
                raise NotSupported('truncated alignment record')
    finally:
        os.close(fd)

    cls_out = cls_out[:n_blk]
    start_out = start_out[:n_blk]
    stop_out = stop_out[:n_blk]
    blocks = []
    for c in range(4):
        sel = cls_out == c
        blocks.append((start_out[sel], stop_out[sel]))
    keys, counts = np.unique(key_out[:n_key], return_counts=True)
    mask = np.uint64(2 ** 31 - 1)
    key_class = (keys >> np.uint64(62)).astype(np.int64)
    introns = []
    for c in range(4):
        sel = key_class == c
        k = keys[sel]
        arr = np.empty((k.shape[0], 3), dtype='uint32')
        arr[:, 0] = (k >> np.uint64(31)) & mask
        arr[:, 1] = k & mask
        arr[:, 2] = counts[sel]
        introns.append(arr)
    return blocks, introns
