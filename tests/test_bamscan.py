"""The compiled BAM reader must give exactly the evidence the pysam loop gives, or hand over to it."""

import itertools
import os
import struct
import zlib

import numpy as np
import pysam
import pytest

from spladder import bamscan, evidence

from . import synthetic


def same(a, b):
    """Complete equality of two ContigEvidence objects (dtypes and shapes included)."""
    for (p1, v1), (p2, v2) in zip(a.cov, b.cov):
        if p1.dtype != p2.dtype or v1.dtype != v2.dtype or not np.array_equal(p1, p2) or not np.array_equal(v1, v2):
            return False
    for key in a.introns:
        for (x1, s1), (x2, s2) in zip(a.introns[key], b.introns[key]):
            if x1.dtype != x2.dtype or x1.shape != x2.shape or not np.array_equal(x1, x2) or s1.dtype != s2.dtype or not np.array_equal(s1, s2):
                return False
    return True


def both(fname, contig, spec, monkeypatch):
    """(evidence of the pysam loop, evidence of the default path); errors are compared as such."""
    results = []
    for disable in ('1', ''):
        monkeypatch.setenv(evidence.ENV_NO_SCAN, disable)
        monkeypatch.setenv(evidence.ENV_DEBUG_SCAN, '1')
        try:
            results.append(evidence.build_contig_evidence(fname, contig, spec))
        except evidence.Unsupported:
            results.append('unsupported')
    return results


SPECS = [evidence.Spec('', None, True, 2, 9, 'NM'), evidence.Spec('', None, False, 0, 20, 'NM'), evidence.Spec('', None, True, None, None, None),
         evidence.Spec('', None, False, 1, None, 'NM'), evidence.Spec('', None, True, None, 12, None)]


@pytest.fixture(scope='module')
def bam(tmp_path_factory):
    return synthetic.make_bam(tmp_path_factory.mktemp('scan') / 'syn.bam', seed=5, n_reads=12000)


@pytest.mark.parametrize('spec', SPECS)
@pytest.mark.parametrize('contig', ['chrA', 'chrB', 'chrE', 'MT'])
def test_scanner_matches_pysam_loop(bam, monkeypatch, spec, contig):
    spec = spec._replace(fname=bam)
    reference, got = both(bam, contig, spec, monkeypatch)
    assert reference != 'unsupported' and got != 'unsupported'
    assert same(reference, got)


@pytest.mark.parametrize('window,capacity', [(1, 4), (700, 64), (5000, 100000)])
def test_records_across_windows_and_growing_outputs(bam, monkeypatch, window, capacity):
    ### the records straddle the windows handed to the compiled function, the output arrays have to grow
    monkeypatch.setattr(bamscan, '_WINDOW', window)
    monkeypatch.setattr(bamscan, '_INITIAL_CAPACITY', capacity)
    for spec, contig in itertools.product(SPECS[:3], ['chrA', 'chrB']):
        spec = spec._replace(fname=bam)
        reference, got = both(bam, contig, spec, monkeypatch)
        assert same(reference, got)


def test_the_scanner_is_used(bam, monkeypatch):
    calls = []
    original = bamscan.scan_contig
    monkeypatch.setattr(bamscan, 'scan_contig', lambda *a, **k: calls.append(1) or original(*a, **k))
    monkeypatch.setenv(evidence.ENV_NO_SCAN, '')
    evidence.build_contig_evidence(bam, 'chrA', SPECS[0]._replace(fname=bam))
    assert calls == [1]
    monkeypatch.setenv(evidence.ENV_NO_SCAN, '1')
    evidence.build_contig_evidence(bam, 'chrA', SPECS[0]._replace(fname=bam))
    assert calls == [1]


def make_tagged_bam(path, nm_type=None, xs_type='A', with_nm=True, seed=3):
    """A small sorted BAM whose NM / XS tags have the requested types."""
    import random
    rng = random.Random(seed)
    header = {'HD': {'VN': '1.4', 'SO': 'coordinate'}, 'SQ': [{'SN': 'c1', 'LN': 5000}, {'SN': 'c2', 'LN': 5000}]}
    reads = []
    for i in range(400):
        tid = rng.choice([0, 1])
        pos = rng.randint(0, 3000)
        spliced = rng.random() < 0.5
        cigar = [(0, rng.randint(8, 40)), (3, rng.randint(50, 500)), (0, rng.randint(8, 40))] if spliced else [(0, rng.randint(30, 90))]
        reads.append((tid, pos, cigar, rng.choice([0, 1, 2, 3, 5, 9]), rng.choice(['+', '-', '.'])))
    reads.sort(key=lambda r: (r[0], r[1]))
    with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
        for i, (tid, pos, cigar, nm, xs) in enumerate(reads):
            a = pysam.AlignedSegment(out.header)
            a.query_name = 'r%d' % i
            a.reference_id = tid
            a.reference_start = pos
            a.mapping_quality = 60
            a.flag = rng.choice([0, 16])
            a.cigartuples = cigar
            a.query_sequence = 'A' * sum(ln for op, ln in cigar if op in (0, 1))
            if with_nm:
                a.set_tag('NM', nm, value_type=nm_type) if nm_type else a.set_tag('NM', nm)
            a.set_tag('NH', 1)
            if spliced_or(cigar):
                if xs_type == 'A':
                    a.set_tag('XS', xs, value_type='A')
                elif xs_type == 'Z':
                    a.set_tag('XS', xs, value_type='Z')
            out.write(a)
    pysam.index(str(path))
    return str(path)


def spliced_or(cigar):
    return any(op == 3 for op, _ in cigar)


@pytest.mark.parametrize('nm_type', [None, 'c', 'C', 's', 'S', 'i', 'I'])
def test_integer_types_of_the_mismatch_tag(tmp_path, monkeypatch, nm_type):
    fname = make_tagged_bam(tmp_path / 't.bam', nm_type=nm_type)
    for contig in ('c1', 'c2'):
        spec = evidence.Spec(fname, None, True, 2, 9, 'NM')
        reference, got = both(fname, contig, spec, monkeypatch)
        assert same(reference, got)


def test_tags_the_scanner_hands_over(tmp_path, monkeypatch):
    ### a string XS tag: the compiled reader gives up, the answer is the one of the pysam loop
    fname = make_tagged_bam(tmp_path / 'xs.bam', xs_type='Z')
    spec = evidence.Spec(fname, None, True, 2, 9, 'NM')
    monkeypatch.setenv(evidence.ENV_NO_SCAN, '')
    with pytest.raises(bamscan.NotSupported):
        bamscan.scan_contig(fname, 0, 2, True, 'NM', 2, 9)
    reference, got = both(fname, 'c1', spec, monkeypatch)
    assert same(reference, got)
    ### no NM tag at all: both fail the same way
    fname = make_tagged_bam(tmp_path / 'nonm.bam', with_nm=False)
    spec = evidence.Spec(fname, None, True, 2, 9, 'NM')
    reference, got = both(fname, 'c1', spec, monkeypatch)
    assert reference == got == 'unsupported'
    ### nothing to check: the missing tag does not matter
    spec = evidence.Spec(fname, None, True, None, 9, None)
    reference, got = both(fname, 'c1', spec, monkeypatch)
    assert reference != 'unsupported' and same(reference, got)


def test_files_the_scanner_does_not_handle(bam, tmp_path, monkeypatch):
    ### no BAI index: the scanner hands over (and so does the general reader, since there is nothing to fetch with)
    with pytest.raises(bamscan.NotSupported):
        bamscan.scan_contig(str(tmp_path / 'missing.bam'), 0, 4, True, 'NM', 2, 9)
    with pytest.raises(bamscan.NotSupported):
        bamscan.scan_contig(bam + '.bai', 0, 4, True, 'NM', 2, 9)
    ### unusual filter values
    with pytest.raises(bamscan.NotSupported):
        bamscan.scan_contig(bam, 0, 4, True, 'NM', 2.5, 9)
    with pytest.raises(bamscan.NotSupported):
        bamscan.scan_contig(bam, 0, 4, True, 'NMX', 2, 9)
    ### wrong number of references
    with pytest.raises(bamscan.NotSupported):
        bamscan.scan_contig(bam, 0, 7, True, 'NM', 2, 9)


def test_index_without_the_pseudo_bin(bam, tmp_path, monkeypatch):
    """An index that does not hold the offsets of the first and last alignment of a reference: they are taken from all the chunks."""
    data = open(bam + '.bai', 'rb').read()
    n_ref = struct.unpack_from('<i', data, 4)[0]
    pos = 8
    out = bytearray(data[:8])
    for _ in range(n_ref):
        n_bin = struct.unpack_from('<i', data, pos)[0]
        pos += 4
        bins = []
        for _ in range(n_bin):
            bin_id, n_chunk = struct.unpack_from('<Ii', data, pos)
            chunk_bytes = data[pos + 8:pos + 8 + 16 * n_chunk]
            pos += 8 + 16 * n_chunk
            if bin_id != bamscan._PSEUDO_BIN:
                bins.append(struct.pack('<Ii', bin_id, n_chunk) + chunk_bytes)
        n_intv = struct.unpack_from('<i', data, pos)[0]
        intv = data[pos + 4:pos + 4 + 8 * n_intv]
        pos += 4 + 8 * n_intv
        out += struct.pack('<i', len(bins)) + b''.join(bins) + struct.pack('<i', n_intv) + intv
    out += data[pos:]
    stripped = str(tmp_path / 'stripped.bam')
    os.symlink(bam, stripped)
    open(stripped + '.bai', 'wb').write(bytes(out))
    for contig in ('chrA', 'chrB', 'chrE'):
        spec = SPECS[0]._replace(fname=stripped)
        reference, got = both(stripped, contig, spec, monkeypatch)
        assert same(reference, got)


def test_truncated_file_is_not_read_partially(bam, tmp_path):
    """A file that ends in the middle of the alignments of a contig is handed over, whatever the cut position is."""
    data = open(bam, 'rb').read()
    n_ref = len(pysam.AlignmentFile(bam).references)
    for fraction in (0.3, 0.55, 0.8, 0.97):
        cut = str(tmp_path / ('cut%d.bam' % int(fraction * 100)))
        open(cut, 'wb').write(data[:int(len(data) * fraction)])
        open(cut + '.bai', 'wb').write(open(bam + '.bai', 'rb').read())
        outcomes = set()
        for tid in range(n_ref):
            try:
                bamscan.scan_contig(cut, tid, n_ref, True, 'NM', 2, 9)
                outcomes.add('ok')
            except bamscan.NotSupported:
                outcomes.add('handed over')
        assert 'handed over' in outcomes


def test_cigar_in_the_cg_tag_is_left_to_the_general_reader(tmp_path, monkeypatch):
    ### an alignment with more than 65535 operations: the record has `<length of read>S<length on reference>N`, the operations are in the CG tag
    import array
    n_pairs = 35000
    ops = array.array('I', [(10 << 4) | 0, (5 << 4) | 3] * n_pairs + [(10 << 4) | 0])
    read_length, span = 10 * (n_pairs + 1), 15 * n_pairs + 10
    header = {'HD': {'VN': '1.4', 'SO': 'coordinate'}, 'SQ': [{'SN': 'c1', 'LN': 2000000}]}
    fname = str(tmp_path / 'long.bam')
    with pysam.AlignmentFile(fname, 'wb', header=header) as out:
        for name, pos, cigar, cg in (('before', 50, [(0, 30), (3, 100), (0, 30)], None), ('long', 1000, [(4, read_length), (3, span)], ops),
                                     ('after', 700000, [(0, 30), (3, 100), (0, 30)], None)):
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.reference_id, a.reference_start, a.mapping_quality, a.flag = name, 0, pos, 60, 0
            a.cigartuples = cigar
            a.query_sequence = 'A' * sum(ln for op, ln in cigar if op in (0, 1, 4))
            a.set_tag('NM', 0)
            a.set_tag('XS', '+', value_type='A')
            if cg is not None:
                a.set_tag('CG', cg)
            out.write(a)
    pysam.index(fname)
    spec = evidence.Spec(fname, None, True, 2, 9, 'NM')
    monkeypatch.setenv(evidence.ENV_NO_SCAN, '')
    with pytest.raises(bamscan.NotSupported):
        bamscan.scan_contig(fname, 0, 1, True, 'NM', 2, 9)
    reference, got = both(fname, 'c1', spec, monkeypatch)
    assert reference != 'unsupported' and same(reference, got)
    ### the evidence has the introns of the long alignment (and not the one intron of the placeholder)
    assert sum(arr.shape[0] for arr, _ in reference.introns[None]) == 2 + n_pairs
    ### without a CG tag the same record is what it looks like: a soft clip and an intron
    plain = str(tmp_path / 'plain.bam')
    with pysam.AlignmentFile(plain, 'wb', header=header) as out:
        a = pysam.AlignedSegment(out.header)
        a.query_name, a.reference_id, a.reference_start, a.mapping_quality, a.flag = 'x', 0, 100, 60, 0
        a.cigartuples = [(4, 20), (3, 500)]
        a.query_sequence = 'A' * 20
        a.set_tag('NM', 0)
        a.set_tag('XS', '+', value_type='A')
        out.write(a)
    pysam.index(plain)
    reference, got = both(plain, 'c1', evidence.Spec(plain, None, True, 2, None, 'NM'), monkeypatch)
    assert reference != 'unsupported' and same(reference, got)
    assert sum(arr.shape[0] for arr, _ in got.introns[None]) == 1


def test_a_damaged_block_is_an_error_as_with_pysam(tmp_path, monkeypatch):
    ### (blocks without compression: the damage is not noticed when inflating, only by the checksum)
    header = {'HD': {'VN': '1.4', 'SO': 'coordinate'}, 'SQ': [{'SN': 'c1', 'LN': 100000}]}
    fname = str(tmp_path / 'stored.bam')
    with pysam.AlignmentFile(fname, 'wb0', header=header) as out:
        for i in range(3000):
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.reference_id, a.reference_start, a.mapping_quality, a.flag = 'r%d' % i, 0, 100 + 20 * i, 60, 0
            a.cigartuples = [(0, 30), (3, 100), (0, 30)]
            a.query_sequence = 'ACGT' * 15
            a.set_tag('NM', 0)
            a.set_tag('XS', '+', value_type='A')
            out.write(a)
    pysam.index(fname)
    original = bytearray(open(fname, 'rb').read())
    blocks, pos = [], 0
    while pos < len(original):
        size = struct.unpack_from('<H', original, pos + 16)[0] + 1
        blocks.append((pos, size))
        pos += size
    assert len(blocks) > 4
    start, size = blocks[len(blocks) // 2]

    def damaged(name, change):
        data = bytearray(original)
        change(data)
        path = str(tmp_path / name)
        open(path, 'wb').write(data)
        open(path + '.bai', 'wb').write(open(fname + '.bai', 'rb').read())
        return path

    flipped = damaged('flipped.bam', lambda d: d.__setitem__(start + 18 + 40, d[start + 18 + 40] ^ 4))
    monkeypatch.setenv(evidence.ENV_NO_SCAN, '')
    with pytest.raises(bamscan.NotSupported):
        bamscan.scan_contig(flipped, 0, 1, True, 'NM', 2, 9)
    ### the general reader gives the error of pysam
    for disable in ('1', ''):
        monkeypatch.setenv(evidence.ENV_NO_SCAN, disable)
        with pytest.raises(OSError):
            evidence.build_contig_evidence(flipped, 'c1', evidence.Spec(flipped, None, True, 2, 9, 'NM'))
    ### a size that cannot be right (pysam does not look at it): handed over, the answer is the one of pysam
    huge = damaged('huge.bam', lambda d: struct.pack_into('<I', d, start + size - 4, 70000))
    monkeypatch.setenv(evidence.ENV_NO_SCAN, '')
    with pytest.raises(bamscan.NotSupported):
        bamscan.scan_contig(huge, 0, 1, True, 'NM', 2, 9)
    reference, got = both(huge, 'c1', evidence.Spec(huge, None, True, 2, 9, 'NM'), monkeypatch)
    assert reference != 'unsupported' and same(reference, got)
    ### and the undamaged file is read by the scanner
    monkeypatch.setenv(evidence.ENV_NO_SCAN, '')
    bamscan.scan_contig(fname, 0, 1, True, 'NM', 2, 9)
    reference, got = both(fname, 'c1', evidence.Spec(fname, None, True, 2, 9, 'NM'), monkeypatch)
    assert reference != 'unsupported' and same(reference, got)


def test_cigar_of_another_length_than_the_sequence_is_an_error_as_with_pysam(tmp_path, monkeypatch):
    ### htslib refuses such an alignment ("error -4 while reading"); the compiled reader has to hand over and not read on
    header = {'HD': {'VN': '1.4', 'SO': 'coordinate'}, 'SQ': [{'SN': 'c1', 'LN': 100000}]}
    fname = str(tmp_path / 'good.bam')
    with pysam.AlignmentFile(fname, 'wb0', header=header) as out:
        for i in range(5):
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.reference_id, a.reference_start, a.mapping_quality, a.flag = 'r%d' % i, 0, 100 + 50 * i, 60, 0
            a.cigartuples = [(0, 30), (3, 100), (0, 30)]
            a.query_sequence = 'ACGT' * 15
            a.set_tag('NM', 0)
            a.set_tag('XS', '+', value_type='A')
            out.write(a)
    pysam.index(fname)
    raw = bytearray(open(fname, 'rb').read())
    blocks, pos = [], 0
    while pos < len(raw):
        size = struct.unpack_from('<H', raw, pos + 16)[0] + 1
        blocks.append((pos, size))
        pos += size
    ### (the header is in the first block, the alignments in the second; the blocks are not compressed)
    start, size = blocks[1]
    isize = struct.unpack_from('<I', raw, start + size - 4)[0]
    first = start + 18 + 5
    data = bytearray(raw[first:first + isize])
    cigar = 4 + 32 + data[4 + 8]
    assert struct.unpack_from('<I', data, cigar + 8)[0] == (30 << 4)
    struct.pack_into('<I', data, cigar + 8, 25 << 4)           # the last aligned block is 5 shorter than the sequence says
    raw[first:first + isize] = data
    struct.pack_into('<I', raw, start + size - 8, zlib.crc32(bytes(data)))
    bad = str(tmp_path / 'bad.bam')
    open(bad, 'wb').write(raw)
    open(bad + '.bai', 'wb').write(open(fname + '.bai', 'rb').read())

    spec = evidence.Spec(bad, None, True, 2, 9, 'NM')
    monkeypatch.setenv(evidence.ENV_NO_SCAN, '')
    with pytest.raises(bamscan.NotSupported):
        bamscan.scan_contig(bad, 0, 1, True, 'NM', 2, 9)
    for disable in ('1', ''):
        monkeypatch.setenv(evidence.ENV_NO_SCAN, disable)
        with pytest.raises(OSError):
            evidence.build_contig_evidence(bad, 'c1', spec)
    ### the file before the change is read by the scanner
    monkeypatch.setenv(evidence.ENV_NO_SCAN, '')
    bamscan.scan_contig(fname, 0, 1, True, 'NM', 2, 9)
    reference, got = both(fname, 'c1', evidence.Spec(fname, None, True, 2, 9, 'NM'), monkeypatch)
    assert reference != 'unsupported' and same(reference, got)
