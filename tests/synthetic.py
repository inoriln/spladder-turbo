"""Deterministic synthetic alignments used by the equivalence tests.

The generated BAM files deliberately contain every kind of read that the read
filters of SplAdder have to deal with (soft clips, insertions, deletions,
spliced reads with several / very short segments, secondary, supplementary and
unmapped records, all flavours of the XS tag, mismatches, ...) so that the
optimized read access paths can be compared against the original
implementation on a non-trivial input.
"""

import random

import pysam

CONTIGS = [('chrA', 20000), ('chrB', 8000), ('MT', 3000), ('chrE', 5000)]

### exon structures [(start, stop), ...] (0-based, half open) and strand per gene
GENES = {
    'chrA': [('+', [(1000, 1300), (2000, 2200), (3200, 3500), (5000, 5400)]),
             ('-', [(7000, 7300), (8500, 8800), (9500, 9900)]),
             ('+', [(12000, 12400), (13500, 13700), (15000, 15500), (17000, 17300)]),
             ### very short introns (30 and 40 bp): window ends close to the intron ends matter here
             ('+', [(18000, 18100), (18130, 18230), (18270, 18400)])],
    'chrB': [('-', [(500, 800), (1500, 1700), (3000, 3300), (4500, 4900), (6000, 6400)])],
}


def _query_len(cigar):
    return sum(ln for op, ln in cigar if op in (0, 1, 4, 7, 8))


def _ref_len(cigar):
    return sum(ln for op, ln in cigar if op in (0, 2, 3, 7, 8))


def _spliced_read(rng, exons):
    """Sample a read crossing one (sometimes two) junctions of an isoform; flanks are often very short."""

    iso = [exons[0]] + [e for e in exons[1:-1] if rng.random() < 0.7] + [exons[-1]]
    j = rng.randrange(len(iso) - 1)
    left = min(rng.choice([rng.randint(2, 12), rng.randint(5, 60)]), iso[j][1] - iso[j][0])
    right_exon = iso[j + 1][1] - iso[j + 1][0]
    right = min(rng.choice([rng.randint(2, 12), rng.randint(5, 60), right_exon]), right_exon)
    pos = iso[j][1] - left
    cigar = [(0, left), (3, iso[j + 1][0] - iso[j][1]), (0, right)]
    ### cross a second junction if the read covers the complete middle exon
    if right == right_exon and j + 2 < len(iso) and rng.random() < 0.5:
        tail = min(rng.randint(2, 40), iso[j + 2][1] - iso[j + 2][0])
        cigar += [(3, iso[j + 2][0] - iso[j + 1][1]), (0, tail)]
    return pos, cigar


def _unspliced_read(rng, exons):
    a, b = rng.choice(exons)
    ln = rng.randint(30, min(100, b - a))
    return rng.randint(a, b - ln), [(0, ln)]


def _decorate(rng, cigar):
    """Add soft clips, deletions and insertions to a CIGAR."""

    cigar = list(cigar)
    if rng.random() < 0.15:
        cigar.insert(0, (4, rng.randint(1, 5)))
    if rng.random() < 0.15:
        cigar.append((4, rng.randint(1, 5)))
    if rng.random() < 0.08:
        ### split the first M block by a deletion
        for i, (op, ln) in enumerate(cigar):
            if op == 0 and ln > 12:
                cigar[i:i + 1] = [(0, ln // 2), (2, rng.randint(1, 3)), (0, ln - ln // 2)]
                break
    if rng.random() < 0.08:
        for i, (op, ln) in enumerate(cigar):
            if op == 0 and ln > 12:
                cigar[i:i + 1] = [(0, ln // 2), (1, rng.randint(1, 3)), (0, ln - ln // 2)]
                break
    return cigar


def make_bam(path, seed=1, n_reads=4000, with_nm=True, nm_missing_region=None, contig_choices=None):
    """Write a sorted + indexed BAM to `path`.

    with_nm             -- write the NM tag to all mapped reads
    nm_missing_region   -- (contig, start, stop): reads starting inside lack the NM tag
    contig_choices      -- the contigs the reads are drawn from (repeating a contig makes it more frequent)
    """

    rng = random.Random(seed)
    header = {'HD': {'VN': '1.4', 'SO': 'unsorted'},
              'SQ': [{'SN': n, 'LN': ln} for n, ln in CONTIGS]}
    lens = dict(CONTIGS)
    tids = {n: i for i, (n, _) in enumerate(CONTIGS)}
    tmp = str(path) + '.unsorted.bam'

    with pysam.AlignmentFile(tmp, 'wb', header=header) as out:
        for i in range(n_reads):
            r = rng.random()
            contig = rng.choice(contig_choices or ['chrA', 'chrA', 'chrA', 'chrB', 'chrB', 'MT', 'chrE'])
            strand = None
            if contig in GENES and r < 0.85:
                strand, exons = rng.choice(GENES[contig])
                if rng.random() < 0.65:
                    pos, cigar = _spliced_read(rng, exons)
                    ### a fraction of spliced reads gets a novel intron
                    if rng.random() < 0.05 and len(cigar) == 3:
                        cigar = [cigar[0], (3, cigar[1][1] + rng.randint(-20, 20)), cigar[2]]
                else:
                    pos, cigar = _unspliced_read(rng, exons)
            else:
                ln = rng.randint(30, 100)
                pos = rng.choice([0, 1, rng.randint(0, lens[contig] - ln), lens[contig] - ln])
                cigar = [(0, ln)]
            cigar = _decorate(rng, cigar)
            if pos + _ref_len(cigar) > lens[contig] or pos < 0:
                continue

            a = pysam.AlignedSegment(out.header)
            a.query_name = 'r%d' % i
            a.reference_id = tids[contig]
            a.reference_start = pos
            a.mapping_quality = 255
            u = rng.random()
            if u < 0.02:
                ### unmapped record placed at a position (as produced by some aligners)
                a.flag = 4
                a.query_sequence = 'A' * 30
                out.write(a)
                continue
            flag = rng.choice([0, 16])
            if u < 0.07:
                flag |= 256
            elif u < 0.09:
                flag |= 2048
            a.flag = flag
            a.cigartuples = cigar
            a.query_sequence = 'A' * _query_len(cigar)
            if with_nm and not (nm_missing_region and contig == nm_missing_region[0] and nm_missing_region[1] <= pos < nm_missing_region[2]):
                a.set_tag('NM', rng.choice([0, 0, 0, 0, 0, 1, 1, 2, 3]))
            spliced = any(op == 3 for op, _ in cigar)
            v = rng.random()
            if strand is None:
                strand = rng.choice(['+', '-'])
            if spliced:
                if v < 0.88:
                    a.set_tag('XS', strand, value_type='A')
                elif v < 0.95:
                    a.set_tag('XS', '+' if strand == '-' else '-', value_type='A')
                elif v < 0.97:
                    a.set_tag('XS', '.', value_type='A')
            elif v < 0.2:
                a.set_tag('XS', strand, value_type='A')
            out.write(a)

    pysam.sort('-o', str(path), tmp)
    pysam.index(str(path))
    return str(path)
