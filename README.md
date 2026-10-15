# spladder-turbo

**An unofficial, performance-optimized fork of [SplAdder](https://github.com/ratschlab/spladder) 3.1.1** for the detection and quantification of alternative splicing events from RNA-Seq alignments.

spladder-turbo is meant as a drop-in replacement: the same command line, the same output files and, as far as we have been able to test, the same results, in a fraction of the time. Like SplAdder it builds the splicing graph of an annotation, augments it with the evidence of the alignments, extracts alternative splicing events from the graph (exon skips, multiple exon skips, intron retentions, alternative 3' and 5' splice sites, mutually exclusive exons) and quantifies them.

> **Status: development version 0.1.0, not yet released on PyPI.**
> spladder-turbo is **not affiliated with or endorsed by** the SplAdder authors or their institutions. The method, the definitions of the events and the file formats are theirs. If you use this software, please cite the SplAdder paper (see [Citation](#citation)).

## Performance

The whole workflow of `spladder build` (graph augmentation, merging, quantification, extraction and verification of events), two human RNA-Seq samples ([GSE272028](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272028), runs SRR29779619 and SRR29779620; aligned with STAR 2.7.7a to hg38, 85 bp reads, about 44 million alignments each), GENCODE v34, whole genome, starting from the GTF file:

| | SplAdder 3.1.1 | spladder-turbo 0.1.0 | |
|---|---|---|---|
| Wall-clock time | 3 h 04 min | 4 min 28 s (`--parallel 8`) | 41 x |
| CPU time | 15,650 s | 1,425 s | 11 x |
| Memory of the largest process | 2.9 GB | 3.0 GB | |

The same workflow on chr20 only, three samples (SRR29779619 to SRR29779621):

| | SplAdder 3.1.1 | spladder-turbo 0.1.0 | |
|---|---|---|---|
| One process | 5 min 04 s | 40 s | 7.5 x |
| `--parallel 8` | (not measured) | 23 s | 13 x, compared with one process of SplAdder |

How to read these numbers:

- Measured on one Linux server with two AMD EPYC 7601 processors (128 threads) and 314 GB of memory, with Python 3.13, numpy 2.3, numba 0.62 and pysam 0.23. Times are from `/usr/bin/time -v`.
- SplAdder was run in the way that is most favourable to it: with single-threaded linear algebra libraries, and for the whole genome with the two per-sample graph steps started as two simultaneous processes. SplAdder's own `--parallel` only distributes the detection and verification of events. We have **not yet measured** SplAdder with `--parallel 8` on the whole genome.
- The speed-up of one process (7.5 x) comes from reading each contig of an alignment file once with a compiled reader, from compiled kernels and from vectorized code. The rest comes from running samples, contigs and event types in parallel processes (Linux only). The CPU time shows how much less work is done in total.
- The first run compiles numba functions, which takes about half a minute (68 s instead of 40 s for the chr20 run above). The compiled code is cached on disk.
- With `--parallel 8` the whole-genome run uses about 9 GB in total, in up to about 24 processes at the same time.
- The scripts and the raw logs of these measurements are not yet part of this repository (see [docs/PERFORMANCE_CHANGES.md](docs/PERFORMANCE_CHANGES.md) for the method and the single measurements).

## Are the results the same?

The goal is that every output of `spladder build` has the same content as with SplAdder 3.1.1. What has been compared so far:

- **Whole genome** (the run above): all 38 output files (text, GFF3, HDF5, event files, the splicing graphs of the samples and the merged graph) and the annotation file. 29 of the 38 files are byte-identical (GFF3, HDF5, event and index files), and so is the annotation file. The six compressed text files (`.txt.gz`) differ only in the time stamp in the gzip header; their text is identical. The three graph files (`genes_graph_*.pickle`) differ in their bytes, but they load to completely equal objects.
- **42 command line configurations** (confidence levels, merge strategies, event types, output formats, sparse BAM files, ...) on the chr20 subset, each run with both programs: identical, or failing in the same way (two configurations fail in SplAdder as well).
- **More than 50 randomly generated data sets** (random annotations and reads with novel exons, retained introns, alternative sites and mutually exclusive exons), run through both programs, serial and parallel: identical.
- **The compiled BAM reader** against pysam on every contig of two real samples (678 combinations of contig and filter) and on 800 random BAM files: identical.
- **250 tests** in this repository compare the optimized functions with verbatim copies of the SplAdder 3.1.1 implementations (`tests/reference.py`) on random and adversarial input.

Known differences:

- spladder-turbo fixes defects of SplAdder 3.1.1 (crashes for annotations without edges or with unknown strands, merging samples with different gene sets, undefined names in rarely used code, stale index files, ...). Where SplAdder worked, the output is the same; where it crashed, spladder-turbo continues. The list is in [docs/PERFORMANCE_CHANGES.md](docs/PERFORMANCE_CHANGES.md) (section 9).
- Peculiarities of SplAdder that change the output are kept, not fixed (for example in `Gene.label_alt`).
- With `--parallel`, the output of the parallel tasks is shown one task after the other, and `--logfile` is written per sample.
- Parallel processing needs Linux; elsewhere everything runs serially (slower, but the same results).

Not validated yet: paired-end data, aligners other than STAR, species other than human, cohorts of more than three samples, CRAM files, `--var-aware`, `--sparse-bam` for a whole genome, and operating systems other than Linux. The code is written to behave like SplAdder in these cases (anything that is not understood falls back to the code of SplAdder), but this has only been tested on small synthetic data. Please compare the outputs of both programs on a small part of your own data before you switch (see [Switching from SplAdder](#switching-from-spladder)), and report differences as issues.

## Installation

spladder-turbo needs Python 3.9 or newer and the packages in [requirements.txt](requirements.txt). We recommend a package manager such as conda to organize the dependencies.

```bash
pip install git+https://github.com/inoriln/spladder-turbo.git
```

or, from a clone of this repository:

```bash
pip install .
```

**Do not install spladder-turbo into an environment that contains SplAdder.** Both provide the Python package `spladder` and the command `spladder`, so one overwrites the files of the other. Use a separate environment. The command `spladder-turbo` is the same program as `spladder`.

## Usage

The command line is that of SplAdder; see the [SplAdder documentation](https://spladder.readthedocs.io) for all options and for the description of the output files. A minimal run:

```bash
spladder build -a annotation.gtf -o results -b sample1.bam,sample2.bam --parallel 8
```

`--parallel N` uses N processes for the graphs of the samples, for the counting, and for the collection and verification of events. The number of processes can briefly exceed N. The results do not depend on N. The alignment files must be sorted and indexed (`.bai`) for the fast reader; other files are read, more slowly, with pysam.

Environment variables:

| Variable | Effect |
|---|---|
| `SPLADDER_EVIDENCE_MAX_GB` | memory limit for the alignment information kept per process, default 6 |
| `SPLADDER_NO_BAMSCAN=1` | read the alignments with pysam instead of the compiled reader |
| `SPLADDER_NO_EVIDENCE=1` | read the alignments gene by gene, as SplAdder does |

### Switching from SplAdder

Run both programs, each in its own environment, on a small part of your data (for example one chromosome of a few samples). Use the same options for both runs, and give each run its own output directory and its own copy of the annotation file (the annotation is converted to a `.pickle` file next to it). Then compare the two output directories with the function that we use ([tests/compare.py](tests/compare.py)): it compares text files, HDF5 datasets and the objects in pickle files, and returns the list of differences. In a clone of this repository, in the environment of spladder-turbo:

```bash
python -c "from tests.compare import compare_dirs; print(compare_dirs('run_spladder', 'run_turbo') or 'identical')"
```

## Tests

```bash
pip install pytest
pytest tests -q
```

takes about three minutes. With numpy 1.x, 6 tests of `tests/test_end_to_end.py` fail because the expected results of the SplAdder test data were written with numpy 2; they fail in the same way with SplAdder.

## Documentation

- [SplAdder documentation](https://spladder.readthedocs.io) for the command line and the output files.
- [docs/PERFORMANCE_CHANGES.md](docs/PERFORMANCE_CHANGES.md): what was changed and why, how the results are kept identical, and the record of the code review (in Chinese for now).
- [ATTRIBUTION.md](ATTRIBUTION.md): credits and licensing.

## Citation

spladder-turbo is a derivative of SplAdder. Please cite the SplAdder paper:

Kahles A, Ong CS, Zhong Y, Rätsch G. SplAdder: identification, quantification and testing of alternative splicing events from RNA-Seq data. *Bioinformatics* 2016;32(12):1840-1847. doi:[10.1093/bioinformatics/btw076](https://doi.org/10.1093/bioinformatics/btw076)

For reproducibility, please also state that you used spladder-turbo and its version (`spladder --help` shows it).

## License and attribution

SplAdder was written by Cheng Soon Ong, Gunnar Raetsch and Andre Kahles (see [AUTHORS](AUTHORS)) and is distributed under the BSD 3-Clause license ([COPYRIGHT](COPYRIGHT)). The modifications of this fork are distributed under the same license ([LICENSE](LICENSE)). The names of the copyright holders and contributors are not used to endorse or promote this fork. [ATTRIBUTION.md](ATTRIBUTION.md) has the details, including a note on one source file that carries a GPL notice upstream.

The software is provided "as is", without warranty of any kind; see the license for the exact terms.
