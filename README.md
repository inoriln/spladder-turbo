# spladder-turbo

[English](README.md) | [简体中文](README.zh-CN.md)

**An unofficial, performance-optimized fork of [SplAdder](https://github.com/ratschlab/spladder) 3.1.1** for the detection and quantification of alternative splicing events from RNA-Seq alignments.

spladder-turbo is meant as a drop-in replacement: the same command line, the same output files and, as far as we have been able to test, the same results, in a fraction of the time. Like SplAdder it builds the splicing graph of an annotation, augments it with the evidence of the alignments, extracts alternative splicing events from the graph (exon skips, multiple exon skips, intron retentions, alternative 3' and 5' splice sites, mutually exclusive exons) and quantifies them.

> **Status: development version 0.1.0, installable from GitHub, not yet released on PyPI.**
> spladder-turbo is **not affiliated with or endorsed by** the SplAdder authors or their institutions. The method, the definitions of the events and the file formats are theirs. If you use this software, please cite the SplAdder paper (see [Citation](#citation)).

**Contents:** [Installation](#installation) · [Usage](#usage) · [Tests](#tests) · [Performance](#performance) · [What was optimized, where, and why](#what-was-optimized-where-and-why) · [Are the results the same?](#are-the-results-the-same) · [Documentation](#documentation) · [Citation](#citation) · [License and attribution](#license-and-attribution)

## Installation

**Recommended: Python 3.11 in a dedicated Conda or micromamba environment.** The package declares Python 3.9 or newer; Python 3.11 has been tested with a fresh installation from GitHub on Linux, including serial and parallel builds. Runtime dependencies in [requirements.txt](requirements.txt) are installed automatically by pip.

**Do not install spladder-turbo into an environment that already contains SplAdder.** Both provide the Python package `spladder` and the command `spladder`, so installing one can overwrite the other. The commands `spladder-turbo` and `spladder` run the same program.

### Conda

```bash
conda create -n spladder-turbo -c conda-forge python=3.11 pip git -y
conda activate spladder-turbo
python -m pip install git+https://github.com/inoriln/spladder-turbo.git
```

### micromamba

```bash
micromamba create -n spladder-turbo -c conda-forge python=3.11 pip git -y
micromamba activate spladder-turbo
python -m pip install git+https://github.com/inoriln/spladder-turbo.git
```

Use either manager to create the environment, then install the project from GitHub with pip. If your shell is not configured for micromamba activation, use `micromamba run -n spladder-turbo python -m pip install git+https://github.com/inoriln/spladder-turbo.git` instead.

After installation, check dependencies and the command entry point in the new environment:

```bash
python -m pip check
spladder-turbo build --help
```

For a local clone, run `python -m pip install .` from its root directory in the same dedicated environment.

## Usage

The command line is the same as SplAdder's; see the [SplAdder documentation](https://spladder.readthedocs.io) for all options and output files. A minimal example:

```bash
spladder build -a annotation.gtf -o results -b sample1.bam,sample2.bam --parallel 8
```

`--parallel N` uses N processes for the sample graphs, counting, and the collection and verification of events. The number of processes may briefly exceed N. Results do not depend on N. The fast reader requires sorted alignment files with a `.bai` index; other files can still be read with pysam, more slowly.

Environment variables:

| Variable | Effect |
|---|---|
| `SPLADDER_EVIDENCE_MAX_GB` | Maximum memory per process for alignment evidence, default 6 |
| `SPLADDER_NO_BAMSCAN=1` | Read alignments with pysam instead of the compiled reader |
| `SPLADDER_NO_EVIDENCE=1` | Read alignments per gene, as SplAdder does |

### Switching from SplAdder

Run both programs in separate environments on a small part of your own data, such as a few samples on one chromosome. Use the same options, separate output directories and separate annotation copies: the annotation is converted to a `.pickle` file next to the input. From a clone of this repository in the spladder-turbo environment, compare the outputs with [tests/compare.py](tests/compare.py), which checks text, HDF5 datasets and deserialized pickle objects:

```bash
python -c "from tests.compare import compare_dirs; print(compare_dirs('run_spladder', 'run_turbo') or 'identical')"
```

## Tests

The repository test suite requires a source checkout. In the dedicated environment created above:

```bash
git clone https://github.com/inoriln/spladder-turbo.git
cd "spladder-turbo"
python -m pip install . pytest
python -m pytest tests -q
```

If you already have a clone, run the last two commands from its root directory. The suite takes about three minutes, depending on the machine. With numpy 1.x, six tests in `tests/test_end_to_end.py` fail because their expected fixtures were written with numpy 2; upstream SplAdder fails in the same way.

A separate installation smoke test on Linux / Python 3.11 installed the GitHub package into a fresh Conda environment without upstream SplAdder, passed `pip check` and all four subcommand help checks, and completed serial and parallel builds on two synthetic samples. All 36 output files matched in content, datasets and objects. This does not establish compatibility with every Python or dependency version.

## Performance

New whole-genome measurement (2026-10-15 to 2026-10-16): two human STAR paired-end samples from [PRJNA312948](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA312948), SRR3184279 and SRR3184280, about 299 million mapped alignments combined, hg38 / GENCODE v34. Both versions run the complete `spladder build` workflow from the GTF, with `--parallel 8 --merge-strat merge_graphs --extract-ase -n 100`:

| | SplAdder 3.1.1 | spladder-turbo 0.1.0 | Speed-up |
|---|---|---|---|
| Wall-clock time | 20,262.3 s (5 h 37 min 42 s) | 320.1 s (5 min 20 s) | 63.3 x |
| CPU time, including child processes | 21,108.3 s | 1,718.3 s | 12.3 x |

Both exit successfully. All 38 output files have equal content, datasets and objects, and the converted annotation objects are equal. Each version ran once on the same shared EPYC 7C13 server, turbo first, without clearing the OS file cache. This is a single measurement; details, exact timings and limits are in [the PRJNA312948 benchmark record](docs/BENCHMARK_PRJNA312948.md).

Historical measurement of the whole workflow of `spladder build` (graph augmentation, merging, quantification, extraction and verification of events), two human RNA-Seq samples ([GSE272028](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272028), runs SRR29779619 and SRR29779620; aligned with STAR 2.7.7a to hg38, 85 bp reads, about 44 million alignments each), GENCODE v34, whole genome, starting from the GTF file:

| | SplAdder 3.1.1 | spladder-turbo 0.1.0 | |
|---|---|---|---|
| Wall-clock time, both with `--parallel 8` | 4 h 12 min 45 s | 4 min 28.5 s and 4 min 33.1 s in two runs | 55.5-56.5 x |
| CPU time | 16,354 s | 1,425 s and 1,605 s | 10.2-11.5 x |

The whole-genome `--parallel 8` baseline and its successful output comparison come from a historical Claude shell transcript supplied by the maintainer. The original temporary files are no longer accessible in the current environment; this is not a new whole-genome measurement. See [the benchmark record](docs/BENCHMARK_PARALLEL8.md) for the evidence and its limits.

The same workflow on chr20 only, three samples (SRR29779619 to SRR29779621). The historical and new measurements used different servers, identified below:

| | SplAdder 3.1.1 | spladder-turbo 0.1.0 | |
|---|---|---|---|
| EPYC 7601, one process (historical) | 5 min 04 s | 40 s | 7.5 x |
| EPYC 7601, turbo `--parallel 8` (historical) | 5 min 04 s (one process) | 23 s | 13 x, compared with one process of SplAdder |
| EPYC 7C13, both `--parallel 8` (2026-10-15) | 253.8 s | 39.3 s | 6.5 x |

The new chr20 measurement runs each version once, from the GTF, without an extra warm-up or forced clearing of compiled caches. Both exit successfully; all 39 output files have identical content, and the converted annotation objects are equal. This compares file contents, datasets and objects, not necessarily raw file bytes. Commands, exact timings and limitations are in [the benchmark record](docs/BENCHMARK_PARALLEL8.md#chr20-rerun-2026-10-15).

How to read these numbers:

- Historical measurements used a Linux server with two AMD EPYC 7601 processors (128 threads) and 314 GB of memory, with Python 3.13, numpy 2.3, numba 0.62 and pysam 0.23, timed with `/usr/bin/time -v`. The new chr20 measurement used an EPYC 7C13 server (256 logical CPUs), Python 3.13.9, numpy 2.3.4, numba 0.62.1 and pysam 0.23.3, timed with `validation/benchmark.py`. Do not compare timings across these servers as a speed-up measurement.
- The historical GSE272028 whole-genome baseline runs SplAdder in one command with `--parallel 8`. An earlier baseline manually ran the two per-sample graph steps simultaneously: 3 h 04 min 15 s, or 40.5-41.2 x slower than turbo. Both baselines limit the linear algebra libraries to one thread. SplAdder's own `--parallel` only distributes the detection and verification of events.
- The speed-up of one process (7.5 x) comes from reading each contig of an alignment file once with a compiled reader, from compiled kernels and from vectorized code. The rest comes from running samples, contigs and event types in parallel processes (Linux only). The CPU time shows how much less work is done in total.
- The first run compiles numba functions, which takes about half a minute (68 s instead of 40 s for the chr20 run above). The compiled code is cached on disk.
- With `--parallel 8` the whole-genome run uses about 9 GB in total, in up to about 24 processes at the same time.
- Reusable timing and comparison scripts are in [validation/](validation/). Historical raw logs are not included; the available transcript evidence is recorded in [docs/BENCHMARK_PARALLEL8.md](docs/BENCHMARK_PARALLEL8.md).

## What was optimized, where, and why

spladder-turbo changes *how* SplAdder computes, not *what* it computes. This section lists the optimizations: where in the code, which problem of SplAdder 3.1.1 each one addresses, what is done instead, and what it brought. The complete record, with every changed function, is in [docs/PERFORMANCE_CHANGES.md](docs/PERFORMANCE_CHANGES.md).

### Where the time went in SplAdder 3.1.1

Sampling profile (py-spy) of SplAdder 3.1.1 on a chr20 subset of human data, and a look at the code, showed five things:

1. **Reading the alignments took about 70% of the run time.** For every gene and in each of five processing steps (introns of the sample, cassette exons, intron retentions, new intron edges, counting), SplAdder fetches the alignments of the gene window with `pysam.fetch` and filters them one by one in a Python loop. An alignment that overlaps a gene is decoded once per step, which reads the data about 5 times.
2. **Almost everything was serial.** `--parallel` only distributed the detection and verification of events. Graphs, merging and counting ran one sample after the other in one process.
3. **Many tiny Python and numpy calls per gene.** Labelling the genes, testing intron retentions, building segment graphs, removing duplicate exons and counting each do dozens of small numpy calls per gene, for 60,000 genes.
4. **Repeated work.** The same 217 MB graph file was loaded about 12 times in a run; `compute_gene_expression` re-read the whole table of gene names inside the loop over the genes; numba functions were compiled again by every process (about 17 s); linear algebra libraries ran hidden threads on small matrices.
5. **Sending data to worker processes.** The pools for the detection and verification of events pickled gene objects for every task. The main process spent 20 s on serialization alone.

### Where the time goes now

Whole genome, two samples, from the GTF file. This stage breakdown uses the earlier SplAdder baseline with manually parallelized graph augmentation (3 h 04 min), rather than the one-command `--parallel 8` baseline above. The times of spladder-turbo are taken from the time stamps of one run (4 min 33 s). Tasks of a parallel stage print their output when they finish, so the boundaries are accurate to a few seconds.

| Step | SplAdder 3.1.1 | spladder-turbo (`--parallel 8`) |
|---|---|---|
| Parsing the annotation (once per annotation) | 2 min 05 s | 42 s |
| Graph augmentation of the two samples | 1 h 34 min (two processes) | 1 min 53 s |
| Merging, counting, gene expression, events | 1 h 28 min | 1 min 58 s: merging 43 s, counting 47 s, gene expression 13 s, verification and writing 15 s (the events are detected during the counting) |
| **Total** | **3 h 04 min** | **4 min 33 s** |

### The optimizations

Paths are relative to `spladder/`. "New" marks files that do not exist in SplAdder.

**A. Reading the alignments**

| Where | Problem in SplAdder 3.1.1 | What spladder-turbo does | Effect |
|---|---|---|---|
| `evidence.py` (new), `reads.py` (`add_reads_from_bam`) | Every gene and every step fetches its window again and decodes the alignments in Python (about 5 times per alignment). | One sequential pass over each contig builds a compact form: the coverage as a step function per read class (by XS tag) and sorted arrays of the introns with their read counts. Every window request is then answered from memory by binary search, with the same values, dtypes and shapes. The memory is limited (`SPLADDER_EVIDENCE_MAX_GB`, LRU). | Removes the read amplification of the dominant cost. |
| `bamscan.py` (new) | Even one pass with pysam costs a Python object and several attribute lookups per alignment, about 4 times the time it takes to decompress the file. | A compiled (numba) reader for indexed BAM files: it finds the contig with the BAI index, inflates the BGZF blocks (CRC32 and sizes are checked) and decodes the records with the same filters in the same order. It gives up, and pysam is used, for CRAM, other indices, unusual tags, records that htslib refuses or expands, and damaged files. | chr1 of a human sample (about 4 million aligned blocks): 11.9 s to 3.2 s. All contigs of two samples: 3.0 times less CPU time than the pysam loop. |
| `reads.py` (`_get_alignment_file`, `_get_reads_collapsed`) | The per-window reading that remains (several BAM files, `--var-aware`) opened the BAM file, with header and index, for every gene, and converted all tags of every alignment. | Open files are cached per process (LRU); the CIGAR is decoded once; only the needed tags are read; the coverage is accumulated in a difference array. | Faster fallback path. |

**B. Annotation and splicing graphs**

| Where | Problem in SplAdder 3.1.1 | What spladder-turbo does | Effect |
|---|---|---|---|
| `init.py` (`init_genes_gtf`), `classes/splicegraph.py` (`from_gene`) | All attributes of all 2.9 million GTF lines were parsed, and the graph of a gene was built with a loop over all vertices for every pair of exons. | Only the attributes of the lines that are used are parsed (genes in the first pass, exons in the second); the graph is built with a dictionary. The annotation file is identical. | GENCODE v34: 2 min 05 s to 42-48 s. |
| `classes/gene.py` (`label_alt`) | A Python loop with many small numpy calls per gene. | A numba kernel for the usual graphs, a numpy implementation for large graphs and unusual dtypes, SplAdder's original code for degenerate graphs. The kernel reproduces two index quirks of the original, so the labels stay the same. | Compared with the original on 2.16 million random graphs and 3 x 60,669 real genes. |
| `classes/gene.py` (`to_sparse`, `from_sparse`, `add_exon`) | Converting the adjacency matrices through scipy's `coo_matrix` for every gene; `np.r_` and `expand_dims` for every exon. | The same entries are computed with numpy (scipy only for unusual shapes); `np.concatenate`. | About 10 times faster conversion. |
| `classes/segmentgraph.py` (`from_gene`), `classes/splicegraph.py` (`uniquify`) | Loops over segments and over all pairs of exons. | Vectorized (`searchsorted`, broadcasting); `uniquify` returns early when there are no duplicates. | |
| `editgraph.py` (`insert_intron_retentions`, `insert_cassette_exons`) | For every pair of connected exons, and for every pair of introns, several numpy calls inside Python loops. | The coverage test for intron retentions is a numba function that follows the order of operations of the numpy code, so the floating point results are the same. The candidate pairs for cassette exons are found with vectorized operations. Unusual inputs use the original loops. | |
| `core/gen_graphs.py`, `core/spladdercore.py` | The editing of the genes ran in one process. | The editing of a gene only depends on the genes of its contig, so the contigs are put together to shards of similar size and edited in different processes. The side effects of the steps on the options are replayed in the main process. The evidence of the contigs is read in parallel. | Part of the graph stage: 1 h 34 min to 1 min 53 s with the other changes. |

**C. Merging and counting**

| Where | Problem in SplAdder 3.1.1 | What spladder-turbo does | Effect |
|---|---|---|---|
| `merge.py` | The merged graph was written, then loaded by the counting step, which built the segment graphs and wrote the 217 MB file a second time. | If the merged graph is quantified next and not changed afterwards (no `--validate-sg`), the segment graphs are made while merging. | One write of the large file less. |
| `count.py` (`_count_gene`, `_count_contig`, `_count_shard`, `_collect_counts`) | For every segment of every gene, the mean coverage was computed with numpy calls; for every edge, the whole intron list was scanned. Samples and contigs were counted one after the other. | The means of all segments of a gene come from cumulative sums (exact, because the coverage values are integers, so the numbers are the same). The intron counts come from a sorted search. Contigs and samples are counted in different processes. The counts of all genes are put together in one step. | |
| `count.py` (`compute_gene_expression`) | The whole `gene_names` table was read from the HDF5 file inside the loop over the genes. | It is read once. | 44.6 s to 9.9 s. |

**D. Events**

| Where | Problem in SplAdder 3.1.1 | What spladder-turbo does | Effect |
|---|---|---|---|
| `alt_splice/detect.py`, `alt_splice/verify.py` | Every task of the process pools received pickled gene objects. For each event, the counts of its gene were read from the HDF5 file again. | The workers inherit the graph by fork; tasks only carry indices. The counts of a gene are read once for all its events; rows are found by binary search. Numba functions are cached on disk. | Verification of exon skips: 36 s to about 10 s. The 20 s of serialization in the main process are gone. |
| `spladder_build.py` | The events were detected after the quantification, and the event types were verified one after the other. | Events depend on the splicing graph, not on the counts, so they are collected in three background jobs while the graph is quantified. The six event types are verified at the same time. | The detection of events no longer appears as a step of its own. |

**E. Orchestration and input/output**

| Where | Problem in SplAdder 3.1.1 | What spladder-turbo does | Effect |
|---|---|---|---|
| `parallel.py` (new) | No parallel processing of the samples. They cannot simply be run side by side, because `make_introns_feasible` changes `options.read_filter` in place, so in a serial run a sample depends on the samples before it. | A fork-based framework. The graphs of the samples are computed speculatively in parallel, each from the state it would have if no sample changed the filter. The results are accepted in sample order and only if that state was right; otherwise the sample is computed again in order. A worker that dies raises an error; workers end with the main process; output and errors are shown by the main process. | The result is exactly that of the serial run. |
| `utils.py` (`load_pickle`, `dump_pickle`, `load_pickle_shared`) | The splicing graphs consist of millions of small objects: the garbage collector made loading about twice as slow. A run loaded the same graph file about 12 times. A killed run could leave a part of a file. | The garbage collector is paused while loading and writing. Files are written to a temporary file and renamed. A file that is unchanged is not loaded again by the same process. | Loading about 2 times faster, writing about one third faster. |
| `utils.py` (`jit_cached`) | Every process compiled the numba functions again (about 17 s). | The compiled code is cached on disk. | The first run compiles (about half a minute), later runs and all worker processes do not. |
| `spladder.py` | The linear algebra libraries (OpenBLAS, OpenMP, MKL) ran threads on small matrices, which only burn CPU, and oversubscribe the cores when several processes run. | One thread, unless the user sets the variables. | |

### How the results are kept the same

- **Reproduce, do not approximate.** Compiled kernels follow the order of the floating point operations of the original code. Where sums are replaced by cumulative sums, the values are integers, so the sums are exact.
- **Fall back to SplAdder's code.** Whatever a fast path does not cover (several alignment files, `--var-aware`, CRAM, a missing NM tag, unusual tags, damaged files, degenerate graphs) is handled by the code of SplAdder.
- **The original stays as the reference.** `tests/reference.py` contains verbatim copies of the SplAdder 3.1.1 implementations of the optimized functions; the tests compare against them with random and adversarial input.
- **Speculation is checked.** Parallel processing of the samples is only accepted if the state each sample started from is the one of the serial run.
- **Peculiarities are kept.** Where SplAdder behaves in a way that looks unintended but changes the output (for example in `Gene.label_alt`), spladder-turbo does the same.

### What was not changed

The method: the rules and thresholds for augmenting the splicing graph, the definitions of the events, the criteria for verifying them, the quantification, the command line options and the output formats.

## Are the results the same?

The goal is that every output of `spladder build` has the same content as with SplAdder 3.1.1. What has been compared so far:

- **Whole genome, paired-end STAR data** (PRJNA312948, SRR3184279 and SRR3184280): both versions with `--parallel 8`, all 38 output files have equal content, HDF5 datasets and pickle objects, including both sample graphs, the merged graph, counts, gene expression and all six event types. The converted annotation objects are also equal. See [the full record](docs/BENCHMARK_PRJNA312948.md) for the tested options and limits.
- **Whole genome** (the historical GSE272028 run): all 38 output files (text, GFF3, HDF5, event files, the splicing graphs of the samples and the merged graph) and the annotation file. 29 of the 38 files are byte-identical (GFF3, HDF5, event and index files), and so is the annotation file. The six compressed text files (`.txt.gz`) differ only in the time stamp in the gzip header; their text is identical. The three graph files (`genes_graph_*.pickle`) differ in their bytes, but they load to completely equal objects.
- **42 command line configurations** (confidence levels, merge strategies, event types, output formats, sparse BAM files, ...) on the chr20 subset, each run with both programs: identical, or failing in the same way (two configurations fail in SplAdder as well).
- **More than 50 randomly generated data sets** (random annotations and reads with novel exons, retained introns, alternative sites and mutually exclusive exons), run through both programs, serial and parallel: identical.
- **The compiled BAM reader** against pysam on every contig of two real samples (678 combinations of contig and filter) and on 800 random BAM files: identical.
- **250 tests** in this repository compare the optimized functions with verbatim copies of the SplAdder 3.1.1 implementations (`tests/reference.py`) on random and adversarial input.

Known differences:

- spladder-turbo fixes defects of SplAdder 3.1.1 (crashes for annotations without edges or with unknown strands, merging samples with different gene sets, undefined names in rarely used code, stale index files, ...). Where SplAdder worked, the output is the same; where it crashed, spladder-turbo continues. The list is in [docs/PERFORMANCE_CHANGES.md](docs/PERFORMANCE_CHANGES.md) (section 9).
- Peculiarities of SplAdder that change the output are kept, not fixed (for example in `Gene.label_alt`).
- With `--parallel`, the output of the parallel tasks is shown one task after the other, and `--logfile` is written per sample.
- Parallel processing needs Linux; elsewhere everything runs serially (slower, but the same results).

Not validated yet: aligners other than STAR, species other than human, cohorts of more than three samples, CRAM files, `--var-aware`, `--sparse-bam` for a whole genome, and operating systems other than Linux. Paired-end validation is limited to the two PRJNA312948 samples and options recorded above. The code is written to behave like SplAdder in these cases (anything that is not understood falls back to the code of SplAdder), but this has only been tested on small synthetic data. Please compare the outputs of both programs on a small part of your own data before you switch (see [Switching from SplAdder](#switching-from-spladder)), and report differences as issues.

## Documentation

- [SplAdder documentation](https://spladder.readthedocs.io) for the command line and the output files.
- [docs/PERFORMANCE_CHANGES.md](docs/PERFORMANCE_CHANGES.md): what was changed and why, how the results are kept identical, and the record of the code review (in Chinese for now).
- [ATTRIBUTION.md](ATTRIBUTION.md): credits and licensing.
- [README.zh-CN.md](README.zh-CN.md): this page in Chinese.

## Citation

spladder-turbo is a derivative of SplAdder. Please cite the SplAdder paper:

Kahles A, Ong CS, Zhong Y, Rätsch G. SplAdder: identification, quantification and testing of alternative splicing events from RNA-Seq data. *Bioinformatics* 2016;32(12):1840-1847. doi:[10.1093/bioinformatics/btw076](https://doi.org/10.1093/bioinformatics/btw076)

For reproducibility, please also state that you used spladder-turbo and its version (`spladder --help` shows it).

## License and attribution

SplAdder was written by Cheng Soon Ong, Gunnar Raetsch and Andre Kahles (see [AUTHORS](AUTHORS)) and is distributed under the BSD 3-Clause license ([COPYRIGHT](COPYRIGHT)). The modifications of this fork are distributed under the same license ([LICENSE](LICENSE)). The names of the copyright holders and contributors are not used to endorse or promote this fork. [ATTRIBUTION.md](ATTRIBUTION.md) has the details, including a note on one source file that carries a GPL notice upstream.

The software is provided "as is", without warranty of any kind; see the license for the exact terms.
