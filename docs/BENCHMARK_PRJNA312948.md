# Whole-Genome Comparison on PRJNA312948 Paired-End Samples

## Scope and Inputs

On 2026-10-15, PRJNA312948 was selected from locally stored public datasets.
The first two samples in accession order, SRR3184279 and SRR3184280, were used
for a whole-genome comparison. No chromosomes were selected, no reads were
subsampled, and no historical graphs, counts or converted annotations were reused.

| Sample | BAM size (bytes) | Mapped alignments | Unmapped alignments |
|---|---|---|---|
| SRR3184279 | 8,492,134,017 | 194,089,948 | 793,100 |
| SRR3184280 | 4,566,068,740 | 104,934,102 | 373,624 |

Alignment counts come from newly built BAM indices. They include multimapping
and individual paired-end records; they are not counts of read pairs.
Both BAM files are coordinate-sorted, with 194 reference sequences in their
headers. The first 1,000 mapped alignments of each file were inspected, and all
had the paired flag set. Of these, 986 and 994 records, respectively, were
100 bp long; the remainder were shorter. The original STAR logs describe
paired-end input with an average input length of 199 bp across both ends,
so both programs used `-n 100`.

The original BAM files are local copies of public data. Storage paths are
represented by placeholders to avoid publishing internal directory layouts:

```text
<local-data-root>/PRJNA312948/aligned/SRR3184279/SRR3184279_Aligned.sortedByCoord.out.bam
<local-data-root>/PRJNA312948/aligned/SRR3184280/SRR3184280_Aligned.sortedByCoord.out.bam
```

The annotation was a local hg38 / GENCODE v34 file. The original STAR alignment
and featureCounts logs both referenced that hg38 reference directory.

Original inputs were read-only. The work directory contained symlinks to the
BAM and GTF files, with new BAI indices built there using samtools/pysam.
Index preparation took 46.6 s and 21.7 s, respectively, and was excluded from
program timings.

## Environment and Method

- Linux 5.14.0-687.15.1.el9_8.x86_64, AMD EPYC 7C13 64-Core Processor, 256 logical CPUs.
- Python 3.13.9, numpy 2.3.4, numba 0.62.1, pysam 0.23.3.
- Upstream SplAdder 3.1.1: `65ceec839b9ff0cf96703c1605ee43667662f410`.
- spladder-turbo 0.1.0: `a16b76ed2a33ea81db41eb9fdaa6a859940ae958`.
- Both source trees were independently exported from their Git commits.
  Actual import paths and versions were checked; neither run used workspace source files.
- Identical options: `--parallel 8 --merge-strat merge_graphs --extract-ase -n 100`.
- `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, `MKL_NUM_THREADS` and `NUMBA_NUM_THREADS` were all 1.
- Each program ran once, sequentially, turbo first and upstream second.
  Each had an independent output directory and converted annotation file.
- Source directories were fresh, without a compilation warm-up.
  Compilation during execution was included in the timing.
- OS file caches were not cleared; index preparation had already read the BAM files.
  The server was shared, without isolation from other workloads.
- While upstream was running, the completed annotation and two sample graphs
  were compared separately. These checks introduced some additional load.
- Timing used `timed_run` from `validation/benchmark.py`: wall-clock time,
  cumulative CPU time including child processes, and largest individual
  process RSS in GiB, not total process-tree memory.
- matplotlib and fontconfig caches were kept inside the work directory
  to avoid writing to user configuration directories.

## Results

Both programs exited successfully with code 0. The comparison ran on
2026-10-15 to 2026-10-16 (Asia/Shanghai):

| Metric | SplAdder 3.1.1 | spladder-turbo 0.1.0 | Upstream / turbo |
|---|---:|---:|---:|
| Wall-clock (s) | 20,262.288435 | 320.140949 | 63.2918 x |
| Wall-clock (h/min/s) | 5 h 37 min 42.3 s | 5 min 20.1 s | |
| Cumulative CPU (s, including children) | 21,108.339352 | 1,718.284282 | 12.2845 x |
| Largest process RSS (GiB) | 2.7644 | 3.3475 | |

RSS is not the total for the process tree and cannot establish total peak memory.
These are single measurements on the same machine with identical options,
not speed-up guarantees for every dataset. Shared server load, OS file caches
and execution order can affect timings. No historical GSE272028 timings were
used to calculate these speed-ups.

`tests.compare.compare_dirs` compared both output directories, excluding
the temporary `tmp/` directory. All **38 output files had equal content**:

| Output category | Files | Comparison |
|---|---:|---|
| Two sample graphs and merged graph | 3 | Recursive comparison of complete deserialized pickle objects |
| Six event types and their confirmed events | 12 | Recursive comparison of complete deserialized pickle objects |
| Counts for six event types | 6 | HDF5 dataset names, values, shapes and dtypes |
| Gene counts and expression | 2 | HDF5 dataset names, values, shapes and dtypes |
| Count indices | 3 | Recursive comparison of deserialized arrays |
| GFF3 | 6 | File content |
| Compressed text | 6 | Decompressed gzip content |

The six event types were exon_skip, mult_exon_skip, intron_retention,
alt_3prime, alt_5prime and mutex_exons. A separate `deep_equal` comparison
found the converted annotation objects equal.
Equality does not require identical gzip header timestamps, pickle
serialization bytes or HDF5 file bytes.

The complete `comparison.json` result was:

```json
{
  "output_files": 38,
  "differences": [],
  "annotation_difference": null
}
```

This validates the complete whole-genome workflow for these two STAR paired-end
samples with the recorded options. It does not establish validation for every
paired-end dataset, aligner or option combination.

## Commands and Artifacts

Each source snapshot ran the following command in its own directory, with
independent annotation and output paths. Replace `<workdir>` with the local
absolute work directory. The full executed commands and environment are recorded
in `status.json`; `run.py` used the existing validation tools for orchestration
and timing:

```bash
python -m spladder.spladder build \
    -a "annotation/gencode.v34.annotation.gtf" \
    -o "out" \
    -b "<workdir>/inputs/SRR3184279.bam,<workdir>/inputs/SRR3184280.bam" \
    --parallel 8 --merge-strat merge_graphs --extract-ase -n 100
```

Local artifacts are under the repository-relative directory
`build/as_pair_20261015/`, excluded by the existing `.gitignore`:

- `run.py`: orchestration script using the existing timing, environment and comparison tools.
- `inputs.json`: input sources, sizes, modification times, index counts and source commits.
- `status.json`: run state, full commands, actual import paths, timestamps and timings.
- `src_{upstream,turbo}/`: the two fixed source snapshots.
- `inputs/`: symlinks to original inputs and newly built BAI indices.
- `{upstream,turbo}/log.txt`: complete logs.
- `{upstream,turbo}/{out,annotation}/`: outputs and converted annotations.
- `comparison.json`: comparison results for all outputs and annotation objects.

These artifacts remain local. BAM files, indices, source snapshots and computed
results are not committed to Git.
