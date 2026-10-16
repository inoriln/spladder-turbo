# `--parallel 8` Benchmark Record

## Whole Genome: Historical Terminal Record

Source: Claude commands and output supplied by the maintainer on 2026-10-15.
The original temporary directories are no longer accessible in the current
environment. The figures below are excerpts from that historical record, not
new measurements or new output comparisons.

The workload used GENCODE v34, the whole human genome, and samples SRR29779619
and SRR29779620 from GSE272028. SplAdder 3.1.1 started from the GTF in a single
`build` command with `--parallel 8`. The full command was truncated by `cut`
in the terminal output, and the original `run.sh` could not be inspected.
The dataset and remaining options follow the existing performance record;
the truncated output does not independently establish every option or thread setting.

| Metric | SplAdder 3.1.1 | Two historical spladder-turbo runs |
|---|---|---|
| Wall-clock | 15,165 s (4:12:45) | 268.5 s, 273.1 s |
| User CPU | 15,505.03 s | Not reported separately in the excerpt |
| System CPU | 849.41 s | Not reported separately in the excerpt |
| Total CPU | 16,354.44 s | 1,425 s, 1,605 s |
| Largest process RSS | 2,955,236 KiB (about 2.82 GiB) | Not reported in the excerpt |
| Exit code | 0 | Existing records report successful completion |

Relative to each turbo run, wall-clock speed-ups are 56.5 x and 55.5 x;
CPU speed-ups are 11.5 x and 10.2 x. These are ratios of individual runs,
not medians of repeated measurements.

The historical output comparison command is shown below with local paths
replaced by placeholders. Replace them with quoted paths before running it:

```bash
cd "<benchmark-workdir>"
python compare_runs.py \
    "<historical-gold>/out" "<historical-baseline>/out"
```

The terminal reported the following (translated into English):

```text
35 comparisons, ALL IDENTICAL
```

This command compared the upstream `--parallel 8` output with the historical
gold standard. It is not a new comparison of all upstream and turbo output files.

An earlier baseline of 11,055 s (3:04:15) manually ran the two sample graph
steps concurrently, giving a 40.5-41.2 x speed-up for turbo. That scheduling
differs from the single-command `--parallel 8` baseline; their stage timings
and speed-up ratios must not be mixed.

## chr20 Rerun: 2026-10-15

### Environment and Inputs

- Linux 5.14.0-687.15.1.el9_8.x86_64, AMD EPYC 7C13 64-Core Processor, 256 logical CPUs.
- Both programs used the same Python environment: Python 3.13.9, numpy 2.3.4,
  numba 0.62.1, pysam 0.23.3, h5py 3.15.1.
- Upstream source was exported with `git archive` from
  `65ceec839b9ff0cf96703c1605ee43667662f410`, package version 3.1.1;
  turbo source was `b00d6c281f6f65e6b45a2bb6724ee0eb9b34dbf3`.
- The annotation was a local hg38 / GENCODE v34 file, retaining the header
  and `chr20` lines. Each program converted the GTF independently, without
  reusing historical annotation pickles.
- All chr20 alignments were extracted from the original sorted BAM files
  using `pysam.view`, and new indices were built:

| Sample | Mapped alignments on chr20 |
|---|---|
| SRR29779619 | 1,028,164 |
| SRR29779620 | 1,108,664 |
| SRR29779621 | 1,004,647 |

The programs ran sequentially to avoid competing with each other. Both used
`--parallel 8 --merge-strat merge_graphs --extract-ase -n 85`.
The validation tool set `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`,
`MKL_NUM_THREADS` and `NUMBA_NUM_THREADS` to 1 in each subprocess environment.
Each version ran once, without `--discard-first` or forced clearing of existing
numba caches. These are not statistics from repeated steady-state measurements.
The server was shared, without isolation from other users' workloads.
Logs contained warnings about unwritable matplotlib/fontconfig cache directories;
both programs completed successfully.

### Results

| Metric | SplAdder 3.1.1 | spladder-turbo 0.1.0 |
|---|---|---|
| Wall-clock | 253.818110 s | 39.335032 s |
| Total CPU | 411.865391 s | 264.587424 s |
| Largest process RSS | 0.372540 GiB | 0.391872 GiB |
| Exit code | 0 | 0 |

The wall-clock speed-up is about 6.5 x; the CPU speed-up is about 1.6 x.
Memory is the largest individual process RSS, not the total for the process tree.
The timing tool's JSON field is named `memory_gb`, but its conversion uses GiB.
CPU time comes from `os.wait4` and includes cumulative child-process usage.

`tests.compare.compare_dirs` compared both output directories, including
decompressed text, GFF3, HDF5 datasets and pickle objects. A separate
`deep_equal` comparison checked the converted annotation objects:

```json
{
  "output_files": 39,
  "output_differences": [],
  "annotation_difference": null
}
```

All 39 output files had equal content, and the separately compared annotation
objects were equal. This does not assert byte-for-byte equality of the files.
This server differs from the historical measurement machine; dividing the new
upstream time by a historical turbo time would not be a valid speed-up measurement.

### Commands and Artifacts

Run the following from the repository root, using the chr20 inputs prepared
for this comparison. Use a fresh `--workdir` when repeating the benchmark,
because the timing script rebuilds its output directory.

```bash
python -B validation/benchmark.py \
    --upstream-src "$PWD/build/parallel8_20261015/src_orig" \
    --annotation "$PWD/build/parallel8_20261015/inputs/gencode.v34.chr20.gtf" \
    --bams "$PWD/build/parallel8_20261015/inputs/SRR29779619.chr20.bam,$PWD/build/parallel8_20261015/inputs/SRR29779620.chr20.bam,$PWD/build/parallel8_20261015/inputs/SRR29779621.chr20.bam" \
    --workdir "$PWD/build/parallel8_20261015/benchmark" \
    --repeat 1 -- --parallel 8 --merge-strat merge_graphs --extract-ase -n 85
```

Compare the generated outputs and annotation objects:

```bash
python -B - <<'PY'
import pickle
from pathlib import Path
from tests.compare import compare_dirs, deep_equal

work = Path('build/parallel8_20261015/benchmark')
differences = compare_dirs(str(work / 'upstream/out'), str(work / 'turbo/out'))
with (work / 'upstream/annotation/gencode.v34.chr20.gtf.pickle').open('rb') as a, \
     (work / 'turbo/annotation/gencode.v34.chr20.gtf.pickle').open('rb') as b:
    annotation_difference = deep_equal(pickle.load(a), pickle.load(b), 'annotation')
print('outputs:', differences or 'IDENTICAL')
print('annotation:', annotation_difference or 'IDENTICAL')
raise SystemExit(bool(differences or annotation_difference))
PY
```

Local artifacts are under the repository-relative directory
`build/parallel8_20261015/`, excluded by the existing `.gitignore`:

- `inputs.json`: input sources, chr20 alignment counts and both source commits.
- `inputs/`: extracted GTF, BAM files and indices.
- `src_orig/`: unmodified SplAdder 3.1.1 source.
- `benchmark/results.json`: complete timings and environment information.
- `benchmark/comparison.json`: output and annotation object comparison results.
- `benchmark/{upstream,turbo}/log.txt`: complete program logs.
- `benchmark/{upstream,turbo}/{out,annotation}/`: outputs and converted annotations.

Large datasets, logs and outputs are not tracked by Git. This record retains
the exact timings, comparison findings and reproduction parameters.
