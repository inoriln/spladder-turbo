# `--parallel 8` 测试记录

## 全基因组：历史终端记录

来源：维护者于 2026-10-15 提供的 Claude 命令及输出。原始目录
`/tmp/baseline_par8`、`/tmp/spladder_bench`、`/tmp/spladder_gold` 在当前环境中
均不可访问；以下数据是历史终端记录的摘录，不是本次重新测量或重新比对的结果。

场景为 GENCODE v34、人类全基因组、GSE272028 的 SRR29779619 和 SRR29779620，
从 GTF 开始，SplAdder 3.1.1 使用一条 `build` 命令及 `--parallel 8`。
终端输出中的完整命令被 `cut` 截断，原始 `run.sh` 未能读取；数据和其余选项沿用
已有性能记录，不能从这份截断输出独立验证全部参数和线程设置。

| 指标 | SplAdder 3.1.1 | spladder-turbo 的历史两次运行 |
|---|---|---|
| Wall-clock | 15,165 秒（4:12:45） | 268.5 秒、273.1 秒 |
| User CPU | 15,505.03 秒 | 未在该终端摘录中单独列出 |
| System CPU | 849.41 秒 | 未在该终端摘录中单独列出 |
| 总 CPU | 16,354.44 秒 | 1,425 秒、1,605 秒 |
| 最大单进程 RSS | 2,955,236 KiB（约 2.82 GiB） | 未在该终端摘录中列出 |
| 退出码 | 0 | 既有记录为成功完成 |

按每次 turbo 运行分别计算，wall-clock 加速为 56.5 倍、55.5 倍，CPU 加速为
11.5 倍、10.2 倍。这些都是单次运行的比值，不是多次重复的中位数。

历史记录中的输出比对命令：

```bash
cd /tmp/spladder_bench
/data/home/zyh/micromamba/envs/spladder/bin/python compare_runs.py \
    /tmp/spladder_gold/out /tmp/baseline_par8/out
```

终端报告：

```text
35 项比对, ALL IDENTICAL
```

该命令比较的是原版 `--parallel 8` 输出与历史金标准；不能将它描述为本次对
原版与 turbo 全部输出文件的重新比对。

此前 11,055 秒（3:04:15）的基线手工把两个样本的建图步骤并行运行，对应 turbo
的 40.5-41.2 倍加速。它与单条命令 `--parallel 8` 的基线是两种不同的调度方式，
不能混用阶段耗时或加速比。

## chr20: 2026-10-15 复测

### 环境与输入

- Linux 5.14.0-687.15.1.el9_8.x86_64，AMD EPYC 7C13 64-Core Processor，256 逻辑 CPU。
- 同一解释器：`/data/home/zyh/micromamba/envs/spladder/bin/python`；Python 3.13.9、
  numpy 2.3.4、numba 0.62.1、pysam 0.23.3、h5py 3.15.1。
- 原版源码通过 `git archive` 从 `65ceec839b9ff0cf96703c1605ee43667662f410`
  导出，包版本为 3.1.1；turbo 源码为 `b00d6c281f6f65e6b45a2bb6724ee0eb9b34dbf3`。
- 注释来自 `/data/data2/reference/hg38_ek12/gencode.v34.annotation.gtf`，
  保留注释头与 `chr20` 行，双方各自转换 GTF，不使用历史注释 pickle。
- 从原始排序 BAM 中用 `pysam.view` 提取全部 chr20 比对，并重新建立索引：

| 样本 | chr20 的 mapped alignments |
|---|---|
| SRR29779619 | 1,028,164 |
| SRR29779620 | 1,108,664 |
| SRR29779621 | 1,004,647 |

两个程序依次运行，避免相互争抢资源。均使用 `--parallel 8 --merge-strat merge_graphs
--extract-ase -n 85`；`OPENBLAS_NUM_THREADS`、`OMP_NUM_THREADS`、`MKL_NUM_THREADS`
和 `NUMBA_NUM_THREADS` 由验证工具在子进程环境中设为 1。每个版本运行一次，没有
`--discard-first`，也没有强制清除已有 numba 缓存，不能视为多次重复的稳定态统计。
这是共享服务器上的观测，没有隔离其他用户负载。日志有 matplotlib/fontconfig
缓存目录不可写的警告，双方均正常完成。

### 结果

| 指标 | SplAdder 3.1.1 | spladder-turbo 0.1.0 |
|---|---|---|
| Wall-clock | 253.818110 秒 | 39.335032 秒 |
| 总 CPU | 411.865391 秒 | 264.587424 秒 |
| 最大单进程 RSS | 0.372540 GiB | 0.391872 GiB |
| 退出码 | 0 | 0 |

Wall-clock 加速约 6.5 倍，CPU 加速约 1.6 倍。内存是最大单个进程的 RSS，
不是整个进程树的总内存；计时工具 JSON 字段名为 `memory_gb`，实际换算单位为 GiB。
CPU 时间由 `os.wait4` 获取，包含子进程累计用量。

`tests.compare.compare_dirs` 比较双方输出目录，包含压缩文本解压后的内容、GFF3、
HDF5 数据集和 pickle 对象；另用 `deep_equal` 比较转换后的注释对象：

```json
{
  "output_files": 39,
  "output_differences": [],
  "annotation_difference": null
}
```

39 个输出文件内容一致，另行比较的注释对象也一致。这不是逐字节一致性的结论。
本次服务器与历史测量机器不同，不能用本次原版耗时除以历史 turbo 耗时计算加速比。

### 命令与产物

以下命令在仓库根目录执行，输入是本次提取的 chr20 数据。重新执行时使用新的
`--workdir`，因为计时脚本会重建其输出目录。

```bash
/data/home/zyh/micromamba/envs/spladder/bin/python -B validation/benchmark.py \
    --upstream-src "$PWD/build/parallel8_20261015/src_orig" \
    --annotation "$PWD/build/parallel8_20261015/inputs/gencode.v34.chr20.gtf" \
    --bams "$PWD/build/parallel8_20261015/inputs/SRR29779619.chr20.bam,$PWD/build/parallel8_20261015/inputs/SRR29779620.chr20.bam,$PWD/build/parallel8_20261015/inputs/SRR29779621.chr20.bam" \
    --workdir "$PWD/build/parallel8_20261015/benchmark" \
    --repeat 1 -- --parallel 8 --merge-strat merge_graphs --extract-ase -n 85
```

比较该次生成的输出和注释对象：

```bash
/data/home/zyh/micromamba/envs/spladder/bin/python -B - <<'PY'
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

本机产物位于仓库内 `build/parallel8_20261015/`，由现有 `.gitignore` 排除：

- `inputs.json`：输入来源、chr20 比对数和两份源码提交。
- `inputs/`：提取后的 GTF、BAM 和索引。
- `src_orig/`：未修改的 SplAdder 3.1.1 源码。
- `benchmark/results.json`：计时工具生成的完整数值和环境。
- `benchmark/comparison.json`：本次输出与注释对象比对结果。
- `benchmark/{upstream,turbo}/log.txt`：各程序完整日志。
- `benchmark/{upstream,turbo}/{out,annotation}/`：各程序输出与转换后的注释。

大型数据、日志和输出未纳入 Git；本记录保留精确计时、比对结论及复现参数。
