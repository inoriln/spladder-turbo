# PRJNA312948 双端样本的全基因组比较

## 范围与输入

2026-10-15，从 `/data/data2/IO_RNA/AS/` 选择 PRJNA312948，使用按编号排序的
前两个样本 SRR3184279、SRR3184280，进行全基因组比较。未截取染色体、未下采样，
未复用历史剪接图、计数结果或转换后的注释。

| 样本 | BAM 大小（字节） | Mapped alignments | Unmapped alignments |
|---|---|---|---|
| SRR3184279 | 8,492,134,017 | 194,089,948 | 793,100 |
| SRR3184280 | 4,566,068,740 | 104,934,102 | 373,624 |

这里的 alignment 数是新建 BAM 索引的记录数，包含多重比对和双端的各条记录，
不是 read pair 数。两份 BAM 均按坐标排序，头部含 194 个参考序列。各抽查前
1,000 条 mapped alignments，均设置 paired 标记；其中长度为 100 bp 的记录分别
有 986 条和 994 条，其余为较短读段。STAR 原始日志显示双端输入，平均输入长度
为 199 bp（两端合计），所以双方均设置 `-n 100`。

原始 BAM 路径为：

```text
/data/data2/IO_RNA/AS/PRJNA312948/aligned/SRR3184279/SRR3184279_Aligned.sortedByCoord.out.bam
/data/data2/IO_RNA/AS/PRJNA312948/aligned/SRR3184280/SRR3184280_Aligned.sortedByCoord.out.bam
```

注释为 `/data/data2/reference/hg38_ek12/gencode.v34.annotation.gtf`。
原始 STAR 比对与 featureCounts 日志均指向该 hg38 参考目录。

原始数据只读使用。测试目录通过软链接引用 BAM 和 GTF，并使用 samtools/pysam
在测试目录重建 BAI；两份索引准备分别耗时 46.6 秒、21.7 秒，不计入程序耗时。

## 环境与方法

- Linux 5.14.0-687.15.1.el9_8.x86_64，AMD EPYC 7C13 64-Core Processor，256 逻辑 CPU。
- Python 3.13.9、numpy 2.3.4、numba 0.62.1、pysam 0.23.3。
- 原版 SplAdder 3.1.1：`65ceec839b9ff0cf96703c1605ee43667662f410`。
- spladder-turbo 0.1.0：`a16b76ed2a33ea81db41eb9fdaa6a859940ae958`。
- 两份源码均从相应 Git 提交独立导出，并核对实际导入路径和版本；运行期间不读取工作区源码。
- 双方参数完全相同：`--parallel 8 --merge-strat merge_graphs --extract-ase -n 100`。
- `OPENBLAS_NUM_THREADS`、`OMP_NUM_THREADS`、`MKL_NUM_THREADS`、`NUMBA_NUM_THREADS` 均为 1。
- 两个程序各运行一次，先 turbo 后原版，依次运行。每次使用独立输出目录和注释转换文件。
- 使用新源码目录，没有提前运行编译预热；运行中产生的编译开销计入耗时。
- 没有清除操作系统文件缓存；准备索引时已读取过 BAM。服务器为共享环境，未隔离其他负载。
- 原版运行期间，另行比较过已生成的注释和两个样本图对象；这些检查也会带来少量额外负载。
- 使用现有 `validation/benchmark.py` 的 `timed_run` 计时：wall-clock、包含子进程的
  累计 CPU、最大单进程 RSS（GiB，不是整个进程树内存合计）。
- matplotlib 和 fontconfig 使用测试目录内的缓存，避免向用户配置目录写入。

## 结果

双方均正常退出（退出码 0）。2026-10-15 至 2026-10-16（Asia/Shanghai）完成本次测试：

| 指标 | SplAdder 3.1.1 | spladder-turbo 0.1.0 | 原版 / turbo |
|---|---:|---:|---:|
| Wall-clock（秒） | 20,262.288435 | 320.140949 | 63.2918 倍 |
| Wall-clock（时分秒） | 5 小时 37 分 42.3 秒 | 5 分 20.1 秒 | |
| 累计 CPU（秒，含子进程） | 21,108.339352 | 1,718.284282 | 12.2845 倍 |
| 最大单进程 RSS（GiB） | 2.7644 | 3.3475 | |

RSS 不是进程树的内存合计，不能据此推断双方的总内存峰值。
这是同机、同参数、各一次的实测，不代表所有数据集的加速幅度；共享服务器负载、
系统文件缓存和执行顺序均可能影响耗时。没有用此前 GSE272028 的历史耗时计算本次加速比。

使用 `tests.compare.compare_dirs` 比较两个输出目录（排除 `tmp/` 临时目录），
全部 **38 个输出文件内容一致**：

| 输出类别 | 文件数 | 比较方式 |
|---|---:|---|
| 两个样本图及合并图 | 3 | pickle 加载后的完整对象递归比较 |
| 六类事件及其 confirmed 事件 | 12 | pickle 加载后的完整对象递归比较 |
| 六类事件计数 | 6 | HDF5 数据集名称、值、shape、dtype |
| 基因计数与基因表达 | 2 | HDF5 数据集名称、值、shape、dtype |
| 计数索引 | 3 | pickle 加载后的数组递归比较 |
| GFF3 | 6 | 文件内容 |
| 压缩文本 | 6 | gzip 解压后的内容 |

六类事件包括 exon_skip、mult_exon_skip、intron_retention、alt_3prime、alt_5prime、
mutex_exons。另用 `deep_equal` 比较转换后的注释对象，结果一致。
这里不要求 gzip 头时间戳、pickle 序列化或 HDF5 文件的原始字节一致。

`comparison.json` 的完整结果为：

```json
{
  "output_files": 38,
  "differences": [],
  "annotation_difference": null
}
```

因此，本次验证覆盖这两个 STAR 双端样本在上述选项下的全基因组完整流程；
不能据此泛化为所有双端数据、比对器或选项组合均已验证。

## 命令与产物

两份源码分别在各自目录运行以下命令，注释和输出路径独立；实际完整命令及环境
记录在 `status.json`，编排和计时由 `run.py` 调用现有验证工具完成：

```bash
python -m spladder.spladder build \
    -a "annotation/gencode.v34.annotation.gtf" \
    -o "out" \
    -b "/data/home/zyh/others/spladder/build/as_pair_20261015/inputs/SRR3184279.bam,/data/home/zyh/others/spladder/build/as_pair_20261015/inputs/SRR3184280.bam" \
    --parallel 8 --merge-strat merge_graphs --extract-ase -n 100
```

本机原始产物位于 `build/as_pair_20261015/`，由现有 `.gitignore` 排除：

- `run.py`：本次运行脚本，复用现有计时、程序环境及输出比对工具。
- `inputs.json`：输入来源、文件大小、修改时间、索引计数和源码提交。
- `status.json`：运行状态、完整命令、实际导入路径、开始结束时间及计时。
- `src_{upstream,turbo}/`：两份固定的源码快照。
- `inputs/`：原始输入的软链接和新建 BAI。
- `{upstream,turbo}/log.txt`：完整日志。
- `{upstream,turbo}/{out,annotation}/`：各程序输出及转换后的注释。
- `comparison.json`：全部输出和注释对象的比较结果。

以上原始产物保留在本机，BAM、索引、源码快照和计算结果不提交到 Git。
