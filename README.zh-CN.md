# spladder-turbo

[English](README.md) | [简体中文](README.zh-CN.md)

**[SplAdder](https://github.com/ratschlab/spladder) 3.1.1 的非官方性能优化分支（fork）**，用于从 RNA-Seq 比对结果中检测并量化可变剪接事件。

spladder-turbo 的目标是可以直接替换 SplAdder：命令行相同，输出文件相同，并且就我们所能测试的范围而言，结果也相同，但耗时只是原来的一小部分。和 SplAdder 一样，它根据注释构建剪接图，用比对数据中的证据扩充剪接图，从图中提取可变剪接事件（外显子跳跃、多外显子跳跃、内含子保留、可变 3' 和 5' 剪接位点、互斥外显子），并对事件进行量化。

> **状态：开发版 0.1.0，可从 GitHub 安装，尚未发布到 PyPI。**
> spladder-turbo **与 SplAdder 的作者及其所属机构没有关联，也未获得他们的认可或背书**。方法、事件的定义和文件格式都是他们的工作。如果你使用了本软件，请引用 SplAdder 的论文（见[引用](#引用)）。

**目录：** [性能](#性能) · [优化了什么、在哪里、为什么](#优化了什么在哪里为什么) · [结果是否相同](#结果是否相同) · [安装](#安装) · [使用](#使用) · [测试](#测试) · [文档](#文档) · [引用](#引用) · [许可与署名](#许可与署名)

## 性能

`spladder build` 的完整流程（图的扩充、合并、量化、事件的提取与验证），两个人类 RNA-Seq 样本（[GSE272028](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272028)，样本 SRR29779619 和 SRR29779620；用 STAR 2.7.7a 比对到 hg38，读长 85 bp，每个样本约 4400 万条比对），GENCODE v34 注释，全基因组，从 GTF 文件开始：

| | SplAdder 3.1.1 | spladder-turbo 0.1.0 | |
|---|---|---|---|
| 耗时（wall-clock） | 3 小时 04 分 | 两次运行分别为 4 分 28 秒和 4 分 33 秒（`--parallel 8`） | 40-41 倍 |
| CPU 时间 | 15,650 秒 | 1,425 秒和 1,605 秒 | 10-11 倍 |
| 最大单个进程的内存 | 2.9 GB | 3.0 GB | |

同一流程只跑 chr20，三个样本（SRR29779619 到 SRR29779621）：

| | SplAdder 3.1.1 | spladder-turbo 0.1.0 | |
|---|---|---|---|
| 单进程 | 5 分 04 秒 | 40 秒 | 7.5 倍 |
| `--parallel 8` | （未测） | 23 秒 | 13 倍（与 SplAdder 单进程相比） |

如何理解这些数字：

- 测试机器是一台 Linux 服务器，两颗 AMD EPYC 7601（128 线程），314 GB 内存；Python 3.13、numpy 2.3、numba 0.62、pysam 0.23。耗时来自 `/usr/bin/time -v`。
- SplAdder 的跑法是对它最有利的：线性代数库限制为单线程；全基因组时把两个样本的建图步骤手工拆成两个同时运行的进程。SplAdder 自带的 `--parallel` 只分担事件的检测与验证。我们**还没有**在全基因组上测过 SplAdder 加 `--parallel 8` 的耗时。
- 单进程的加速（7.5 倍）来自：每个 contig 只用编译版读取器读一遍比对文件、编译内核、向量化代码。其余来自让样本、contig 和事件类型在多个进程里并行（仅 Linux）。CPU 时间反映总共少做了多少工作。
- 第一次运行要编译 numba 函数，大约半分钟（上面 chr20 的运行是 68 秒，之后是 40 秒）。编译结果会缓存在磁盘上。
- `--parallel 8` 时，全基因组的运行总共占用约 9 GB 内存，同时最多约 24 个进程。
- 这些测量用的脚本和原始日志还没有放进本仓库（方法和单次测量见 [docs/PERFORMANCE_CHANGES.md](docs/PERFORMANCE_CHANGES.md)）。

## 优化了什么、在哪里、为什么

spladder-turbo 改变的是 SplAdder *怎样*计算，而不是*计算什么*。本节列出各项优化：在代码的哪个位置、解决 SplAdder 3.1.1 的什么问题、改成了什么做法、带来了什么效果。每个被改动的函数的完整记录见 [docs/PERFORMANCE_CHANGES.md](docs/PERFORMANCE_CHANGES.md)。

### SplAdder 3.1.1 的时间花在哪里

对 SplAdder 3.1.1 在人类数据 chr20 子集上做采样分析（py-spy），再结合代码，得到五点结论：

1. **读取比对约占运行时间的 70%。** 对每个基因，在五个处理步骤（样本的内含子、cassette 外显子、内含子保留、新内含子边、计数）里，SplAdder 都用 `pysam.fetch` 取出基因窗口里的比对，再用 Python 循环逐条过滤。与基因重叠的一条比对，每个步骤都要解码一次，数据相当于被读了约 5 遍。
2. **几乎全是串行。** `--parallel` 只分担事件的检测和验证。建图、合并、计数在一个进程里，样本一个接一个地处理。
3. **每个基因有大量很小的 Python 和 numpy 调用。** 基因的标注、内含子保留的判断、分段图的构建、重复外显子的去除、计数，每个基因都要做几十次小的 numpy 调用，而基因有 6 万个。
4. **重复劳动。** 同一个 217 MB 的图文件在一次运行里被加载约 12 次；`compute_gene_expression` 在遍历基因的循环里反复读取整张基因名表；每个进程都重新编译一遍 numba 函数（约 17 秒）；线性代数库对很小的矩阵也开着隐藏的线程。
5. **向工作进程传送数据。** 事件检测和验证的进程池，每个任务都要 pickle 基因对象。仅主进程里的序列化就用了 20 秒。

### 现在的时间花在哪里

全基因组，两个样本，从 GTF 文件开始（与上面的表一致）。spladder-turbo 的各步耗时取自一次运行（4 分 33 秒）日志行上的时间戳。并行阶段的任务在结束时才输出日志，所以各步的边界精确到几秒。

| 步骤 | SplAdder 3.1.1 | spladder-turbo（`--parallel 8`） |
|---|---|---|
| 解析注释（每个注释只做一次） | 2 分 05 秒 | 42 秒 |
| 两个样本的图扩充 | 1 小时 34 分（两个进程） | 1 分 53 秒 |
| 合并、计数、基因表达、事件 | 1 小时 28 分 | 1 分 58 秒：合并 43 秒，计数 47 秒，基因表达 13 秒，验证与写出 15 秒（事件在计数期间检测） |
| **合计** | **3 小时 04 分** | **4 分 33 秒** |

### 各项优化

路径相对于 `spladder/`。“新”表示 SplAdder 里没有的文件。

**A. 读取比对**

| 位置 | SplAdder 3.1.1 的问题 | spladder-turbo 的做法 | 效果 |
|---|---|---|---|
| `evidence.py`（新）、`reads.py`（`add_reads_from_bam`） | 每个基因、每个步骤都重新取一遍窗口，并用 Python 解码比对（每条比对约 5 次）。 | 对每个 contig 顺序读一遍，建立紧凑的结构：按读段类别（XS 标签）分的覆盖度阶梯函数，以及带读段计数的有序内含子数组。之后每次窗口查询都在内存里用二分查找回答，值、dtype 和 shape 与原来一致。内存有上限（`SPLADDER_EVIDENCE_MAX_GB`，LRU 淘汰）。 | 消除了主要开销里的读放大。 |
| `bamscan.py`（新） | 即使只读一遍，pysam 每条比对也要创建 Python 对象并做多次属性查找，约为解压文件所需时间的 4 倍。 | 针对有索引的 BAM 文件的编译（numba）读取器：用 BAI 索引定位 contig，解压 BGZF 块（校验 CRC32 和大小），按相同的顺序、用相同的过滤条件解码记录。遇到 CRAM、其他索引、异常的标签、htslib 会拒绝或展开的记录、损坏的文件时放弃，改用 pysam。 | 人类样本的 chr1（约 400 万个比对块）：11.9 秒降到 3.2 秒。两个样本的全部 contig：CPU 时间是 pysam 循环的三分之一。 |
| `reads.py`（`_get_alignment_file`、`_get_reads_collapsed`） | 仍然逐窗口读取的情形（多个 BAM 文件、`--var-aware`）对每个基因都重新打开 BAM 文件（含头部和索引），并转换每条比对的全部标签。 | 打开的文件按进程缓存（LRU）；CIGAR 只解码一次；只读需要的标签；覆盖度用差分数组累加。 | 回退路径更快。 |

**B. 注释与剪接图**

| 位置 | SplAdder 3.1.1 的问题 | spladder-turbo 的做法 | 效果 |
|---|---|---|---|
| `init.py`（`init_genes_gtf`）、`classes/splicegraph.py`（`from_gene`） | 解析了 290 万行 GTF 的全部属性；构建基因的图时，对每一对外显子都循环扫描所有顶点。 | 只解析用到的行的属性（第一遍只看基因行，第二遍只看 exon 行）；用字典构建图。生成的注释文件与原来相同。 | GENCODE v34：2 分 05 秒降到 42-48 秒。 |
| `classes/gene.py`（`label_alt`） | 每个基因一个 Python 循环，里面是很多小的 numpy 调用。 | 常见的图用 numba 内核，大图和少见的 dtype 用 numpy 实现，退化的图用 SplAdder 的原始代码。内核复现了原版的两个下标怪癖，所以标注结果不变。 | 已与原版对比：216 万个随机图，以及 3 x 60,669 个真实基因。 |
| `classes/gene.py`（`to_sparse`、`from_sparse`、`add_exon`） | 每个基因的邻接矩阵都经过 scipy 的 `coo_matrix` 转换；每个外显子用 `np.r_` 和 `expand_dims`。 | 用 numpy 算出同样的条目（少见的形状仍走 scipy）；改用 `np.concatenate`。 | 转换约快 10 倍。 |
| `classes/segmentgraph.py`（`from_gene`）、`classes/splicegraph.py`（`uniquify`） | 对分段、对所有外显子对做循环。 | 向量化（`searchsorted`、广播）；没有重复时 `uniquify` 提前返回。 | |
| `editgraph.py`（`insert_intron_retentions`、`insert_cassette_exons`） | 对每一对相连的外显子、每一对内含子，在 Python 循环里调用多次 numpy。 | 内含子保留的覆盖度判断改为 numba 函数，运算顺序与原来的 numpy 代码一致，所以浮点结果相同。cassette 外显子的候选对用向量化运算找出。输入不常见时仍用原来的循环。 | |
| `core/gen_graphs.py`、`core/spladdercore.py` | 基因的编辑在一个进程里进行。 | 一个基因的编辑只依赖同一 contig 的基因，所以把 contig 凑成大小相近的分片，放在不同进程里编辑；这些步骤对 options 的副作用由主进程重放。各 contig 的证据并行读取。 | 属于建图阶段：连同其他改动，1 小时 34 分降到 1 分 53 秒。 |

**C. 合并与计数**

| 位置 | SplAdder 3.1.1 的问题 | spladder-turbo 的做法 | 效果 |
|---|---|---|---|
| `merge.py` | 合并后的图先写出，再由计数步骤加载并构建分段图，然后把 217 MB 的文件再写一遍。 | 如果合并后的图接下来要做量化且之后不再改动（没有 `--validate-sg`），就在合并时建好分段图。 | 少写一遍大文件。 |
| `count.py`（`_count_gene`、`_count_contig`、`_count_shard`、`_collect_counts`） | 每个基因的每个分段都用 numpy 调用计算平均覆盖度；每条边都扫描一遍完整的内含子列表；样本和 contig 一个接一个地计数。 | 一个基因所有分段的平均值由累加和得出（精确，因为覆盖度是整数，数值相同）；内含子计数用有序查找；不同的 contig 和样本在不同进程里计数；所有基因的计数一次性拼装。 | |
| `count.py`（`compute_gene_expression`） | 在遍历基因的循环里，每次都从 HDF5 文件读取整张 `gene_names` 表。 | 只读一次。 | 44.6 秒降到 9.9 秒。 |

**D. 事件**

| 位置 | SplAdder 3.1.1 的问题 | spladder-turbo 的做法 | 效果 |
|---|---|---|---|
| `alt_splice/detect.py`、`alt_splice/verify.py` | 进程池的每个任务都接收 pickle 过的基因对象；每个事件都要重新从 HDF5 文件读取它所在基因的计数。 | 工作进程通过 fork 继承图，任务里只传下标；一个基因的计数对它所有事件只读一次，行用二分查找定位；numba 函数缓存在磁盘上。 | 外显子跳跃的验证：36 秒降到约 10 秒。主进程里 20 秒的序列化没有了。 |
| `spladder_build.py` | 事件在量化之后才检测，各事件类型依次验证。 | 事件只依赖剪接图而不依赖计数，所以在图量化的同时用三个后台作业收集；六种事件类型同时验证。 | 事件检测不再作为单独的一步出现。 |

**E. 调度与输入输出**

| 位置 | SplAdder 3.1.1 的问题 | spladder-turbo 的做法 | 效果 |
|---|---|---|---|
| `parallel.py`（新） | 样本不能并行处理。不能直接让样本并排运行，因为 `make_introns_feasible` 会就地修改 `options.read_filter`，串行运行时后一个样本依赖前面的样本。 | 基于 fork 的并行框架。各样本的图推测式地并行计算，每个样本从“没有任何样本改动过过滤器”时应有的状态出发；结果按样本顺序验收，只有该状态猜对时才采用，否则按顺序重算该样本。worker 意外死亡会报错；worker 随主进程一起结束；输出和错误由主进程显示。 | 结果与串行运行完全一样。 |
| `utils.py`（`load_pickle`、`dump_pickle`、`load_pickle_shared`） | 剪接图由数百万个小对象组成：垃圾回收使加载慢了约一倍；一次运行把同一个图文件加载约 12 次；被杀掉的运行可能留下写了一半的文件。 | 加载和写入时暂停垃圾回收；文件先写成临时文件再改名；同一进程不重复加载没有变化的文件。 | 加载约快 2 倍，写出约快三分之一。 |
| `utils.py`（`jit_cached`） | 每个进程都重新编译 numba 函数（约 17 秒）。 | 编译结果缓存在磁盘上。 | 第一次运行要编译（约半分钟），之后的运行和所有工作进程不再编译。 |
| `spladder.py` | 线性代数库（OpenBLAS、OpenMP、MKL）对很小的矩阵也开线程，只是白白占用 CPU，并且多个进程同时运行时会让核数超额。 | 默认一个线程，用户自己设置了环境变量则不覆盖。 | |

### 怎样保证结果不变

- **复现，而不是近似。** 编译内核按原来代码的浮点运算顺序计算。用累加和代替求和的地方，被加的值都是整数，所以和是精确的。
- **做不到的退回 SplAdder 的代码。** 快速路径没有覆盖的情形（多个比对文件、`--var-aware`、CRAM、缺少 NM 标签、异常的标签、损坏的文件、退化的图）都由 SplAdder 的代码处理。
- **原版保留作参照。** `tests/reference.py` 里是被优化函数的 SplAdder 3.1.1 实现的逐字拷贝；测试用随机和对抗性的输入与它们比对。
- **推测要经过检查。** 样本的并行处理，只有每个样本的起始状态与串行运行一致时才被采用。
- **保留原版的特性。** SplAdder 里看起来不是有意为之、但会改变输出的行为（例如 `Gene.label_alt` 里的），spladder-turbo 也照样保留。

### 没有改的

方法本身：图扩充的规则与阈值、事件的定义、验证事件的标准、量化、命令行选项和输出格式。

## 结果是否相同

目标是 `spladder build` 的每个输出，内容都与 SplAdder 3.1.1 相同。到目前为止做过的比对：

- **全基因组**（上面的运行）：全部 38 个输出文件（文本、GFF3、HDF5、事件文件、各样本的剪接图和合并后的图）以及注释文件。38 个文件中有 29 个逐字节相同（GFF3、HDF5、事件文件和索引文件），注释文件也是。6 个压缩文本文件（`.txt.gz`）只有 gzip 头里的时间戳不同，文本内容相同。3 个图文件（`genes_graph_*.pickle`）的字节不同，但加载后的对象完全相等。
- **42 种命令行配置**（置信水平、合并策略、事件类型、输出格式、稀疏 BAM 文件等），在 chr20 子集上用两个程序各跑一遍：结果相同，或以同样的方式失败（有两种配置在 SplAdder 里同样失败）。
- **50 多个随机生成的数据集**（随机的注释和读段，含新外显子、保留的内含子、可变位点和互斥外显子），两个程序都跑一遍，串行和并行都有：结果相同。
- **编译版 BAM 读取器**与 pysam 比对：两个真实样本的每个 contig（678 种 contig 与过滤条件的组合）以及 800 个随机 BAM 文件：结果相同。
- 本仓库的 **250 个测试**，用随机和对抗性的输入，把被优化的函数与 SplAdder 3.1.1 实现的逐字拷贝（`tests/reference.py`）作比对。

已知的差异：

- spladder-turbo 修复了 SplAdder 3.1.1 的一些缺陷（注释没有任何边或链方向未知时崩溃、合并基因集合不同的样本、少用代码里未定义的名字、过期的索引文件等）。SplAdder 能正常运行的地方，输出相同；SplAdder 崩溃的地方，spladder-turbo 可以继续运行。清单见 [docs/PERFORMANCE_CHANGES.md](docs/PERFORMANCE_CHANGES.md)（第 9 部分）。
- SplAdder 里会改变输出的特性被保留，没有修复（例如 `Gene.label_alt`）。
- 使用 `--parallel` 时，并行任务的输出一个任务接一个任务地显示，`--logfile` 按样本分别写出。
- 并行处理需要 Linux；其他系统上一切串行运行（更慢，但结果相同）。

尚未验证：双端数据、STAR 以外的比对器、人类以外的物种、三个以上样本的队列、CRAM 文件、`--var-aware`、全基因组的 `--sparse-bam`，以及 Linux 以外的操作系统。代码的设计是在这些情形下与 SplAdder 保持一致（凡是不能确认的都退回 SplAdder 的代码），但目前只在小的合成数据上测试过。切换之前，请先在你自己数据的一小部分上比较两个程序的输出（见[从 SplAdder 切换](#从-spladder-切换)），发现差异请提交 issue。

## 安装

spladder-turbo 需要 Python 3.9 或更新版本，以及 [requirements.txt](requirements.txt) 里的包。建议用 conda 之类的包管理器管理依赖。

```bash
pip install git+https://github.com/inoriln/spladder-turbo.git
```

或者在本仓库的克隆里：

```bash
pip install .
```

**不要把 spladder-turbo 装进已经有 SplAdder 的环境。** 两者都提供 Python 包 `spladder` 和命令 `spladder`，后装的会覆盖先装的文件。请使用单独的环境。命令 `spladder-turbo` 与 `spladder` 是同一个程序。

## 使用

命令行与 SplAdder 相同；全部选项和输出文件的说明见 [SplAdder 文档](https://spladder.readthedocs.io)。最小的例子：

```bash
spladder build -a annotation.gtf -o results -b sample1.bam,sample2.bam --parallel 8
```

`--parallel N` 用 N 个进程处理各样本的图、计数，以及事件的收集和验证。进程数短时间内可能超过 N。结果不依赖于 N。使用快速读取器时，比对文件必须已排序并建好索引（`.bai`）；其他文件仍可读取，但用 pysam，较慢。

环境变量：

| 变量 | 作用 |
|---|---|
| `SPLADDER_EVIDENCE_MAX_GB` | 每个进程保存比对信息的内存上限，默认 6 |
| `SPLADDER_NO_BAMSCAN=1` | 用 pysam 读取比对，不用编译版读取器 |
| `SPLADDER_NO_EVIDENCE=1` | 像 SplAdder 那样逐基因读取比对 |

### 从 SplAdder 切换

在各自的环境里运行两个程序，用你数据的一小部分（例如几个样本的一条染色体）。两次运行使用相同的选项，各自使用独立的输出目录和独立的注释文件副本（注释会被转换成它旁边的一个 `.pickle` 文件）。然后用我们使用的函数（[tests/compare.py](tests/compare.py)）比较两个输出目录：它比较文本文件、HDF5 数据集和 pickle 文件里的对象，并返回差异的列表。在本仓库的克隆里，在 spladder-turbo 的环境中运行：

```bash
python -c "from tests.compare import compare_dirs; print(compare_dirs('run_spladder', 'run_turbo') or 'identical')"
```

## 测试

```bash
pip install pytest
pytest tests -q
```

大约需要三分钟。在 numpy 1.x 下，`tests/test_end_to_end.py` 里有 6 个测试会失败，因为 SplAdder 测试数据的期望结果是用 numpy 2 写出的；SplAdder 本身也同样失败。

## 文档

- [SplAdder 文档](https://spladder.readthedocs.io)：命令行与输出文件。
- [docs/PERFORMANCE_CHANGES.md](docs/PERFORMANCE_CHANGES.md)：改了什么、为什么改、怎样保证结果相同，以及代码审查的记录（目前只有中文）。
- [ATTRIBUTION.md](ATTRIBUTION.md)：署名与许可。
- [README.md](README.md)：本页的英文版。

## 引用

spladder-turbo 是 SplAdder 的衍生作品。请引用 SplAdder 的论文：

Kahles A, Ong CS, Zhong Y, Rätsch G. SplAdder: identification, quantification and testing of alternative splicing events from RNA-Seq data. *Bioinformatics* 2016;32(12):1840-1847. doi:[10.1093/bioinformatics/btw076](https://doi.org/10.1093/bioinformatics/btw076)

为了可复现，请同时说明你使用了 spladder-turbo 及其版本（`spladder --help` 会显示版本）。

## 许可与署名

SplAdder 由 Cheng Soon Ong、Gunnar Raetsch 和 Andre Kahles 编写（见 [AUTHORS](AUTHORS)），以 BSD 3-Clause 许可发布（[COPYRIGHT](COPYRIGHT)）。本分支的修改以相同的许可发布（[LICENSE](LICENSE)）。版权所有者和贡献者的名字没有被用来为本分支背书或做宣传。详情见 [ATTRIBUTION.md](ATTRIBUTION.md)，其中包括一个在上游带有 GPL 声明的源文件的说明。

本软件按“原样”提供，不附带任何形式的担保；具体条款见许可证。
