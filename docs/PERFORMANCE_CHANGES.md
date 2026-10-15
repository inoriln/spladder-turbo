# SplAdder 性能优化改动记录

本文记录对 SplAdder（基线：`65ceec8 version bump`，3.1.1）所做的性能优化：目标、设计原则、
逐模块改动、保证“结果不变”的方法与证据、行为变化、复现验证的方法，以及后续的代码审查与
Bug 修复（第 9 部分）。

所有改动基于上游 3.1.1（提交 `65ceec8`），在本仓库的历史里按主题分成多个提交。

---

## 1. 目标与约束

- 目标：大幅缩短 `spladder build` 的运行时间（原版比 rMATS / MAJIQ / SUPPA2 慢很多）。
- 硬约束：**计算结果不能变**。输出文件（文本、gff3、hdf5、事件 pickle）和中间文件
  （样本剪接图 pickle、合并图 pickle、注释 pickle）的内容必须与原版一致。
- 次要约束：不改变命令行接口；不引入新的必需依赖（numba 本来就是依赖）。

### 总体结果

测试数据：GSE272028 的 2 个真实样本（人，STAR 比对），GENCODE v34 全基因组，
命令 `spladder build --merge-strat merge_graphs --extract-ase -n 85`。

| 场景 | 原版 | 优化后 | 提升 |
|---|---|---|---|
| 全基因组 · 2 样本 · `--parallel 8` · 含 GTF 解析（从零开始） | 3:04:15 | **4:28–4:44**（2 次运行；最终代码 4:28） | 39–41× |
| 同上，注释 pickle 已存在（同一注释的第二次起） | 3:02:10 | 3:39–3:56（4 次运行，最终代码） | 46–50× |
| 从零开始的 CPU 时间 | 15,650 s | 1,425–1,723 s | 9–11× |
| GENCODE v34 注释解析（每个注释只做一次） | 2:05 | 48 s | 2.6× |
| chr20 · 3 样本 · 单进程 | 5:04 | 40 s | 7.5× |
| chr20 · 3 样本 · `--parallel 8` | 5:04（单进程） | 23 s | 13× |

说明：

- 原版基线已经是对原版最有利的跑法（BLAS 线程设为 1，两个样本手工拆成 2 个进程并行建图）。
- 最大单进程 RSS 约 2.6–2.9 GB（原版同一流程的峰值是 2.95 GB）；`--parallel 8` 峰值时有 19–24 个进程
  同时在跑，PSS 合计约 9 GB。
- 各阶段耗时（`--parallel 8`）：注释解析 48 s（首次）、建图约 105 s、合并约 50 s、量化约
  60 s（含基因表达）、事件验证与写出约 20 s。事件收集与量化重叠进行。

---

## 2. 瓶颈分析（改动的依据）

对原版做 py-spy 采样得到的主要结论：

1. **读取 BAM 的方式**：对每个基因、每个处理步骤都做一次随机访问 `pysam.fetch`，
   逐条 read 用纯 Python 循环过滤，读放大约 5.1 倍，占运行时间约 70%。
2. **串行**：`--parallel` 原来只覆盖事件检测与验证；建图、计数、合并都是串行的，样本逐个处理。
3. **逐基因的 Python 小循环**：`Gene.label_alt`、内含子保留判断、`Segmentgraph.from_gene`、
   `Splicegraph.uniquify`、计数等，每个基因要做几十次很小的 numpy 调用。
4. **重复工作**：同一个 217 MB 的图 pickle 被重复加载约 12 次；`compute_gene_expression` 在基因
   循环里每次重读整张 `gene_names` 表；numba 函数每个进程重新编译（约 17 s）；隐藏的 BLAS 多线程
   空转 CPU。
5. **进程池的数据传输**：事件验证与检测的进程池每个任务都 pickle 一份基因对象，主进程光序列化
   就占 20 s。

---

## 3. 设计原则

1. **单遍读取 + 内存里回答窗口查询**（`evidence.py`）：每个 contig 顺序读一遍 BAM，建立覆盖度
   阶梯函数和内含子数组，之后所有基因窗口的查询都在内存里回答，答案与逐窗口读取完全一致
   （包括 dtype 和 shape）。
2. **做不了就退回原路径**：任何不能保证等价的情况（多 BAM 合并、`--var-aware`、缺 NM 标签、
   没有索引、CRAM、异常标签……）自动退回原来的逐窗口实现，行为保持原样。
3. **原版逐字保留为参照**：每个被优化的函数，原版实现逐字放在 `tests/reference.py`
   （个别在包内作为兜底），随机/穷举/边界测试都以它为标准。
4. **并行不改变结果**：并行路径要么是互相独立的任务（按顺序合并结果），要么是对原版跨样本状态泄漏
   的精确复现（见 4.2）。
5. **能用真实数据全量比对就全量比对**：除单元测试外，用 6 万个基因的真实图做差分验证。

---

## 4. 改动清单

### 4.1 读取比对

| 文件 | 改动 |
|---|---|
| `spladder/evidence.py`（新） | `ContigEvidence`：每个 contig 一遍扫描得到各 XS 类别（无 / + / − / 其他）的覆盖度阶梯函数和内含子数组。`EvidenceScope` / `evidence_scope()`：持有已建好的 contig 证据，带内存上限（`SPLADDER_EVIDENCE_MAX_GB`，默认 6）和 LRU 淘汰。`prefetch_contigs`：用多进程并行预读 contig。`add_reads_from_evidence`：按原版 `add_reads_from_bam` 的语义回答窗口查询（包含判定 `start > win_start & stop < win_stop`、`mincount` 过滤、`unstranded` 时正负链合并排序），不支持时返回 `NotImplemented` 退回原路径。 |
| `spladder/bamscan.py`（新） | 编译版 BAM 读取器：用 BAI 索引（优先用伪 bin 37450，缺失时取全部 chunk 的包络）定位 contig，自己解压 BGZF 块（`zlib` 原始 deflate，并校验每块的 CRC32 与长度），用 numba 函数 `_scan` 直接解码比对记录，过滤顺序与条件与 pysam 循环逐条对应。对 chr1（约 400 万个比对块）从 11.9 s 降到 3.2 s（含块校验）。遇到 CRAM、非 BAI 索引、非整数 NM、非 `A` 型 XS、缺 NM、截断或损坏的文件、CIGAR 放在 `CG` 标签里的超长比对等就抛 `NotSupported` 由调用方退回 pysam（htslib 在读出记录之后做的那些处理，读取器要么复现，要么放弃，见 9.2）。 |
| `spladder/reads.py` | `_get_reads_reference`（原版逐字保留）；`_get_reads_collapsed`（逐窗口快速路径：CIGAR 只解码一次、只读需要的标签、差分数组累加覆盖度）；`_get_alignment_file`（按进程 id 区分的打开文件 LRU 缓存，避免每个基因重新打开 BAM 并加载索引）；`add_reads_from_bam` 增加 `use_cache` 参数并先询问 `evidence`；`get_intron_list` 增加 `cache` 参数。 |
| `spladder/helpers.py` | `MAX_FEASIBLE_INTRONS = 200` 提为常量（便于测试）；`make_introns_feasible` 的重新查询使用 `cache=False`。 |

### 4.2 并行框架（`spladder/parallel.py`，新）

- `ordered_map(function, tasks, workers, daemon)`：每个任务一个 fork 出来的进程，按任务顺序产出结果。
  与 `multiprocessing.Pool` 的区别是：某个 worker 被杀（例如内存不足）会报错而不是永远阻塞；
  fork 前 `gc.freeze()`，避免子进程的垃圾回收遍历并复制父进程的全部对象。
- **按样本并行建图（推测执行）**：原版 `helpers.make_introns_feasible` 会就地修改
  `options.read_filter`，导致后一个样本的建图依赖前一个样本；而且 `options.read_filter` 与
  `options.intron_retention['read_filter']` 只在第一个样本里是同一个字典对象（IR 步骤的
  `.copy()` 把它们拆开）。并行时每个样本从“没有任何样本改动过状态”的猜测状态出发，结果按样本
  顺序验收，只有当样本实际起始状态与串行运行一致时才采用；否则在主进程按串行状态重算该样本。
  猜错只影响速度，不影响结果。
- `forked_pool` / `release_forked_pool` / `shared()`：进程池的 worker 通过 fork 继承数据，
  任务里只传下标，不再 pickle 基因对象。
- `balanced_shards` / `run_shards`：把 contig 按基因数装箱成大小接近的分片，每个分片一个 fork
  进程，按分片顺序返回结果，输出由主进程回放。
- `Background`：在 fork 出的后台进程里运行一个函数，`join()` 时回放输出并抛出错误；
  `abort_background()` 在出错时结束尚未 join 的作业。
- **进程生命周期**：每个 worker（`ordered_map`、`forked_pool`、`Background`）启动时先把继承来的 SIGTERM
  处理恢复为默认动作，再用 `prctl(PR_SET_PDEATHSIG, SIGTERM)` 请内核在父进程结束时（无论怎样结束，包括
  `kill -9`）终止自己；嵌套的 worker 各自设置，形成链式终止。主进程在 `spladder()` 期间把 SIGTERM
  转成 `SystemExit`（`sigterm_unwinds`，调用者已有处理函数时不覆盖），使 `finally` 能结束后台作业、
  清理临时文件。Pool 的 worker 忽略 SIGINT，由主进程统一处理（9.2 的 B1）。
- 仅 Linux 使用 fork 并行；其他平台自动串行。

### 4.3 建图阶段

| 文件 | 改动 |
|---|---|
| `spladder/core/gen_graphs.py` | 把“插入 cassette 外显子 / 内含子保留 / 新内含子边 / 重新标注”抽成 `_augment_genes`（串行路径与原版语义一致）。`_augment_genes_in_shards`：基因的编辑只依赖同一 contig 的基因，因此按 contig 分片并行；并在主进程精确重放这些步骤对 `options` 的副作用（`read_filter` 的换入换出、`insert_intron_iterations`、`debug`、`chrm_lookup`）。`remove_se`、`sparse_bam`、单 contig 时走串行路径。 |
| `spladder/editgraph.py` | `_all_introns`、`_cassette_candidate_pairs`：cassette 外显子候选对的向量化筛选。`_retentions_of_gene`（numba）：内含子保留的覆盖度判定，浮点运算顺序、`max`/`min` 对 NaN 的行为、中位数的计算方式都与原 numpy 实现一致；输入不是常见形态（稀疏 BAM、非 float64 的 tracks、越界的区间……）时退回原循环。 |
| `spladder/classes/gene.py` | `label_alt`：四层实现——numba 内核（常见情形）→ 通用 numpy 实现（大图、非常见 dtype）→ 原版逐字实现（含长度为负的外显子等退化图）；内核复现了原版的两个怪癖（`j` 是同端点列表中的位置而非外显子编号、`idx_next` 偏移少 1）。`to_sparse` / `from_sparse`：用 numpy 替换 scipy `coo_matrix`（约 10 倍快），非常见形态仍走 scipy。`add_exon`：`np.concatenate` 替换 `np.r_` + `expand_dims`。 |
| `spladder/classes/segmentgraph.py` | `from_gene`：向量化（`searchsorted` + 广播）。 |
| `spladder/classes/splicegraph.py` | `uniquify`：向量化，没有重复外显子时提前返回。`from_gene`：用字典 + 列表代替对每对外显子扫描全部顶点的 O(V²) 循环（保持“最后一个同坐标外显子”的语义）。 |

### 4.4 合并、量化

| 文件 | 改动 |
|---|---|
| `spladder/merge.py` | 加载用 `load_pickle`；写出用 `dump_pickle`；量化且不用 `--validate-sg` 时，最终合并图在合并时就带分段图（原来由 count 阶段补算并再写一遍 217 MB 的 pickle）。 |
| `spladder/count.py` | `_count_gene`：用前缀和一次算出所有分段的均值与覆盖位置数（覆盖度是整数，所以和与原来的逐段求和精确相等），内含子计数用 int64 键 + `unique` / `bincount` / `searchsorted`。`_count_contig` / `_count_shard`：按 contig 分片并行计数；分片进程自己读自己的 contig，不再预读。`_collect_counts`：整体拼装各基因的计数矩阵（异常形态退回逐基因拼装）。`compute_gene_expression`：`gene_names` 只读一次（44.6 s → 9.9 s）。样本级并行 `parallel.count_samples`。 |

### 4.5 事件阶段

| 文件 | 改动 |
|---|---|
| `spladder/alt_splice/detect.py` | numba 函数 `cache=True`（`utils.jit_cached`）；进程池改为 fork 继承基因数据（`DETECT_BINSIZE = 25`，原来每任务 3 个基因并 pickle 基因对象）。 |
| `spladder/alt_splice/verify.py` | `verify_wrapper`：同一基因的计数只从 hdf5 读一次，行定位用二分查找；进程池改为 fork 继承（验证 exon_skip 从 36 s 降到约 10 s）；`graph_filenames`、`load_count_index`（快速索引文件，原子写出；比计数文件旧的索引视为过期并重写，见 9.2 的 B6）。 |
| `spladder/alt_splice/analyze.py` / `collect.py` | 图 pickle 用 `load_pickle` / `load_pickle_shared`，避免每个事件类型重复加载。 |
| `spladder/spladder_build.py` | 事件收集拆成 3 个后台作业（`_EVENT_GROUPS`），在图量化的同时进行（事件检测只依赖剪接图，不依赖计数）；6 种事件类型并发验证（`_analyze_events`，预加载只读取剪接图，不产生文件）；`_spladder` 外面包一层，出错时结束后台作业；`todo` 中同名样本只算一次；`_prep_workdir` 用 `exist_ok=True`。 |

### 4.6 I/O、注释解析、杂项

| 文件 | 改动 |
|---|---|
| `spladder/utils.py` | `load_pickle`（加载大 pickle 时暂停 GC，约 2 倍快）；`dump_pickle`（暂停 GC，写临时文件后改名，跟随符号链接）；`load_pickle_shared`（同一进程内未改变的文件不重复加载）；`frozen_heap`；`jit_cached`；`unique_rows` / `ismember` 向量化。 |
| `spladder/init.py` | GTF 解析：只对需要的行解析属性（第一遍只解析基因行，第二遍只解析 exon 行）；注释 pickle 用 `dump_pickle`。 |
| `spladder/spladder.py` | 命令行启动时把 `OPENBLAS_NUM_THREADS` / `OMP_NUM_THREADS` / `MKL_NUM_THREADS` 默认设为 1（用户已设置则不覆盖）。 |
| `spladder/settings.py` | 创建输出目录用 `exist_ok=True`（并行启动时不再有“先检查后创建”的竞态）。 |
| `spladder/spladder_prep.py` | 注释加载用 `load_pickle`。 |
| `docs/source/general.rst` | 新增 “Run time, parallel processing and memory” 一节。 |

### 4.7 新增环境变量

| 变量 | 作用 |
|---|---|
| `SPLADDER_EVIDENCE_MAX_GB` | 保存读证据的内存上限（默认 6）；超限后按 LRU 淘汰，用到时重读。 |
| `SPLADDER_NO_EVIDENCE=1` | 关闭单遍读证据，回到逐基因读取。 |
| `SPLADDER_NO_BAMSCAN=1` | 不用编译版 BAM 读取器，改用 pysam。 |
| `SPLADDER_DEBUG_BAMSCAN=1` | 编译版读取器出任何错误时直接抛出（默认静默退回 pysam）。仅用于调试。 |

---

## 5. 保证“结果不变”的方法与证据

### 5.1 方法

1. **参照实现**：被改动的函数的原版实现逐字保存（`tests/reference.py`，`label_alt` 和
   `get_reads` 的原版也留在包内作为兜底）。
2. **差分测试**：随机 / 穷举 / 边界输入，比较输出的值、dtype、shape（`tests/compare.py`
   的 `deep_equal` 对整个对象递归比较）。
3. **变异测试**：对关键的向量化代码注入错误，检查测试能否抓住；据此补强了测试生成器，
   剩余存活的变异体均已证明与原代码等价。
4. **真实数据回放**：把新旧实现分别用在真实剪接图（每个文件 60,669 个基因）上逐基因比较。
5. **端到端比对**：原版与新版各跑一遍，逐项比较全部输出。

### 5.2 证据（均在最终代码上）

| 验证 | 结果 |
|---|---|
| 测试套件 | Python 3.13：250 个测试全部通过。最旧的受支持环境（Python 3.9 / numpy 1.26 / numba 0.60 / pysam 0.23）：244 个通过，另 6 个端到端测试失败——仓库里的期望结果 pickle 由 numpy 2 写成，numpy 1.26 读不了，原版代码在同一环境同样失败 |
| 全基因组（从 GTF 开始）与原版金标准比对 | 35 项（文本 / gff3 / hdf5 / 事件 pickle / 剪接图）全部一致；两个样本的建图 pickle 整个对象一致；合并图 pickle（含分段图）整个对象一致；24 个事件输出文件一致；注释 pickle 一致；输出文件列表相同。最终代码、从零开始（含 GTF 解析，注释 pickle 与原版逐字节相同）：4:28，CPU 1,425 s |
| 逐字节比较（全基因组，38 个输出文件） | 29 个逐字节相同（6 个 gff3、8 个 hdf5、12 个事件 pickle、3 个 `.quick_*` 索引）；6 个 `.txt.gz` 只有 gzip 头里的时间戳（第 5–8 字节）不同，解压后的文本逐字节相同；3 个图 pickle（两个样本图、合并图）字节不同，但加载后对象完全相等（值、dtype、shape），合并图上两边的对象数量（248 万个）和可变对象的共享情况也一致（都没有被多处引用的列表、字典或数组）——差别只在 pickle 记录“相等的不可变对象（如 numpy dtype）是同一个还是几个”，文件大 0.7%。chr20 的配置矩阵里字节差异的范围更清楚：样本图在串行和 `--parallel 4` 下都字节相同（只有一个 contig，不分片）；合并图在量化的配置里不同（合并时就建分段图），在不量化的工作流（cohort、chunked 合并）和 `merge_bams` 里字节相同。也就是说，只有经过 contig 分片（多个 contig 且 `--parallel` 大于 1 的样本图）或合并时建分段图的 pickle 字节会变 |
| 配置矩阵（原版 vs 新版） | 42 组：20 种参数配置 + 6 种工作流（cohort、chunked-merge、并行、sparse 等）+ 16 组 `--parallel 4` 变体；40 组逐文件一致，2 组（一次传 3 个 BAM 的 `single` 与 `qmode single`）在原版中就以同样方式失败 |
| 端到端随机流水线 | 随机注释 + 随机 reads（新外显子、保留内含子、可变位点、跳过两个外显子的新接头、互斥外显子），原版 vs 新版的串行与 `--parallel 3`：Python 3.13 上 24 个种子（共确认 6 种事件类型约 1,000 个事件）、Python 3.9 环境上 8 个种子，全部逐文件一致 |
| 边界输入差分 | 14 种非常规输入 × 3 种运行方式 = 42 组：26 组逐文件一致，6 组（缺 NM 标签、没有索引）两边以同样方式失败，10 组原版崩溃而新版完成（B7、B8、B10），没有任何输出不同 |
| `label_alt` 差分 | 216 万个随机图与原版一致；真实数据 3 × 60,669 个基因一致 |
| BAM 读取器 | 两个真实样本的所有 contig × 3 种过滤设置共 678 组合，与 pysam 循环逐一一致；800 个随机 BAM（每条记录随机 CIGAR 操作、flag、标签类型）共 12,396 次比较无差异 |
| 进程生命周期 | 主进程（带嵌套 worker）被 SIGTERM / SIGINT / SIGKILL 结束后，残留进程数为 0（修复前 SIGTERM 后残留 8–10 个） |
| chr20 单进程 | 39 个输出文件与原版一致 |

---

## 6. 行为变化与兼容性

1. `--parallel N` 的含义扩展：并行样本级建图、同一样本内按 contig 并行编辑与计数、事件收集与量化
   重叠、事件类型并发验证。**瞬时进程数可能超过 N**（事件收集后台作业和并发验证会叠加）。
   结果不依赖进程数。
2. 日志：并行时各任务的输出被捕获并按任务顺序回放，不再实时交织；指定 `--logfile` 时，
   并行样本写入 `<logfile>.<样本名>`，分片写入 `<logfile>.<样本名>.shard<k>`。
3. 合并图 pickle 在合并时就带分段图（量化且不用 `--validate-sg`）。最终文件内容与原版相同。
4. 大文件（图 pickle、注释 pickle、快速索引文件）改为“写临时文件 + 改名”，读它的其他进程总能读到完整的
   旧版或新版；临时文件名为 `<目标>.tmp<随机标识>`。
5. 首次运行要编译 numba 函数，约 30 s（chr20 单进程实测：首次 68 s，之后 40 s），之后读缓存（缓存在 `__pycache__`，不可写时 numba 退回用户
   缓存目录）。
6. 仅 Linux 使用 fork 并行；Windows / macOS 自动串行（读证据与编译内核不受影响）。
7. 环境里安装的 `spladder` 是独立副本，不会自动带上这些改动；使用前需 `pip install .` 或设置
   `PYTHONPATH`。
8. 图 pickle 的字节与原版不同：样本图和合并图的 pickle 加载后与原版完全相等，但文件的字节不同（见 5.2）。如果有流程用 md5 / sha 校验这些文件，校验值会变；文本、gff3、hdf5 和事件 pickle 的字节不变。
9. 发行名是 `spladder-turbo`，Python 包名和命令仍是 `spladder`（输出和 pickle 与上游互通），另外提供命令 `spladder-turbo`。不要与上游 SplAdder 装在同一环境，二者提供同名的包和命令。
10. 验证范围：真实数据端到端验证用的是 GENCODE v34 + GSE272028（人，STAR）。CRAM 输入走 pysam
   路径，`--var-aware` 走逐窗口路径，二者没有做全基因组端到端验证。

---

## 7. 复现验证

- 单元 / 集成测试：`pytest tests -q`（约 2–3 分钟）。包含：与原版逐字参照实现的差分测试、BAM 读取器与 pysam 的
  对比（`tests/test_bamscan.py`）、进程生命周期（`tests/test_parallel.py::test_workers_end_with_the_main_process`，
  需要 Linux 和 `pgrep`）、第 9 部分每个缺陷的回归测试，以及 `tests/test_static.py`（所有函数用到的全局名都有定义）。
- 参照实现与差分工具：`tests/reference.py`、`tests/compare.py`（`compare_dirs`、`deep_equal`）、
  `tests/synthetic.py`（合成 BAM 生成器）。
- 最旧的受支持环境（Python 3.9、numpy 1.26、numba 0.60、pysam 0.23）：同样用 `pytest tests -q`。其中
  `tests/test_end_to_end.py` 的 6 个测试会失败，因为 `tests/testcase_*/results_*` 里的期望结果 pickle 是 numpy 2
  写的，numpy 1.26 读不了；原版代码在同一环境下同样失败。
- 端到端验证用到的脚本（配置矩阵、与金标准逐项比对、随机流水线模糊测试、边界输入差分、`label_alt` 模糊测试、
  BAM 读取器全 contig 比对与随机 BAM 模糊测试等）放在开发机的临时目录里，还没有整理进仓库；
  方法见 5.1 和 9.1。

---

## 8. 已知限制与后续方向

- 建图阶段剩余的主要成本：分片内的 `deepcopy` / `isequal`（约占分片 CPU 的 1/3）、父进程接收分片
  结果和写出样本 pickle（各约 10 s）、内含子列表仍在父进程里串行计算。
- 合并仍是串行的（每增加一个样本约 +20 s）；样本数远大于进程数的队列，合并和每样本约 550 CPU-s
  是主要成本，适合做分层合并。
- 图阶段与计数阶段对同一个 BAM 各扫描一遍（过滤条件不同），可以合并成一遍。
- 事件对象的构造与后处理（`Event`、`post_process_event_struct`）仍是逐事件的 Python 循环。
- 运行被强行中断（`kill -9`、断电）时，可能留下 `*.tmp<12 位十六进制>` 临时文件（`dump_pickle` 写到一半的内容）。
  它们不影响下一次运行，可以直接删除；正式文件要么是旧的完整版本，要么是新的完整版本。
- 内存：读证据按 contig 保存，每个进程默认上限 6 GB（`SPLADDER_EVIDENCE_MAX_GB`），构造时的瞬时峰值约为证据
  大小的 3 到 4 倍；进程数多而样本很深（几亿条比对）时，需要调小这个上限。

---

## 9. 代码审查与严重 Bug 修复记录

审查范围：本次全部改动（21 个修改的源文件，3 个新模块 `parallel.py` / `evidence.py` / `bamscan.py`）
和它们依赖的原有代码路径。共发现并修复 **16 个缺陷**：5 个是本次优化引入的，11 个是原版就有的。
**每个缺陷都有回归测试，并且这些测试在修复前的代码上会失败**（逐个验证过）。
另外有几处会改变输出结果的原版问题，按“结果不变”的约束只记录、不修改（9.3）。

### 9.1 审查方法

| 方法 | 做法 | 发现 |
|---|---|---|
| 静态检查 | pyflakes；用 `dis` 遍历所有函数加载的全局名，找出模块里没有定义的名字（星号导入后面的也能查到，pyflakes 做不到）；`ast.parse(feature_version=(3, 9))` 检查语法兼容；以“警告即错误”编译全部文件 | B12–B16 |
| 通读改动 | 逐文件读 diff。对每个涉及进程、信号、临时文件的函数问：被杀、磁盘满、同名并发、输入异常时会怎样。对 BAM 读取器逐项核对 htslib 在“读出记录”之后还做了什么 | B1–B4 |
| 边界输入差分 | 14 种非常规输入（空 BAM、无任何边的注释、链方向为 `.` 的基因、单外显子基因、完全相同的转录本、BAM 缺 contig、没有索引、缺 NM / XS 标签、`MT` contig、重叠基因、大量小 contig……）× 3 种运行方式（串行 / `--parallel 3` / `--parallel 3 -c 1`），原版与新版比较退出码与全部输出 | B7、B8、B10 |
| 进程生命周期实验 | 起一个带嵌套 worker 的主进程，分别用 SIGTERM / SIGKILL / SIGINT 结束它，数残留进程 | B1 |
| 模糊测试 | ① BAM 读取器 vs pysam 循环：800 个随机 BAM 文件（每条记录随机 CIGAR 操作、flag、标签类型）、12,396 次比较；② 端到端：随机注释 + 随机 reads（新外显子、保留内含子、可变位点、跳过两个外显子的新接头、互斥外显子），原版 vs 新版（串行与 `--parallel 3`），最终代码上 32 个种子（Python 3.13 上 24 个、Python 3.9 环境上 8 个），修改过程中另跑过约 100 个 | 全部一致（B3、B4 不是这样发现的，是通读时发现的，见下） |
| 最旧的受支持环境 | Python 3.9 + numpy 1.26 + numba 0.60 + pysam 0.23（`requirements.txt` 声明的下限）上跑全部测试和端到端模糊种子 | 无新问题（说明见 9.3 的第 6 条） |

模糊测试没有发现 B3 和 B4：随机 BAM 里没有超过 65535 个 CIGAR 操作的记录，也没有被改过的压缩块。这两个是
通过列出“pysam 在读出记录之后做的所有事”并逐条对照发现的。教训是：与参照实现做差分，只能发现输入
生成器覆盖到的差异，所以除差分之外还要逐条核对参照实现的隐含行为。

### 9.2 缺陷清单

严重程度：**高** = 资源失控或可能静默给出错误结果；**中** = 特定输入下崩溃或结果出错；**低** = 只在不常用路径出错。

| 编号 | 严重程度 | 来源 | 位置 | 一句话 |
|---|---|---|---|---|
| B1 | 高 | 优化引入 | `parallel.py`、`spladder_build.py` | 主进程被杀后 worker 变成孤儿进程，继续占用 CPU 和内存 |
| B2 | 中 | 优化引入 | `utils.py` | `dump_pickle` 的临时文件名只含 PID，同名并发写会写出损坏的文件 |
| B3 | 高（概率极低） | 优化引入 | `bamscan.py` | 忽略 `CG` 标签：超长 CIGAR 的比对被静默读错 |
| B4 | 中 | 优化引入 | `bamscan.py` | 不校验 BGZF 块：损坏的 BAM 被静默读下去，pysam 则报错 |
| B5 | 低 | 优化引入 | `spladder_build.py` | 并发验证的预加载有副作用：无事件时多生成 `.quick_*` 文件，坏文件时比串行更早报错 |
| B6 | 中 | 原版 | `alt_splice/verify.py` | 计数文件被重写后，旧的 `.quick_*` 索引仍被使用，结果静默出错 |
| B7 | 中 | 原版 | `classes/gene.py` | 基因链方向不是 `+`/`-`（如 GTF 里的 `.`）时 `TypeError` |
| B8 | 中 | 原版 | `count.py` | 注释里没有任何一条边（全是单外显子基因）时写计数文件 `AttributeError` |
| B9 | 中 | 原版 | `merge.py` | 各样本基因集合不同时合并 `IndexError` / `ValueError` |
| B10 | 低 | 原版 | `helpers.py` | `log_progress` 在 total 为 0 时 `ZeroDivisionError` |
| B11 | 低 | 原版 | `alt_splice/analyze.py` | `options.spladderfile` 路径括号错位，`TypeError` |
| B12 | 中 | 原版 | `spladder_prep.py` | `prep_sparse_bam_*` 用了没有导入的 `pysam` 等名字，`NameError` |
| B13 | 低 | 原版 | `viz/graph.py` | `plot_graph(highlight=...)` 用了没有导入的 `patches`，`NameError` |
| B14 | 低 | 原版 | `spladder_viz.py` | 缺 `compute_gene_expression` 的导入，`NameError` |
| B15 | 低 | 原版 | `alt_splice/quantify.py` | 不支持的事件类型要抛的 `Error` 类不存在，得到 `NameError` |
| B16 | 低 | 原版 | `count.py` | `get_size_factors(kind='uq')` 用的 `scoreatpercentile` 在去掉 scipy 时没有替换 |

#### 本次优化引入的缺陷

**B1 · 主进程被杀后 worker 变成孤儿进程**（高）

- 现象：用 `kill`（SIGTERM）或作业调度器结束 `spladder build` 的主进程后，所有 worker 以及 worker 自己
  启动的子进程继续跑到各自任务结束（PPID 变成 1）。实测 SIGTERM 之后残留 8–10 个进程，每个 0.5–2 GB，
  继续写 `*.part<pid>` 临时文件。在集群上这会占着节点的内存和 CPU，拖垮同节点的下一个作业。
- 原因：fork 出来的 worker 不会随父进程退出；`multiprocessing` 的 daemon 机制只在父进程“正常退出”时
  才会终止子进程；主进程对 SIGTERM 用默认动作，立刻消失，`finally` 里的清理一行也不执行。
- 修复：见 4.2“进程生命周期”。每个 worker 启动时恢复默认的 SIGTERM 动作并设置 `PR_SET_PDEATHSIG`
  （设置后再检查父进程是否还在，堵住竞态）；主进程在运行期间把 SIGTERM 转成 `SystemExit`，让
  `finally` 结束后台作业并清理临时文件。
- 补充：中止后的清理还漏了一种文件。worker 在写样本图 pickle 的途中被结束，会留下 `dump_pickle` 的临时文件
  `<图>.part<主进程 pid>.tmp<随机>`，原来的清理只认 `.part<pid>` 结尾。这是机器很忙时 `test_failing_sample_aborts_like_the_serial_run`
  偶发失败才发现的（Python 3.9 的那次测试运行里遇到）；现在连同这种文件一起删除。
- 同一类问题还有事件 pickle：事件收集现在与量化重叠，主进程出错或被中断时，后台作业可能正在写
  `merge_graphs_<类型>_C3.pickle`。这个文件存在就表示“该事件类型已完成”，被截断的文件会让下一次运行读到
  坏数据。这些文件现在也用 `dump_pickle` 原子写出（`test_event_files_are_written_in_one_piece`）。
- 验证：`test_workers_end_with_the_main_process`（SIGTERM / SIGINT / SIGKILL 各一）：主进程 + 2 个 worker +
  4 个嵌套 worker，发信号后 15 s 内一个进程都不应残留；SIGTERM 和 SIGINT 时 `finally` 的清理要执行。
  修复前 SIGTERM 实测残留 8–10 个进程。`test_temporary_files_of_workers_that_were_ended_are_removed`：确定性地
  检查上面那种遗留文件被删除，别的运行的临时文件不动。

**B2 · `dump_pickle` 的临时文件重名**（中）

- 现象：两个写同一目标文件的进程如果 PID 相同（共享目录上的不同主机、容器的 PID 命名空间），会互相截断
  对方的临时文件：轻则 `os.replace` 抛 `FileNotFoundError`，重则把两个进程交错写成的损坏 pickle
  改名成正式文件。
- 原因：临时文件名是 `<目标>.tmp<pid>`。
- 修复：`<目标>.tmp<随机 12 位十六进制>`。同时删掉 `verify.py` 里重复实现的 `_dump_atomically`，
  统一用 `dump_pickle`。
- 验证：`test_dump_pickle_of_concurrent_writers_is_never_seen_partially`：3 个写进程（强制 `os.getpid()`
  相同）各写 25 轮 30 万元素的数组，主进程同时不停地读，读到的必须是完整文件，结束后没有临时文件残留。

**B3 · BAM 读取器忽略 `CG` 标签**（高，概率极低）

- 现象：BAM 规范规定，CIGAR 操作数超过 65535 时，记录里的 CIGAR 字段写成占位的 `<读长>S<参考长度>N`，
  真实的 CIGAR 放进 `CG:B:I` 标签。htslib（pysam）读出记录时会自动还原，编译版读取器直接解码记录里
  的字段，把这条读段当成“一个软剪切 + 一个横跨整个参考长度的内含子”：覆盖度缺失，凭空多出一个巨大
  的内含子。不报错，结果静默出错。只有超长读段（长读长测序的比对）会遇到。
- 原因：写读取器时按“记录里写了什么”解码，没有核对 htslib 在读出之后还做了什么。
- 修复：遇到“恰好两个操作，依次是 `S` 和 `N`，且带 `CG` 标签（或标签区格式不对）”的记录，读取器放弃，
  该 contig 退回 pysam 路径。没有 `CG` 标签的同样形状的记录照常解码（htslib 也不会改它）。
- 验证：`test_cigar_in_the_cg_tag_is_left_to_the_general_reader`：含一条 70,001 个操作（35,000 个内含子）的
  读段；修复前 `scan_contig` 不抛 `NotSupported`，读出的是占位的那一个巨大内含子。

**B4 · BAM 读取器不校验 BGZF 块和记录**（中）

- 现象：BAM 文件损坏（磁盘位翻转、传输损坏）时，htslib 校验每个块的 CRC32，报 `CRC32 checksum mismatch`，
  pysam 抛 `OSError`。编译版读取器不校验，读到被改过的字节照样往下解码，可能静默得到错误结果。另外块尾
  声明的解压后大小也没有检查，被改成很大的值时会按它去分配内存。
- 修复：每个块检查 `isize ≤ 65536`（BGZF 上限）、解压后长度等于 `isize`、`crc32` 与块尾记录一致；不一致
  就抛 `NotSupported`，由 pysam 路径给出原版的行为。代价：读 BAM 约多 19% 的时间（chr1：2.7 s → 3.2 s），
  整个流程 CPU 约多 2%。
- 验证：`test_a_damaged_block_is_an_error_as_with_pysam`（用不压缩的块做出“解压不报错、只有校验和不符”的损坏
  文件；以及声明大小不合理的块）。修复前读取器对损坏文件不报错。
- 同类的还有记录本身：CIGAR 的查询长度与序列长度不一致（序列长度不为 0）的比对，htslib 读到会报错（`error -4`），
  读取器原来照常解码。现在遇到这样的记录就放弃，由 pysam 给出同样的错误。验证：
  `test_cigar_of_another_length_than_the_sequence_is_an_error_as_with_pysam`。

**B5 · 并发验证的预加载有副作用**（低）

- 现象：6 种事件类型并发验证之前，主进程先加载共享的剪接图并生成计数索引（`.quick_*`）。这有两个后果：
  ① 图文件或计数文件缺失、图文件损坏时，异常比串行运行更早抛出（串行运行只有真有事件要验证时才读这些文件）；
  ② 没有任何事件时，串行运行不会生成 `.quick_*` 文件，并行运行却生成了，`--parallel` 改变了输出目录里的
  文件集合（边界输入差分里比较 `many_small_contigs`、`identical_transcripts`、`single_exon_genes` 的串行与并行
  输出目录时发现）。
- 修复：预加载只读取已经存在的剪接图（读失败就留给任务去失败），不生成任何文件；计数索引由真有事件要验证的
  任务生成（原子写出，几个任务同时生成也是安全的）。
- 验证：`test_preload_of_the_verification_tasks_only_reads_what_exists` / `_makes_no_files` /
  `_survives_a_graph_that_cannot_be_read`。

#### 原版就有的缺陷

**B6 · 过期的快速索引文件**（中）

- 现象：验证事件时，每个 `.count.hdf5` 旁边会缓存三个 `.quick_*` 索引文件（基因编号等）。计数文件被重写
  （换了样本、重新量化）而旧的 `.quick_*` 还在时，旧索引被继续使用，事件验证读到错位的计数，结果静默出错。
- 修复：`load_count_index` 只使用不比计数文件旧的索引，否则从计数文件重新读取并重写。
- 验证：`test_count_index_follows_the_count_file`（修复前读到的是旧索引）。

**B7 · 基因链方向不是 `+` / `-`**（中）

- 现象：GTF 中链方向为 `.` 的基因，`Gene.__init__` 抛 `TypeError: unsupported operand type(s) for %`。
- 原因：`warnings.warn('...%s...') % str(strand)` —— 括号位置错了，`%` 作用在 `warn` 的返回值 `None` 上。
- 修复：括号挪对；警告照原意发出，链方向按 `+` 处理。
- 验证：`test_gene_on_an_unknown_strand_is_read_as_plus_strand`；边界输入差分中的 `unstranded_dot_gene`
  （原版崩溃，新版完成）。

**B8 · 注释里没有任何一条边**（中）

- 现象：注释只有单外显子基因时，写计数 hdf5 抛 `AttributeError: 'list' object has no attribute 'dtype'`。
- 原因：没有任何基因有边时 `counts['edges']` 仍是空的 Python 列表，写文件时取了它的 `dtype`。
- 修复：先 `np.asarray`。
- 验证：`test_annotation_without_any_edge`；边界输入差分中的 `single_exon_genes`。

**B9 · 各样本基因集合不同时合并崩溃**（中）

- 现象：`merge_genes_by_splicegraph` 合并基因集合不同的样本时（例如用不同注释建的图，或某个样本多出基因），
  遇到前面样本没有的基因：新基因排在所有已有基因之后时 `IndexError`，排在中间时 `ValueError`。
- 原因：两处循环条件 `g_idx <= genes2.shape[0]` 差一，可以访问 `genes2[len]`；追加新基因用的
  `np.c_[genes2, gene]` 把一维对象数组和单个对象拼成二维，形状不符。
- 修复：条件改为 `<`，追加改用 `np.append`（下一个样本开始前本来就会重新按名字排序）。
- 验证：`tests/test_merge.py` 的 24 个测试（随机的基因集合 × 2–4 个样本，结果必须是所有样本基因的并集，
  边的计数等于持有该边的样本数）。

**B10 · `log_progress` 在 total 为 0 时除零**（低）

- 现象：`--parallel N -v` 且没有任何事件需要验证时，验证阶段最后一次进度输出 `ZeroDivisionError`。
- 修复：total 不大于 0 时直接返回。
- 验证：`test_progress_of_nothing_is_not_shown`、`test_verbose_run_with_processes_and_nothing_to_detect`。

**B11 · `analyze.py` 括号错位**（低）

- 现象：`options.spladderfile` 路径下 `TypeError`：`pickle.load(open(f), 'rb')` —— `'rb'` 传给了
  `pickle.load` 而不是 `open`。
- 修复：括号挪对。验证：`test_gene_information_of_the_events_from_a_given_graph_file`。

**B12 · `spladder_prep.py` 缺少导入**（中）

- 现象：`prep_sparse_bam_filtered` / `prep_sparse_bam_full` 在需要读 BAM 头（还没有 `chrm_lookup`）时
  `NameError: name 'pysam' is not defined`（`append_chrms`、`parse_header` 也一样）。
- 修复：补上三个导入。验证：`test_preparing_sparse_alignment_files`（用合成 BAM 跑两个准备函数）。

**B13–B16 · 不常用路径上的未定义名字**（低）

- B13 `viz/graph.py`：`plot_graph(..., highlight=...)` 用的 `patches` 没有导入（文件里导入的是 `mpatches`）。
  测试：`test_graph_plot_with_a_highlighted_range`。
- B14 `spladder_viz.py`：`_load_size_factors` 在基因表达文件缺失时调用 `compute_gene_expression`，但没有导入。
  测试：`tests/test_static.py`。
- B15 `alt_splice/quantify.py`：`quantify_from_counted_events` 和 `get_event_ids` 对不支持的事件类型要
  `raise Error(...)`，`Error` 不存在，得到的是 `NameError` 而不是说明原因的错误。改为 `ValueError`。
  测试：`test_event_ids_of_an_unknown_event_type_are_refused_with_an_error`。
- B16 `count.py`：`get_size_factors(kind='uq')` 用的 `scoreatpercentile` 在“去掉 scipy 依赖”时没有替换。
  改用 `np.percentile`（线性插值，与 scipy 的默认方式一致）。测试：`test_size_factors_by_upper_quartile`。

`tests/test_static.py` 把这一类问题堵住：它检查包里所有函数加载的全局名都有定义，今后再出现笔误会直接失败。
当前只放行 9 处原版就未完成的代码（见 9.3 的第 2 条），这个名单本身也被检查，代码一旦被修好或删除，
名单里多出来的条目会让测试失败。

### 9.3 发现但没有改动的问题

1. **`Gene.label_alt` 的两个下标笔误（会改变输出，所以保留）。** 一是 `j` 在代码里是“同端点外显子列表中的
   位置”，却被当作外显子编号用；二是 `idx_next` 的偏移少 1。优化后的内核逐位复现了这两个行为。
   在真实数据的 60,669 个基因里，按“代码注释所述的本意”修正后，有 727 个基因（1.2%）的
   `is_alt_spliced` 会不同；`is_alt` 和所有事件输出不受影响。要修需要维护者决定，并接受输出变化。
2. **原版未完成的代码（修不了，也没有入口）。** `--re-infer-sg` 调用的 `infer_splice_graph` 在整个代码库
   里不存在（只有注释掉的桩），选项一用就 `NameError`；`write_events_txt(anno_fn=...)` 里的
   `load` / `sort` / `genes` 是没有移植完的 MATLAB 风格代码；`reads.get_all_data_uncollapsed`、
   `editgraph.infer_splice_graph_caller` 里是笔误和未完成的函数体；`quantify.quantify_from_graph` 一族是
   Python 2 时代的遗留（文本模式读写 pickle、`np.ravel_multi_index` 的 `dims` 传了整数），没有任何地方
   调用它。这些在 `tests/test_static.py` 的名单里。
3. **`get_reads(collapse=False)`** 在窗口越过 contig 末尾时崩溃（原版就有）。
4. **文档与代码不符：** `--chunked-merge` 的层级说明与实际行为不一致。
5. **一次传入多个 BAM 的 `single` 与 `qmode single`** 在原版中就会崩溃，新版以同样方式崩溃（配置矩阵里
   这 2 组两边都是退出码 1）。
6. **旧环境下的测试：** 仓库里 `tests/testcase_*/results_*` 的期望结果 pickle 是 numpy 2 写的，numpy 1.26
   读不了（`No module named 'numpy._core.numeric'`），所以 6 个端到端测试在 numpy 1.26 下失败。原版代码
   在同一环境下同样失败，与本次改动无关。
7. **`--sparse-bam`** 的覆盖度与普通路径略有不同（原版如此）。
