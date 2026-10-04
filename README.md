# Cell Flow

单细胞组学分析管线：从表达矩阵完成质量控制、归一化、降维聚类与差异表达分析，输出可复现的分析结果与图表数据。

## 范围

本仓库从零开始实现上述方向的可用工具，不依赖外部同类实现。

## 状态

初始基线：只有本说明，尚无实现。

## 输入

入口仍为 `cell-flow analyze --input <路径> --output-dir <结果目录>`，
新增可选 `--input-format`（取值 `tsv` 或 `mtx`，默认 `tsv`）、
`--gene-sets` 与 `--replicate-metadata`。`--key value` 与 `--key=value`
两种写法均支持。跨样本批次校正见下文 `--cell-metadata`。

### `tsv`（默认，UTF-8 制表符基因计数矩阵）

首列是唯一基因 ID，其余列名是唯一细胞 ID，取值为非负整数 UMI 计数。
`--input` 可直接指向未压缩文本，也可指向其 gzip 压缩文件：是否压缩按
文件内容（gzip 魔数）自动识别，与文件名或后缀无关。gzip 仅允许单成员
且流后无尾随数据；gzip 头、CRC、长度校验失败或流被截断一律报输入错误
（退出码 2）。压缩只改变承载方式，解压后的矩阵与行列顺序不变。

### `mtx`（标准 10x MatrixMarket 稀疏目录）

`--input` 指向一个目录，其中须含三个文件。每个文件可保持原名，也可
整体改用 gzip 压缩名，但三个文件**必须同为未压缩或同为 gzip**，混用
报输入错误（退出码 2）：

- `matrix.mtx` 或 `matrix.mtx.gz`：仅接受 MatrixMarket 头
  `%%MatrixMarket matrix coordinate integer general` 或
  `coordinate real general`；行是基因、列是细胞、索引从 1 开始；
  零值可省略，显式值必须能无损解析为有限非负整数（`real` 字段同样
  要求恰为整数值，超大整数按精确十进制解析，不做 float64 舍入）。
- `barcodes.tsv` 或 `barcodes.tsv.gz`：每个非空行是一个细胞 ID，
  按文件行序进入分析。
- `features.tsv` 或 `features.tsv.gz`：每个非空数据行取第一列作为
  基因 ID，按文件行序进入分析。

同一基名的压缩与未压缩文件不得并存。gzip 仅允许单成员且无尾随数据；
gzip 头、CRC、长度校验失败或流被截断一律报输入错误（退出码 2）。

`matrix.mtx` 声明的行/列维度必须分别等于基因/细胞 ID 数量。
缺文件、格式头不支持、索引越界、重复坐标、空或重复 ID、声明元素数不符、
非法数值一律报输入错误（退出码 2）；`--input-format` 取其他值报配置错误
（退出码 3）。任何输入失败都不会创建或修改结果目录。

两种格式给出的等价矩阵产生完全一致的分析结果（排序、数值、图表数据）。
`tsv` 运行的 `run.json` 字段与基线一致（gzip 载体的 `input.sha256`
按压缩文件的原始字节计算）；`mtx` 运行仅在其 `input` 结构中
额外记录 `"input_format": "mtx"` 以及三个输入文件的实际文件名
（压缩时为 `*.gz`）和各自按原始字节计算的 SHA-256。除 `input`
来源字段外，矩阵等价的两种 TSV 或两种 MTX 来源，每个结果文件
逐字节一致。

### `--metadata`（可选，细胞分组差异表达）

`--metadata <路径>`（`--key value` 与 `--key=value` 均可）指向一个
UTF-8 制表符文本，表头恰为 `cell_id` 和 `group` 两列；`cell_id`
唯一且与表达矩阵的全部细胞一一对应（不多不少），`group` 非空。
文件可为纯文本或单成员 gzip（按内容识别，与文件名无关）。内容不合法
或 gzip 非单成员流报输入错误（退出码 2），且不改动任何已有结果；
参数缺值、未知参数或调用不合法报配置错误（退出码 3）。质控后非空
分组少于两个报数据错误（退出码 4）。

提供元数据时在既有结果之外新增两个文件（其余结果与无元数据运行
逐字节一致）：

- `group_markers.tsv`：只用质控后保留细胞与基因的 log 归一化表达。
  每个分组先做 one-vs-rest，再按 group 升序对每对分组做两两检验。
  前三列为 `comparison_type`、`group_a`、`group_b`（one-vs-rest 的
  `group_b` 为空），其余列沿用成对 marker 的口径：`gene_id`、
  `mean_in_a`、`mean_in_b`、`log_fc_a_vs_b`（均值差）、`t_stat`、
  `p_value`（Welch t 检验）、`p_value_adj`（每个比较内 BH 校正）。
  每个比较内按校正 P 值升序、`log_fc_a_vs_b` 降序、`gene_id` 升序排序。
- `group_marker_chart.tsv`：每个比较取前 20 个基因，含 `rank` 列。

`run.json` 的 `input.metadata` 记录元数据文件名、原始字节 SHA-256
与各分组细胞数；既有 `input` 字段不变。计算完成后事务性发布；目标
非空报 `OutputPathError`（退出码 5），同输入同版本的新增结果逐字节
一致。

### `--batch-metadata`（可选，批次均值中心化校正）

`--batch-metadata <路径>`（`--key value` 与 `--key=value` 均可）指向
一个 UTF-8 制表符文本，表头恰为 `cell_id` 和 `batch` 两列；
`cell_id` 唯一且与表达矩阵的全部细胞一一对应（不多不少），`batch`
非空。文件可为纯文本或单成员 gzip（按内容识别，与文件名无关）。
内容不合法或 gzip 非单成员流报输入错误（退出码 2），且不改动任何
已有结果；参数缺值、未知参数或调用不合法报配置错误（退出码 3）。
质控后保留细胞覆盖的批次少于两个报数据错误（退出码 4）。不提供
该参数时，全部行为与结果文件与基线逐字节一致。

校正按基因进行：质控与 log 归一化沿用既有规则，随后每个归一化值
减去对应批次保留细胞的基因均值，再加回全部保留细胞的基因总均值
（batch mean centering）。校正值进入高变基因选择、PCA、聚类、
markers、成对 markers、`--metadata` 分组差异表达及图表数据；
`normalized_expression.tsv` 仍写未校正值。提供批次元数据时在既有
结果之外新增两个文件（行列顺序与既有表达矩阵一致）：

- `batch_corrected_expression.tsv`：校正后的表达矩阵，行（保留基因）
  列（保留细胞）顺序与 `normalized_expression.tsv` 相同。
- `batch_summary.tsv`：按 `batch_id` 升序，每行给出 `batch_id`、
  `n_cells`（批次内保留细胞数）、`total_counts`、`detected_genes`、
  `mitochondrial_fraction`，后三项为批次内保留细胞的均值。

`run.json` 的 `input.batch_metadata` 记录批次元数据文件名、原始字节
SHA-256 与质控前后各批次细胞数；`parameters` 增加
`"batch_mean_centering": true`；既有字段不变。`--batch-metadata`
可与 `--metadata` 并用；TSV、MTX、gzip 各输入承载方式下等价矩阵的
结果一致性保持不变。目标非空或暂存、写出、发布失败报
`OutputPathError`（退出码 5），同输入同版本的新增结果逐字节一致。

### `--cell-metadata`（可选，跨样本批次校正）

`--cell-metadata <路径>`（`--key value` 与 `--key=value` 均可）指向一个
UTF-8 制表符文本：首列必须是细胞条码列 `cell_id`，并须包含样本标识列
（默认 `sample_id`，可用 `--sample-column <列名>` 指定）与批次标签列
（默认 `batch`，可用 `--batch-column <列名>` 指定）；其余列任意，
原样保留并透传到公开细胞结果，不覆盖、不重命名任何既有字段。
`cell_id` 唯一且与表达矩阵的全部细胞一一对应（不多不少）；样本标识与
批次标签均非空——批次标签缺失或为空同样报错，不静默归入未知批次。
文件可为纯文本或单成员 gzip（按 gzip 魔数识别，与文件名无关）。
`--cell-metadata` 与 `--batch-metadata` 不能同时使用；`--batch-column`
与 `--sample-column` 不能为空、不能为 `cell_id`、也不能彼此相同，
违反报配置错误（退出码 3）。

表达矩阵含重复细胞条码、元数据含重复细胞条码、元数据缺少样本标识列或
批次标签列、表达矩阵与元数据的细胞集合不一致、批次标签缺失或为空，
一律抛出 `ValueError`（命令行表现为输入错误，退出码 2），消息指出具体
冲突类型，且在启动降维聚类之前失败，不生成任何部分成功结果。

校正发生在质量控制之后、降维聚类之前：质控与 log 归一化沿用既有规则，
随后每个归一化值减去对应批次保留细胞的基因均值，再加回全部保留细胞的
基因总均值（batch mean centering，与 `--batch-metadata` 同一口径）。
校正值进入高变基因选择、PCA、聚类、markers、成对 markers、
`--metadata` 分组差异表达、`--gene-sets` 评分及图表数据；
`normalized_expression.tsv` 仍写未校正值。质控后仅一个批次时不施加
任何扰动（校正即恒等），全部既有结果文件与无批次运行逐字节一致。
不提供该参数时，全部行为与结果文件与基线逐字节一致；已有的质量控制
阈值、聚类参数、差异表达口径、结果结构和图表数据字段均不变。

提供该参数时（无论一个还是多个批次）在既有结果之外新增以下文件：

- `cell_metadata.tsv`：保留细胞的元数据透传，表头与字段保持输入原样
  （含全部额外列），行序同 `clusters.tsv`；每行给出该细胞的原始批次。
- `batch_corrected_expression.tsv`：校正后的表达矩阵，行列顺序与
  `normalized_expression.tsv` 相同；单批次时与未校正值逐字节一致
  （即后续聚类和差异表达实际使用数据的等价公开产物）。
- `batch_summary.tsv`：按 `batch_id` 升序，每行给出 `batch_id`、
  `n_cells`（批次内保留细胞数）、`total_counts`、`detected_genes`、
  `mitochondrial_fraction`，后三项为批次内保留细胞的均值。
- `batch_pca_scatter.tsv`：按批次着色的降维坐标，列为
  `cell_id`、`batch`、`PC1`、`PC2`，取自实际用于聚类的 PCA 空间。
- `batch_mixing.tsv`：校正前后各一个批次混合分数（`stage` 为
  `before`/`after`）。分数按每个细胞的 15 个最近邻（细胞数不足时取
  细胞数减一）中来自其他批次的细胞占比计算，再对全部细胞求平均；
  距离为 PCA 坐标上的欧氏距离，并列按细胞顺序确定。`before` 在未校正的
  归一化 → 高变基因 → PCA 坐标上计算，`after` 在实际用于聚类的（校正后）
  PCA 坐标上计算；单批次时两者均为 0。

`run.json` 的 `input.cell_metadata` 记录元数据文件名、原始字节 SHA-256、
批次/样本列名与质控前后各批次细胞数；`parameters` 增加
`"batch_column"`、`"sample_column"`，实际施加校正时另增
`"batch_mean_centering": true`；新增顶层 `batch_correction` 记录校正
方法（`batch_mean_centering` 或 `none_single_batch`）、批次数与校正
前后混合分数；既有字段不变。相同输入、参数与随机种子结果逐字节一致，
未显式指定种子时使用固定默认种子，结果同样确定。目标非空或暂存、
写出、发布失败报 `OutputPathError`（退出码 5）。

### `--gene-sets`（可选，基因集评分）

`--gene-sets <路径>`（`--key value` 与 `--key=value` 均可）指向一个
UTF-8 制表符文本，表头恰为 `set_id`、`gene_id` 两列；每行一个
成员关系，字段非空、`(set_id, gene_id)` 组合唯一，且至少一条数据行。
文件可为纯文本或单成员 gzip（按 gzip 魔数识别，与文件名无关）。基因 ID
与表达矩阵首列**精确匹配**，不做大小写、别名、前缀转换；矩阵外基因
不计分但计入集合总基因数。文件不是 UTF-8、表头不符、字段为空、成员
重复、没有数据行，或 gzip 多成员、尾随数据、截断、CRC/长度错误，
一律报输入错误（退出码 2）且不改动结果。`--gene-sets` 缺值、
`--gene-sets=` 空路径或出现未知参数报配置错误（退出码 3）。

评分只用质控后保留的细胞与基因：有 `--batch-metadata` 时取批次均值
中心化值，否则取 log 归一化值。每个集合的 score 是该集合与保留基因
交集内表达值的算术平均；任一集合交集为空即报数据错误（退出码 4），
不产出结果。提供该参数时在既有结果之外新增两个文件（其余结果与不提供时
逐字节一致，`run.json` 除外）：

- `gene_set_scores.tsv`：列为 `set_id`、`n_genes_total`、
  `n_genes_used`、`cell_id`、`score`，按 `set_id` 升序、集合内
  沿用保留细胞原顺序排列。
- `gene_set_score_chart.tsv`：以 `cell_id`、`cluster` 开头，各集合列
  按 `set_id` 升序排列，细胞顺序沿用 `clusters.tsv`。

`run.json` 的 `input.gene_sets` 记录基因集文件名、原始字节 SHA-256
与各集合总基因数（`set_sizes`）和实际使用基因数（`set_sizes_used`，
均按 `set_id` 升序）；`parameters` 增加 `"gene_set_scoring": true`；
既有字段不变。未提供 `--gene-sets` 时全部行为与基线逐字节一致。
浮点格式沿用既有结果口径（最短往返表示），重复运行逐字节一致；目标非空
或暂存、写出、发布失败报 `OutputPathError`（退出码 5）。

### `--replicate-metadata`（可选，按生物学重复的 pseudobulk 差异表达）

`--replicate-metadata <路径>`（`--key value` 与 `--key=value` 均可）指向
一个 UTF-8 制表符文本，表头恰为 `cell_id`、`sample_id`、`group` 三列；
每个输入细胞恰好一行，`cell_id` 唯一且与表达矩阵的全部细胞一一对应
（不多不少），`sample_id` 与 `group` 非空，同一 `sample_id` 只能归属一个
`group`。文件可为纯文本或单成员 gzip（按 gzip 魔数识别，与文件名无关）。
表头不符、行列数不一致、ID 重复、覆盖不全或含矩阵之外的细胞、空样本或空
分组、同一样本归属多个分组，或 gzip 多成员、尾随数据、截断、CRC/长度错误，
一律报输入错误（退出码 2）且不改动结果目录；`--replicate-metadata` 缺值、
空路径或出现未知参数报配置错误（退出码 3）。质控后某样本无保留细胞、某
group 有效重复少于两个或非空 group 不足两个报数据错误（退出码 4）。

分析只用既有质控保留的细胞与基因，输入、QC、聚类、批次校正、基因集评分与
其余图表行为均不变；未提供该参数时，全部结果与当前版本逐字节一致。
pseudobulk 按 `sample_id` 对**原始计数**求和形成，样本文库归一到 10000 后
取 log1p；样本文库总计数按该样本在保留基因上的计数和计算。差异表达以样本
为观测单位：每个 group 先做 one-vs-rest，再按 group 升序两两比较，
执行 Welch t 检验、双侧 P 值，并在每个比较内跨基因做 BH 校正。提供该参数
时在既有结果之外新增三个文件：

- `pseudobulk_expression.tsv`：保留基因（沿用保留顺序）× 样本的 pseudobulk
  log 归一化表达；样本按其输入细胞在矩阵列序中的首次出现排列，质控后无保留
  细胞的样本不出现。
- `pseudobulk_group_markers.tsv`：前三列为 `comparison_type`、`group_a`、
  `group_b`（one-vs-rest 的 `group_b` 为空），其余列沿用现有 group marker
  统计列：`gene_id`、`mean_in_a`、`mean_in_b`、`log_fc_a_vs_b`（均值差）、
  `t_stat`、`p_value`、`p_value_adj`。每个比较内按校正 P 值升序、
  `log_fc_a_vs_b` 降序、`gene_id` 升序排列。
- `pseudobulk_marker_chart.tsv`：前三列同上，随后为 `rank`、`gene_id`、
  `log_fc_a_vs_b`、`p_value`、`p_value_adj`、`neg_log10_p_adj`；每个比较
  取前 20 个基因，沿用同一排序，`rank` 自 1 起在每个比较内单独编号。

`run.json` 的 `input.replicate_metadata` 记录元数据文件名、原始字节
SHA-256、质控前各样体细胞数（`sample_sizes`，按输入细胞首次出现顺序）、
质控后各样体保留细胞数（`sample_sizes_after_qc`，按样本 ID 升序）、
质控前各分组细胞数（`group_sizes`，按 group 升序）、质控前后重复数
（`replicates_before_qc`/`replicates_after_qc`）与质控后各 group 的有效
重复数（`group_replicates_after_qc`，按 group 升序）；`parameters` 增加
`"pseudobulk_de": true`；`stage_counts` 增加 `pseudobulk_samples`、
`pseudobulk_groups`、`pseudobulk_comparisons` 与
`pseudobulk_marker_tests`（比较数 × 保留基因数）。其余字段不变。该参数可与
`--metadata`、`--batch-metadata`、`--cell-metadata`、`--gene-sets` 并用；
新增文件沿用事务性发布与最短往返浮点格式，同输入同版本逐字节一致；
目标目录问题仍由 `OutputPathError`（退出码 5）报告。

### `--detect-doublets`（可选，双细胞识别与过滤）

无值开关 `--detect-doublets` 启用双细胞识别；数值参数
`--expected-doublet-rate F` 给定期望双细胞比例，只接受有限的 `[0, 1)`
值，缺省 `0.08`。只给 `--expected-doublet-rate` 而未启用
`--detect-doublets`、参数缺值、取值非法（含 `nan`/`inf`/区间外/非数值）、
给无值开关传值或出现未知参数，一律报配置错误（退出码 3）。不启用时
输入、全部结果文件与退出码与基线逐字节一致。

启用后以现有 QC 候选细胞（质控通过、尚未做双细胞过滤的细胞）的原始计数
子矩阵（全部输入基因）计算确定性 `doublet_score`：用固定 `--seed` 的
内置随机数抽取等量候选细胞两两配对，把两个细胞的计数谱相加，再按候选
细胞平均文库深度做无放回降采样（超几何抽样），得到“两个候选细胞计数
叠加”的合成谱；观测谱与合成谱都按各自文库大小归一到 10000 后取
`ln(x + 1)`。每个观测细胞在全部基因的 log 归一化空间中取最近的
`min(15, 2*候选数-1)` 个谱（欧氏距离平方，并列按原列序、观测谱优先于
合成谱），其中合成谱所占比例即该细胞的分数。评分区分两个候选细胞计数
叠加谱与单细胞谱，同输入、同 seed、同参数下跨运行逐字节复现。

目标标记数按 `floor(候选数 * rate)` 计算，并保证至少留一个细胞
（rate 为 0 时不标记任何细胞）。所有候选细胞按分数降序排名，同分按原
矩阵列序，`doublet_rank` 自 1 起；前目标数个细胞 `doublet_flag` 为真
并被过滤。过滤后按 `--min-cells` 在最终细胞上重新计算基因保留，后续
归一化、HVG、PCA、聚类、差异表达、批次校正、基因集评分、pseudobulk 与
图表数据一律只使用最终细胞与最终基因；被过滤细胞在 `cells.tsv` 的
`retained` 记为 `false`。

启用时在既有结果之外新增两个文件（其余结果沿用最终细胞与基因，与不启用
时同口径）：

- `doublet_scores.tsv`：列 QC 候选细胞（原矩阵列序），表头为
  `cell_id`、`doublet_score`、`doublet_rank`、`doublet_flag`、
  `retained_after_doublet_filter`；`doublet_rank` 自 1 起，
  `doublet_flag` 为是否被过滤，浮点用最短往返表示。
- `doublet_score_chart.tsv`：把分数划分为 20 个等宽箱，列为
  `bin_start`、`bin_end`、`cell_count`；空箱保留，边界确定
  （分数全相同时宽度为 0，全部细胞计入首箱）。

`run.json` 的 `parameters` 增加 `"detect_doublets": true` 与
`"expected-doublet-rate"`；`stage_counts` 增加 `doublet_candidates`
（QC 候选细胞数）、`doublets_flagged`（标记过滤数）、
`cells_after_doublet_filter`（最终细胞数）与
`genes_after_doublet_filter`（最终基因数）；既有 `cells_after_qc` 等
字段仍记录双细胞过滤前的 QC 结果。QC 候选细胞或过滤后细胞少于两个、
过滤后无保留基因报数据错误（退出码 4）。新增结果与既有结果一起事务性
发布：目录非空、暂存或发布失败报 `OutputPathError`（退出码 5），分析
阶段失败不创建或改动结果目录。该功能可与 `--metadata`、
`--batch-metadata`、`--cell-metadata`、`--gene-sets`、
`--replicate-metadata` 并用，TSV、MTX、gzip 各承载方式下等价矩阵结果
一致。

## 约定

- 公开行为以 README 与源码为准。
- 后续需求在此基线上增量实现。
