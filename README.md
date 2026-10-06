# Cell Flow

单细胞组学分析管线：从表达矩阵完成质量控制、归一化、降维聚类与差异表达分析，输出可复现的分析结果与图表数据。

## 范围

本仓库从零开始实现上述方向的可用工具，不依赖外部同类实现。

## 状态

初始基线：只有本说明，尚无实现。

## 输入

入口仍为 `cell-flow analyze --input <路径> --output-dir <结果目录>`，
新增可选 `--input-format`（取值 `tsv` 或 `mtx`，默认 `tsv`）、
`--gene-sets`、`--replicate-metadata` 与 `--cell-type-reference`；
marker 基因集富集通过 `--enrich-markers` 开关（仅在 `--gene-sets` 基线上
生效）与 `--enrichment-alpha`、`--enrichment-min-log-fc` 控制；
pseudobulk 基因集分组差异通过 `--pseudobulk-gene-set-de` 无值开关启用
（仅在同时提供 `--gene-sets` 与 `--replicate-metadata` 时生效）；
簇级样本差异丰度通过 `--differential-abundance` 无值开关启用（仅在提供
`--replicate-metadata` 时生效）。
`--key value` 与 `--key=value` 两种写法均支持。跨样本批次校正见下文
`--cell-metadata`。

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

### `--pseudobulk-gene-set-de`（可选，pseudobulk 基因集分组差异）

`--pseudobulk-gene-set-de` 是无值开关，只在同时提供 `--gene-sets` 与
`--replicate-metadata` 时可用：缺任一基线参数即启用开关，或把开关写成
带值形式（如 `--pseudobulk-gene-set-de=true`），一律报配置错误
（退出码 3）。开关本身无取值写法，`--key value`/`--key=value` 的取值
参数写法不受影响。未启用时，全部结果文件与 `run.json` 与基线逐字节一致。

评分只使用质控后保留的基因与细胞汇总出的 pseudobulk log1p 值
（与 `pseudobulk_expression.tsv` 同口径）：样本按其输入细胞在矩阵列序中的
首次出现排列，即沿用 pseudobulk 样本次序。`set_id` 只取与保留基因 ID 的
交集，score 为交集基因 pseudobulk 值的算术平均，并记录 `n_genes_total`
（文件中的成员总数，含矩阵外与被 QC 剔除成员）与 `n_genes_used`
（交集大小）。任一集合与保留基因交集为空即报数据错误（退出码 4），
不产出结果。启用时在既有结果之外新增四个文件（其余文件与不启用时
逐字节一致，`run.json` 除外）：

- `pseudobulk_gene_set_scores.tsv`：沿用 `gene_set_scores.tsv` 的评分布局，
  以 `sample_id` 取代 `cell_id`；列为 `set_id`、`n_genes_total`、
  `n_genes_used`、`sample_id`、`score`，按 `set_id` 升序、集合内按
  pseudobulk 样本顺序排列。
- `pseudobulk_gene_set_score_chart.tsv`：沿用基因集评分图表的宽表布局，
  以 `sample_id`、`group` 开头，各集合列按 `set_id` 升序，样本顺序
  沿用 pseudobulk 样本次序。
- `pseudobulk_gene_set_de.tsv`：沿用 marker 分组差异布局，以 `set_id`
  取代 `gene_id`、`score_difference_a_vs_b` 取代 `log_fc_a_vs_b`；
  列为 `comparison_type`、`group_a`、`group_b`、`set_id`、`mean_in_a`、
  `mean_in_b`、`score_difference_a_vs_b`、`t_stat`、`p_value`、
  `p_value_adj`。每个集合按 group 升序先做 one-vs-rest（`group_b` 为空，
  差值为 mean_in_a 减 mean_in_b，rest 为该组之外的全部样本），再按
  group 升序两两比较；检验为 Welch t 双侧，BH 校正以单个比较内的全部
  集合为一个校正家族。比较先 one-vs-rest 后 pairwise；每个比较内按
  校正 P 值升序、`score_difference_a_vs_b` 降序、`set_id` 升序排列。
- `pseudobulk_gene_set_chart.tsv`：每个比较取前 20 个集合，前三列同上，
  在 `group_b` 之后加入自 1 起、每比较内单独编号的 `rank`，列为
  `comparison_type`、`group_a`、`group_b`、`rank`、`set_id`、
  `score_difference_a_vs_b`、`p_value`、`p_value_adj`、`neg_log10_p_adj`
  （末列）。

`run.json` 的 `parameters` 增加 `"pseudobulk_gene_set_de": true`；
新增顶层 `gene_set_de` 汇总，含 `set_count`、`sample_count`、
`comparison_count`、`test_count`（比较数 × 集合数）与
`min_p_value_adj`（全部检验中的最小校正 P 值）。基因集或重复元数据本身的
错误仍按各自基线报输入错误（退出码 2），集合与保留基因空交集等数据问题
报数据错误（退出码 4）；目录非空或暂存、写出、发布失败报
`OutputPathError`（退出码 5），任何失败都不创建或改动结果目录。
浮点格式沿用既有最短往返表示；相同输入、配置与种子下新增输出逐字节一致。

### `--differential-abundance`（可选，簇级样本差异丰度）

`--differential-abundance` 是无值开关，只在提供 `--replicate-metadata` 时
可用：未提供该基线参数即启用开关，或把开关写成带值形式（如
`--differential-abundance=true`），或取值参数缺值、出现未知参数，一律报
配置错误（退出码 3）。`--replicate-metadata` 的文件格式、gzip 识别、
覆盖与分组约束完全沿用其当前语义；未启用本开关时，全部结果文件与
`run.json` 与基线逐字节一致。

分析只用质控后保留的细胞（启用 `--detect-doublets` 时为双细胞过滤后的
最终细胞）、最终簇标签与重复元数据给出的样本分组，不重新聚类：

- 对每个样本统计其保留细胞在每个最终簇上的细胞数 `n_cells`（样本中缺簇
  计 0）、该样本保留细胞总数 `total_cells`，以及比例
  `proportion = n_cells / total_cells`，覆盖样本 × 簇的全部组合。
- 以样本为观测单位，对每个 group 先做 one-vs-rest（rest 为该组之外的
  全部样本），再按 group 升序两两比较；对每个簇的样本比例执行双侧
  Welch t 检验，差异值 `proportion_difference_a_vs_b` 为 a 组比例均值减
  b 组比例均值，BH 校正以单个比较内的全部簇为一个校正家族。两组方差均为
  0 且均值相同时 `t_stat = 0`、`p_value = 1`（均值不同时为 `inf`/0）。

启用时在既有结果之外新增三个文件（其余文件与不启用时逐字节一致，
`run.json` 除外）：

- `cluster_abundance.tsv`：样本 × 簇全组合，列为 `sample_id`、`group`、
  `cluster`、`n_cells`、`total_cells`、`proportion`，按 `sample_id`、
  `cluster` 升序排列。
- `cluster_differential_abundance.tsv`：列为 `comparison_type`、
  `group_a`、`group_b`（one-vs-rest 时为空）、`cluster`、`mean_in_a`、
  `mean_in_b`、`proportion_difference_a_vs_b`、`t_stat`、`p_value`、
  `p_value_adj`；比较先全部 one-vs-rest（group 升序）再全部 pairwise
  （group 升序），每个比较内按 `p_value_adj` 升序、差异值降序、
  `cluster` 升序排列。
- `cluster_abundance_chart.tsv`：按簇升序，每簇按比例降序、`sample_id`
  升序取前 20 个样本；列为 `cluster`、`rank`（自 1 起、每簇内单独编号）、
  `sample_id`、`group`、`n_cells`、`total_cells`、`proportion`。

`run.json` 的 `parameters` 增加 `"differential_abundance": true`；新增顶层
`differential_abundance` 汇总，含 `sample_count`、`cluster_count`、
`comparison_count` 与 `min_p_value_adj`（全部簇检验中的最小校正 P 值）。
质控后某样本无保留细胞、非空 group 不足两个或任一 group 有效重复不足两个，
报数据错误（退出码 4）；目录非空或暂存、写出、发布失败报
`OutputPathError`（退出码 5），任何失败都不创建或改动结果目录。浮点格式
沿用既有最短往返表示；相同输入、配置与种子下三个新增 TSV 与 `run.json`
逐字节一致。

### `--detect-doublets` / `--expected-doublet-rate`（可选，双细胞识别与过滤）

`--detect-doublets` 是无值开关；`--expected-doublet-rate` 取 `[0, 1)` 内的
有限数值，缺省 `0.08`，仅在启用开关时有效。只给 rate 而不给开关、rate
缺值或非法（非数值、非有限、越界）、未知参数一律报配置错误（退出码 3）。

启用后以 QC 候选细胞 × 候选基因的原始计数子矩阵计算确定性的
`doublet_score`：用 `--seed` 派生的确定性伪随机合成“两个候选细胞计数
叠加”的模拟双细胞，候选细胞与其最相似模拟双细胞的余弦相似度均值即评分，
用以区分双细胞叠加谱与单细胞谱；同输入、同 seed、同参数下逐比特可复现。
目标标记数按候选细胞数与 rate 计算，分数降序、同分按原列序依次标记；
rate 为 0 不标记，且任何情况下至少保留一个细胞。过滤后按 `--min-cells`
在最终细胞上重算保留基因，归一化、HVG、PCA、聚类、差异表达等下游阶段
只用最终细胞与基因。候选或过滤后细胞不足两个、过滤后无保留基因报数据
错误（退出码 4）；目录非空、暂存或发布失败仍报 `OutputPathError`
（退出码 5），分析阶段失败不创建或改动结果目录。

新增两个结果文件（随既有结果事务性发布，未启用时文件集与基线一致）：

- `doublet_scores.tsv`：列 QC 候选细胞原序，含 `cell_id`、`doublet_score`、
  `doublet_rank`、`doublet_flag`、`retained_after_doublet_filter`；
  `doublet_rank` 自 1 起，`doublet_flag` 为是否被过滤，浮点最短往返。
- `doublet_score_chart.tsv`：按分数分 20 个等宽箱，列 `bin_start`、
  `bin_end`、`cell_count`；空箱保留，边界由分数最小/最大值确定。

`run.json` 的 `parameters` 增加 `"detect_doublets": true` 与
`"expected-doublet-rate"`；`stage_counts` 增加 `doublet_candidates`、
`doublets_flagged`、`cells_after_doublet_filter`、
`genes_after_doublet_filter`、`final_cells` 与 `final_genes`。
其余字段不变；未启用时输入、结果与退出码与基线完全一致。

### `--stability-analysis`（可选，聚类稳定性分析）

`--stability-analysis` 是无值开关；仅当显式启用时才执行附加计算，未启用时
主分析的输出内容、文件位置与默认执行顺序与当前基线逐字节一致。配套三个
取值参数，只在开关启用时有效（`--key value` 与 `--key=value` 均可）：

- `--stability-n-samples`：抽样次数，整数且 `>= 2`，缺省 `100`。
- `--stability-sample-fraction`：抽样比例，`(0, 1)` 开区间内的有限数值，
  缺省 `0.8`。比例是每次参与分析的细胞占**通过质量控制细胞数**（启用
  `--detect-doublets` 时为双细胞过滤后的最终细胞数）的比例；每次抽样规模
  为 `floor(比例 × 可用细胞数)`，抽样不放回。
- `--stability-seed`：抽样随机种子，整数，缺省 `20240617`。

只给上述任一取值参数而不给 `--stability-analysis`、取值缺值或非法
（非整数、非数值、非有限、越界），一律报配置错误（退出码 3）。
抽样次数小于二、抽样比例不在零到一之间、随机种子不是整数，在库调用层
统一以 `ValueError` 终止并指出对应配置字段；输入矩阵缺少细胞标识或特征
标识、矩阵为空，沿用既有输入错误（退出码 2，亦为 `ValueError`）；
通过质量控制后少于两个细胞、或抽样规模不足以重新聚类，报数据错误
（退出码 4，亦为 `ValueError`），消息说明可用细胞数。任何失败都在触碰
输出目录之前终止，不生成任何稳定性文件或半成品。

分析以**完整数据的一次降维聚类结果作为参照**，随后逐次抽样：每次对通过
质控的细胞做不放回随机抽样，在抽样子矩阵上沿用与主分析完全相同的
质量控制、`--min-cells` 基因保留、log 归一化、高变基因选择、PCA、聚类
配置（固定 `--n-clusters` 用同一 k，`auto` 用同一自动选择流程；有批次
校正时沿用同一批次均值中心化）生成标签，聚类随机仍只来自主分析
`--seed`，抽样随机只来自 `--stability-seed`，两者相互独立。每次只在
抽样细胞与参照标签都存在的交集上计算调整兰德指数（Adjusted Rand Index）。

启用时在既有结果之外新增三个文件（未启用时文件集与基线一致）：

- `stability_scores.csv`：逐次稳定性分数，表头
  `sample,n_sampled_cells,n_intersection_cells,adjusted_rand_index`，
  一行一次抽样，按抽样次序排列，浮点最短往返。
- `summary.json`：含 `n_samples`、`sample_fraction`、`seed`（以及
  `n_available_cells`、`sample_size`）与 `mean`、`median`、`min`、`max`。
- `stability_score_distribution.json`：供绘图使用的分数分布数据，含
  逐次 `scores` 序列与 20 个等宽箱（边界由分数最小/最大值确定，空箱保留）。

`run.json` 的 `parameters` 增加 `"stability_analysis": true` 与三个稳定性
参数；新增顶层 `stability_analysis` 记录抽样次数、比例、种子、可用细胞数、
每次抽样规模与均值/中位数/最小/最大调整兰德指数；既有字段不变。

输出目录中已有的其他分析产物必须保留：若目标目录已是含 `run.json` 的
既有结果目录且本次启用稳定性分析，则增量覆盖上述三个稳定性文件与
`run.json`（各自先写隐藏临时文件再原子改名，run.json 最后替换），
其余文件原样保留、字节不动；全新目录仍随主分析事务性一次性发布。
本次稳定性分析生成的同名文件在再次运行时被覆盖。目录非普通目录、
父目录不存在或暂存、写出、发布失败仍报 `OutputPathError`（退出码 5）。

相同输入、配置与随机种子重复执行时，逐次抽样选择、稳定性分数、
`summary.json` 与图表数据完全一致；改变 `--stability-seed` 只改变抽样
选择与稳定性统计，绝不改变完整数据主分析的降维坐标、聚类标签或差异
表达结果。该分析可与 `--metadata`、`--batch-metadata`、`--cell-metadata`、
`--gene-sets`、`--replicate-metadata`、`--detect-doublets` 并用；
TSV、MTX、gzip 各输入承载方式下等价矩阵的结果一致性保持不变。

### `--cell-type-reference`（可选，自动细胞类型注释）

`--cell-type-reference <路径>`（`--key value` 与 `--key=value` 均可）指向
一个 UTF-8 制表符文本，表头恰为 `cell_type`、`gene_id` 两列；每行一个
标记关系，两列均非空、`(cell_type, gene_id)` 组合唯一，且至少一条数据行。
文件可为纯文本或单成员 gzip（按 gzip 魔数识别，与文件名无关）。基因 ID
与质控后保留基因**精确匹配**，不做大小写、别名、前缀转换；矩阵外或被
QC 剔除的标记只计入该类型总标记数，不参与评分。文件不存在/不可读、不是
普通文件、不是 UTF-8、表头不符、字段为空、标记关系重复、没有数据行，或
gzip 多成员、尾随数据、截断、CRC/长度错误，一律报输入错误（退出码 2）
且不触碰输出目录；参数缺值、`--cell-type-reference=` 空路径或出现未知
参数报配置错误（退出码 3）。全部标记与保留基因无交集报数据错误
（退出码 4），同样不产出任何结果。

注释只用质控后保留的细胞与基因，取值与实际聚类批次校正口径一致：有
`--batch-metadata`/`--cell-metadata` 的实际批次均值中心化时取校正值，
否则取 log 归一化值。对每个最终簇与每个细胞类型，仅使用同时列于该类型
标记表且命中保留基因的标记，逐基因求“该簇均值减其余簇均值”，再对这些
标记取算术平均作为该（簇, 类型）候选分数。没有任何命中保留基因标记的
类型在该簇不参与竞争；允许多个簇注释为同一类型。每簇取分数最高的类型，
`cell_type` 并列时按 Unicode 码点升序取第一；最高分不大于 0 时
`annotation_status` 为 `unassigned`，但仍保留最高分候选类型与其分数。

启用时在既有结果之外新增两个文件（其余结果与不提供时逐字节一致，
`run.json` 除外）：

- `cluster_annotations.tsv`：按 `cluster` 升序，列为 `cluster`、`n_cells`、
  `cell_type`、`annotation_status`、`score`、`n_markers_total`、
  `n_markers_used`；`score` 为获胜候选分数，`n_markers_total` 为获胜类型
  参考表标记总数，`n_markers_used` 为命中保留基因并实际参与平均的标记数。
- `cluster_annotation_chart.tsv`：细胞顺序沿用 `pca_scatter.tsv`，仅含最终
  细胞，在其 `cell_id`、`PC1`、`PC2`、`cluster` 四列上追加 `cell_type`、
  `annotation_status`、`score`；浮点用最短往返表示。

`run.json` 的 `input.cell_type_reference` 记录参考文件名、原始字节
SHA-256（gzip 载体按压缩字节计算）与各类型总标记数（`marker_counts`，
按 `cell_type` 升序）；`parameters` 增加 `"cell_type_annotation": true`；
新增顶层 `annotation`，按类型汇总其作为获胜类型被注释
（`assigned_clusters`）与最高分不大于 0 而未注释
（`unassigned_clusters`）的簇数。未提供该参数时全部行为与基线逐字节
一致；该参数可与其余全部可选功能并用，不改变它们的结果与口径。浮点格式
沿用既有结果口径（最短往返表示），同输入、配置与种子下注释、排序与图表
数据逐字节一致；目录非空或暂存、写出、发布失败报 `OutputPathError`
（退出码 5）。

### `--enrich-markers` / `--enrichment-alpha` / `--enrichment-min-log-fc`（可选，marker 基因集富集）

`--enrich-markers` 是无值开关，只在 `--gene-sets` 评分基线上可用：未提供
`--gene-sets` 时启用开关，或只给两个取值参数而不启用开关，一律报配置错误
（退出码 3）。配套取值参数（`--key value` 与 `--key=value` 均可）：

- `--enrichment-alpha`：命中与显著性阈值，缺省 `0.05`，必须为 `(0, 1)`
  开区间内的有限数值（0 与 1 均不合法）。
- `--enrichment-min-log-fc`：命中的最小对数倍数变化，缺省 `0`，必须为有限
  数值（允许负数；`inf`/`-inf`/`nan` 不合法）。

开关写成带值形式（如 `--enrich-markers=true`）、取值参数缺值或取值非法
（非数值、非有限、越界）、未知参数，一律报配置错误（退出码 3）。

富集沿用最终细胞、保留基因、最终簇与 `markers.tsv` 的 one-versus-rest
Welch t 检验结果，不重新检验。每簇命中基因（marker hits）为该簇
`p_value_adj <= alpha` 且 `log_fc >= --enrichment-min-log-fc` 的基因；
背景为全部质控后保留基因。对每个 `set_id` 取集合与背景的交集
（`n_set_used`），以超几何分布单侧生存函数求“命中数不少于观测值”的 P 值
（总体为背景、成功为交集基因、抽取为该簇命中），并在每个簇内跨集合做
Benjamini-Hochberg 校正。集合在文件中的成员总数（含矩阵外、被 QC 剔除的
成员）记为 `n_set_total`。

`expected_overlap = n_markers * n_set_used / n_background`；
`fold_enrichment = n_overlap / expected_overlap`，除零取 0；
`odds_ratio` 按“命中是否属于集合”的 2×2 表（命中且在集合、命中不在集合、
非命中且在集合、非命中不在集合）四格统一加 0.5 计算。无命中属于集合
（`n_overlap == 0`）的集合输出 `p_value = 1`、`fold_enrichment = 0`；
某簇无任何命中基因时该簇全部集合同样按此口径输出。

启用时在既有结果（含 marker、基因集评分、批次校正、注释等全部既有产物）
之外新增两个文件，其余结果文件与不启用时逐字节一致（`run.json` 除外）：

- `marker_gene_set_enrichment.tsv`：列为 `cluster`、`set_id`、
  `n_set_total`、`n_set_used`、`n_markers`（命中数）、`n_overlap`、
  `expected_overlap`、`fold_enrichment`、`odds_ratio`、`p_value`、
  `p_value_adj`；全部簇 × 全部集合，按 `cluster`、`p_value_adj`、`set_id`
  升序，`odds_ratio` 降序排列。
- `marker_gene_set_enrichment_chart.tsv`：每簇取上述排序前 20 行，在
  `cluster` 后插入自 1 起、每簇内单独编号的 `rank` 列，其余列同上。

`run.json` 的 `parameters` 增加 `"enrich_markers": true`、
`"enrichment_alpha"` 与 `"enrichment_min_log_fc"`；新增顶层 `enrichment`
记录全部簇命中基因数之和（`marker_hits`）与 `p_value_adj <= alpha` 的
显著（簇, 集合）组合数（`significant_sets`）；既有字段不变。基因集文件
本身的错误仍报输入错误（退出码 2），基因集与保留基因交集为空等数据问题
报数据错误（退出码 4）；目录非空或暂存、写出、发布失败报
`OutputPathError`（退出码 5）。未启用 `--enrich-markers` 时全部既有行为、
TSV/MTX 读取与 gzip 处理与基线逐字节一致。

## 约定
- 公开行为以 README 与源码为准。
- 后续需求在此基线上增量实现。
