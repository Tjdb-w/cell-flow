# Cell Flow

单细胞组学分析管线：从表达矩阵完成质量控制、归一化、降维聚类与差异表达分析，输出可复现的分析结果与图表数据。

## 范围

本仓库从零开始实现上述方向的可用工具，不依赖外部同类实现。

## 状态

初始基线：只有本说明，尚无实现。

## 输入

入口仍为 `cell-flow analyze --input <路径> --output-dir <结果目录>`，
新增可选 `--input-format`（取值 `tsv` 或 `mtx`，默认 `tsv`）与
`--gene-sets`。`--key value` 与 `--key=value` 两种写法均支持。

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

### `--batch-metadata`（可选，跨样本批次校正）

`--batch-metadata <路径>`（`--key value` 与 `--key=value` 均可）指向
一个 UTF-8 制表符文本，必须同时包含三个带名列：固定的细胞条码列
`cell_id`、样本标识列与批次标签列。样本列默认名为 `sample_id`，批次列
默认名为 `batch`，可分别通过 `--sample-column <列名>` 与
`--batch-column <列名>` 指定其他列名；两列名不得相同，也不得使用
`cell_id`。表头中的其他列原样保留并透传到逐细胞公开结果，不覆盖也不
重命名既有字段。文件可为纯文本或单成员 gzip（按内容识别，与文件名无关）。

`cell_id` 唯一且与表达矩阵的全部细胞条码一一对应（不多不少），样本标识
与批次标签均非空——批次标签缺失或为空一律报错，不会静默归入未知批次。
下列冲突统一抛 `ValueError`（经命令行运行时为退出码 2），消息指出具体
冲突类型，且都在质控之后、降维聚类启动之前失败，不生成任何部分结果：

- 表达矩阵包含重复细胞条码；
- 批次元数据包含重复细胞条码；
- 元数据缺少样本标识列或批次标签列（或指定列不存在）；
- 样本标识或批次标签存在空值；
- 元数据与表达矩阵的细胞集合不一致（缺少或多出细胞）。

**单批次数据按原流程继续分析**：不做校正、不引入任何数值扰动，质控、
归一化、高变基因、PCA、聚类、差异表达等既有结果与不提供该参数时逐字节
一致（批次功能新增的文件除外），已有的单样本使用方式不变。质控后仍保留
两个及以上不同批次时，才在 log 归一化之后按基因做批次均值中心化：每个
归一化值减去对应批次保留细胞的基因均值，再加回全部保留细胞的基因总均值
（batch mean centering）。校正值进入高变基因选择、PCA、聚类、markers、
成对 markers、`--metadata` 分组差异表达及图表数据；
`normalized_expression.tsv` 始终写未校正值。

提供批次元数据时（单批次与多批次均如此）在既有结果之外新增三个文件：

- `cell_metadata.tsv`：每个质控后保留细胞一行，列为 `cell_id`、
  `sample_id`、`batch`（语义列名固定，实际来源列名见 `run.json`），
  后随元数据中的其他列（原列名、原值、原顺序）；细胞顺序与
  `clusters.tsv` 一致。
- `batch_pca_scatter.tsv`：按批次与样本着色的降维坐标，列为
  `cell_id`、`PC1`、`PC2`、`cluster`、`batch`、`sample_id`；
  坐标即后续聚类实际使用的校正后低维表示（单批次时与 `pca.tsv` 相同）。
- `batch_mixing.tsv`：校正前、校正后各一行，列为 `stage`、
  `n_neighbors`、`mixing_score`。混合分数按每个细胞的 15 个最近邻
  （细胞不足时为全部其他细胞）中来自其他批次的细胞占比计算，再对全部
  细胞求平均；最近邻在降维坐标上按欧氏距离确定，距离并列以细胞顺序
  打破，不使用随机数。校正前分数在未校正归一化的同一管线 PCA 上计算
  （单批次未发生校正，两行分数相同）。

多批次实际执行校正时另产出 `batch_corrected_expression.tsv`：校正后的
表达矩阵，行（保留基因）列（保留细胞）顺序与 `normalized_expression.tsv`
相同；单批次不校正故不产出该文件。批次功能下还产出
`batch_summary.tsv`：按 `batch_id` 升序，每行给出 `batch_id`、
`n_cells`（批次内保留细胞数）、`total_counts`、`detected_genes`、
`mitochondrial_fraction`，后三项为批次内保留细胞的均值。

`run.json` 的 `input.batch_metadata` 记录批次元数据文件名、原始字节
SHA-256、实际样本/批次列名、透传列名以及质控前后各批次与各样本的细胞数；
`parameters` 在批次功能下增加 `batch_correction`（含方法、是否实际应用、
批次数量与混合分数近邻数），实际执行中心化时另含
`"batch_mean_centering": true`；顶层 `batch_mixing` 记录校正前后混合
分数；`stage_counts` 增加批次与样本计数。不提供 `--batch-metadata` 时
`parameters` 与 `stage_counts` 与基线逐字节一致。相同输入、参数与随机
种子（含不显式提供种子时的固定默认种子）给出逐字节相同的结果。
`--batch-metadata` 可与 `--metadata`、`--gene-sets` 并用；TSV、MTX、
gzip 各输入承载方式下等价矩阵的结果一致性保持不变。目标非空或暂存、
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

## 约定

- 公开行为以 README 与源码为准。
- 后续需求在此基线上增量实现。
