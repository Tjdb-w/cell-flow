# Cell Flow

单细胞组学分析管线：从表达矩阵完成质量控制、归一化、降维聚类与差异表达分析，输出可复现的分析结果与图表数据。

## 范围

本仓库从零开始实现上述方向的可用工具，不依赖外部同类实现。

## 状态

初始基线：只有本说明，尚无实现。

## 输入

入口仍为 `cell-flow analyze --input <路径> --output-dir <结果目录>`，
新增可选 `--input-format`（取值 `tsv` 或 `mtx`，默认 `tsv`）与
`--metadata <分组文件>`（用于细胞分组差异表达，缺省时不做分组分析）。
`--key value` 与 `--key=value` 两种写法均支持。

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

`--metadata` 指向 UTF-8 TSV，表头恰为 `cell_id` 与 `group` 两列；
每个数据行给出一个细胞 ID 与其非空分组名，`cell_id` 唯一且恰好覆盖
矩阵的全部细胞（不缺、不多）。承载方式与表达矩阵一致：未压缩文本或
单成员 gzip，按 gzip 魔数自动识别，与文件名或后缀无关；gzip 多成员、
尾随数据、CRC、长度校验失败或截断均按内容判定报输入错误。

内容不合法（缺文件、不可读、非 UTF-8、表头不符、行列数不符、空或
重复 `cell_id`、空 `group`、未覆盖全部细胞）或 gzip 不是单成员流，
一律报输入错误（退出码 2）；`--metadata` 缺值或出现未知参数报配置
错误（退出码 3）。任何输入失败都不会创建或修改结果目录。

提供该参数后，仅用质控后保留细胞与保留基因的 log 归一化表达做分组
差异表达：先对每个分组按 `group` 升序做 one-vs-rest，再按 `group`
升序对每对分组（a < b）做两两检验。每个比较内跨全部保留基因做
Benjamini-Hochberg 校正，记录按校正 P 值升序、`log_fc` 降序、
`gene_id` 升序排序。新增两个结果文件：

- `group_markers.tsv`：前三列为 `comparison_type`（`one_vs_rest` 或
  `pairwise`）、`group_a`、`group_b`（one-vs-rest 时 `group_b` 为空列），
  其余列为 `gene_id`、`mean_in_a`、`mean_in_b`、`log_fc`、`t_stat`、
  `p_value`、`p_value_adj`，口径与 pairwise marker 一致
  （均值差即 `log_fc`，统计量与 P 值用 Welch t 检验）。
- `group_marker_chart.tsv`：取每个比较排序后的前 20 个基因，前三列同上，
  另含比较内 `rank`（从 1 起）以及 `log_fc`、`p_value`、`p_value_adj`、
  `neg_log10_p_adj`。

若质控后非空分组少于两个，报数据错误（退出码 4）。`run.json` 的
`input` 结构在提供该参数时额外记录 `metadata`：文件名 `name`、按原始
字节计算的 SHA-256 `sha256`，以及质控后各分组的细胞数
`group_cell_counts`（按分组名升序）。不提供 `--metadata` 时 `input`
结构与基线逐字一致，也不产出上述两个文件。

## 约定

- 公开行为以 README 与源码为准。
- 后续需求在此基线上增量实现。
