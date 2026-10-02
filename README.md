# Cell Flow

单细胞组学分析管线：从表达矩阵完成质量控制、归一化、降维聚类与差异表达分析，输出可复现的分析结果与图表数据。

## 范围

本仓库从零开始实现上述方向的可用工具，不依赖外部同类实现。

## 状态

初始基线：只有本说明，尚无实现。

## 输入

入口仍为 `cell-flow analyze --input <路径> --output-dir <结果目录>`，
新增可选 `--input-format`（取值 `tsv` 或 `mtx`，默认 `tsv`）。
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

## 约定

- 公开行为以 README 与源码为准。
- 后续需求在此基线上增量实现。
