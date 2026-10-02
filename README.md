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
文件可以是 UTF-8 文本，也可以是同一文本的 gzip 压缩（按 gzip 魔数识别，
文件名不限）；gzip 只接受单成员且无尾随数据，头部、CRC、长度校验失败
或数据截断均报输入错误（退出码 2）。

### `mtx`（标准 10x MatrixMarket 稀疏目录）

`--input` 指向一个目录，其中须含三个文件，三者必须同为未压缩或同为
gzip（`matrix.mtx.gz`、`barcodes.tsv.gz`、`features.tsv.gz`），混用报
输入错误（退出码 2）：

- `matrix.mtx`：仅接受 MatrixMarket 头
  `%%MatrixMarket matrix coordinate integer general` 或
  `coordinate real general`；行是基因、列是细胞、索引从 1 开始；
  零值可省略，显式值必须能无损解析为有限非负整数（`real` 字段同样
  要求恰为整数值，超大整数按精确十进制解析，不做 float64 舍入）。
- `barcodes.tsv`：每个非空行是一个细胞 ID，按文件行序进入分析。
- `features.tsv`：每个非空数据行取第一列作为基因 ID，按文件行序进入分析。

`matrix.mtx` 声明的行/列维度必须分别等于基因/细胞 ID 数量。
缺文件、压缩与未压缩混用、gzip 无法无损还原、格式头不支持、索引越界、
重复坐标、空或重复 ID、声明元素数不符、非法数值一律报输入错误
（退出码 2）；`--input-format` 取其他值报配置错误（退出码 3）。
任何输入失败都不会创建或修改结果目录。

压缩只改变承载方式：矩阵等价的两种 TSV（文本 / gzip）或两种 MTX
（原名 / `.gz` 名）来源，除 `run.json` 的 `input` 来源字段外，每个结果
文件完全一致。`input.sha256` 与 `input.files` 中各 SHA-256 均按实际文件
原始字节（压缩态）计算，`files` 依次记录 matrix、barcodes、features 的
实际文件名。

两种格式给出的等价矩阵产生完全一致的分析结果（排序、数值、图表数据）。
`tsv` 运行的 `run.json` 字段与基线一致；`mtx` 运行仅在其 `input` 结构中
额外记录 `"input_format": "mtx"` 以及三个输入文件名和各自的 SHA-256。

## 约定

- 公开行为以 README 与源码为准。
- 后续需求在此基线上增量实现。
