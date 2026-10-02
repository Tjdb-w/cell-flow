# Cell Flow

单细胞组学分析管线：从表达矩阵完成质量控制、归一化、降维聚类与差异表达分析，输出可复现的分析结果与图表数据。

## 入口

```
cell-flow analyze --input <输入路径> --output-dir <结果目录> \
    [--input-format tsv|mtx] [其余参数同基线]
```

`--input-format` 默认 `tsv`，`--key value` 与 `--key=value` 两种写法均支持。

### 输入格式

- `tsv`（默认）：单个 UTF-8 制表符计数矩阵，首列是唯一基因 ID，表头其余列名是唯一细胞 ID，取值为非负整数 UMI 计数。
- `mtx`：标准 10x MatrixMarket **目录**，须含 `matrix.mtx`、`barcodes.tsv`、`features.tsv` 三个文件：
  - `matrix.mtx` 仅接受 MatrixMarket 头 `coordinate integer general` 或 `coordinate real general`；行＝基因、列＝细胞、索引从 1 开始；零值可省略；显式值须无损解析为有限非负整数（real 头下 `1.0`、`5e2`、`1.2e3` 等整数值写法可接受，小数、负数、nan/inf 拒绝）；
  - `barcodes.tsv`：每个非空行是一个细胞 ID，按文件行序进入分析；
  - `features.tsv`：每个非空数据行取第一列作为基因 ID，按文件行序进入分析；
  - 声明维度须与 ID 数量一致；坐标不得重复；ID 不得为空或重复。

两种格式读入后展开为同一基因行序、细胞列序、逐格整数计数的矩阵，下游质控、归一化、PCA、聚类、差异表达与全部结果文件的字段、行序、数值完全一致。

### 错误语义

- 输入文件/目录缺失、格式头不支持、越界索引、重复坐标、空或重复 ID、非法数值：`CellFlowInputError`（退出码 2）；
- `--input-format` 取值非法：`CellFlowConfigError`（退出码 3）；
- 输入失败发生在任何输出动作之前，不会创建或修改结果目录。

### 来源记录（run.json）

- `tsv` 运行的 `run.json` 与基线逐字段一致，不增删或改字段；
- `mtx` 运行仅在 `input` 结构内追加 `input_format: "mtx"` 与 `files`：`matrix.mtx`、`barcodes.tsv`、`features.tsv` 三个文件名及各自原始字节的 SHA-256。

### 结果目录

只接受不存在或空目录；全部计算成功后原子写出，且不覆盖任何已有文件。

## 范围

本仓库从零开始实现上述方向的可用工具，不依赖外部同类实现。

## 约定

- 公开行为以 README 与源码为准。
- 后续需求在此基线上增量实现。
