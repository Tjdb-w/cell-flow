"""Cell Flow 的异常体系与对应退出码。"""


class CellFlowError(Exception):
    """所有 Cell Flow 受控错误的基类，携带确定退出码。"""

    exit_code = 1


class CellFlowInputError(CellFlowError):
    """输入矩阵不合法（空路径、不可读、表头缺失、ID 重复、空矩阵、非法计数）。"""

    exit_code = 2


class CellFlowBatchError(CellFlowError, ValueError):
    """批次元数据与表达矩阵细胞条码冲突，或批次/样本标签不合法。

    同时是 :class:`ValueError`：调用方可直接 ``except ValueError`` 捕获；
    经命令行运行时按受控输入错误以退出码 2 报告。所有此类失败都发生在
    质控之后、降维聚类之前，且在触碰输出目录之前，不产生任何部分结果。
    """

    exit_code = 2


class CellFlowConfigError(CellFlowError):
    """命令行参数冲突或取值非法。"""

    exit_code = 3


class CellFlowDataError(CellFlowError):
    """数据无法支撑分析（质控后细胞不足、无可用基因、PCA/聚类无法成立）。"""

    exit_code = 4


class OutputPathError(CellFlowError):
    """结果目录已存在且非空。"""

    exit_code = 5
