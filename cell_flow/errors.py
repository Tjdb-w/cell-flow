"""Cell Flow 的异常体系与对应退出码。"""


class CellFlowError(Exception):
    """所有 Cell Flow 受控错误的基类，携带确定退出码。"""

    exit_code = 1


class CellFlowInputError(CellFlowError, ValueError):
    """输入矩阵不合法（空路径、不可读、表头缺失、ID 重复、空矩阵、非法计数）。

    同时是 :class:`ValueError`：输入数据冲突（重复条码、缺少样本标识或
    批次标签、细胞集合不一致、空批次等）对调用方统一表现为 ValueError，
    命令行退出码仍保持 2。
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
