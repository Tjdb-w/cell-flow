"""Typed errors with stable exit codes."""


class CellFlowError(Exception):
    """Base class for all cell-flow errors."""

    exit_code = 1


class CellFlowInputError(CellFlowError):
    """Invalid or unreadable input matrix."""

    exit_code = 2


class CellFlowConfigError(CellFlowError):
    """Invalid or conflicting parameters."""

    exit_code = 3


class CellFlowDataError(CellFlowError):
    """Data cannot support the requested analysis."""

    exit_code = 4


class OutputPathError(CellFlowError):
    """Output directory cannot be used."""

    exit_code = 5
