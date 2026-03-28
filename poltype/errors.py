"""
poltype.errors – custom exception hierarchy.

All Poltype-specific exceptions inherit from :class:`PoltypeError` so
callers can catch a single base class when they want a broad handler.

This module has **no dependency on PoltypeModules** and can be
unit-tested in isolation.
"""

from __future__ import annotations


class PoltypeError(Exception):
    """Base exception for all Poltype errors."""


class ConfigError(PoltypeError):
    """Raised when configuration loading or validation fails.

    Examples include a missing ``poltype.ini`` file, an unrecognised
    keyword, or a value that cannot be coerced to the expected type.
    """


class StageError(PoltypeError):
    """Raised when a pipeline stage encounters an unrecoverable error.

    Attributes
    ----------
    stage_name:
        Name of the stage that failed.
    """

    def __init__(self, message: str, stage_name: str = "") -> None:
        self.stage_name = stage_name
        super().__init__(message)


class BackendError(PoltypeError):
    """Raised when a QM backend operation fails.

    Examples include a missing executable, a failed QM calculation, or
    an unsupported method/basis-set combination.

    Attributes
    ----------
    backend_name:
        Name of the backend that raised the error.
    """

    def __init__(self, message: str, backend_name: str = "") -> None:
        self.backend_name = backend_name
        super().__init__(message)
