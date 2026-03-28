"""
poltype.logging_config – structured logging for the pipeline.

Configures Python's stdlib :mod:`logging` with a consistent format and
sensible defaults.  Call :func:`setup_logging` once at application
start-up (typically from ``__main__``) before running the pipeline.

Usage::

    from poltype.logging_config import setup_logging

    setup_logging(verbose=True)

This module has **no dependency on PoltypeModules** and can be
unit-tested in isolation.
"""

from __future__ import annotations

import logging
import sys

#: Default format for log messages.
LOG_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"

#: Date format used in log messages.
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

#: The root logger name for the poltype package.
LOGGER_NAME = "poltype"


def setup_logging(
    *,
    verbose: bool = False,
    quiet: bool = False,
    log_file: str | None = None,
) -> logging.Logger:
    """Configure the ``poltype`` logger hierarchy.

    Parameters
    ----------
    verbose:
        If ``True``, set the log level to ``DEBUG``.
    quiet:
        If ``True``, set the log level to ``WARNING``.  ``quiet``
        takes precedence over ``verbose`` if both are set.
    log_file:
        Optional path to a file where log output is also written.

    Returns
    -------
    logging.Logger
        The configured root ``poltype`` logger.
    """
    if quiet:
        level = logging.WARNING
    elif verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(level)

    # Avoid duplicate handlers when called more than once (e.g. in tests).
    if not logger.handlers:
        formatter = logging.Formatter(LOG_FORMAT, datefmt=DATE_FORMAT)

        console = logging.StreamHandler(sys.stderr)
        console.setLevel(level)
        console.setFormatter(formatter)
        logger.addHandler(console)

        if log_file is not None:
            fh = logging.FileHandler(log_file)
            fh.setLevel(level)
            fh.setFormatter(formatter)
            logger.addHandler(fh)
    else:
        # Update existing handler levels.
        for handler in logger.handlers:
            handler.setLevel(level)
        logger.setLevel(level)

    return logger
