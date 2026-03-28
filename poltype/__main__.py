"""
poltype.__main__ – CLI entry point for ``python -m poltype``.

Parses command-line arguments, loads configuration, builds the default
pipeline, and runs it.  This wires together every layer introduced in
Phases 1–3.

Usage::

    python -m poltype --config poltype.ini --work-dir ./output
    python -m poltype --verbose
    python -m poltype --version
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from poltype.config.loader import load_config
from poltype.errors import ConfigError, PoltypeError
from poltype.logging_config import setup_logging
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.factory import build_default_pipeline
from poltype.pipeline.stage import StageStatus
from poltype.qm.factory import select_backend

logger = logging.getLogger("poltype.cli")


def _build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser."""
    parser = argparse.ArgumentParser(
        prog="poltype2",
        description=(
            "Poltype 2 — automated AMOEBA force-field parameterisation "
            "for small molecules."
        ),
    )
    parser.add_argument(
        "-c",
        "--config",
        default="poltype.ini",
        help="Path to the poltype.ini configuration file (default: poltype.ini)",
    )
    parser.add_argument(
        "-w",
        "--work-dir",
        default=None,
        help="Working directory for file I/O (default: current directory)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose (DEBUG-level) logging",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress all output except warnings and errors",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="poltype2 0.0.0-dev",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the Poltype 2 pipeline.

    Parameters
    ----------
    argv:
        Command-line arguments.  Defaults to ``sys.argv[1:]``.

    Returns
    -------
    int
        Exit code: ``0`` for success, ``1`` for pipeline failure,
        ``2`` for configuration or usage errors.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)

    # -- set up logging --
    setup_logging(verbose=args.verbose, quiet=args.quiet)

    # -- load configuration --
    try:
        config = load_config(args.config)
    except FileNotFoundError:
        logger.error("Configuration file not found: %s", args.config)
        return 2
    except ConfigError as exc:
        logger.error("Configuration error: %s", exc)
        return 2

    logger.info("Configuration loaded from %s", args.config)

    # -- resolve work directory --
    work_dir = Path(args.work_dir).resolve() if args.work_dir else Path.cwd()
    logger.info("Working directory: %s", work_dir)

    # -- resolve relative structure path against work directory --
    if config.structure is not None and not config.structure.is_absolute():
        config.structure = (work_dir / config.structure).resolve()

    # -- select QM backend --
    backend = select_backend(config, work_dir=work_dir)
    logger.info("QM backend: %s", backend.name)

    # -- build and run pipeline --
    runner = build_default_pipeline()
    context = PipelineContext(
        config=config,
        backend=backend,
        work_dir=work_dir,
    )

    logger.info("Starting pipeline (%s)", runner)

    try:
        context = runner.run(context)
    except PoltypeError as exc:
        logger.error("Pipeline error: %s", exc)
        return 1
    except Exception:
        logger.exception("Unexpected error during pipeline execution")
        return 1

    # -- report outcome --
    failed = [
        name
        for name, result in context.stage_results.items()
        if result.status is StageStatus.FAILED
    ]
    if failed:
        logger.error("Pipeline finished with failures in: %s", ", ".join(failed))
        return 1

    logger.info("Pipeline completed successfully")
    return 0


if __name__ == "__main__":
    sys.exit(main())
