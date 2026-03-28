"""
poltype.__main__ – CLI entry point for ``python -m poltype``.

Parses command-line arguments, loads configuration, builds the default
pipeline, and runs it.  This wires together every layer introduced in
Phases 1–5.

Usage::

    python -m poltype --config poltype.ini --work-dir ./output
    python -m poltype --verbose
    python -m poltype --resume
    python -m poltype --version
"""

from __future__ import annotations

import argparse
import logging
import sys
from importlib.metadata import version as pkg_version
from pathlib import Path

from poltype.config.loader import load_config
from poltype.errors import ConfigError, PoltypeError
from poltype.logging_config import setup_logging
from poltype.pipeline.checkpoint import CheckpointManager
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.events import EventBus, EventData, PipelineEvent
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
        "--resume",
        action="store_true",
        help="Resume from the last checkpoint (skip already-completed stages)",
    )
    parser.add_argument(
        "--no-checkpoint",
        action="store_true",
        help="Disable checkpoint saving during this run",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate inputs and report what would be done without running QM",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"poltype2 {pkg_version('poltype2')}",
    )
    return parser


def _logging_hook(data: EventData) -> None:
    """Default event hook that logs every pipeline event."""
    if data.event is PipelineEvent.STAGE_STARTED:
        logger.info("[event] Stage '%s' started", data.stage_name)
    elif data.event is PipelineEvent.STAGE_COMPLETED:
        elapsed = data.extra.get("elapsed_seconds", 0.0)
        logger.info(
            "[event] Stage '%s' completed (%.2fs)",
            data.stage_name,
            elapsed,
        )
    elif data.event is PipelineEvent.STAGE_SKIPPED:
        logger.info("[event] Stage '%s' skipped", data.stage_name)
    elif data.event is PipelineEvent.STAGE_FAILED:
        logger.error(
            "[event] Stage '%s' failed: %s", data.stage_name, data.message,
        )
    elif data.event in (
        PipelineEvent.PIPELINE_STARTED,
        PipelineEvent.PIPELINE_COMPLETED,
    ):
        logger.info("[event] %s", data.message)


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

    # -- apply dry-run flag --
    if args.dry_run:
        config.dry_run = True
        logger.info("Dry-run mode: QM stages will be skipped")

    # -- select QM backend --
    backend = select_backend(config, work_dir=work_dir)
    logger.info("QM backend: %s", backend.name)

    # -- set up event bus --
    event_bus = EventBus()
    event_bus.register(_logging_hook)

    # -- set up checkpoint manager --
    checkpoint: CheckpointManager | None = None
    if not args.no_checkpoint:
        checkpoint = CheckpointManager(checkpoint_dir=work_dir)
        if args.resume:
            try:
                checkpoint.load()
                logger.info("Resuming from checkpoint")
            except FileNotFoundError:
                logger.warning(
                    "No checkpoint file found — starting from scratch"
                )

    # -- build and run pipeline --
    runner = build_default_pipeline()
    runner.event_bus = event_bus
    runner.checkpoint_manager = checkpoint

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
