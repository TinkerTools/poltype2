"""
poltype.pipeline.stages.validation – parameter validation stage.

Runs one or more :class:`Validator` instances against the current
pipeline context to check that generated parameters meet quality
thresholds.  The :class:`ValidationReport` is stored as the
``"validation_report"`` artifact for consumption by the finalization
stage.

This stage is skipped when ``config.skip_validation`` is ``True`` or
the pipeline is running in dry-run mode.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, List, Optional

from poltype.pipeline.stage import Stage, StageResult, StageStatus
from poltype.validation.validation_result import ValidationReport, ValidationSeverity
from poltype.validation.validator import EnergyValidator, Validator

if TYPE_CHECKING:
    from poltype.pipeline.context import PipelineContext

logger = logging.getLogger(__name__)


class ValidationStage(Stage):
    """Validate generated parameters before final output.

    Parameters
    ----------
    validators:
        A list of :class:`Validator` instances.  If ``None``, a
        default :class:`EnergyValidator` is used.
    """

    def __init__(self, validators: Optional[List[Validator]] = None) -> None:
        super().__init__(name="validation")
        self._validators: List[Validator] = (
            list(validators) if validators is not None else []
        )

    @property
    def validators(self) -> List[Validator]:
        """Configured validators."""
        return list(self._validators)

    def execute(self, context: "PipelineContext") -> StageResult:
        """Run all validators and aggregate results.

        Returns
        -------
        StageResult
            ``COMPLETED`` with a ``validation_report`` artifact on
            success, ``FAILED`` if validation encounters errors that
            cannot be recovered from.
        """
        if context.molecule is None:
            return StageResult(
                status=StageStatus.FAILED,
                message="No molecule in context",
            )

        # Build validators lazily from config if not injected.
        validators = self._validators
        if not validators:
            validators = [EnergyValidator()]

        all_records = []
        for validator in validators:
            logger.info("Running validator: %s", validator.name)
            report = validator.validate(context)
            all_records.extend(report.records)

        combined = ValidationReport(
            records=all_records,
            molecule_name=context.molecule.name,
            validation_method=", ".join(v.name for v in validators),
            metadata={"num_validators": len(validators)},
        )
        context.set_artifact("validation_report", combined)

        # Log results
        for record in all_records:
            if not record.passed:
                if record.severity is ValidationSeverity.ERROR:
                    logger.error(
                        "Validation FAILED [%s]: %s",
                        record.check_name,
                        record.message,
                    )
                else:
                    logger.warning(
                        "Validation WARNING [%s]: %s",
                        record.check_name,
                        record.message,
                    )
            else:
                logger.info(
                    "Validation PASSED [%s]: %s",
                    record.check_name,
                    record.message,
                )

        msg = (
            f"{combined.num_passed} passed, {combined.num_failed} failed "
            f"({combined.num_records} total)"
        )
        logger.info("Validation complete: %s", msg)

        return StageResult(
            status=StageStatus.COMPLETED,
            message=msg,
            artifacts={"validation_report": combined},
        )

    def should_skip(self, context: "PipelineContext") -> bool:
        """Skip when validation is disabled or dry-run mode is active."""
        return context.config.skip_validation or context.config.dry_run
