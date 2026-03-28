"""
poltype.pipeline.stages.input_prep – input preparation stage.

Validates that the configuration specifies a structure file, creates a
:class:`Molecule` from it, and stores the molecule in the pipeline
context.
"""

from __future__ import annotations

from poltype.molecule.molecule import Molecule
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus


class InputPreparationStage(Stage):
    """Prepare the input molecule for the pipeline.

    If ``context.molecule`` is already set (e.g. the caller injected one
    directly), this stage does nothing and returns ``COMPLETED``.
    Otherwise it reads the structure file from ``context.config.structure``.
    """

    def __init__(self) -> None:
        super().__init__(name="input_preparation")

    def execute(self, context: PipelineContext) -> StageResult:
        """Load the molecule from the configured structure file.

        Returns
        -------
        StageResult
            ``COMPLETED`` with a ``molecule`` artifact.
        """
        if context.molecule is not None:
            return StageResult(
                status=StageStatus.COMPLETED,
                message="Molecule already present in context",
                artifacts={"molecule": context.molecule},
            )

        structure = context.config.structure
        if structure is None:
            return StageResult(
                status=StageStatus.FAILED,
                message="No structure path in config",
            )

        mol = Molecule.from_file(structure, work_dir=context.work_dir)
        context.update_molecule(mol)
        return StageResult(
            status=StageStatus.COMPLETED,
            message=f"Loaded molecule from {structure}",
            artifacts={"molecule": mol},
        )

    def should_skip(self, context: PipelineContext) -> bool:
        """Skip when only checking input validity."""
        return context.config.check_input_only
