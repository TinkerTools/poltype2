"""
poltype.pipeline.factory – default pipeline construction.

The :func:`build_default_pipeline` function assembles a
:class:`PipelineRunner` with the standard ten-stage sequence that
mirrors the monolithic ``GenerateParameters()`` method.

Usage::

    from poltype.pipeline.factory import build_default_pipeline

    runner = build_default_pipeline()
    ctx = runner.run(context)

This module has **no dependency on PoltypeModules** and can be
unit-tested in isolation.
"""

from __future__ import annotations

from poltype.pipeline.runner import PipelineRunner
from poltype.pipeline.stages.atom_typing import AtomTypingStage
from poltype.pipeline.stages.database_match import DatabaseMatchStage
from poltype.pipeline.stages.esp import ESPFittingStage
from poltype.pipeline.stages.finalization import FinalizationStage
from poltype.pipeline.stages.fragmentation import FragmentationStage
from poltype.pipeline.stages.geometry_opt import GeometryOptimizationStage
from poltype.pipeline.stages.input_prep import InputPreparationStage
from poltype.pipeline.stages.multipole import MultipoleStage
from poltype.pipeline.stages.torsion import TorsionFittingStage
from poltype.pipeline.stages.validation import ValidationStage


def build_default_pipeline() -> PipelineRunner:
    """Create a :class:`PipelineRunner` with the standard stage sequence.

    The returned runner contains, in order:

    1. :class:`InputPreparationStage`
    2. :class:`GeometryOptimizationStage`
    3. :class:`ESPFittingStage`
    4. :class:`MultipoleStage`
    5. :class:`AtomTypingStage`
    6. :class:`DatabaseMatchStage`
    7. :class:`FragmentationStage`
    8. :class:`TorsionFittingStage`
    9. :class:`ValidationStage`
    10. :class:`FinalizationStage` (with default output writers)

    Returns
    -------
    PipelineRunner
        A fully assembled runner ready to execute.
    """
    # Lazy imports to avoid circular dependency (output → pipeline → factory → output)
    from poltype.output.key_writer import KeyFileWriter
    from poltype.output.summary import SummaryWriter
    from poltype.output.xyz_writer import XYZFileWriter

    writers = [
        XYZFileWriter(),
        KeyFileWriter(),
        SummaryWriter(),
    ]

    return (
        PipelineRunner()
        .add_stage(InputPreparationStage())
        .add_stage(GeometryOptimizationStage())
        .add_stage(ESPFittingStage())
        .add_stage(MultipoleStage())
        .add_stage(AtomTypingStage())
        .add_stage(DatabaseMatchStage())
        .add_stage(FragmentationStage())
        .add_stage(TorsionFittingStage())
        .add_stage(ValidationStage())
        .add_stage(FinalizationStage(writers=writers))
    )
