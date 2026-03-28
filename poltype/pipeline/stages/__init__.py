"""
poltype.pipeline.stages – concrete pipeline stage implementations.

Each stage maps to a major phase of the monolithic
``GenerateParameters()`` method in the legacy code.  Stages delegate
to the :class:`~poltype.qm.backend.QMBackend` interface for QM
computations.

Re-exports all concrete stages for convenience::

    from poltype.pipeline.stages import InputPreparationStage
"""

from __future__ import annotations

from poltype.pipeline.stages.esp import ESPFittingStage
from poltype.pipeline.stages.finalization import FinalizationStage
from poltype.pipeline.stages.fragmentation import FragmentationStage
from poltype.pipeline.stages.geometry_opt import GeometryOptimizationStage
from poltype.pipeline.stages.input_prep import InputPreparationStage
from poltype.pipeline.stages.multipole import MultipoleStage
from poltype.pipeline.stages.torsion import TorsionFittingStage
from poltype.pipeline.stages.validation import ValidationStage

__all__ = [
    "InputPreparationStage",
    "GeometryOptimizationStage",
    "ESPFittingStage",
    "MultipoleStage",
    "FragmentationStage",
    "TorsionFittingStage",
    "ValidationStage",
    "FinalizationStage",
]
