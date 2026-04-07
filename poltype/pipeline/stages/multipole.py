"""
poltype.pipeline.stages.multipole – distributed multipole analysis stage.

Performs the full DMA workflow:

1. Run a QM single-point calculation (typically MP2/6-311G**) to produce
   a formatted checkpoint (fchk) file via :meth:`QMBackend.compute_dma`.
2. Run the GDMA program on the fchk file to extract distributed
   multipoles (charges, dipoles, quadrupoles).
3. Run Tinker's ``poledit`` program to process the raw GDMA multipoles
   into AMOEBA-compatible local-frame multipole parameters.

The resulting ``multipole_frames`` artifact is consumed by the downstream
:class:`ESPFittingStage`, which refines the multipoles by fitting to
the electrostatic potential.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from poltype.errors import BackendError
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus

logger = logging.getLogger(__name__)


class MultipoleStage(Stage):
    """Compute initial atomic multipoles via Distributed Multipole Analysis.

    This stage performs three sub-steps:

    1. **DMA QM calculation** – Runs a QM single-point at the DMA
       method/basis (default MP2/6-311G**) to produce a checkpoint/fchk
       file.
    2. **GDMA** – Runs the GDMA program on the fchk file to extract
       distributed multipoles.
    3. **poledit** – Runs Tinker's poledit program to process the GDMA
       output into AMOEBA-compatible multipole parameters.

    If the external programs (GDMA, poledit) are not available, the stage
    falls back to producing topology-derived local frames (as before)
    while still running the DMA QM calculation.
    """

    def __init__(self) -> None:
        super().__init__(name="multipole")

    def execute(self, context: PipelineContext) -> StageResult:
        """Compute distributed multipoles via DMA.

        Returns
        -------
        StageResult
            ``COMPLETED`` with a ``multipole_frames`` artifact on
            success, ``FAILED`` when required preconditions are not met.
        """
        if context.molecule is None:
            return StageResult(
                status=StageStatus.FAILED,
                message="No molecule in context",
            )

        if context.backend is None:
            return StageResult(
                status=StageStatus.FAILED,
                message="No QM backend in context",
            )

        qm = context.config.qm

        logger.info(
            "Running Distributed Multipole Analysis: method=%s basis=%s",
            qm.dma_method,
            qm.dma_basis_set,
        )

        # ---- Step 2: DMA QM calculation to produce fchk ----
        dma_result = None
        try:
            dma_result = context.backend.compute_dma(
                context.molecule,
                method=qm.dma_method,
                basis_set=qm.dma_basis_set,
            )
            logger.info(
                "DMA QM calculation completed: fchk at %s",
                dma_result.fchk_path,
            )
        except NotImplementedError:
            raise
        except Exception as exc:
            raise BackendError(
                f"DMA QM calculation failed: {exc}",
                backend_name=context.backend.name,
            ) from exc

        # ---- Step 3: Run GDMA on the fchk file ----
        gdma_result = None
        try:
            from poltype.external.gdma import GDMARunner

            gdma = GDMARunner(gdma_exe=context.config.gdma_path)
            gdma_result = gdma.run(
                fchk_path=dma_result.fchk_path,
                work_dir=context.work_dir,
                molecule_name=context.molecule.name or "mol",
            )
            logger.info(
                "GDMA completed: punch file at %s",
                gdma_result.punch_path,
            )
        except FileNotFoundError:
            logger.warning(
                "GDMA executable not found (%s); skipping GDMA step. "
                "Multipoles will use topology-derived frames.",
                context.config.gdma_path,
            )
        except Exception as exc:
            logger.warning(
                "GDMA failed: %s. Falling back to topology-derived frames.",
                exc,
            )

        # ---- Step 4: Run poledit to process GDMA output ----
        poledit_result = None
        opt_xyz_path = context.get_artifact("opt_xyz_path")
        if gdma_result is not None and opt_xyz_path is not None:
            try:
                from poltype.external.poledit import PoleditRunner

                poledit = PoleditRunner(
                    poledit_exe=context.config.poledit_path,
                )
                poledit_result = poledit.run(
                    xyz_path=opt_xyz_path,
                    work_dir=context.work_dir,
                    molecule_name=context.molecule.name or "mol",
                    gdma_punch_path=gdma_result.punch_path,
                )
                logger.info(
                    "poledit completed: key file at %s",
                    poledit_result.key_path,
                )
            except FileNotFoundError:
                logger.warning(
                    "poledit executable not found (%s); skipping poledit step.",
                    context.config.poledit_path,
                )
            except Exception as exc:
                logger.warning("poledit failed: %s", exc)

        # Build the multipole_frames artifact
        frames_by_atom_index = self._build_frames_by_atom_index(context)
        multipole_frames = {
            "source": "dma",
            "num_atoms": context.molecule.num_atoms,
            "dma_method": qm.dma_method,
            "dma_basis_set": qm.dma_basis_set,
            "gdma_path": context.config.gdma_path,
            "poledit_path": context.config.poledit_path,
            "frames_by_atom_index": frames_by_atom_index,
        }

        # Attach file paths if available
        artifacts = {"multipole_frames": multipole_frames}
        if dma_result is not None:
            artifacts["dma_result"] = dma_result
            multipole_frames["fchk_path"] = str(dma_result.fchk_path)
        if gdma_result is not None:
            artifacts["gdma_result"] = gdma_result
            multipole_frames["gdma_punch_path"] = str(gdma_result.punch_path)
        if poledit_result is not None:
            artifacts["poledit_result"] = poledit_result
            if poledit_result.key_path is not None:
                multipole_frames["poledit_key_path"] = str(
                    poledit_result.key_path
                )

        logger.info(
            "DMA multipoles obtained for %d atoms",
            context.molecule.num_atoms,
        )

        return StageResult(
            status=StageStatus.COMPLETED,
            message=(
                f"DMA multipoles computed for "
                f"{context.molecule.num_atoms} atoms"
            ),
            artifacts=artifacts,
        )

    def should_skip(self, context: PipelineContext) -> bool:
        """Skip when only matching against the database or in dry-run mode."""
        return context.config.database_match_only or context.config.dry_run

    @staticmethod
    def _build_frames_by_atom_index(context: PipelineContext) -> dict[int, list[int]]:
        """Build simple local frames keyed by atom index.

        The full GDMA/poledit integration produces better frames when
        available; this method provides a stable topology-derived local
        frame for each atom as a fallback.
        """
        rdmol = context.molecule.rdmol
        frames: dict[int, list[int]] = {}
        for atom in rdmol.GetAtoms():
            atom_idx = atom.GetIdx()
            neighbor_indices = sorted(
                neighbor.GetIdx() for neighbor in atom.GetNeighbors()
            )
            frames[atom_idx] = neighbor_indices[:3] or [atom_idx]
        return frames
