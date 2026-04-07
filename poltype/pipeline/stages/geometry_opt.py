"""
poltype.pipeline.stages.geometry_opt – QM geometry optimisation stage.

Delegates to :meth:`QMBackend.optimize_geometry` to run a QM geometry
optimisation and updates the molecule in the pipeline context with the
optimised coordinates.  On success, writes the optimised geometry to
``{fname}_opt.xyz`` in the working directory.
"""

from __future__ import annotations

import logging
from pathlib import Path

from poltype.errors import BackendError, StageError
from poltype.pipeline.context import PipelineContext
from poltype.pipeline.stage import Stage, StageResult, StageStatus

logger = logging.getLogger(__name__)


class GeometryOptimizationStage(Stage):
    """Run QM geometry optimisation on the current molecule.

    Calls :meth:`QMBackend.optimize_geometry` with the method and
    basis set from ``context.config.qm``, then replaces the molecule
    in the context with one carrying the optimised coordinates and
    writes the result to ``{fname}_opt.xyz``.
    """

    def __init__(self) -> None:
        super().__init__(name="geometry_optimization")

    def execute(self, context: PipelineContext) -> StageResult:
        """Optimise geometry via the configured QM backend.

        Returns
        -------
        StageResult
            ``COMPLETED`` with ``opt_result`` and ``molecule`` artifacts
            on success, ``FAILED`` on error.
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
        molecule = context.molecule

        logger.info(
            "Running geometry optimisation: method=%s basis=%s",
            qm.opt_method,
            qm.opt_basis_set,
        )

        try:
            opt_result = context.backend.optimize_geometry(
                molecule,
                method=qm.opt_method,
                basis_set=qm.opt_basis_set,
            )
        except NotImplementedError:
            raise
        except Exception as exc:
            raise BackendError(
                f"Geometry optimisation failed: {exc}",
                backend_name=context.backend.name,
            ) from exc

        if not opt_result.converged:
            return StageResult(
                status=StageStatus.FAILED,
                message="Geometry optimisation did not converge",
                artifacts={"opt_result": opt_result},
            )

        new_mol = molecule.with_updated_coordinates(opt_result.coordinates)
        context.update_molecule(new_mol)

        # Write optimised geometry to {fname}_opt.xyz
        opt_xyz_path = self._write_opt_xyz(new_mol, context.work_dir)

        logger.info(
            "Geometry optimisation converged: energy=%.6f Hartree",
            opt_result.energy,
        )

        artifacts = {
            "opt_result": opt_result,
            "molecule": new_mol,
        }
        if opt_xyz_path is not None:
            artifacts["opt_xyz_path"] = opt_xyz_path

        return StageResult(
            status=StageStatus.COMPLETED,
            message=f"Optimised geometry (energy={opt_result.energy:.6f} Ha)",
            artifacts=artifacts,
        )

    def should_skip(self, context: PipelineContext) -> bool:
        """Skip when only matching against the database or in dry-run mode."""
        return context.config.database_match_only or context.config.dry_run

    @staticmethod
    def _write_opt_xyz(molecule, work_dir: Path) -> Path:
        """Write the optimised geometry to a Tinker-style XYZ file.

        The file is named ``{molecule.name}_opt.xyz`` and placed in the
        working directory.

        Parameters
        ----------
        molecule:
            Molecule with optimised coordinates.
        work_dir:
            Directory to write the file into.

        Returns
        -------
        Path
            Absolute path to the written XYZ file.
        """
        stem = molecule.name or "mol"
        xyz_path = Path(work_dir) / f"{stem}_opt.xyz"

        coords = molecule.coordinates
        rdmol = molecule.rdmol
        n_atoms = molecule.num_atoms

        lines = [f"    {n_atoms}  {stem} optimised geometry"]
        for i in range(n_atoms):
            atom = rdmol.GetAtomWithIdx(i)
            sym = atom.GetSymbol()
            x, y, z = coords[i]
            # 1-based atom index, symbol, coords, atom type (placeholder 0),
            # then bonded atom indices (1-based)
            neighbors = sorted(
                nb.GetIdx() + 1 for nb in atom.GetNeighbors()
            )
            nb_str = "  ".join(str(n) for n in neighbors)
            lines.append(
                f"  {i + 1:>4}  {sym:<2}  {x:>12.6f}  {y:>12.6f}  "
                f"{z:>12.6f}     0  {nb_str}"
            )

        xyz_path.write_text("\n".join(lines) + "\n")
        logger.info("Wrote optimised geometry to %s", xyz_path)
        return xyz_path
