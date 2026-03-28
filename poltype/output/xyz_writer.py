"""
poltype.output.xyz_writer – Tinker-format XYZ file writer.

Writes the current molecule coordinates in the Tinker ``.xyz`` format::

    N  molecule_name
    1  C     0.000000    0.000000    0.000000    1    2    3
    2  H     1.089000    0.000000    0.000000    5    1
    ...

The writer uses the molecule present in the pipeline context (which may
be the optimised geometry if geometry optimisation was performed).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from poltype.pipeline.context import PipelineContext

from poltype.output.writer import OutputWriter

logger = logging.getLogger(__name__)


class XYZFileWriter(OutputWriter):
    """Write a Tinker-format ``.xyz`` file from the pipeline context.

    Parameters
    ----------
    output_dir:
        Directory for the output file.  Falls back to ``context.work_dir``.
    filename:
        Name of the output file.  Defaults to ``"final.xyz"``.
    """

    def __init__(
        self,
        output_dir: Path | None = None,
        filename: str = "final.xyz",
    ) -> None:
        super().__init__(output_dir=output_dir)
        self._filename = filename

    @property
    def name(self) -> str:
        return "xyz_writer"

    def write(self, context: "PipelineContext") -> Path:
        """Write the Tinker XYZ file.

        Returns
        -------
        Path
            Absolute path to the written ``final.xyz``.

        Raises
        ------
        ValueError
            If no molecule is available in the context.
        """
        mol = context.molecule
        if mol is None:
            raise ValueError("No molecule in context; cannot write XYZ file")

        out_dir = self._resolve_dir(context)
        out_path = out_dir / self._filename

        coords = mol.coordinates
        rdmol = mol.rdmol  # RDKit Mol for atom info
        num_atoms = mol.num_atoms

        lines = [f" {num_atoms}  {mol.name}"]

        for i in range(num_atoms):
            atom = rdmol.GetAtomWithIdx(i)
            symbol = atom.GetSymbol()
            x, y, z = coords[i]
            # Tinker uses 1-based indexing
            tinker_idx = i + 1
            # Collect bonded neighbours (1-based)
            neighbours = sorted(
                n.GetIdx() + 1 for n in atom.GetNeighbors()
            )
            neighbour_str = "  ".join(str(n) for n in neighbours)
            # Type placeholder: use atom index as type until real typing
            atom_type = tinker_idx
            lines.append(
                f" {tinker_idx:>5d}  {symbol:<2s}"
                f"  {x:12.6f}{y:12.6f}{z:12.6f}"
                f"  {atom_type:>5d}"
                f"  {neighbour_str}"
            )

        out_path.write_text("\n".join(lines) + "\n")
        logger.info("Wrote Tinker XYZ file: %s", out_path)
        return out_path
