"""
poltype.external.gdma – wrapper for the GDMA program.

GDMA (Gaussian Distributed Multipole Analysis) reads a formatted
checkpoint (fchk) file and produces atomic multipoles up to a
specified rank.
"""

from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


@dataclass
class GDMAResult:
    """Outcome of a GDMA run.

    Attributes
    ----------
    punch_path:
        Path to the GDMA punch file containing the distributed
        multipoles.
    log_path:
        Path to the GDMA log / output file.
    """

    punch_path: Path
    log_path: Optional[Path] = None


class GDMARunner:
    """Run the GDMA program to extract distributed multipoles.

    Parameters
    ----------
    gdma_exe:
        Path to the GDMA executable.
    """

    def __init__(self, gdma_exe: str = "gdma") -> None:
        self.gdma_exe = gdma_exe

    def run(
        self,
        fchk_path: Path,
        work_dir: Path,
        molecule_name: str = "mol",
        multipole_rank: int = 2,
    ) -> GDMAResult:
        """Execute GDMA on the given fchk file.

        Parameters
        ----------
        fchk_path:
            Path to the formatted checkpoint file.
        work_dir:
            Directory for GDMA input/output files.
        molecule_name:
            Base name for generated files.
        multipole_rank:
            Maximum multipole rank (0=charge, 1=dipole, 2=quadrupole).

        Returns
        -------
        GDMAResult
            Contains the path to the punch file with multipoles.
        """
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        fchk_path = Path(fchk_path).resolve()
        if not fchk_path.exists():
            raise FileNotFoundError(
                f"fchk file does not exist: {fchk_path}"
            )
        if fchk_path.suffix != ".fchk":
            raise ValueError(
                f"GDMA requires a .fchk file, but got: {fchk_path}. "
                f"Ensure the QM backend produces a formatted checkpoint file."
            )

        gdma_input_path = work_dir / f"{molecule_name}_gdma.inp"
        gdma_output_path = work_dir / f"{molecule_name}_gdma.out"
        punch_path = work_dir / f"{molecule_name}_gdma.punch"

        # Write GDMA input file
        gdma_input = self._build_input(
            fchk_path=fchk_path,
            punch_path=punch_path,
            multipole_rank=multipole_rank,
        )
        gdma_input_path.write_text(gdma_input)

        logger.info(
            "Running GDMA: %s < %s > %s",
            self.gdma_exe,
            gdma_input_path,
            gdma_output_path,
        )

        with open(gdma_input_path, "r") as fin, \
             open(gdma_output_path, "w") as fout:
            result = subprocess.run(
                [self.gdma_exe],
                stdin=fin,
                stdout=fout,
                stderr=subprocess.PIPE,
                cwd=str(work_dir),
                text=True,
            )

        if result.returncode != 0:
            raise RuntimeError(
                f"GDMA failed (exit code {result.returncode}).\n"
                f"stderr: {result.stderr[:2000]}"
            )

        if not punch_path.exists():
            raise RuntimeError(
                f"GDMA did not produce punch file: {punch_path}"
            )

        logger.info("GDMA completed: punch file at %s", punch_path)
        return GDMAResult(punch_path=punch_path, log_path=gdma_output_path)

    @staticmethod
    def _build_input(
        fchk_path: Path,
        punch_path: Path,
        multipole_rank: int = 2,
    ) -> str:
        """Build the GDMA input file content.

        The input instructs GDMA to read the fchk file, analyse the
        density, and write the multipoles to a punch file.
        """
        return (
            f"Title \"DMA analysis\"\n"
            f"\n"
            f"File {fchk_path}\n"
            f"Multipoles\n"
            f"  Limit {multipole_rank}\n"
            f"  Punch {punch_path}\n"
            f"Start\n"
            f"\n"
            f"Finish\n"
        )
