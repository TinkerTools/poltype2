"""
poltype.external.poledit – wrapper for Tinker's poledit program.

Poledit processes the raw distributed multipoles from GDMA and
prepares them for use in the AMOEBA force field (e.g. rotating
multipoles into local frames, assigning polarization groups).
"""

from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


@dataclass
class PoleditResult:
    """Outcome of a poledit run.

    Attributes
    ----------
    key_path:
        Path to the key file with processed multipole parameters.
    log_path:
        Path to the poledit log / output file.
    """

    key_path: Optional[Path] = None
    log_path: Optional[Path] = None


class PoleditRunner:
    """Run Tinker's poledit program to process multipoles.

    Parameters
    ----------
    poledit_exe:
        Path to the poledit executable.
    """

    def __init__(self, poledit_exe: str = "poledit") -> None:
        self.poledit_exe = poledit_exe

    def run(
        self,
        xyz_path: Path,
        work_dir: Path,
        molecule_name: str = "mol",
        gdma_punch_path: Optional[Path] = None,
    ) -> PoleditResult:
        """Execute poledit to process multipoles from GDMA output.

        Parameters
        ----------
        xyz_path:
            Path to the Tinker XYZ file for the molecule.
        work_dir:
            Directory for poledit input/output files.
        molecule_name:
            Base name for generated files.
        gdma_punch_path:
            Optional path to the GDMA punch file.  If provided, poledit
            reads multipoles from this file.

        Returns
        -------
        PoleditResult
            Contains the path to the generated key file.
        """
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        log_path = work_dir / f"{molecule_name}_poledit.out"

        # Build the poledit command.  Poledit typically reads from stdin
        # for interactive prompts; here we pipe the expected responses.
        cmd = [self.poledit_exe]

        # Poledit expects option number followed by the xyz file
        # Option 6 = "Process GDMA Multipole File"
        input_text = f"6\n{xyz_path}\n"
        if gdma_punch_path is not None:
            input_text += f"{gdma_punch_path}\n"
        input_text += "\n"

        logger.info("Running poledit: %s", " ".join(cmd))

        with open(log_path, "w") as fout:
            result = subprocess.run(
                cmd,
                input=input_text,
                stdout=fout,
                stderr=subprocess.PIPE,
                cwd=str(work_dir),
                text=True,
            )

        if result.returncode != 0:
            raise RuntimeError(
                f"poledit failed (exit code {result.returncode}).\n"
                f"stderr: {result.stderr[:2000]}"
            )

        # Look for the generated key file
        key_path = xyz_path.with_suffix(".key")
        if not key_path.exists():
            key_path = work_dir / f"{molecule_name}.key"
            if not key_path.exists():
                key_path = None

        logger.info("poledit completed: key file at %s", key_path)
        return PoleditResult(key_path=key_path, log_path=log_path)
