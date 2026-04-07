"""
poltype.external.potential – wrapper for Tinker's potential program.

The potential program fits atomic multipoles to reproduce the
quantum-mechanical electrostatic potential computed on a grid around
the molecule.
"""

from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


@dataclass
class PotentialResult:
    """Outcome of a Tinker potential fitting run.

    Attributes
    ----------
    key_path:
        Path to the key file with fitted multipole parameters.
    log_path:
        Path to the potential program log / output file.
    """

    key_path: Optional[Path] = None
    log_path: Optional[Path] = None


class PotentialRunner:
    """Run Tinker's potential program to fit multipoles to ESP.

    Parameters
    ----------
    potential_exe:
        Path to the potential executable.
    """

    def __init__(self, potential_exe: str = "potential") -> None:
        self.potential_exe = potential_exe

    def run(
        self,
        xyz_path: Path,
        work_dir: Path,
        molecule_name: str = "mol",
        esp_grid_path: Optional[Path] = None,
    ) -> PotentialResult:
        """Execute the potential program for ESP multipole fitting.

        Parameters
        ----------
        xyz_path:
            Path to the Tinker XYZ file for the molecule.
        work_dir:
            Directory for potential input/output files.
        molecule_name:
            Base name for generated files.
        esp_grid_path:
            Optional path to the ESP grid file (potential reads this).

        Returns
        -------
        PotentialResult
            Contains the path to the fitted key file.
        """
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)

        log_path = work_dir / f"{molecule_name}_potential.out"

        cmd = [self.potential_exe]

        # Potential program option 6 = "Fit Electrostatic Parameters to Map"
        input_text = f"6\n{xyz_path}\nY\n"
        if esp_grid_path is not None:
            input_text += f"{esp_grid_path}\n"
        input_text += "\n"

        logger.info("Running potential: %s", " ".join(cmd))

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
                f"potential failed (exit code {result.returncode}).\n"
                f"stderr: {result.stderr[:2000]}"
            )

        # Look for the generated/updated key file
        key_path = xyz_path.with_suffix(".key_2")
        if not key_path.exists():
            key_path = xyz_path.with_suffix(".key")
            if not key_path.exists():
                key_path = work_dir / f"{molecule_name}.key"
                if not key_path.exists():
                    key_path = None

        logger.info("potential completed: key file at %s", key_path)
        return PotentialResult(key_path=key_path, log_path=log_path)
