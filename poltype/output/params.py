"""
poltype.output.params – typed parameter records for AMOEBA key file output.

These dataclasses are the canonical in-memory representation of fitted
AMOEBA parameters.  Pipeline stages that perform fitting (ESP multipole
fitting, polarization fitting, torsion fitting) store lists of these
records as pipeline artifacts so the :class:`~poltype.output.key_writer.KeyFileWriter`
can write them directly to the ``final.key`` file.

Artifact keys
-------------
``multipoles``
    ``list[MultipoleParam]`` – one per unique atom type.
``polarization_params``
    ``list[PolarizeParam]`` – one per unique atom type.
``torsion_params``
    ``list[TorsionParam]`` – one per fitted torsion.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple


@dataclass
class MultipoleParam:
    """AMOEBA multipole parameters for one atom type.

    Attributes
    ----------
    atom_type:
        Tinker atom type integer (1-based).
    frame:
        List of 1–3 frame-defining atom types.  Negative values signal a
        bisector frame (AMOEBA convention).
    charge:
        Monopole (charge) in elementary-charge units.
    dipole:
        Dipole components ``(dx, dy, dz)`` in Debye-Å (AMOEBA units).
    quadrupole:
        Upper-triangular traceless quadrupole tensor components in the
        order ``(Qxx, Qyx, Qyy, Qzx, Qzy, Qzz)``, in Debye-Å².
    """

    atom_type: int
    frame: List[int]
    charge: float
    dipole: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    quadrupole: Tuple[float, float, float, float, float, float] = (
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    )


@dataclass
class PolarizeParam:
    """AMOEBA polarization parameters for one atom type.

    Attributes
    ----------
    atom_type:
        Tinker atom type integer (1-based).
    alpha:
        Isotropic polarizability in Å³.
    thole:
        Thole damping factor (dimensionless).
    neighbors:
        List of directly bonded atom types used for group polarization.
    """

    atom_type: int
    alpha: float
    thole: float = 0.3900
    neighbors: List[int] = field(default_factory=list)


@dataclass
class TorsionFold:
    """One cosine term in an AMOEBA torsion potential.

    Attributes
    ----------
    amplitude:
        Force constant in kcal mol⁻¹.
    phase:
        Phase offset in degrees (typically 0.0 or 180.0).
    periodicity:
        Fold number (1, 2, 3, …).
    """

    amplitude: float
    phase: float
    periodicity: int


@dataclass
class TorsionParam:
    """AMOEBA torsion parameters for one dihedral type.

    Attributes
    ----------
    atom_types:
        Four Tinker atom types ``(a, b, c, d)`` defining the dihedral.
    folds:
        Cosine terms; AMOEBA typically uses up to three folds.
    """

    atom_types: Tuple[int, int, int, int]
    folds: List[TorsionFold] = field(default_factory=list)
