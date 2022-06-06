# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

from binana._utils.shim import _set_default
from binana.interactions.default_params import (
    PI_PADDING_DIST,
    PI_PI_INTERACTING_DIST_CUTOFF,
    PI_STACKING_ANGLE_TOLERANCE,
    T_STACKING_ANGLE_TOLERANCE,
    T_STACKING_CLOSEST_DIST_CUTOFF,
)
import binana
from binana.load_ligand_receptor import (
    _get_ligand_receptor_aromatic_dists,
)

# _get_ligand_receptor_dists,
from binana._utils.utils import hashtable_entry_add_one, list_alphebetize_and_combine
from binana._structure.mol import Mol
from binana._utils._math_functions import project_point_onto_plane

# __pragma__ ('skip')
# Python
from math import fabs

# __pragma__ ('noskip')

"""?
# Transcrypt
from binana._utils.shim import fabs
?"""


def _t_stacking(
    ligand,
    receptor,
    ligand_aromatic,
    receptor_aromatic,
    dist,
    angle_between_planes,
    t_stacking_angle_tol,
    t_stacking_closest_dist_cutoff,
    pi_padding,
    pi_pi_interactions,
    pdb_pi_t,
    t_stacking_labels,
):

    if (
        min(fabs(angle_between_planes - 90), fabs(angle_between_planes - 270))
        < t_stacking_angle_tol
    ):
        # so they're more or less perpendicular, it's probably a pi-edge
        # interaction

        # having looked at many structures, I noticed the algorithm was
        # identifying T-pi reactions when the two rings were in fact quite
        # distant, often with other atoms in between. Eye-balling it, requiring
        # that at their closest they be at least 5 A apart seems to separate the
        # good T's from the bad
        min_dist = 100.0
        for ligand_ind in ligand_aromatic.indices:
            ligand_at = ligand.all_atoms[ligand_ind]
            for receptor_ind in receptor_aromatic.indices:
                receptor_at = receptor.all_atoms[receptor_ind]
                dist = ligand_at.coordinates.dist_to(receptor_at.coordinates)
                if dist < min_dist:
                    min_dist = dist

        if min_dist <= t_stacking_closest_dist_cutoff:
            # so at their closest points, the two rings come within 5 A of each
            # other.

            # okay, is the ligand pi pointing into the receptor pi, or the other
            # way around? first, project the center of the ligand pi onto the
            # plane of the receptor pi, and vs. versa

            # This could be directional somehow, like a hydrogen bond.

            pt_on_receptor_plane = project_point_onto_plane(
                ligand_aromatic.center, receptor_aromatic.plane_coeff
            )
            pt_on_lignad_plane = project_point_onto_plane(
                receptor_aromatic.center, ligand_aromatic.plane_coeff
            )

            # now, if it's a true pi-T interaction, this projected point should
            # fall within the ring whose plane it's been projected into.
            if (
                pt_on_receptor_plane.dist_to(receptor_aromatic.center)
                <= receptor_aromatic.radius + pi_padding
            ) or (
                pt_on_lignad_plane.dist_to(ligand_aromatic.center)
                <= ligand_aromatic.radius + pi_padding
            ):
                # so it is in the ring on the projected plane.
                structure = receptor.all_atoms[receptor_aromatic.indices[0]].structure
                if structure == "":
                    # since it could be interacting with a cofactor or something
                    structure = "OTHER"

                key = "T-SHAPED_" + structure

                for index in ligand_aromatic.indices:
                    pdb_pi_t.add_new_atom(ligand.all_atoms[index].copy_of())
                for index in receptor_aromatic.indices:
                    pdb_pi_t.add_new_atom(receptor.all_atoms[index].copy_of())

                hashtable_entry_add_one(pi_pi_interactions, key)

                t_stacking_labels.append(
                    _make_pi_pi_interaction_label(
                        ligand,
                        ligand_aromatic,
                        receptor,
                        receptor_aromatic,
                        {
                            "distance": dist,
                            "angle": min(
                                fabs(angle_between_planes - 0),
                                fabs(angle_between_planes - 180),
                            ),
                        },
                    )
                )

    return (
        pi_pi_interactions,
        pdb_pi_t,
        t_stacking_labels,
    )


def _make_pi_pi_interaction_label(
    ligand, ligand_aromatic, receptor, receptor_aromatic, metric
):
    return (
        "["
        + " / ".join(
            [ligand.all_atoms[index].string_id() for index in ligand_aromatic.indices]
        )
        + "]",
        "["
        + " / ".join(
            [
                receptor.all_atoms[index].string_id()
                for index in receptor_aromatic.indices
            ]
        )
        + "]",
        metric,
    )


def _pi_pi_detect_by_projecting_all_ring_atoms(
    mol1, mol1_aromatic, mol2_aromatic, pi_padding
):
    for mol1_ring_index in mol1_aromatic.indices:
        # project the mol1 atom onto the plane of the mol2 ring
        pt_on_mol2_plane = project_point_onto_plane(
            mol1.all_atoms[mol1_ring_index].coordinates,
            mol2_aromatic.plane_coeff,
        )
        if (
            pt_on_mol2_plane.dist_to(mol2_aromatic.center)
            <= mol2_aromatic.radius + pi_padding
        ):
            # Detected
            return True
            # pi_pi = True
            # break
    return False


def _pi_stacking(
    ligand,
    receptor,
    ligand_aromatic,
    receptor_aromatic,
    dist,
    angle_between_planes,
    pi_stacking_angle_tol,
    pi_padding,
    pi_pi_interactions,
    pdb_pistack,
    pi_stacking_labels,
):
    angle = min(fabs(angle_between_planes - 0), fabs(angle_between_planes - 180))
    if angle < pi_stacking_angle_tol:
        # so they're more or less parallel, it's probably pi-pi
        # stackingoutput_dir now, pi-pi are not usually right on top of each
        # other. They're often staggared. So I don't want to just look at the
        # centers of the rings and compare. Let's look at each of the atoms. do
        # atom of the atoms of one ring, when projected onto the plane of the
        # other, fall within that other ring?

        # Check the ligand atoms projected onto the receptor ring.
        pi_pi = _pi_pi_detect_by_projecting_all_ring_atoms(
            ligand, ligand_aromatic, receptor_aromatic, pi_padding
        )

        # if you've already determined it's a pi-pi stacking interaction, no
        # need to keep trying.
        if not pi_pi:
            # Check the receptor atoms projected onto the ligand ring.
            pi_pi = _pi_pi_detect_by_projecting_all_ring_atoms(
                receptor, receptor_aromatic, ligand_aromatic, pi_padding
            )

        if pi_pi:
            # It's a pi-pi stacking interaction.
            structure = receptor.all_atoms[receptor_aromatic.indices[0]].structure

            # since it could be interacting with a cofactor or something
            structure = "OTHER" if structure == "" else structure
            key = "STACKING_" + structure

            for index in ligand_aromatic.indices:
                pdb_pistack.add_new_atom(ligand.all_atoms[index].copy_of())
            for index in receptor_aromatic.indices:
                pdb_pistack.add_new_atom(receptor.all_atoms[index].copy_of())

            hashtable_entry_add_one(pi_pi_interactions, key)

            pi_stacking_labels.append(
                _make_pi_pi_interaction_label(
                    ligand,
                    ligand_aromatic,
                    receptor,
                    receptor_aromatic,
                    {"distance": dist, "angle": angle},
                )
            )
        pi_stacking_detected = True
    else:
        pi_stacking_detected = False

    return (pi_pi_interactions, pdb_pistack, pi_stacking_labels, pi_stacking_detected)


# Be sure to update the corresponding function in
# binana.interactions.__init__.py as well!


def get_pi_pi(
    ligand,
    receptor,
    pi_pi_general_dist_cutoff=None,
    pi_stacking_angle_tol=None,
    t_stacking_angle_tol=None,
    t_stacking_closest_dist_cutoff=None,
    pi_padding=None,
):
    """Identifies and counts the number of pi-pi stacking and T-shaped
    interactions between the protein and ligand. Output is formatted like
    this::

        {
            'labels': {
                'T_stacking': [
                    ('[A:CHT(1):C6(4) / A:CHT(1):C7(5) / A:CHT(1):C8(6) / A:CHT(1):C9(7) / A:CHT(1):O2(8)]', '[A:PHE(233):CG(657) / A:PHE(233):CD1(658) / A:PHE(233):CE1(660) / A:PHE(233):CZ(662) / A:PHE(233):CE2(661) / A:PHE(233):CD2(659)]'),
                    ('[A:CHT(1):C2(17) / A:CHT(1):O1(18) / A:CHT(1):C5(19) / A:CHT(1):C4(20) / A:CHT(1):C3(21)]', '[A:TRP(43):CG(28) / A:TRP(43):CD1(29) / A:TRP(43):NE1(31) / A:TRP(43):CE2(32) / A:TRP(43):CD2(30)]')
                ],
                'pi_stacking': [
                    ('[A:CHT(1):C6(4) / A:CHT(1):C7(5) / A:CHT(1):C8(6) / A:CHT(1):C9(7) / A:CHT(1):O2(8)]', '[A:TRP(90):CG(100) / A:TRP(90):CD1(101) / A:TRP(90):NE1(103) / A:TRP(90):CE2(104) / A:TRP(90):CD2(102)]'),
                    ('[A:CHT(1):C6(4) / A:CHT(1):C7(5) / A:CHT(1):C8(6) / A:CHT(1):C9(7) / A:CHT(1):O2(8)]', '[A:TRP(90):CE2(104) / A:TRP(90):CD2(102) / A:TRP(90):CE3(105) / A:TRP(90):CZ3(107) / A:TRP(90):CH2(108) / A:TRP(90):CZ2(106)]')
                ]
            },
            'counts': {
                'STACKING_BETA': 2,
                'T-SHAPED_OTHER': 3
            },
            'mols': {
                'T_stacking': <binana._structure.mol.Mol instance at 0x7feb20478fc8>,
                'pi_stacking': <binana._structure.mol.Mol instance at 0x7feb20478f80>
            }
        }

    Args:
        ligand (binana._structure.mol.Mol): The ligand molecule to analyze.
        receptor (binana._structure.mol.Mol): The receptor molecule to analyze.
        pi_pi_general_dist_cutoff (float, optional): The distance cutoff used
            for all pi-pi interactions (stacking and T-shaped). Defaults to
            PI_PI_INTERACTING_DIST_CUTOFF.
        pi_stacking_angle_tol (float, optional): The angle tolerance for the
            pi-pi stacking interactions. Defaults to
            PI_STACKING_ANGLE_TOLERANCE.
        t_stacking_angle_tol (float, optional): The angle tolerance for the
            T-shaped interactions. Defaults to T_STACKING_ANGLE_TOLERANCE.
        t_stacking_closest_dist_cutoff (float, optional): The distance cutoff
            for T-shaped interactions specifically. Defaults to
            T_STACKING_CLOSEST_DIST_CUTOFF.
        pi_padding (float, optional): The amount by which the radius of each pi
            ring should be artificially expanded, to be sure to catch the
            interactions. Defaults to PI_PADDING_DIST.

    Returns:
        dict: Contains the atom tallies ("counts"), the
        binana._structure.mol.Mol objects with the participating atoms
        ("mols"), and the labels to use in the log file ("labels").
    """

    pi_pi_general_dist_cutoff = _set_default(
        pi_pi_general_dist_cutoff, PI_PI_INTERACTING_DIST_CUTOFF
    )
    pi_stacking_angle_tol = _set_default(
        pi_stacking_angle_tol, PI_STACKING_ANGLE_TOLERANCE
    )
    t_stacking_angle_tol = _set_default(
        t_stacking_angle_tol, T_STACKING_ANGLE_TOLERANCE
    )
    t_stacking_closest_dist_cutoff = _set_default(
        t_stacking_closest_dist_cutoff, T_STACKING_CLOSEST_DIST_CUTOFF
    )
    pi_padding = _set_default(pi_padding, PI_PADDING_DIST)

    # Calculate the distances.
    ligand_receptor_aromatic_dists = _get_ligand_receptor_aromatic_dists(
        ligand, receptor, pi_pi_general_dist_cutoff
    )

    pi_interactions = {}
    pdb_pistack = Mol()
    pdb_pi_t = Mol()
    pi_stacking_labels = []
    t_stacking_labels = []

    # "PI-Stacking Interactions ALIVE AND WELL IN PROTEINS" says distance of 7.5
    # A is good cutoff. This seems really big to me, except that pi-pi
    # interactions (parallel) are actuall usually off centered. Interesting
    # paper. Note that adenine and tryptophan count as two aromatic rings. So,
    # for example, an interaction between these two, if positioned correctly,
    # could count for 4 pi-pi interactions.
    for (
        ligand_aromatic,
        receptor_aromatic,
        dist,
        angle_between_planes,
    ) in ligand_receptor_aromatic_dists:

        (
            pi_interactions,
            pdb_pistack,
            pi_stacking_labels,
            pi_stacking_detected,
        ) = _pi_stacking(
            ligand,
            receptor,
            ligand_aromatic,
            receptor_aromatic,
            dist,
            angle_between_planes,
            pi_stacking_angle_tol,
            pi_padding,
            pi_interactions,
            pdb_pistack,
            pi_stacking_labels,
        )

        if not pi_stacking_detected:
            pi_interactions, pdb_pi_t, t_stacking_labels = _t_stacking(
                ligand,
                receptor,
                ligand_aromatic,
                receptor_aromatic,
                dist,
                angle_between_planes,
                t_stacking_angle_tol,
                t_stacking_closest_dist_cutoff,
                pi_padding,
                pi_interactions,
                pdb_pi_t,
                t_stacking_labels,
            )

    return {
        "counts": pi_interactions,
        "mols": {"pi_stacking": pdb_pistack, "T_stacking": pdb_pi_t},
        "labels": {"pi_stacking": pi_stacking_labels, "T_stacking": t_stacking_labels},
    }
