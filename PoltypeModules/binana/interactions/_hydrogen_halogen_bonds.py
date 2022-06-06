# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

from binana._utils.shim import _set_default
from binana.interactions.default_params import (
    HYDROGEN_HALOGEN_BOND_ANGLE_CUTOFF,
    HYDROGEN_BOND_DIST_CUTOFF,
    HALOGEN_BOND_DIST_CUTOFF,
)
import binana
from binana.load_ligand_receptor import _get_ligand_receptor_dists
from binana._utils.utils import hashtable_entry_add_one
from binana._structure.mol import Mol
from binana._utils._math_functions import angle_between_three_points
from binana._structure.consts import to_deg

import __future__

import math

# __pragma__ ('skip')
# Python
from math import fabs

# __pragma__ ('noskip')

"""?
# Transcrypt
from binana._utils.shim import fabs
?"""

# Be sure to update the corresponding function in
# binana.interactions.__init__.py as well!


def _get_potential_donors_acceptors(ligand, receptor, dist_cutoff, hydrogen_bond=True):
    # Any that are close to each other (not considering orientation yet).

    donors_and_acceptors = ["O", "N", "S"]

    # The donor can be a carbon if halogen bond. See
    # https://macmillan.princeton.edu/wp-content/uploads/HalogenBonding-min.pdf
    if not hydrogen_bond:
        donors_and_acceptors.append("C")

    # Calculate the distances.
    ligand_receptor_dists = _get_ligand_receptor_dists(
        ligand, receptor, dist_cutoff, donors_and_acceptors
    )

    return [
        [ligand_atom, receptor_atom, dist]
        for ligand_atom, receptor_atom, dist in ligand_receptor_dists
    ]


def _update_mol_and_data(
    pdb_hbonds,
    hbonds,
    hbonds_labels,
    lig_donor_or_accept,
    receptor_atom,
    ligand_atom,
    center_atom,
    dist,
    angle,
):
    comment = "RECEPTOR" if lig_donor_or_accept == "ACCEPTOR" else "LIGAND"

    hbonds_key = (
        "HDONOR_"
        + comment
        + "_"
        + receptor_atom.side_chain_or_backbone()
        + "_"
        + receptor_atom.structure
    )

    pdb_hbonds.add_new_atom(ligand_atom.copy_of())
    pdb_hbonds.add_new_atom(center_atom.copy_of())
    pdb_hbonds.add_new_atom(receptor_atom.copy_of())
    hashtable_entry_add_one(hbonds, hbonds_key)

    hbonds_labels.append(
        (
            ligand_atom.string_id(),
            center_atom.string_id(),
            receptor_atom.string_id(),
            comment,
            {
                "distance": dist,
                
                # Because if no hydrogens, angle might be None
                "angle": angle,
            },
        )
    )


def _product(lst1, lst2):
    # Below because I had problems compiling itertools.product using
    # transcrypt.
    combos = []
    for l1 in lst1:
        for l2 in lst2:
            combos.append([l1, l2])
    # combos = product(lig_atm_hbond_infs, recep_atm_hbond_infs)
    return combos


def _collect_bonds(bonds_organized_by_donor, pdb_hbonds, hbonds, hbonds_labels):
    # You need to sort by the index. This is just for backwards compatibility.
    unwrapped_bond_infos = []
    for donor_key in bonds_organized_by_donor.keys():
        for bond_inf in bonds_organized_by_donor[donor_key]:
            unwrapped_bond_infos.append(bond_inf)
    unwrapped_bond_infos = sorted(unwrapped_bond_infos, key=lambda i: i[0])

    # Now update pdb_hbonds, hbonds, and hbonds_labels (in place)
    for (
        idx,
        lig_donor_or_accept,
        receptor_atom,
        ligand_atom,
        center_atom,
        dist,
        angle,
    ) in unwrapped_bond_infos:
        _update_mol_and_data(
            pdb_hbonds,
            hbonds,
            hbonds_labels,
            lig_donor_or_accept,
            receptor_atom,
            ligand_atom,
            center_atom,
            dist,
            angle,
        )


def _score_angle_deviation_from_sp3_sp2(angle, donor_has_sp3_geometry):
    # Should be 109 if sp3, or 120 if sp2.
    # if angle < 79 or angle > 150:
    #     # catastrophically bad
    #     return 10000
    # else:

    if donor_has_sp3_geometry is None:
        diff = min(fabs(109 - angle), fabs(120 - angle))
        min_angle = 89
        max_angle = 150
    elif donor_has_sp3_geometry == True:
        diff = fabs(109 - angle)
        min_angle = 89
        max_angle = 129
    else:
        diff = fabs(120 - angle)
        min_angle = 100
        max_angle = 150

    return (diff, angle < min_angle or angle > max_angle)

    # return fabs(109 - angle) if donor_has_sp3_geometry else fabs(120 - angle)

    # if angle < 70 or angle > 160:
    #     return 2
    # elif angle < 90 or angle > 140:
    #     return 1
    # return 0


def _select_acceptor(lig, recep, lig_donor_or_accept):
    return lig if lig_donor_or_accept == "ACCEPTOR" else recep


def _remove_extra_noh_hydrogen_bonds(
    ligand, receptor, acceptor_donor_atoms, bonds_organized_by_donor
):
    # This only used if no hydrogens specified.
    for donor_key in bonds_organized_by_donor.keys():
        bond_infos = bonds_organized_by_donor[donor_key]
        donor_atom = acceptor_donor_atoms[donor_key]
        donor_coor = donor_atom.coordinates

        lig_donor_or_accept = bond_infos[0][1]  # Just pick first (0th) item.
        donor_mol = ligand if lig_donor_or_accept == "DONOR" else receptor

        num_neighbors = donor_atom.number_of_neighbors()
        if num_neighbors == 0:
            donor_mol.create_bond_by_distance(donor_atom)
            num_neighbors = donor_atom.number_of_neighbors()

        # if donor_atom.residue[-3:] == "ASN" and donor_atom.atom_name.strip() == "ND2":
        # if donor_atom.residue[-3:] == "TRP" and donor_atom.atom_name.strip() == "NE1":
        # if donor_atom.residue[-3:] == "ARG" and donor_atom.atom_name.strip() == "NH2":
        # Should be sp2, but will return sp3...
        # import pdb; pdb.set_trace()

        if donor_atom.belongs_to_protein():
            donor_has_sp3_geometry = donor_atom.has_sp3_geometry(donor_mol)
        else:
            donor_has_sp3_geometry = (
                None if num_neighbors == 1 else donor_atom.has_sp3_geometry(donor_mol)
            )

        donor_neighbor_coors = [
            donor_mol.all_atoms[i].coordinates
            for i in donor_atom.indecies_of_atoms_connecting
        ]

        max_hydrogen_atoms = 0
        if donor_atom.element in ["O", "S"]:
            max_hydrogen_atoms = 2 - num_neighbors
        if donor_atom.element == "N":
            # Note using 3 here, not 4 (neutral).
            max_hydrogen_atoms = 3 - num_neighbors

        # If the number of hydrogen bonds is less than the max number of
        # hydrogens, no need to change anything.
        # if len(bond_infos) <= max_hydrogen_atoms:
        # continue

        # First, score based on the angle between the donor neighbor, donor, and
        # acceptor.
        scores_and_bond_infs = []
        for bond_inf in bond_infos:
            lig_donor_or_accept = bond_inf[1]
            receptor_atom = bond_inf[2]
            ligand_atom = bond_inf[3]
            bad_score = 0

            acceptor_coor = _select_acceptor(
                ligand_atom.coordinates, receptor_atom.coordinates, lig_donor_or_accept
            )

            neighbor_donor_acceptor_angles = [
                angle_between_three_points(
                    donor_neighbor_coor, donor_coor, acceptor_coor
                )
                * to_deg
                for donor_neighbor_coor in donor_neighbor_coors
            ]

            bad_scores = [
                _score_angle_deviation_from_sp3_sp2(
                    neighbor_donor_acceptor_angle, donor_has_sp3_geometry
                )
                for neighbor_donor_acceptor_angle in neighbor_donor_acceptor_angles
            ]

            if True in [a[1] for a in bad_scores]:
                # There's a catastrophic problem
                bad_score += 10000
            else:
                # Nothing too catastrophic (use worst score based on angle).
                bad_score += max(s[0] for s in bad_scores)

            scores_and_bond_infs.append([bad_score, bond_inf])

        # Now go through all donor-acceptor-donor angles. Flag pairs of hydrogen
        # bonds that are too small. Should be 109 if neighbor is sp3, or 120 if
        # sp2.
        for idx1 in range(len(scores_and_bond_infs) - 1):
            score1, bond_inf1 = scores_and_bond_infs[idx1]
            acceptor1 = _select_acceptor(
                bond_inf1[3].coordinates,  # ligand
                bond_inf1[2].coordinates,  #  receptor
                bond_inf1[1],  # lig_donor_or_accept
            )
            dist1 = bond_inf1[5]
            for idx2 in range(idx1 + 1, len(scores_and_bond_infs)):
                score2, bond_inf2 = scores_and_bond_infs[idx2]
                acceptor2 = _select_acceptor(
                    bond_inf2[3].coordinates,  # ligand
                    bond_inf2[2].coordinates,  #  receptor
                    bond_inf2[1],  # lig_donor_or_accept
                )
                dist2 = bond_inf2[5]

                acceptor_donor_acceptor_angle = (
                    angle_between_three_points(acceptor1, donor_coor, acceptor2)
                    * to_deg
                )

                bad_score, catastrophic = _score_angle_deviation_from_sp3_sp2(
                    acceptor_donor_acceptor_angle, donor_has_sp3_geometry
                )

                # If catastrophic, penalize only one. Otherwise, penalize both.
                # Use distance to decide.
                if dist1 > dist2:
                    scores_and_bond_infs[idx1][0] = score1 + (
                        10000 if catastrophic else bad_score
                    )
                    scores_and_bond_infs[idx2][0] = score2 + bad_score
                else:
                    scores_and_bond_infs[idx2][0] = score2 + (
                        10000 if catastrophic else bad_score
                    )
                    scores_and_bond_infs[idx1][0] = score1 + bad_score

        # Remove any with catastrophic scores.
        scores_and_bond_infs = [
            sbinf for sbinf in scores_and_bond_infs if sbinf[0] < 10000
        ]

        scores_and_bond_infs = sorted(scores_and_bond_infs, key=lambda i: i[0])
        scores_and_bond_infs = scores_and_bond_infs[:max_hydrogen_atoms]

        # Update, but remove score
        bonds_organized_by_donor[donor_key] = [s[1] for s in scores_and_bond_infs]


def _get_hydrogen_or_halogen_bonds(
    ligand, receptor, dist_cutoff=None, angle_cutoff=None, hydrogen_bond=True
):
    """Identifies and counts the number of hydrogen or halogen bonds between
    the protein and ligand. Output is formatted like this::

        {
            'counts': {
                'HDONOR_RECEPTOR_SIDECHAIN_OTHER': 1,
                'HDONOR_LIGAND_SIDECHAIN_OTHER': 2
            },
            'labels': [
                ('A:CHT(1):N1(14)', 'A:CHT(1):H1(16)', 'A:ASP(157):OD2(285)', 'LIGAND'),
                ('A:CHT(1):O6(22)', 'A:ASN(156):2HD2(276)', 'A:ASN(156):ND2(274)', 'RECEPTOR'),
                ('A:CHT(1):O6(22)', 'A:CHT(1):HO6(23)', 'A:ASP(157):OD1(284)', 'LIGAND')
            ],
            'mol': <binana._structure.mol.Mol instance at 0x7feb20478518>
        }

    Args:
        ligand (binana._structure.mol.Mol): The ligand molecule to analyze.
        receptor (binana._structure.mol.Mol): The receptor molecule to analyze.
        dist_cutoff (float, optional): The distance cutoff. Defaults to
            HYDROGEN_BOND_DIST_CUTOFF or HALOGEN_BOND_DIST_CUTOFF.
        angle_cutoff (float, optional): The angle cutoff. Defaults to
            HYDROGEN_HALOGEN_BOND_ANGLE_CUTOFF.
        hydrogen_bond (boolean, optional): If True, calculates hydrogen bonds.
            Otherwise, calculates halogen bonds. Defaults to True.

    Returns:
        dict: Contains the atom tallies ("counts"), a binana._structure.mol.Mol
        object with the participating atoms ("mol"), and the labels to use in
        the log file ("labels").
    """

    hbonds = {}
    pdb_hbonds = Mol()
    hbonds_labels = []

    dist_cutoff = _set_default(
        dist_cutoff,
        HYDROGEN_BOND_DIST_CUTOFF if hydrogen_bond else HALOGEN_BOND_DIST_CUTOFF,
    )
    angle_cutoff = _set_default(angle_cutoff, HYDROGEN_HALOGEN_BOND_ANGLE_CUTOFF)

    # Check if hydrogen atoms added.
    lig_and_recep_have_hydrogens = ligand.has_hydrogens and receptor.has_hydrogens

    # Get all donor-acceptor pairs that are near each other.
    close_donors_acceptors = _get_potential_donors_acceptors(
        ligand, receptor, dist_cutoff, hydrogen_bond
    )

    acceptor_donor_atoms = {}
    bonds_organized_by_donor = {}
    idx = 0

    # Go through those pairs and find ones with complementary receptor/ligand
    # labels.
    for ligand_atom, receptor_atom, dist in close_donors_acceptors:
        # hbond_detected = False
        lig_atm_hbond_infs = ligand.is_hbond_donor_acceptor(ligand_atom, hydrogen_bond)
        recep_atm_hbond_infs = receptor.is_hbond_donor_acceptor(
            receptor_atom, hydrogen_bond
        )

        combos = _product(lig_atm_hbond_infs, recep_atm_hbond_infs)
        for lig_atm_hbond_inf, recep_atm_hbond_inf in combos:
            lig_donor_or_accept, lig_center_atom = lig_atm_hbond_inf
            recep_donor_or_accept, accept_center_atom = recep_atm_hbond_inf

            if lig_donor_or_accept == recep_donor_or_accept:
                # Both acceptors or both donors. Doesn't work.
                continue

            center_atom = (
                lig_center_atom
                if lig_donor_or_accept == "DONOR"
                else accept_center_atom
            )
            # center_atom_key = center_atom.string_id()

            # print(hydrogen_bond, center_atom.element)

            # Now that you've got the atoms, check the angles if appropriate.
            angle = None
            if lig_and_recep_have_hydrogens or not hydrogen_bond:
                # Hydrogens present and you're detecting hydrogen bonds, or
                # you're detecting halogen bonds (so hydrogens don't matter).

                angle = fabs(
                    180
                    - angle_between_three_points(
                        ligand_atom.coordinates,
                        center_atom.coordinates,
                        receptor_atom.coordinates,
                    )
                    * to_deg
                )

                if angle > angle_cutoff:
                    # Angle is too big.
                    continue

            donor_atom = (
                receptor_atom if lig_donor_or_accept == "ACCEPTOR" else ligand_atom
            )
            donor_key = donor_atom.string_id()
            acceptor_donor_atoms[donor_key] = donor_atom

            if donor_key not in bonds_organized_by_donor:
                bonds_organized_by_donor[donor_key] = []
            bonds_organized_by_donor[donor_key].append(
                (
                    idx,
                    lig_donor_or_accept,
                    receptor_atom,
                    ligand_atom,
                    center_atom,
                    dist,
                    angle,
                )
            )
            idx += 1

            # If you get here, it's identified a hydrogen bond. No need to keep
            # checking for a hydrogen bond between these two atoms.
            break

    if not lig_and_recep_have_hydrogens and hydrogen_bond:
        # If no hydrogens and trying to predict hydrogen bonds, further process
        # the detected bonds to account for additional angles. I found that
        # without this additional step, too many hydrogen bond donors are
        # detected.
        _remove_extra_noh_hydrogen_bonds(
            ligand, receptor, acceptor_donor_atoms, bonds_organized_by_donor
        )

    _collect_bonds(bonds_organized_by_donor, pdb_hbonds, hbonds, hbonds_labels)

    return {
        "counts": hbonds,
        "mol": pdb_hbonds,
        "labels": hbonds_labels,
    }


def get_hydrogen_bonds(ligand, receptor, dist_cutoff=None, angle_cutoff=None):
    """Identifies and counts the number of hydrogen bonds between the protein
    and ligand. Output is formatted like this::

        {
            'counts': {
                'HDONOR_RECEPTOR_SIDECHAIN_OTHER': 1,
                'HDONOR_LIGAND_SIDECHAIN_OTHER': 2
            },
            'labels': [
                ('A:CHT(1):N1(14)', 'A:CHT(1):H1(16)', 'A:ASP(157):OD2(285)', 'LIGAND'),
                ('A:CHT(1):O6(22)', 'A:ASN(156):2HD2(276)', 'A:ASN(156):ND2(274)', 'RECEPTOR'),
                ('A:CHT(1):O6(22)', 'A:CHT(1):HO6(23)', 'A:ASP(157):OD1(284)', 'LIGAND')
            ],
            'mol': <binana._structure.mol.Mol instance at 0x7feb20478518>
        }

    Args:
        ligand (binana._structure.mol.Mol): The ligand molecule to analyze.
        receptor (binana._structure.mol.Mol): The receptor molecule to analyze.
        dist_cutoff (float, optional): The distance cutoff. Defaults to
            HYDROGEN_BOND_DIST_CUTOFF.
        angle_cutoff (float, optional): The angle cutoff. Defaults to
            HYDROGEN_HALOGEN_BOND_ANGLE_CUTOFF.

    Returns:
        dict: Contains the atom tallies ("counts"), a binana._structure.mol.Mol
        object with the participating atoms ("mol"), and the labels to use in
        the log file ("labels").
    """

    # True means hydrogen bonds instead of halogen bonds.
    return _get_hydrogen_or_halogen_bonds(
        ligand, receptor, dist_cutoff, angle_cutoff, True
    )


def get_halogen_bonds(ligand, receptor, dist_cutoff=None, angle_cutoff=None):
    """Identifies and counts the number of halogen bonds between the protein
    and ligand. Output is formatted like this::

        {
            'counts': {
                'HDONOR_RECEPTOR_SIDECHAIN_OTHER': 1,
                'HDONOR_LIGAND_SIDECHAIN_OTHER': 2
            },
            'labels': [
                ('A:CHT(1):N1(14)', 'A:CHT(1):H1(16)', 'A:ASP(157):OD2(285)', 'LIGAND'),
                ('A:CHT(1):O6(22)', 'A:ASN(156):2HD2(276)', 'A:ASN(156):ND2(274)', 'RECEPTOR'),
                ('A:CHT(1):O6(22)', 'A:CHT(1):HO6(23)', 'A:ASP(157):OD1(284)', 'LIGAND')
            ],
            'mol': <binana._structure.mol.Mol instance at 0x7feb20478518>
        }

    Args:
        ligand (binana._structure.mol.Mol): The ligand molecule to analyze.
        receptor (binana._structure.mol.Mol): The receptor molecule to analyze.
        dist_cutoff (float, optional): The distance cutoff. Defaults to
            HALOGEN_BOND_DIST_CUTOFF.
        angle_cutoff (float, optional): The angle cutoff. Defaults to
            HYDROGEN_HALOGEN_BOND_ANGLE_CUTOFF.

    Returns:
        dict: Contains the atom tallies ("counts"), a binana._structure.mol.Mol
        object with the participating atoms ("mol"), and the labels to use in
        the log file ("labels").
    """

    # False means halogen bonds instead of hydrogen bonds.
    return _get_hydrogen_or_halogen_bonds(
        ligand, receptor, dist_cutoff, angle_cutoff, False
    )