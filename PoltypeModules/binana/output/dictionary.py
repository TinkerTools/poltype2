# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

import re

# __pragma__ ('skip')
# Python
import json
from os import sep
from os.path import basename, dirname

# __pragma__ ('noskip')

"""?
# Transcrypt
import binana
import binana._utils.shim as json
from binana._utils.shim import sep
from binana._utils.shim import basename
from binana._utils.shim import dirname
?"""

# __pragma__ ('skip')
# Python, just alias open
_openFile = open
# __pragma__ ('noskip')

"""?
# Transcrypt
from binana._utils.shim import OpenFile
_openFile = OpenFile
?"""


def _atom_details_to_dict(details):
    return {
        "chain": details[0].strip(),
        "resID": int(details[2]),
        "resName": details[1],
        "atomName": details[3],
        "atomIndex": int(details[4]),
    }


# takes care of pairwise interactions
# interaction_name: string
# interaction_labels: list containind raw data to be parsed into the dictionary
# returns list
def _get_close_atom_list(interaction_labels):
    interaction_list = []
    for i, atom_pairs in enumerate(interaction_labels):
        # add new dictionary to interaction list
        # json_output[interaction_name].append({})
        interaction_list.append({})
        # parse atom_pairs
        ligand_atom_details = re.split(r"[():]", atom_pairs[0])
        receptor_atom_details = re.split(r"[():]", atom_pairs[1])
        # remove whitespace
        for detail in ligand_atom_details:
            if detail == "":
                ligand_atom_details.remove(detail)
        for detail in receptor_atom_details:
            if detail == "":
                receptor_atom_details.remove(detail)

        # add each detail to the appropriate key in ligand_atoms and receptor_atoms
        # json_output[interaction_name][i] = {

        interaction_list[i] = {
            "ligandAtoms": [_atom_details_to_dict(ligand_atom_details)],
            "receptorAtoms": [_atom_details_to_dict(receptor_atom_details)],
            "metrics": atom_pairs[2],
        }

    return interaction_list


def _collect_hydrogen_halogen_bonds(hydrogen_bonds, json_output, hydrogen_bond=True):
    dict_key = "hydrogenBonds" if hydrogen_bond else "halogenBonds"

    # hydrogen bonds
    for i, atom_pairs in enumerate(hydrogen_bonds):
        # add new dictionary to hydrogen_bonds list
        json_output[dict_key].append({})
        # parse atom_trios
        ligand_and_receptor = [
            re.split(r"[():]", atom_pairs[0]),  # Ligand
            re.split(r"[():]", atom_pairs[1]),  # hydrogen (Ligand or receptor)
            re.split(r"[():]", atom_pairs[2]),  # receptor
        ]

        # put atoms in appropriate group
        ligand_atom_details = [ligand_and_receptor[0]]
        receptor_atom_details = [ligand_and_receptor[2]]
        if atom_pairs[3] == "RECEPTOR":
            receptor_atom_details.append(ligand_and_receptor[1])
        else:
            ligand_atom_details.append(ligand_and_receptor[1])

        # for atom in ligand_and_receptor:
        #     if len(atom[0]) > 1: # ligand ("ATP")
        #         ligand_atom_details.append(atom)
        #     else: # receptor ("A")
        #         receptor_atom_details.append(atom)

        # remove whitespace
        for atom in ligand_atom_details:
            for detail in atom:
                if detail == "":
                    atom.remove(detail)
        for atom in receptor_atom_details:
            for detail in atom:
                if detail == "":
                    atom.remove(detail)
        # add each detail to the appropriate key in ligand_atoms and receptor_atoms
        # add "ligandAtoms" and "receptorAtoms" keys to dictionary
        json_output[dict_key][i] = {
            "ligandAtoms": [],
            "receptorAtoms": [],
            "metrics": atom_pairs[4],
        }
        for detail in ligand_atom_details:
            json_output[dict_key][i]["ligandAtoms"].append(
                _atom_details_to_dict(detail)
            )
        for detail in receptor_atom_details:
            json_output[dict_key][i]["receptorAtoms"].append(
                _atom_details_to_dict(detail)
            )


def _collect_pi_pi(pi_stacking_interactions, json_output):
    # pi-pi stacking interactions
    for i, atom_pair in enumerate(pi_stacking_interactions):
        # add new dictionary to pi_stacking list
        json_output["piPiStackingInteractions"].append({})
        # parse atom_pairs into individual atoms
        individual_ligand_atoms = atom_pair[0].split("/")
        individual_receptor_atoms = atom_pair[1].split("/")
        # parse individual atoms into details
        individual_ligand_atoms_details = [
            re.split(r"[\[\]():]", atom)
            for atom in individual_ligand_atoms
            if atom != ""
        ]

        individual_receptor_atoms_details = [
            re.split(r"[\[\]():]", atom)
            for atom in individual_receptor_atoms
            if atom != ""
        ]

        # remove whitespace
        for detail_list in individual_ligand_atoms_details:
            for detail in detail_list:
                if detail == "":
                    detail_list.remove(detail)
        for detail_list in individual_receptor_atoms_details:
            for detail in detail_list:
                if detail == "":
                    detail_list.remove(detail)
        # add each detail to the appropriate key in ligand_atoms and receptor_atoms
        # add "ligandAtoms" and "receptorAtoms" keys to dictionary
        json_output["piPiStackingInteractions"][i] = {
            "ligandAtoms": [],
            "receptorAtoms": [],
            "metrics": atom_pair[2]
        }
        for detail in individual_ligand_atoms_details:
            json_output["piPiStackingInteractions"][i]["ligandAtoms"].append(
                _atom_details_to_dict(detail)
            )
        for detail in individual_receptor_atoms_details:
            json_output["piPiStackingInteractions"][i]["receptorAtoms"].append(
                _atom_details_to_dict(detail)
            )


def _collect_t_stacking(t_stacking_interactions, json_output):
    # t stacking interactions
    for i, atom_pair in enumerate(t_stacking_interactions):
        # add new dictionary to t_stacking list
        json_output["tStackingInteractions"].append({})
        # parse atom_pairs into individual atoms
        individual_ligand_atoms = atom_pair[0].split("/")
        individual_receptor_atoms = atom_pair[1].split("/")
        # parse individual atoms into details
        individual_ligand_atoms_details = [
            re.split(r"[\[\]():]", atom)
            for atom in individual_ligand_atoms
            if atom != ""
        ]

        individual_receptor_atoms_details = [
            re.split(r"[\[\]():]", atom)
            for atom in individual_receptor_atoms
            if atom != ""
        ]

        # remove whitespace
        for detail_list in individual_ligand_atoms_details:
            for detail in detail_list:
                if detail == "":
                    detail_list.remove(detail)
        for detail_list in individual_receptor_atoms_details:
            for detail in detail_list:
                if detail == "":
                    detail_list.remove(detail)
        # add each detail to the appropriate key in ligand_atoms and receptor_atoms
        # add "ligandAtoms" and "receptorAtoms" keys to dictionary
        json_output["tStackingInteractions"][i] = {
            "ligandAtoms": [],
            "receptorAtoms": [],
            "metrics": atom_pair[2]
        }

        for detail in individual_ligand_atoms_details:
            json_output["tStackingInteractions"][i]["ligandAtoms"].append(
                _atom_details_to_dict(detail)
            )
        for detail in individual_receptor_atoms_details:
            json_output["tStackingInteractions"][i]["receptorAtoms"].append(
                _atom_details_to_dict(detail)
            )


def _collect_cat_pi(cat_pi_interactions, json_output):
    # cat-pi stacking interactions
    for i, atom_pair in enumerate(cat_pi_interactions):
        # add new dictionary to cation-pi_stacking list
        json_output["cationPiInteractions"].append({})
        # parse atom_pairs into individual atoms
        individual_ligand_atoms = atom_pair[0].split("/")
        individual_receptor_atoms = atom_pair[1].split("/")
        # parse individual atoms into details
        individual_ligand_atoms_details = [
            re.split(r"[\[\]():]", atom)
            for atom in individual_ligand_atoms
            if atom != ""
        ]

        individual_receptor_atoms_details = [
            re.split(r"[\[\]():]", atom)
            for atom in individual_receptor_atoms
            if atom != ""
        ]

        # remove whitespace
        for detail_list in individual_ligand_atoms_details:
            for detail in detail_list:
                if detail == "":
                    detail_list.remove(detail)
        for detail_list in individual_receptor_atoms_details:
            for detail in detail_list:
                if detail == "":
                    detail_list.remove(detail)
        # add each detail to the appropriate key in ligand_atoms and receptor_atoms
        # add "ligandAtoms" and "receptorAtoms" keys to dictionary
        json_output["cationPiInteractions"][i] = {
            "ligandAtoms": [],
            "receptorAtoms": [],
            "metrics": atom_pair[2],
        }
        for detail in individual_ligand_atoms_details:
            json_output["cationPiInteractions"][i]["ligandAtoms"].append(
                _atom_details_to_dict(detail)
            )
        for detail in individual_receptor_atoms_details:
            json_output["cationPiInteractions"][i]["receptorAtoms"].append(
                _atom_details_to_dict(detail)
            )


def _collect_salt_bridge(salt_bridge_interactions, json_output):
    # salt bridge interactions
    for i, atom_pair in enumerate(salt_bridge_interactions):
        # add new dictionary to salt_bridges list
        json_output["saltBridges"].append({})
        # parse atom_pairs into individual atoms
        individual_ligand_atoms = atom_pair[0].split("/")
        individual_receptor_atoms = atom_pair[1].split("/")
        # parse individual atoms into details
        individual_ligand_atoms_details = [
            re.split(r"[\[\]():]", atom)
            for atom in individual_ligand_atoms
            if atom != ""
        ]

        individual_receptor_atoms_details = [
            re.split(r"[\[\]():]", atom)
            for atom in individual_receptor_atoms
            if atom != ""
        ]

        # remove whitespace
        for detail_list in individual_ligand_atoms_details:
            for detail in detail_list:
                if detail == "":
                    detail_list.remove(detail)
        for detail_list in individual_receptor_atoms_details:
            for detail in detail_list:
                if detail == "":
                    detail_list.remove(detail)
        # add each detail to the appropriate key in ligand_atoms and receptor_atoms
        # add "ligandAtoms" and "receptorAtoms" keys to dictionary
        json_output["saltBridges"][i] = {"ligandAtoms": [], "receptorAtoms": [], "metrics": atom_pair[2]}
        for detail in individual_ligand_atoms_details:
            json_output["saltBridges"][i]["ligandAtoms"].append(
                _atom_details_to_dict(detail)
            )
        for detail in individual_receptor_atoms_details:
            json_output["saltBridges"][i]["receptorAtoms"].append(
                _atom_details_to_dict(detail)
            )


def _collect_metal_coordinations(metal_coordinations_interactions, json_output):
    # metal coordination interactions
    for i, metal_coord_atoms in enumerate(metal_coordinations_interactions):
        # add new dictionary to metalCoordinations list
        json_output["metalCoordinations"].append({})

        # parse atom_pairs into individual atoms
        ligand_atom = metal_coord_atoms[0]
        ligand_atom_details = re.split(r"[\[\]():]", ligand_atom)
        ligand_atom_details = [d for d in ligand_atom_details if d != ""]

        receptor_atom = metal_coord_atoms[1]
        receptor_atom_details = re.split(r"[\[\]():]", receptor_atom)
        receptor_atom_details = [d for d in receptor_atom_details if d != ""]

        json_output["metalCoordinations"][i] = {
            "ligandAtoms": [_atom_details_to_dict(ligand_atom_details)],
            "receptorAtoms": [_atom_details_to_dict(receptor_atom_details)],
            "metrics": metal_coord_atoms[2]
        }


# json output
def collect(
    closest=None,
    close=None,
    hydrophobics=None,
    hydrogen_bonds=None,
    halogen_bonds=None,
    salt_bridges=None,
    metal_coordinations=None,
    pi_pi=None,
    cat_pi=None,
    electrostatic_energies=None,
    active_site_flexibility=None,
    ligand_atom_types=None,
    ligand_rotatable_bonds=None,
):
    """Collects all the characterized interactions between the protein and
    ligand into one dict object, suitable for conversion to JSON.

    Args:
        closest (dict, optional): A dictionary containing information about the
            closest protein/ligand interactions. Defaults to None.
        close (dict, optional): A dictionary containing information about the
            close protein/ligand interactions. Defaults to None.
        hydrophobics (dict, optional): A dictionary containing information
            about the hydrophobic protein/ligand interactions. Defaults to
            None.
        hydrogen_bonds (dict, optional): A dictionary containing information
            about the hydrogen bonds between the protein and ligand. Defaults
            to None.
        halogen_bonds (dict, optional): A dictionary containing information
            about the halogen bonds between the protein and ligand. Defaults
            to None.
        salt_bridges (dict, optional): A dictionary containing information
            about the salt-bridges protein/ligand interactions. Defaults to
            None.
        metal_coordinations (dict, optional): A dictionary containing 
            information about the metal-coordination protein/ligand 
            interactions. Defaults to None.
        pi_pi (dict, optional): A dictionary containing information about the
            pi-pi (stacking and T-shaped) protein/ligand interactions. Defaults
            to None.
        cat_pi (dict, optional): A dictionary containing information about the
            pi-cation protein/ligand interactions. Defaults to None.
        electrostatic_energies (dict, optional): A dictionary containing
            information about the electrostatic energies between protein and
            ligand atoms. Defaults to None.
        active_site_flexibility (dict, optional): A dictionary containing
            information about the flexibility of ligand-adjacent protein atoms. Defaults to None.
        ligand_atom_types (dict, optional): A dictionary containing information
            about the ligand atom types. Defaults to None.
        ligand_rotatable_bonds (int, optional): The number of ligand rotatable
            bonds. Defaults to None.

    Returns:
        dict: A dictionary describing all the detected interactions, suitable
        for conversion to JSON.
    """

    json_output = {}

    # populate the lists of proximity metrics. use helper function to populate
    # lists for pairwise interactions
    if closest is not None:
        json_output["closestContacts"] = _get_close_atom_list(closest["labels"])
    if close is not None:
        json_output["closeContacts"] = _get_close_atom_list(close["labels"])
    if hydrophobics is not None:
        json_output["hydrophobicContacts"] = _get_close_atom_list(
            hydrophobics["labels"]
        )

    # Add in the other metrics that are more difficult to calculate.
    if hydrogen_bonds is not None:
        json_output["hydrogenBonds"] = []
        _collect_hydrogen_halogen_bonds(hydrogen_bonds["labels"], json_output, True)
    if halogen_bonds is not None:
        json_output["halogenBonds"] = []
        _collect_hydrogen_halogen_bonds(halogen_bonds["labels"], json_output, False)
    if pi_pi is not None:
        json_output["piPiStackingInteractions"] = []
        _collect_pi_pi(pi_pi["labels"]["pi_stacking"], json_output)
        json_output["tStackingInteractions"] = []
        _collect_t_stacking(pi_pi["labels"]["T_stacking"], json_output)
    if cat_pi is not None:
        json_output["cationPiInteractions"] = []
        _collect_cat_pi(cat_pi["labels"], json_output)
    if salt_bridges is not None:
        json_output["saltBridges"] = []
        _collect_salt_bridge(salt_bridges["labels"], json_output)
    if metal_coordinations is not None:
        json_output["metalCoordinations"] = []
        _collect_metal_coordinations(metal_coordinations["labels"], json_output)

    # For flexibility and electrostatics, just return the counts
    if active_site_flexibility is not None:
        json_output["activeSiteFlexibility"] = active_site_flexibility["counts"]
    if electrostatic_energies is not None:
        json_output["electrostaticEnergies"] = electrostatic_energies["counts"]
    if ligand_atom_types is not None:
        json_output["ligandAtomTypes"] = ligand_atom_types["counts"]

    if ligand_rotatable_bonds is not None:
        json_output["ligandRotatableBonds"] = ligand_rotatable_bonds

    return json_output


def collect_all(all_interactions):
    """Collects all the characterized interactions between the protein and
    ligand into one dict object, suitable for conversion to JSON. This function
    simply unpacks the contents of `all_interactions` and passes them to
    :py:func:`~binana.output.dictionary.collect`.

    Args:
        all_interactions (dict): A single dictionary containing information
            about all the protein/ligand interactions. The output of
            :py:func:`~binana.interactions.get_all_interactions`

    Returns:
        dict: A dictionary describing all the detected interactions, suitable
        for conversion to JSON.
    """

    # __pragma__ ('skip')
    return collect(
        closest=all_interactions["closest"],
        close=all_interactions["close"],
        hydrophobics=all_interactions["hydrophobics"],
        hydrogen_bonds=all_interactions["hydrogen_bonds"],
        halogen_bonds=all_interactions["halogen_bonds"],
        salt_bridges=all_interactions["salt_bridges"],
        metal_coordinations=all_interactions["metal_coordinations"],
        pi_pi=all_interactions["pi_pi"],
        cat_pi=all_interactions["cat_pi"],
        electrostatic_energies=all_interactions["electrostatic_energies"],
        active_site_flexibility=all_interactions["active_site_flexibility"],
        ligand_atom_types=all_interactions["ligand_atom_types"],
        ligand_rotatable_bonds=all_interactions["ligand_rotatable_bonds"],
    )
    # __pragma__ ('noskip')

    """?
    # Odd that transcrypt requires it this other way.
    return collect(
        all_interactions["closest"],
        all_interactions["close"],
        all_interactions["hydrophobics"],
        all_interactions["hydrogen_bonds"],
        all_interactions["halogen_bonds"],
        all_interactions["salt_bridges"],
        all_interactions["metal_coordinations"],
        all_interactions["pi_pi"],
        all_interactions["cat_pi"],
        all_interactions["electrostatic_energies"],
        all_interactions["active_site_flexibility"],
        all_interactions["ligand_atom_types"],
        all_interactions["ligand_rotatable_bonds"],
    )
    ?"""