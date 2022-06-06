# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

from binana._utils.shim import _set_default
from binana.interactions.default_params import CLOSE_CONTACTS_DIST1_CUTOFF
import binana
from binana.load_ligand_receptor import _get_ligand_receptor_dists
from binana._utils.utils import hashtable_entry_add_one, list_alphebetize_and_combine
from binana._structure.mol import Mol

# Be sure to update the corresponding function in
# binana.interactions.__init__.py as well!


def get_closest(ligand, receptor, cutoff=None):
    """Identifies and counts the number of closest (very close) protein/ligand
    contacts. Output is formatted like this::

        {
            'counts': {
                'HD_OA': 8,
                'A_OA': 3
            },
            'labels': [
                ('A:CHT(1):C9(7)', 'A:TRP(205):CB(467)'),
                ('A:CHT(1):O2(8)', 'A:TRP(205):CG(468)'),
            'mol': <binana._structure.mol.Mol instance at 0x7feb20290908>
        }

    Args:
        ligand (binana._structure.mol.Mol): The ligand molecule to analyze.
        receptor (binana._structure.mol.Mol): The receptor molecule to analyze.
        cutoff (float, optional): The distance cutoff. Defaults to
            CLOSE_CONTACTS_DIST1_CUTOFF.

    Returns:
        dict: Contains the atom tallies ("counts"), a binana._structure.mol.Mol
        object with the participating atoms ("mol"), and the labels to use in
        the log file ("labels").
    """

    cutoff = _set_default(cutoff, CLOSE_CONTACTS_DIST1_CUTOFF)

    ligand_receptor_atom_type_pairs_closest = {}
    pdb_closest_contacts = Mol()
    closest_contacts_labels = []

    # Calculate the distances.
    ligand_receptor_dists = _get_ligand_receptor_dists(ligand, receptor, cutoff)

    # Identify closest contacts
    for ligand_atom, receptor_atom, dist in ligand_receptor_dists:
        # if dist < cutoff:
        # less than 2.5 A
        list_ligand_atom = [ligand_atom.atom_type, receptor_atom.atom_type]
        hashtable_entry_add_one(
            ligand_receptor_atom_type_pairs_closest,
            list_alphebetize_and_combine(list_ligand_atom),
        )
        pdb_closest_contacts.add_new_atom(ligand_atom.copy_of())
        pdb_closest_contacts.add_new_atom(receptor_atom.copy_of())

        closest_contacts_labels.append(
            (
                ligand_atom.string_id(),
                receptor_atom.string_id(),
                {"distance": dist},
            )
        )

    return {
        "counts": ligand_receptor_atom_type_pairs_closest,
        "mol": pdb_closest_contacts,
        "labels": closest_contacts_labels,
    }
