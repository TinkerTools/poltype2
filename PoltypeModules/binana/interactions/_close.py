# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

from binana._utils.shim import _set_default
from binana.interactions.default_params import CLOSE_CONTACTS_DIST2_CUTOFF
import binana
from binana.load_ligand_receptor import _get_ligand_receptor_dists
from binana._utils.utils import hashtable_entry_add_one, list_alphebetize_and_combine
from binana._structure.mol import Mol

# Be sure to update the corresponding function in
# binana.interactions.__init__.py as well!


def get_close(ligand, receptor, cutoff=None):
    """Identifies and counts the number of close protein/ligand contacts.
    Output is formatted like this::

        {
            'counts': {
                'C_C': 5,
                'A_OA': 29,
                'HD_NA': 1,
                'HD_N': 6
            },
            'labels': [
                ('A:CHT(1):C5(1)', 'A:TRP(43):CD2(30)'),
                ('A:CHT(1):C5(1)', 'A:TRP(43):CE2(32)'),
                ('A:CHT(1):C5(1)', 'A:TRP(43):CE3(33)')
            ],
            'mol': <binana._structure.mol.Mol instance at 0x7feb203ce3f8>
        }

    Args:
        ligand (binana._structure.mol.Mol): The ligand molecule to analyze.
        receptor (binana._structure.mol.Mol): The receptor molecule to analyze.
        cutoff (float, optional): The distance cutoff. Defaults to
            CLOSE_CONTACTS_DIST2_CUTOFF.

    Returns:
        dict: Contains the atom tallies ("counts"), a binana._structure.mol.Mol
        object with the participating atoms ("mol"), and the labels to use in
        the log file ("labels").
    """

    cutoff = _set_default(cutoff, CLOSE_CONTACTS_DIST2_CUTOFF)

    ligand_receptor_atom_type_pairs_close = {}
    pdb_close_contacts = Mol()
    close_contacts_labels = []

    # Calculate the distances.
    ligand_receptor_dists = _get_ligand_receptor_dists(ligand, receptor, cutoff)

    for i, (ligand_atom, receptor_atom, dist) in enumerate(ligand_receptor_dists):
        # if dist < cutoff:
        # less than 4 A
        list_ligand_atom = [ligand_atom.atom_type, receptor_atom.atom_type]
        hashtable_entry_add_one(
            ligand_receptor_atom_type_pairs_close,
            list_alphebetize_and_combine(list_ligand_atom),
        )
        pdb_close_contacts.add_new_atom(ligand_atom.copy_of())
        pdb_close_contacts.add_new_atom(receptor_atom.copy_of())

        close_contacts_labels.append(
            (
                ligand_atom.string_id(),
                receptor_atom.string_id(),
                {"distance": dist},
            )
        )

    return {
        "counts": ligand_receptor_atom_type_pairs_close,
        "mol": pdb_close_contacts,
        "labels": close_contacts_labels,
    }
