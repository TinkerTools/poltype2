# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

import binana
# from binana.load_ligand_receptor import _get_ligand_receptor_dists
from binana._utils.utils import hashtable_entry_add_one # , list_alphebetize_and_combine

# Be sure to update the corresponding function in
# binana.interactions.__init__.py as well!


def get_ligand_atom_types(ligand):
    """Tallies the ligand atoms by atom type. Output is formatted like this::

        {
            'counts': {
                'A': 8,
                'C': 5,
                'HD': 3,
                'OA': 5,
                'N': 2
            }
        }

    Args:
        ligand (binana._structure.mol.Mol): The ligand molecule to analyze.

    Returns:
        dict: Contains the atom tallies ("counts").
    """

    ligand_atom_types = {}

    # Get the total number of each atom type in the ligand
    for ligand_atom_index in ligand.all_atoms.keys():
        ligand_atom = ligand.all_atoms[ligand_atom_index]
        hashtable_entry_add_one(ligand_atom_types, ligand_atom.atom_type)

    return {
        "counts": ligand_atom_types,
    }
