# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

from binana._utils.shim import _set_default
from binana.interactions.default_params import ELECTROSTATIC_DIST_CUTOFF
import binana
from binana.load_ligand_receptor import _get_ligand_receptor_dists
from binana._utils.utils import hashtable_entry_add_one, list_alphebetize_and_combine

# Be sure to update the corresponding function in
# binana.interactions.__init__.py as well!


def get_electrostatic_energies(ligand, receptor, cutoff=None):
    """Calculates and tallies the electrostatic energies between receptor and
    ligand atoms that come within a given distance of each other. Output is
    formatted like this::

        {
            'counts': {
                'C_C': 49372.61585423234,
                'A_OA': -311243.9243779809
            }
        }

    Args:
        ligand (binana._structure.mol.Mol): The ligand molecule to analyze.
        receptor (binana._structure.mol.Mol): The receptor molecule to analyze.
        cutoff (float, optional): The distance cutoff. Defaults to
            ELECTROSTATIC_DIST_CUTOFF.

    Returns:
        dict: Contains the tallies ("counts") of the energies by atom-type
        pair.
    """

    cutoff = _set_default(cutoff, ELECTROSTATIC_DIST_CUTOFF)

    ligand_receptor_atom_type_pairs_electrostatic = {}
    # pdb_close_contacts = binana._structure.mol.Mol()
    # close_contacts_labels = []

    # Calculate the distances.
    ligand_receptor_dists = _get_ligand_receptor_dists(ligand, receptor, cutoff)

    # calculate electrostatic energies for all less than 4 A
    for ligand_atom, receptor_atom, dist in ligand_receptor_dists:
        # if dist < cutoff:
        # calculate electrostatic energies for all less than 4 A
        ligand_charge = ligand_atom.charge
        receptor_charge = receptor_atom.charge
        # to convert into J/mol # might be nice to double check this
        coulomb_energy = (ligand_charge * receptor_charge / dist) * 138.94238460104697e4
        list_ligand_atom = [ligand_atom.atom_type, receptor_atom.atom_type]
        hashtable_entry_add_one(
            ligand_receptor_atom_type_pairs_electrostatic,
            list_alphebetize_and_combine(list_ligand_atom),
            coulomb_energy,
        )

    return {
        "counts": ligand_receptor_atom_type_pairs_electrostatic,
    }
