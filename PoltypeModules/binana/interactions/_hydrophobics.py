# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

from binana._utils.shim import _set_default
from binana.interactions.default_params import HYDROPHOBIC_DIST_CUTOFF
import binana
from binana.load_ligand_receptor import _get_ligand_receptor_dists
from binana._utils.utils import hashtable_entry_add_one
from binana._structure.mol import Mol

# Be sure to update the corresponding function in
# binana.interactions.__init__.py as well!


def get_hydrophobics(ligand, receptor, cutoff=None):
    """Identifies and counts the number of hydrophobic (C-C) interactions
    between the protein and ligand. Output is formatted like this::

        {
            'counts': {
                'SIDECHAIN_OTHER': 43,
                'SIDECHAIN_BETA': 29,
                'BACKBONE_OTHER': 2
            },
            'labels': [
                ('A:CHT(1):C5(1)', 'A:TRP(43):CD2(30)'),
                ('A:CHT(1):C5(1)', 'A:TRP(43):CE2(32)'),
                ('A:CHT(1):C5(1)', 'A:TRP(43):CE3(33)')
            ],
            'mol': <binana._structure.mol.Mol instance at 0x7feb000acc68>
        }

    Args:
        ligand (binana._structure.mol.Mol): The ligand molecule to analyze.
        receptor (binana._structure.mol.Mol): The receptor molecule to analyze.
        cutoff (float, optional): The distance cutoff. Defaults to
            HYDROPHOBIC_DIST_CUTOFF.

    Returns:
        dict: Contains the atom tallies ("counts"), a binana._structure.mol.Mol
        object with the participating atoms ("mol"), and the labels to use in
        the log file ("labels").
    """

    cutoff = _set_default(cutoff, HYDROPHOBIC_DIST_CUTOFF)

    hydrophobics = {}
    pdb_hydrophobic = Mol()
    hydrophobic_labels = []

    # Calculate the distances.
    ligand_receptor_dists = _get_ligand_receptor_dists(ligand, receptor, cutoff, ["C"])

    # Now see if there's hydrophobic contacts (C-C contacts)
    for ligand_atom, receptor_atom, dist in ligand_receptor_dists:
        hydrophobic_key = (
            receptor_atom.side_chain_or_backbone() + "_" + receptor_atom.structure
        )
        pdb_hydrophobic.add_new_atom(ligand_atom.copy_of())
        pdb_hydrophobic.add_new_atom(receptor_atom.copy_of())

        hashtable_entry_add_one(hydrophobics, hydrophobic_key)

        hydrophobic_labels.append(
            (
                ligand_atom.string_id(),
                receptor_atom.string_id(),
                {"distance": dist},
            )
        )

    return {
        "counts": hydrophobics,
        "mol": pdb_hydrophobic,
        "labels": hydrophobic_labels,
    }
