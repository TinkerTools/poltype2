# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

from binana._utils.shim import _set_default
from binana.interactions.default_params import ACTIVE_SITE_FLEXIBILITY_DIST_CUTOFF
import binana
from binana.load_ligand_receptor import _get_ligand_receptor_dists
from binana._utils.utils import hashtable_entry_add_one, list_alphebetize_and_combine
from binana._structure.mol import Mol

# Be sure to update the corresponding function in
# binana.interactions.__init__.py as well!


def get_flexibility(ligand, receptor, cutoff=None):
    """Categorizes ligand-adjacent receptor atoms as belonging to a sidechain
    or backbone, as well as an alpha helix, beta sheet, or other secondary
    structure. Output is formatted like this::

        {

            'counts': {
                'SIDECHAIN_OTHER': 136,
                'SIDECHAIN_BETA': 72,
                'BACKBONE_OTHER': 7,
                'BACKBONE_BETA': 3,
                'SIDECHAIN_ALPHA': 18
            },
            'mols': {
                'alpha_helix': <binana._structure.mol.Mol instance at 0x7feb20438170>,
                'beta_sheet': <binana._structure.mol.Mol instance at 0x7feb204381b8>,
                'side_chain': <binana._structure.mol.Mol instance at 0x7feb20438368>,
                'other_2nd_structure': <binana._structure.mol.Mol instance at 0x7feb20438248>,
                'back_bone': <binana._structure.mol.Mol instance at 0x7feb20438320>
            }
        }

    Args:
        ligand (binana._structure.mol.Mol): The ligand molecule to analyze.
        receptor (binana._structure.mol.Mol): The receptor molecule to analyze.
        cutoff (float, optional): The distance cutoff. Defaults to
            ACTIVE_SITE_FLEXIBILITY_DIST_CUTOFF.

    Returns:
        dict: Contains the atom tallies ("counts"), as well as a list of
        binana._structure.mol.Mol objects ("mols"), each with the participating
        atoms that belong to alpha helixes, beta sheets, and other,
        respectively.
    """

    cutoff = _set_default(cutoff, ACTIVE_SITE_FLEXIBILITY_DIST_CUTOFF)

    active_site_flexibility = {}
    pdb_contacts_alpha_helix = Mol()
    pdb_contacts_beta_sheet = Mol()
    pdb_contacts_other_2nd_structure = Mol()
    pdb_back_bone = Mol()
    pdb_side_chain = Mol()

    # close_contacts_labels = []

    # Calculate the distances.
    ligand_receptor_dists = _get_ligand_receptor_dists(ligand, receptor, cutoff)

    # Now get statistics to judge active-site flexibility
    for ligand_atom, receptor_atom, dist in ligand_receptor_dists:
        # if dist < cutoff:
        # first can be sidechain or backbone, second back be alpha,
        # beta, or other, so six catagories
        flexibility_key = (
            receptor_atom.side_chain_or_backbone() + "_" + receptor_atom.structure
        )
        if receptor_atom.structure == "ALPHA":
            pdb_contacts_alpha_helix.add_new_atom(receptor_atom.copy_of())
        elif receptor_atom.structure == "BETA":
            pdb_contacts_beta_sheet.add_new_atom(receptor_atom.copy_of())
        elif receptor_atom.structure == "OTHER":
            pdb_contacts_other_2nd_structure.add_new_atom(receptor_atom.copy_of())

        if receptor_atom.side_chain_or_backbone() == "BACKBONE":
            pdb_back_bone.add_new_atom(receptor_atom.copy_of())
        elif receptor_atom.side_chain_or_backbone() == "SIDECHAIN":
            pdb_side_chain.add_new_atom(receptor_atom.copy_of())

        hashtable_entry_add_one(active_site_flexibility, flexibility_key)

    return {
        "counts": active_site_flexibility,
        "mols": {
            "alpha_helix": pdb_contacts_alpha_helix,
            "beta_sheet": pdb_contacts_beta_sheet,
            "other_2nd_structure": pdb_contacts_other_2nd_structure,
            "back_bone": pdb_back_bone,
            "side_chain": pdb_side_chain,
        },
        # "labels": close_contacts_labels,
    }
