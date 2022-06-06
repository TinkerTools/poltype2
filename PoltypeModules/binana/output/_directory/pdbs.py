# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

def output_dir_pdbs(
    pdb_closest_contacts,
    parameters,
    pdb_close_contacts,
    pdb_contacts_alpha_helix,
    pdb_contacts_beta_sheet,
    pdb_contacts_other_2nd_structure,
    pdb_back_bone,
    pdb_side_chain,
    pdb_hydrophobic,
    pdb_hbonds,
    pdb_halbonds,
    pdb_pistack,
    pdb_pi_T,
    pdb_pi_cat,
    pdb_salt_bridges,
    pdb_metal_coordinations,
    ligand,
    receptor,
):

    # so an output directory has been specified. Write the pdb files
    # out separately

    pdb_closest_contacts.save_pdb(
        parameters.params["output_dir"] + "/close_contacts.pdb"
    )
    pdb_close_contacts.save_pdb(parameters.params["output_dir"] + "/contacts.pdb")
    pdb_contacts_alpha_helix.save_pdb(
        parameters.params["output_dir"] + "/contacts_alpha_helix.pdb"
    )
    pdb_contacts_beta_sheet.save_pdb(
        parameters.params["output_dir"] + "/contacts_beta_sheet.pdb"
    )
    pdb_contacts_other_2nd_structure.save_pdb(
        parameters.params["output_dir"] + "/contacts_other_secondary_structure.pdb"
    )
    pdb_back_bone.save_pdb(parameters.params["output_dir"] + "/back_bone.pdb")
    pdb_side_chain.save_pdb(parameters.params["output_dir"] + "/side_chain.pdb")
    pdb_hydrophobic.save_pdb(parameters.params["output_dir"] + "/hydrophobic.pdb")
    pdb_hbonds.save_pdb(parameters.params["output_dir"] + "/hydrogen_bonds.pdb")
    pdb_halbonds.save_pdb(parameters.params["output_dir"] + "/halogen_bonds.pdb")
    pdb_pistack.save_pdb(parameters.params["output_dir"] + "/pi_pi_stacking.pdb")
    pdb_pi_T.save_pdb(parameters.params["output_dir"] + "/T_stacking.pdb")
    pdb_pi_cat.save_pdb(parameters.params["output_dir"] + "/cat_pi.pdb")
    pdb_salt_bridges.save_pdb(parameters.params["output_dir"] + "/salt_bridges.pdb")
    pdb_metal_coordinations.save_pdb(parameters.params["output_dir"] + "/metal_coordinations.pdb")
    ligand.save_pdb(parameters.params["output_dir"] + "/ligand.pdb")
    receptor.save_pdb(parameters.params["output_dir"] + "/receptor.pdb")
