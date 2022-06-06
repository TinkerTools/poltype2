# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

from binana._utils.shim import _set_default

# __pragma__ ('skip')
# Python
from textwrap import wrap as _wrap

# __pragma__ ('noskip')

"""?
# Transcrypt
import binana
import binana._utils
from binana._utils.shim import wrap as _wrap
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


def write(
    ligand,
    receptor,
    closest=None,
    close=None,
    hydrophobics=None,
    hydrogen_bonds=None,
    halogen_bonds=None,
    salt_bridges=None,
    metal_coordinations=None,
    pi_pi=None,
    cat_pi=None,
    active_site_flexibility=None,
    log_output=None,
    as_str=None,
    pdb_filename=None,
):
    """Writes a single PDB file containing the ligand, receptor, and atoms that
    participate in various interactions (with distinct resnames).

    Args:
        ligand (binana._structure.mol.Mol): The ligand molecule object.
        receptor (binana._structure.mol.Mol): The receptor molecule object.
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
        active_site_flexibility (dict, optional): A dictionary containing
            information about the flexibility of ligand-adjacent protein atoms.
            Defaults to None.
        log_output (str, optional): The log text, returned from
            :py:func:`~binana.output.log.collect`. Defaults to ``""``.
        as_str (bool, optional): Whether to save the file to the disk (or fake
            disk in case of JavaScript), or to return the contents as a string.
            Defaults to False.
        pdb_filename (str, optional): The name of the file where the pdb should
            be saved, assuming as_str is False. Defaults to "results.pdb".

    Returns:
        str: The contents of the PDB file if ``as_str`` is ``True``. Otherwise,
        ``""``.
    """

    # so it's writing to a single file.

    log_output = _set_default(log_output, "")
    as_str = _set_default(as_str, False)
    pdb_filename = _set_default(pdb_filename, "results.pdb")

    # first, make an explaination.

    explain = (
        'The residue named "CCN" contains the closest contacts between the protein and receptor. '
        + '"CON" indicates close contacts. '
        + '"ALP", "BET", and "OTH" indicate receptor contacts whose respective protein residues have the alpha-helix, beta-sheet, or "other" secondary structure. '
        + '"BAC" and "SID" indicate receptor contacts that are part of the protein backbone and sidechain, respectively. '
        + '"HYD" indicates hydrophobic contacts between the protein and ligand. '
        + '"HBN" indicates hydrogen bonds. "HAL" indicates halogen bonds. "SAL" indicates salt bridges. '
        + '"PIS" indicates pi-pi stacking interactions, "PIT" indicates T-stacking interactions, and "PIC" indicates cation-pi interactions. '
        + '"MTL" indicates metal-coordination interactions. '
        + 'Protein residue names are unchanged, but the ligand residue is now named "LIG".'
    )

    # Old text:
    # explain = (
    #     'The residue named "CCN" illustrates close contacts where the protein and ligand atoms come within '
    #     + str(parameters.params["close_contacts_dist1_cutoff"])
    #     + ' of each other. "CON" illustrates close contacts where the protein and ligand atoms come within '
    #     + str(parameters.params["close_contacts_dist2_cutoff"])
    #     + ' of each other. "ALP", "BET", and "OTH" illustrates receptor contacts whose respective protein residues have the alpha-helix, beta-sheet, or "other" secondary structure. "BAC" and "SID" illustrate receptor contacts that are part of the protein backbone and sidechain, respectively. "HYD" illustrates hydrophobic contacts between the protein and ligand. "HBN" illustrates hydrogen bonds. "SAL" illustrates salt bridges. "PIS" illustrates pi-pi stacking interactions, "PIT" illustrates T-stacking interactions, and "PIC" illustrates cation-pi interactions. Protein residue names are unchanged, but the ligand residue is now named "LIG".'
    # )

    log_output = log_output + "REMARK\n"

    lines = _wrap(explain, 71)
    for line in lines:
        log_output = log_output + "REMARK " + line + "\n"

    log_output = log_output + "REMARK\n"

    ligand.set_resname("LIG")
    log_output = (
        log_output
        + receptor.save_pdb_string()
        + "TER\n"
        + ligand.save_pdb_string()
        + "TER\n"
    )

    if closest is not None:
        closest["mol"].set_resname("CCN")
        log_output = log_output + closest["mol"].save_pdb_string() + "TER\n"

    if close is not None:
        close["mol"].set_resname("CON")
        log_output = log_output + close["mol"].save_pdb_string() + "TER\n"

    if active_site_flexibility is not None:
        active_site_flexibility["mols"]["alpha_helix"].set_resname("ALP")
        active_site_flexibility["mols"]["beta_sheet"].set_resname("BET")
        active_site_flexibility["mols"]["other_2nd_structure"].set_resname("OTH")
        active_site_flexibility["mols"]["back_bone"].set_resname("BAC")
        active_site_flexibility["mols"]["side_chain"].set_resname("SID")
        log_output = (
            log_output
            + active_site_flexibility["mols"]["alpha_helix"].save_pdb_string()
            + "TER\n"
            + active_site_flexibility["mols"]["beta_sheet"].save_pdb_string()
            + "TER\n"
            + active_site_flexibility["mols"]["other_2nd_structure"].save_pdb_string()
            + "TER\n"
            + active_site_flexibility["mols"]["back_bone"].save_pdb_string()
            + "TER\n"
            + active_site_flexibility["mols"]["side_chain"].save_pdb_string()
            + "TER\n"
        )

    if hydrophobics is not None:
        hydrophobics["mol"].set_resname("HYD")
        log_output = log_output + hydrophobics["mol"].save_pdb_string() + "TER\n"

    if hydrogen_bonds is not None:
        hydrogen_bonds["mol"].set_resname("HBN")
        log_output = log_output + hydrogen_bonds["mol"].save_pdb_string() + "TER\n"

    if halogen_bonds is not None:
        halogen_bonds["mol"].set_resname("HAL")
        log_output = log_output + halogen_bonds["mol"].save_pdb_string() + "TER\n"

    if pi_pi is not None:
        pi_pi["mols"]["pi_stacking"].set_resname("PIS")
        pi_pi["mols"]["T_stacking"].set_resname("PIT")
        log_output = (
            log_output
            + pi_pi["mols"]["pi_stacking"].save_pdb_string()
            + "TER\n"
            + pi_pi["mols"]["T_stacking"].save_pdb_string()
            + "TER\n"
        )

    if cat_pi is not None:
        cat_pi["mol"].set_resname("PIC")
        log_output = log_output + cat_pi["mol"].save_pdb_string() + "TER\n"

    if salt_bridges is not None:
        salt_bridges["mol"].set_resname("SAL")
        log_output = log_output + salt_bridges["mol"].save_pdb_string() + "TER\n"

    if metal_coordinations is not None:
        metal_coordinations["mol"].set_resname("MTL")
        log_output = log_output + metal_coordinations["mol"].save_pdb_string() + "TER\n"

    if as_str:
        return log_output

    f = _openFile(pdb_filename, "w")
    f.write(log_output)
    f.close()

    return ""


def write_all(
    ligand,
    receptor,
    all_interactions,
    log_output=None,
    as_str=None,
    pdb_filename=None,
):
    """Writes a single PDB file containing the ligand, receptor, and atoms that
    participate in various interactions (with distinct resnames). This function
    simply unpacks the contents of `all_interactions` and passes them to
    :py:func:`~binana.output.pdb_file.write`.

    Args:
        ligand (binana._structure.mol.Mol): The ligand molecule object.
        receptor (binana._structure.mol.Mol): The receptor molecule object.
        all_interactions (dict): A single dictionary containing information
            about all the protein/ligand interactions. The output of
            :py:func:`~binana.interactions.get_all_interactions`
        log_output (str, optional): The log text, returned from
            :py:func:`~binana.output.log.collect`. Defaults to ``""``.
        as_str (bool, optional): Whether to save the file to the disk (or fake
            disk in case of JavaScript), or to return the contents as a string.
            Defaults to False.
        pdb_filename (str, optional): The name of the file where the pdb should
            be saved, assuming as_str is False. Defaults to "results.pdb".

    Returns:
        str: The contents of the PDB file if ``as_str`` is ``True``. Otherwise,
        ``""``.
    """

    log_output = _set_default(log_output, "")
    as_str = _set_default(as_str, False)
    pdb_filename = _set_default(pdb_filename, "results.pdb")

    return write(
        ligand,
        receptor,
        all_interactions["closest"],
        all_interactions["close"],
        all_interactions["hydrophobics"],
        all_interactions["hydrogen_bonds"],
        all_interactions["halogen_bonds"],
        all_interactions["salt_bridges"],
        all_interactions["metal_coordinations"],
        all_interactions["pi_pi"],
        all_interactions["cat_pi"],
        all_interactions["active_site_flexibility"],
        log_output,
        as_str,
        pdb_filename,
    )
