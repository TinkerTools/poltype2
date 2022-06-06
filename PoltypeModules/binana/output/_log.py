# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

preface = "REMARK "

# __pragma__ ('skip')
# Python, just alias open
openFile = open
import json

# __pragma__ ('noskip')

"""?
# Transcrypt
import binana
from binana._utils.shim import OpenFile
import binana._utils.shim as json
openFile = OpenFile
?"""


def _center(string, length):
    while len(string) < length:
        string = " " + string
        if len(string) < length:
            string = string + " "
    return string


def _get_parameters(parameters, output):
    # restate the parameters
    output = output + preface + "Command-line parameters used:" + "\n"
    output = (
        output
        + preface
        + "                 Parameter              |            Value           "
        + "\n"
    )
    output = (
        output
        + preface
        + "   -------------------------------------|----------------------------"
        + "\n"
    )

    for key in sorted(list(parameters.params.keys())):
        value = str(parameters.params[key])
        output = (
            output
            + preface
            + "   "
            + _center(key, 37)
            + "| "
            + _center(value, 27)
            + "\n"
        )

    return output


def _get_close_contacts_dist1_cutoff(
    parameters,
    closest,
    output,
):
    output = output + preface + "\n"
    output = (
        output
        + preface
        + "Atom-type pair counts within "
        + str(parameters.params["close_contacts_dist1_cutoff"])
        + " angstroms:"
        + "\n"
    )
    output = output + preface + "    Atom Type | Atom Type | Count" + "\n"
    output = output + preface + "   -----------|-----------|-------" + "\n"
    for key in sorted(closest["counts"].keys()):
        value = closest["counts"][key]
        key = key.split("_")
        output = (
            output
            + preface
            + "   "
            + _center(key[0], 11)
            + "|"
            + _center(key[1], 11)
            + "|"
            + _center(str(value), 7)
            + "\n"
        )

    output = output + preface + "\n" + preface + "Raw data:\n"
    for atom_pairs in closest["labels"]:
        output = (
            output + preface + "     " + atom_pairs[0] + " - " + atom_pairs[1] + "\n"
        )
    return output


def _get_close_contacts_dist2_cutoff(parameters, close, output):
    output = output + preface + "\n"
    # output = output + preface + "\n"
    output = (
        output
        + preface
        + "Atom-type pair counts within "
        + str(parameters.params["close_contacts_dist2_cutoff"])
        + " angstroms:"
        + "\n"
    )
    output = output + preface + "    Atom Type | Atom Type | Count" + "\n"
    output = output + preface + "   -----------|-----------|-------" + "\n"
    for key in sorted(close["counts"].keys()):
        value = close["counts"][key]
        key = key.split("_")
        output = (
            output
            + preface
            + "   "
            + _center(key[0], 11)
            + "|"
            + _center(key[1], 11)
            + "|"
            + _center(str(value), 7)
            + "\n"
        )

    output = output + preface + "\n" + preface + "Raw data:\n"
    for atom_pairs in close["labels"]:
        output = (
            output + preface + "     " + atom_pairs[0] + " - " + atom_pairs[1] + "\n"
        )

    return output


def _get_ligand_atom_types(ligand_atom_types, electrostatic_energies, output):
    output = output + preface + "" + "\n"
    output = output + preface + "Ligand atom types:" + "\n"
    output = output + preface + "    Atom Type " + "\n"
    output = output + preface + "   -----------" + "\n"
    for key in sorted(ligand_atom_types.keys()):
        output = output + preface + "   " + _center(key, 11) + "\n"

    output = output + preface + "" + "\n"
    output = (
        output
        + preface
        + "Summed electrostatic energy by atom-type pair, in J/mol:"
        + "\n"
    )
    output = output + preface + "    Atom Type | Atom Type | Energy (J/mol)" + "\n"
    output = output + preface + "   -----------|-----------|----------------" + "\n"
    for key in sorted(electrostatic_energies["counts"].keys()):
        value = electrostatic_energies["counts"][key]
        key = key.split("_")

        value2 = str(value)[:13] if value > 0 else str(value)[:14]
        num_decimals = len(value2.split(".")[1])
        value3 = round(value, num_decimals)
        value4 = str(value3)[:13] if value3 > 0 else str(value3)[:14]

        output = (
            output
            + preface
            + "   "
            + _center(key[0], 11)
            + "|"
            + _center(key[1], 11)
            + "|"
            + _center(value4, 16)
            + "\n"
        )
    return output


def _get_rotateable_bonds_count(ligand, output):
    output = output + preface + "" + "\n"
    output = (
        output
        + preface
        + "Number of rotatable bonds in the ligand: "
        + str(ligand.rotatable_bonds_count)
        + "\n"
    )
    return output


def _get_active_site_flexibility(active_site_flexibility, output):
    output = output + preface + "" + "\n"
    output = output + preface + "Active-site flexibility:" + "\n"
    output = (
        output
        + preface
        + "    Sidechain/Backbone | Secondary Structure | Count "
        + "\n"
    )
    output = (
        output
        + preface
        + "   --------------------|---------------------|-------"
        + "\n"
    )
    for key in sorted(active_site_flexibility["counts"].keys()):
        value = active_site_flexibility["counts"][key]
        key = key.split("_")
        output = (
            output
            + preface
            + "   "
            + _center(key[0], 20)
            + "|"
            + _center(key[1], 21)
            + "|"
            + _center(str(value), 7)
            + "\n"
        )

    return output


def _get_hbonds(hbonds, output, hydrogen_bond=True):

    name = "Hydrogen" if hydrogen_bond else "Halogen"

    output = output + preface + "" + "\n"
    output = output + preface + name + " bonds:" + "\n"
    output = (
        output
        + preface
        + "    Location of Donor | Sidechain/Backbone | Secondary Structure | Count "
        + "\n"
    )
    output = (
        output
        + preface
        + "   -------------------|--------------------|---------------------|-------"
        + "\n"
    )
    for key in sorted(hbonds["counts"].keys()):
        value = hbonds["counts"][key]
        key = key.split("_")
        output = (
            output
            + preface
            + "   "
            + _center(key[1], 19)
            + "|"
            + _center(key[2], 20)
            + "|"
            + _center(key[3], 21)
            + "|"
            + _center(str(value), 7)
            + "\n"
        )

    output = output + preface + "\n" + preface + "Raw data:\n"
    for atom_pairs in hbonds["labels"]:
        output = (
            output
            + preface
            + "     "
            + atom_pairs[0]
            + " - "
            + atom_pairs[1]
            + " - "
            + atom_pairs[2]
            + "\n"
        )
    return output


def _get_hydrophobics(hydrophobics, output):
    output = output + preface + "" + "\n"
    output = output + preface + "Hydrophobic contacts (C-C):" + "\n"
    output = (
        output
        + preface
        + "    Sidechain/Backbone | Secondary Structure | Count "
        + "\n"
    )
    output = (
        output
        + preface
        + "   --------------------|---------------------|-------"
        + "\n"
    )
    for key in sorted(hydrophobics["counts"].keys()):
        value = hydrophobics["counts"][key]
        key = key.split("_")
        output = (
            output
            + preface
            + "   "
            + _center(key[0], 20)
            + "|"
            + _center(key[1], 21)
            + "|"
            + _center(str(value), 7)
            + "\n"
        )

    output = output + preface + "\n" + preface + "Raw data:\n"
    for atom_pairs in hydrophobics["labels"]:
        output = (
            output + preface + "     " + atom_pairs[0] + " - " + atom_pairs[1] + "\n"
        )
    return output


def _get_pi_stacking(pi_pi, output):
    stacking = []
    for key in sorted(pi_pi["counts"].keys()):
        value = pi_pi["counts"][key]
        together = key + "_" + str(value)
        if "STACKING" in together:
            stacking.append(together)

    output = output + preface + "" + "\n"
    output = output + preface + "pi-pi stacking interactions:" + "\n"
    output = output + preface + "    Secondary Structure | Count " + "\n"
    output = output + preface + "   ---------------------|-------" + "\n"
    for item in stacking:
        item = item.split("_")
        output = (
            output
            + preface
            + "   "
            + _center(item[1], 21)
            + "|"
            + _center(item[2], 7)
            + "\n"
        )

    output = output + preface + "\n" + preface + "Raw data:\n"
    for atom_pairs in pi_pi["labels"]["pi_stacking"]:
        output = (
            output + preface + "     " + atom_pairs[0] + " - " + atom_pairs[1] + "\n"
        )

    return output


def _get_t_stacking(pi_pi, output):
    t_shaped = []
    for key in sorted(pi_pi["counts"].keys()):
        value = pi_pi["counts"][key]
        together = key + "_" + str(value)
        if "SHAPED" in together:
            t_shaped.append(together)

    output = output + preface + "" + "\n"
    output = output + preface + "T-stacking (face-to-edge) interactions:" + "\n"
    output = output + preface + "    Secondary Structure | Count " + "\n"
    output = output + preface + "   ---------------------|-------" + "\n"
    for item in t_shaped:
        # need to check
        item = item.split("_")
        output = (
            output
            + preface
            + "   "
            + _center(item[1], 21)
            + "|"
            + _center(item[2], 7)
            + "\n"
        )

    output = output + preface + "\n" + preface + "Raw data:\n"
    for atom_pairs in pi_pi["labels"]["T_stacking"]:
        output = (
            output + preface + "     " + atom_pairs[0] + " - " + atom_pairs[1] + "\n"
        )

    return output


def _get_pi_cation(pi_pi, cat_pi, output):
    pi_cation = []
    for key in sorted(pi_pi["counts"].keys()):
        value = pi_pi["counts"][key]
        together = key + "_" + str(value)
        if "CATION" in together:
            pi_cation.append(together)

    output = output + preface + "" + "\n"
    output = output + preface + "Cation-pi interactions:" + "\n"
    output = (
        output
        + preface
        + "    Which residue is charged? | Secondary Structure | Count "
        + "\n"
    )
    output = (
        output
        + preface
        + "   ---------------------------|---------------------|-------"
        + "\n"
    )
    for item in pi_cation:
        # need to check
        item = item.split("_")
        item2 = item[1].split("-")
        output = (
            output
            + preface
            + "   "
            + _center(item2[0], 27)
            + "|"
            + _center(item[2], 21)
            + "|"
            + _center(item[3], 7)
            + "\n"
        )

    output = output + preface + "\n" + preface + "Raw data:\n"
    for atom_pairs in cat_pi["labels"]:
        output = (
            output + preface + "     " + atom_pairs[0] + " - " + atom_pairs[1] + "\n"
        )

    return output


def _get_salt_bridges(salt_bridges, output):
    output = output + preface + "" + "\n"
    output = output + preface + "Salt Bridges:" + "\n"
    output = output + preface + "    Secondary Structure | Count " + "\n"
    output = output + preface + "   ---------------------|-------" + "\n"
    for key in sorted(salt_bridges["counts"].keys()):
        value = salt_bridges["counts"][key]
        key = key.split("_")
        output = (
            output
            + preface
            + "   "
            + _center(key[1], 21)
            + "|"
            + _center(str(value), 7)
            + "\n"
        )

    output = output + preface + "\n" + preface + "Raw data:\n"
    for atom_pairs in salt_bridges["labels"]:
        output = (
            output + preface + "     " + atom_pairs[0] + " - " + atom_pairs[1] + "\n"
        )
    return output

def _get_metal_coordinations(metal_coordinations, output):
    output = output + preface + "" + "\n"
    output = output + preface + "Metal Coordinations:" + "\n"
    # TODO: Perhaps in some future version you could report something in the
    # table below.

    # output = output + preface + "    Secondary Structure | Count " + "\n"
    # output = output + preface + "   ---------------------|-------" + "\n"
    # for key in sorted(metal_coordinations["counts"].keys()):
    #     value = metal_coordinations["counts"][key]
    #     key = key.split("_")
    #     output = (
    #         output
    #         + preface
    #         + "   "
    #         + _center(key[1], 21)
    #         + "|"
    #         + _center(str(value), 7)
    #         + "\n"
    #     )

    output = output + preface + "\n" + preface + "Raw data:\n"
    for atom_pairs in metal_coordinations["labels"]:
        output = (
            output + preface + "     " + atom_pairs[0] + " - " + atom_pairs[1] + "\n"
        )

    return output


def collect(
    parameters,
    ligand,
    closest,
    close,
    hydrophobics,
    hydrogen_bonds,
    halogen_bonds,
    salt_bridges,
    metal_coordinations,
    pi_pi,
    cat_pi,
    electrostatic_energies,
    active_site_flexibility,
    ligand_atom_types,
    json_output,
):
    """Collects all the characterized interactions between the protein and
    ligand (as well as a few other metrics) into a single string (text block).

    Args:
        parameters (binana._cli_params.get_params.CommandLineParameters): An
            object containing the user-specified parameters. See
            :py:func:`~binana.run`.
        ligand (binana._structure.mol.Mol): The ligand object. Used for
            counting the number of rotatable bonds (if PDBQT formatted).
        closest (dict): A dictionary containing information about the closest
            protein/ligand interactions.
        close (dict): A dictionary containing information about the close
            protein/ligand interactions.
        hydrophobics (dict): A dictionary containing information about the
            hydrophobic protein/ligand interactions.
        hydrogen_bonds (dict): A dictionary containing information about the
            hydrogen bonds between the protein and ligand.
        halogen_bonds (dict): A dictionary containing information about the
            halogen bonds between the protein and ligand.
        salt_bridges (dict): A dictionary containing information about the
            salt-bridges protein/ligand interactions.
        metal_coordinations (dict): A dictionary containing information about
            the metal-coordination protein/ligand interactions.
        pi_pi (dict): A dictionary containing information about the pi-pi
            (stacking and T-shaped) protein/ligand interactions.
        cat_pi (dict): A dictionary containing information about the pi-cation
            protein/ligand interactions.
        electrostatic_energies (dict): A dictionary containing information
            about the electrostatic energies between protein and ligand atoms.
        active_site_flexibility (dict): A dictionary containing information
            about the flexibility of ligand-adjacent protein atoms.
        ligand_atom_types (dict): A dictionary containing information about
            the ligand atom types.
        json_output (dict): A dictionary describing the interactions, from
            :py:func:`~binana.output.dictionary.collect`

    Returns:
        str: The contents of the entire log.
    """

    output = ""

    output = _get_parameters(parameters, output)

    # a description of the analysis
    output = _get_close_contacts_dist1_cutoff(
        parameters,
        closest,
        output,
    )
    output = _get_close_contacts_dist2_cutoff(
        parameters,
        close,
        output,
    )
    output = _get_ligand_atom_types(
        ligand_atom_types["counts"], electrostatic_energies, output
    )
    output = _get_rotateable_bonds_count(ligand, output)
    output = _get_active_site_flexibility(active_site_flexibility, output)

    output = _get_hbonds(hydrogen_bonds, output, True)   # hydrogen bonds
    output = _get_hbonds(halogen_bonds, output, False)  # halogen bonds

    output = _get_hydrophobics(hydrophobics, output)

    output = _get_pi_stacking(pi_pi, output)
    output = _get_t_stacking(pi_pi, output)
    output = _get_pi_cation(pi_pi, cat_pi, output)
    output = _get_salt_bridges(salt_bridges, output)
    output = _get_metal_coordinations(metal_coordinations, output)

    # Append the JSON file to the end of the log as well (always).
    output = output + preface + "\n"
    output = output + preface + "JSON Output:\n"
    output = output + preface + "\n"
    output = output + preface + "\n"
    output = (
        output
        + preface
        + json.dumps(
            json_output, indent=2, sort_keys=True, separators=(",", ": ")
        ).replace("\n", "\n" + preface)
        + "\n"
    )

    # Output some files/to the screen.
    if parameters.params["output_dir"] != "":
        f = openFile(parameters.params["output_dir"] + "log.txt", "w")
        f.write(output.replace("REMARK ", ""))
        f.close()

    if parameters.params["output_file"] == "" and parameters.params["output_dir"] == "":
        # so you're not outputing to either a file or a directory. Output to the
        # screen.
        to_print = output.replace("REMARK ", "")
        print(to_print.split("JSON Output:")[0].strip())

    return output
