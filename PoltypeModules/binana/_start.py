# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

import __future__
from binana.output import _write_main
from binana.interactions import get_all_interactions
from binana.load_ligand_receptor import from_files

import math
import binana

# __pragma__ ('skip')
# Python
import textwrap
from math import fabs

# __pragma__ ('noskip')

"""?
# Transcrypt
from binana._utils import shim
textwrap = shim
from binana._utils.shim import fabs
?"""

VERSION = "2.1"


def _get_all_interactions(parameters):
    """Gets all the interactions between the specified ligand and receptor
    files.

    Args:
        parameters (binana._cli_params.get_params.CommandLineParameters):
            The BINANA parameters to use.
    """

    max_cutoff = (
        max(
            i[1]
            for i in parameters.params.items()
            if "dist" in i[0] and "cutoff" in i[0]
        )
        + 15
    )

    ligand, receptor = from_files(
        parameters.params["ligand"], parameters.params["receptor"], max_cutoff
    )

    # This is perhaps controversial. I noticed that often a pi-cation
    # interaction or other pi interaction was only slightly off, but looking at
    # the structure, it was clearly supposed to be a pi-cation interaction. I've
    # decided then to artificially expand the radius of each pi ring. Think of
    # this as adding in a VDW radius, or accounting for poor crystal-structure
    # resolution, or whatever you want to justify it.
    pi_padding = parameters.params["pi_padding_dist"]

    all_interacts = get_all_interactions(
        ligand,
        receptor,
        parameters.params["close_contacts_dist1_cutoff"],
        parameters.params["close_contacts_dist2_cutoff"],
        parameters.params["electrostatic_dist_cutoff"],
        parameters.params["active_site_flexibility_dist_cutoff"],
        parameters.params["hydrophobic_dist_cutoff"],
        parameters.params["hydrogen_bond_dist_cutoff"],
        parameters.params["hydrogen_halogen_bond_angle_cutoff"],
        parameters.params["halogen_bond_dist_cutoff"],
        parameters.params["pi_pi_interacting_dist_cutoff"],
        parameters.params["pi_stacking_angle_tolerance"],
        parameters.params["T_stacking_angle_tolerance"],
        parameters.params["T_stacking_closest_dist_cutoff"],
        parameters.params["cation_pi_dist_cutoff"],
        parameters.params["salt_bridge_dist_cutoff"],
        parameters.params["metal_coordination_dist_cutoff"],
        pi_padding,
    )

    # The original implementation merged all pi-related interactions into
    # one. Do that here too for backwards compatibility.
    for key in all_interacts["cat_pi"]["counts"].keys():
        all_interacts["pi_pi"]["counts"][key] = all_interacts["cat_pi"]["counts"][key]

    # Now save the files
    _write_main(
        parameters,
        ligand,
        receptor,
        all_interacts["closest"],
        all_interacts["close"],
        all_interacts["hydrophobics"],
        all_interacts["hydrogen_bonds"],
        all_interacts["halogen_bonds"],
        all_interacts["salt_bridges"],
        all_interacts["metal_coordinations"],
        all_interacts["pi_pi"],
        all_interacts["cat_pi"],
        all_interacts["electrostatic_energies"],
        all_interacts["active_site_flexibility"],
        all_interacts["ligand_atom_types"],
    )


def _intro():
    print("# BINANA " + VERSION + "\n")

    # __pragma__ ('skip')
    import os

    dir_path = os.path.abspath(
        os.path.dirname(os.path.realpath(__file__)) + "/../COMMAND_LINE_USE.md"
    )

    if os.path.exists(dir_path):
        with open(dir_path) as f:
            print("\n" + f.read().rstrip())

    import sys
    import time
    if sys.version_info[0] < 3:
        print("\nWARNING: Only Python 3 is officially supported!\n")
        time.sleep(2)

    # __pragma__ ('noskip')

    print("\n                            [- END INTRO -]\n")
