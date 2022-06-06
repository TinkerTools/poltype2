# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

# this file contains the CommandLineParameters class
# for binana.py

from binana.interactions.default_params import (
    ACTIVE_SITE_FLEXIBILITY_DIST_CUTOFF,
    CATION_PI_DIST_CUTOFF,
    CLOSE_CONTACTS_DIST1_CUTOFF,
    CLOSE_CONTACTS_DIST2_CUTOFF,
    ELECTROSTATIC_DIST_CUTOFF,
    HYDROGEN_HALOGEN_BOND_ANGLE_CUTOFF,
    HYDROGEN_BOND_DIST_CUTOFF,
    HALOGEN_BOND_DIST_CUTOFF,
    HYDROPHOBIC_DIST_CUTOFF,
    LIGAND,
    OUTPUT_DIR,
    OUTPUT_FILE,
    OUTPUT_JSON,
    OUTPUT_CSV,
    PI_PADDING_DIST,
    PI_PI_INTERACTING_DIST_CUTOFF,
    PI_STACKING_ANGLE_TOLERANCE,
    METAL_COORDINATION_DIST_CUTOFF,
    RECEPTOR,
    SALT_BRIDGE_DIST_CUTOFF,
    T_STACKING_ANGLE_TOLERANCE,
    T_STACKING_CLOSEST_DIST_CUTOFF,
    TEST,
)

# __pragma__ ('skip')
# Python
from os import sep
# __pragma__ ('noskip')

"""?
# Transcrypt
import binana
sep = "/"
?"""

"""
Class CommmandLineParameters
"""


class CommandLineParameters:
    params = {}

    def is_num(self, num):
        try:
            return float(num)
        except ValueError:
            return num

    def __init__(self, parameters):

        # first, set defaults
        self.params["close_contacts_dist1_cutoff"] = CLOSE_CONTACTS_DIST1_CUTOFF
        self.params["close_contacts_dist2_cutoff"] = CLOSE_CONTACTS_DIST2_CUTOFF
        self.params["electrostatic_dist_cutoff"] = ELECTROSTATIC_DIST_CUTOFF
        self.params[
            "active_site_flexibility_dist_cutoff"
        ] = ACTIVE_SITE_FLEXIBILITY_DIST_CUTOFF
        self.params["hydrophobic_dist_cutoff"] = HYDROPHOBIC_DIST_CUTOFF
        self.params["hydrogen_bond_dist_cutoff"] = HYDROGEN_BOND_DIST_CUTOFF
        self.params["hydrogen_halogen_bond_angle_cutoff"] = HYDROGEN_HALOGEN_BOND_ANGLE_CUTOFF
        self.params["halogen_bond_dist_cutoff"] = HALOGEN_BOND_DIST_CUTOFF
        self.params["pi_padding_dist"] = PI_PADDING_DIST
        self.params["pi_pi_interacting_dist_cutoff"] = PI_PI_INTERACTING_DIST_CUTOFF
        self.params["pi_stacking_angle_tolerance"] = PI_STACKING_ANGLE_TOLERANCE
        self.params["T_stacking_angle_tolerance"] = T_STACKING_ANGLE_TOLERANCE
        self.params["T_stacking_closest_dist_cutoff"] = T_STACKING_CLOSEST_DIST_CUTOFF
        self.params["cation_pi_dist_cutoff"] = CATION_PI_DIST_CUTOFF
        self.params["salt_bridge_dist_cutoff"] = SALT_BRIDGE_DIST_CUTOFF
        self.params["metal_coordination_dist_cutoff"] = METAL_COORDINATION_DIST_CUTOFF
        self.params["receptor"] = RECEPTOR
        self.params["ligand"] = LIGAND
        self.params["output_dir"] = OUTPUT_DIR
        self.params["output_file"] = OUTPUT_FILE
        self.params["output_json"] = OUTPUT_JSON
        self.params["output_csv"] = OUTPUT_CSV
        self.params["test"] = TEST

        # now get user inputed values

        for index in range(len(parameters)):
            item = parameters[index]
            if len(item) > 0 and item[0] == "-":
                # so it's a parameter key value
                key = item.replace("-", "")
                value = self.is_num(parameters[index + 1])
                if key in list(self.params.keys()):
                    self.params[key] = value
                    parameters[index] = ""
                    parameters[index + 1] = ""

        # make a list of all the command-line parameters not used
        self.error = ""
        for index in range(1, len(parameters)):
            item = parameters[index]
            if item != "":
                self.error = self.error + item + " "

        # Make sure the output directory, if specified, ends in a /
        if self.params["output_dir"] != "" and self.params["output_dir"][-1:] != sep:
            self.params["output_dir"] = self.params["output_dir"] + sep

        # If an output directory is specified but a log file isn't, set a
        # default logfile
        single_output_files = [("output_file", "pdb"), ("output_json", "json"), ("output_csv", "csv")]
        for single_output_file, ext in single_output_files:
            if (
                self.params["output_dir"] != ""
                and self.params[single_output_file] == ""
            ):
                self.params[single_output_file] = (
                    self.params["output_dir"] + "output." + ext
                )

    def okay_to_proceed(self):
        # at the very least, you need the ligand and the receptor
        return self.params["receptor"] != "" and self.params["ligand"] != ""
