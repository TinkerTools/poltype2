# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

import math

# The maximum distance between a hydrogen- or halogen-bond donor and the middle
# (central atom), but it H or X. O-H distance is 0.96 A, N-H is 1.01 A. See
# http://www.science.uwaterloo.ca/~cchieh/cact/c120/bondel.html
_max_donor_X_dist = {
    "H": 1.3,
    "I": 2.04 * 1.4,  # O-I: 2.04 per avogadro
    "BR": 1.86 * 1.4,  # O-Br: 1.86
    "Br": 1.86 * 1.4,
    "CL": 1.71 * 1.4,  # O-Cl: 1.71
    "Cl": 1.71 * 1.4,
    "F": 1.33 * 1.4,  # O-F: 1.33
}

_alternate_protein_resname = {
    "LYS": ["LYS", "LYN"],
    "HIS": ["HIS", "HID", "HIE", "HIP"],
    "GLU": ["GLU", "GLH", "GLX"],
    "ASP": ["ASP", "ASH", "ASX"],
}

# Information about protein atoms that can be hydrogen-bond donors. Useful for
# when hydrogen atoms haven't been added to the protein model.
_protein_hydro_bond_donors = [
    [["ARG"], ["NE", "NH1", "NH2"]],
    [_alternate_protein_resname["HIS"], ["NE2", "ND1"]],
    [_alternate_protein_resname["LYS"], ["NZ"]],
    [["SER"], ["OG"]],
    [["THR"], ["OG1"]],
    [["ASN"], ["ND2"]],
    [["GLN"], ["NE2"]],
    [["TYR"], ["OH"]],
    [["TRP"], ["NE1"]],
    [["CYS"], ["SG"]]
]

protein_resnames = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "ASH",
    "ASX",
    "CYS",
    "CYM",
    "CYX",
    "GLN",
    "GLU",
    "GLH",
    "GLX",
    "GLY",
    "HIS",
    "HID",
    "HIE",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "LYN",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
]


# Specifies names that must be present in protein residues. Will throw a warning
# otherwise. Important for detecting charged groups and hydrogen bonds in some
# cases.
_required_protein_atom_names = [
    [_alternate_protein_resname["GLU"], ["OE1", "OE2"]],
    [_alternate_protein_resname["ASP"], ["OD1", "OD2"]],
    [["ARG"], ["NH1", "NH2", "NE"]],
    [_alternate_protein_resname["HIS"], ["NE2", "ND1"]],
    [["PHE"], ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]],
    [["TYR"], ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"]],
    [["TRP"], ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"]],
    [_alternate_protein_resname["HIS"], ["CG", "ND1", "CD2", "CE1", "NE2"]],
    [_alternate_protein_resname["LYS"], ["NZ"]],
]

to_deg = 180.0 / math.pi

two_leter_atom_names = [
    "AC", "AG", "AL", "AM", "AU", "BA", "BE", "BK", "CA", "CO", "CU", "DB",
    "DY", "ER", "ES", "EU", "FE", "GA", "GD", "GE", "LA", "LR", "LU", "MD",
    "MG", "MN", "MO", "NI", "PB", "RA", "RE", "RF", "RH", "RU", "TA", "TB",
    "TC", "TH", "TI", "TL", "TM", "YB", "ZN", "ZR", "BR", "CL", "BI", "AS",
    "LI",

    # Rare ones no longer using because of potential for confusion with other
    # elements (not worth it).
    # "CD", "CE",
    # "CF", "CM", "CR", "CS", "FM", "FR", 
    # "HF", "HG", "HO", "IN", "IR", 
    # "NO", "NP", "NB", "ND", 
    # "OS", "PA", "PD", "PM", "PO", "PR", "PT", "PU", 
    # "SB", "SC", "SG", "SM", "SN", "SR", 

]