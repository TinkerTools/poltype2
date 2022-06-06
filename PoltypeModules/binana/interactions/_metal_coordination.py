# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

import math
from binana._utils.shim import _set_default
from binana.interactions.default_params import METAL_COORDINATION_DIST_CUTOFF
import binana
from binana.load_ligand_receptor import _get_ligand_receptor_dists
from binana._utils.utils import hashtable_entry_add_one, list_alphebetize_and_combine
from binana._structure.mol import Mol
from binana._utils._math_functions import angle_between_three_points

# Be sure to update the corresponding function in
# binana.interactions.__init__.py as well!


# TODO: The below represents initial efforts to account for angles when
# detecting metal coordination. I have decided to ignore angles for now because:
#
# 1. There are many different coordinating geometries.
# 2. Even if a coordinating geometry were known, atomic positions in practice
#    often deviate from the ideal angles.
# 3. Coordination positions can be unoccupied (vacancy), further complicating
#    detection.
#
# The good news is that when key atoms are near metal cations, they almost
# always seem to participate in metal coordination bonds to one degree or
# another, so I think detecting such bonds by distance is reasonable. If you
# ever need to find examples of such bonds for further testing, see
# https://metalpdb.cerm.unifi.it/

# def _approve_metal_coordination_geometry(metal_coordinations):
#     metal_atom = metal_coordinations[0]
#     metal_atom_coor = metal_atom.coordinates
#     coord_atoms = metal_coordinations[1:]

#     to_deg = 180.0 / math.pi

#     angles = []
#     for i1 in range(len(coord_atoms) - 1):
#         coord_atom1 = coord_atoms[i1]
#         coord_atom1_coors = coord_atom1.coordinates
#         for i2 in range(i1 + 1, len(coord_atoms)):
#             coord_atom2 = coord_atoms[i2]
#             coord_atom2_coors = coord_atom2.coordinates
#             angles.append(
#                 angle_between_three_points(
#                     coord_atom1_coors, metal_atom_coor, coord_atom2_coors
#                 )
#                 * to_deg
#             )
#             # print(coord_atom1, coord_atom2)
#     angles.sort()
#     print(angles)
#     # See
#     # https://chem.libretexts.org/Bookshelves/General_Chemistry/Chemistry_(OpenSTAX)/19%3A_Transition_Metals_and_Coordination_Chemistry/19.2%3A_Coordination_Chemistry_of_Transition_Metals
#     # 1. Linear (2): 180 degrees
#     # 2. Trigonal planar (3): 60 degrees
#     # 3. Tetrahedral (4): 6 109 degrees
#     # 4. Square planar (4): 4 90 degrees
#     # 5. Trigonal bipyramidal (5): 6 90 degrees, 3 120 degrees
#     # 6. Square pyramidal (5): 8 90 degrees
#     # 7. Octahedral (6): 12 90 degrees
#     # 8. Pentagonal bipyramid (7): 10 90 degrees, 5 72 degrees
#     # 9. Square antiprism (8): Challenging to calculate...
#     #    https://en.wikipedia.org/wiki/Square_antiprismatic_molecular_geometry
#     # 10. Dodecahedron: 41.8 degrees. See
#     #     https://math.stackexchange.com/questions/451943/angles-between-two-vertices-on-a-dodecahedron


def get_metal_coordination(ligand, receptor, cutoff=None):
    """Identifies and counts the number of metal-coordination protein/ligand
    contacts. Output is formatted like this::

        {
            'counts': {
                'N_ZN': 3,
                'O_ZN': 2
            },
            'labels': [
                (
                    'A:ZN(201):ZN(3059)',
                    'A:HIS(97):ND1(1426)',
                    {'distance': 1.974986835399159}
                ),
                (
                    'A:ZN(201):ZN(3059)',
                    'A:HIS(100):NE2(1470)',
                    {'distance': 2.0332422383965976}
                )
            ],
            'mol': <binana._structure.mol.Mol instance at 0x7feb20290908>
        }

    Args:
        ligand (binana._structure.mol.Mol): The ligand molecule to analyze.
        receptor (binana._structure.mol.Mol): The receptor molecule to analyze.
        cutoff (float, optional): The distance cutoff. Defaults to
            METAL_COORDINATION_DIST_CUTOFF.

    Returns:
        dict: Contains the atom tallies ("counts"), a binana._structure.mol.Mol
        object with the participating atoms ("mol"), and the labels to use in
        the log file ("labels").
    """

    cutoff = _set_default(cutoff, METAL_COORDINATION_DIST_CUTOFF)

    # Calculate the distances. See
    # https://chem.libretexts.org/Bookshelves/General_Chemistry/Chemistry_(OpenSTAX)/19%3A_Transition_Metals_and_Coordination_Chemistry/19.2%3A_Coordination_Chemistry_of_Transition_Metals
    # metals = ["Co", "Pt", "Ag", "Cr", "Ni", "Sn", "Mn", "Fe", "Cu", "Zn"]
    metals = [
        "Ac",
        "Ag",
        "Al",
        "Am",
        "Au",
        "Ba",
        "Be",
        "Bi",
        "Bk",
        # "Cd",  # Keep C- ones to min to not confusion with carbon
        # "Ce",
        # "Cf",
        # "Cm",
        # "Cr",
        # "Cs",
        "Ca",
        "Co",
        "Cu",
        "Db",
        "Dy",
        "Er",
        "Es",
        "Eu",
        # "Fm",  # To avoid confusion with F
        # "Fr",
        "Fe",
        "Ga",
        "Gd",
        "Ge",
        # "Hf",  # I worry about confusing with a hydrogen
        # "Hg",
        # "Ho",
        # "In",  # Too much potential for confusion with I
        # "Ir",
        "La",
        "Lr",
        "Lu",
        "Md",
        "Mg",
        "Mn",
        "Mo",
        # "No",  # Keep to a min because of confusion with N.
        # "Np",
        # "Nb",
        # "Nd",
        "Ni",
        # "Os",  # Confusion with O
        # "Pa",  # Minimize confusion with P
        # "Pd",
        # "Pm",
        # "Po",
        # "Pr",
        # "Pt",
        # "Pu",
        "Pb",
        "Ra",
        "Re",
        "Rf",
        "Rh",
        "Ru",
        # "Sb",  # Confusion with S
        # "Sc",
        # "Sg",
        # "Sm",
        # "Sn",
        # "Sr",
        "Ta",
        "Tb",
        "Tc",
        "Th",
        "Ti",
        "Tl",
        "Tm",
        "Yb",
        "Zn",
        "Zr",
    ]

    for m in metals[:]:
        metals.append(m.upper())

    coord_lig_atoms = ["N", "O", "Cl", "F", "Br", "I", "CL", "BR", "S"]

    atoms_together = metals[:]
    atoms_together.extend(coord_lig_atoms)

    metal_coord_dists = _get_ligand_receptor_dists(
        ligand, receptor, cutoff, atoms_together
    )

    atom_type_counts = {}
    pdb_metal_coordinations = Mol()
    metal_coordinations_labels = []

    for ligand_atom, receptor_atom, dist in metal_coord_dists:
        if (ligand_atom.element in metals and receptor_atom.element not in metals) or (
            ligand_atom.element not in metals and receptor_atom.element in metals
        ):
            pdb_metal_coordinations.add_new_atom(ligand_atom.copy_of())
            pdb_metal_coordinations.add_new_atom(receptor_atom.copy_of())

            list_metal_atom = [ligand_atom.atom_type, receptor_atom.atom_type]
            hashtable_entry_add_one(
                atom_type_counts,
                list_alphebetize_and_combine(list_metal_atom),
            )

            metal_coordinations_labels.append(
                (
                    ligand_atom.string_id(),
                    receptor_atom.string_id(),
                    {"distance": dist},
                )
            )

    # for metal_id in metal_coordinations.keys():

    #     # _approve_metal_coordination_geometry(metal_coordinations[metal_id])

    #     atoms = metal_coordinations[metal_id]
    #     metal_atom = atoms[0]
    #     coord_atoms = atoms[1:]

    #     pdb_metal_coordinations.add_new_atom(metal_atom.copy_of())

    #     coord_atoms_labels = []

    #     for coord_atom in coord_atoms:
    #         list_metal_atom = [metal_atom.atom_type, coord_atom.atom_type]
    #         hashtable_entry_add_one(
    #             atom_type_counts,
    #             list_alphebetize_and_combine(list_metal_atom),
    #         )
    #         pdb_metal_coordinations.add_new_atom(coord_atom.copy_of())

    #         coord_atoms_labels.append(coord_atom.string_id())

    #     metal_coordinations_labels.append((metal_atom.string_id(), coord_atoms_labels))

    return {
        "counts": atom_type_counts,
        "mol": pdb_metal_coordinations,
        "labels": metal_coordinations_labels,
    }
