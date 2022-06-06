# This file is part of BINANA, released under the Apache 2.0 License. See
# LICENSE.md or go to https://opensource.org/licenses/Apache-2.0 for full
# details. Copyright 2021 Jacob D. Durrant.

# this file contains the Point class for binana.py

import math
import binana
from binana._utils.shim import r_just, round_to_thousandths_to_str


class Point:
    x = 99999.0
    y = 99999.0
    z = 99999.0

    # Initialize nitialize a point
    # Param x (float): x coordinate
    # Param y (float): y coordinate
    # Param z (float): z coordinate
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    # Returns a copy of a point
    # Param self (Point): point to be copied
    def copy_of(self):
        return Point(self.x, self.y, self.z)

    # Print the coordinates of a point
    # Param self (Point)
    def print_coors(self):
        print((str(self.x) + "\t" + str(self.y) + "\t" + str(self.z)))

    def snap(self, reso):  # snap the point to a grid
        self.x = round(self.x / reso) * reso
        self.y = round(self.y / reso) * reso
        self.z = round(self.z / reso) * reso

    # Returns the distance between two points
    # Param self (Point): this point
    # Param a_point (Point): the other point
    def dist_to(self, apoint):
        return math.sqrt(
            math.pow(self.x - apoint.x, 2)
            + math.pow(self.y - apoint.y, 2)
            + math.pow(self.z - apoint.z, 2)
        )

    # Returns a the coordinates of a point
    # Param self (Point)
    def description(self):
        return str(self.x) + " " + str(self.y) + " " + str(self.z)

    # Returns the magnitude of a point (distance from origin)
    # Param self (Point)
    def magnitude(self):
        return self.dist_to(Point(0, 0, 0))

    # Returns a PDB line for the point
    # Param self (Point)
    # Param index (integer): index of the point
    def create_pdb_line(self, index):
        output = "ATOM "
        output = output + r_just(str(index), 6) + r_just("X", 5) + r_just("XXX", 4)
        output = output + r_just(round_to_thousandths_to_str(self.x), 18)
        output = output + r_just(round_to_thousandths_to_str(self.y), 8)
        output = output + r_just(round_to_thousandths_to_str(self.z), 8)
        output = output + r_just("X", 24)
        return output
