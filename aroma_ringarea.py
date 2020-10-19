#! /usr/bin/env python
# Author : Anuja
# 23.02.2014
# Last Updated : 23.02.2014

# Utility for AROMA to Determine the Area of the Ring

import os
import sys
import string

import aroma_constants
from aroma_constants import *

import aroma_util
from aroma_util import *

import aroma_parser
from aroma_parser import *

import aroma_molecule
from aroma_molecule import *

def ring_area(geomfl, idx):

   geomflext = geomfl[geomfl.rindex(".")+1:len(geomfl)+1]

   for extension in EXTENSIONS_FOR_GAUSSIAN_FILES:
      if (EXTENSIONS_FOR_GAUSSIAN_FILES[extension].count(geomflext) == 1): exttype = extension

   # Read and store molecular geometry
   theParser = ReaderFunctCall[exttype](geomfl)
   geom, hashLine, title, charge, mult = theParser.getInpData()
   conn_mat, Conn = genConnectivityMatrix(geom)

   ring_atoms = getOrderedRing(Conn, idx)
   reoriented_geom, new_points, new_normal = reorient(geom, Conn, ring_atoms)

   vertices = []
   for i in range (0, len(ring_atoms)):
       j = ring_atoms[i]
       vertices.append([reoriented_geom[j][1], reoriented_geom[j][2], reoriented_geom[j][3]])
   area = abs(round(areaOfPolygon(vertices, 'Z'), 3))

   return area


if __name__ == "__main__":
   geomfl = sys.argv[1]
   idx = []
   for i in range (2, len(sys.argv)):
      idx.append(string.atoi(sys.argv[i]))

   area = ring_area(geomfl, idx)
   print("\n\n       Area of Ring is " + repr(round(abs(area),3)) + " sq. unit\n\n")
