#! /usr/bin/env python
# Author : Anuja
# 24.12.2013
# Last Updated : 24.12.2013

# This is an post-aroma run, external utility for filtering CMO-data for user-defined pi-MOs
# It takes the ".log" file as input and generates ".picmo" file with the same name.

import os
import sys
import string
import fpformat

import aroma_constants
from aroma_constants import *

import aroma_util
from aroma_util import *

def grepPiCMO(piMOs, outfl, GaussOutExt):

#  Now plane is always 'XY'
   Plane = 'XY'

   olines = readFile(outfl + GaussOutExt)

#  Find out No. of atoms and Bqs
   for o in olines:
      if ( (o.find("NAT") >= 0) and (o.find("NAtoms") >= 0) ):
         nat = int(o.split()[2])
         natoms = int(o.split()[4])
         nghost = natoms - nat
         break;

# Find out No. of basis funations, No. of electrons and No. of orbitals
   for o in range(0, len(olines)):
      if ( (olines[o].find("NBasis") >= 0) and (olines[o].find("NBE") >= 0) ):
         nbasis = int(olines[o].split()[1])
         nocc = int(olines[o].split()[3])
#         nelec = nocc*2
#         norb = int(olines[o+1].split()[1])
         break;


   lable_string = "pi-MO #  "
   for i in range (0, len(piMOs)):
      lable_string += repr(piMOs[i]) + "    "
   f = open(outfl + ".picmo", "w")

   f.write(lable_string + "   " + "-Sum \n")
   for i in range (nat+1, nat+nghost+1):
      for j in range (0, len(olines)):
         if ((olines[j].find("Full Cartesian NMR shielding tensor (ppm) for atom gh(" + repr(i).rjust(2) + ")") >= 0) and (olines[j+1].find("Canonical MO contributions") >= 0)):
            break;
         elif ((olines[j].find("Full Cartesian NMR shielding tensor (ppm) for atom gh(" + repr(i).rjust(3) + ")") >= 0) and (olines[j+1].find("Canonical MO contributions") >= 0)):
            break;

      data_string = '      '
      prvsmo = 0
      sumval = 0.0
      for k in range (0, nocc):
         if (olines[j+5+k+1].find("Total") >= 0):
            break
         else:
            words = olines[j+5+k].split()
            for l in range (prvsmo, len(piMOs)):
               if (int(float(words[0])) == int(piMOs[l])):
                  data_string += "  " + words[9]
                  sumval += float(words[9])
                  prvsmo = l+1
                  break
            continue

      prvs = j+5+k
      f.write(data_string + "  " + repr(-round(sumval,2)) + "\n")

   f.write("\n")
   f.close()


def main(Prfx):

   piMOs = []
   if (len(sys.argv) > 2):
      for i in range (2, len(sys.argv)):
         piMOs.append(int(sys.argv[i]))
   else: print "Error: MOs are not specified .. Aborting .. "; sys.exit(10)

   outfl = sys.argv[1]
   grepPiCMO(piMOs, outfl, GaussOutExt)

if __name__ == "__main__":
   main (sys.argv[1:len(sys.argv)])
