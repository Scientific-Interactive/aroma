#! /usr/bin/env python
# Author : Anuja
# @ 12.12.2012
# Last Updated : 27-08-2014

# Identifying Pi-MOs

import os
import sys
import time
import string
import struct
import fpformat

import aroma_constants
from aroma_constants import *

import aroma_util
from aroma_util import *

def identifyPiMOs(outfl):

   olines = readFile(outfl)
#  Find out No. of basis funations, No. of electrons and No. of orbitals
   for o in range(0, len(olines)):
      if (olines[o].find("NActive") >= 0): nat = int(olines[o].split()[1])
      if ( (olines[o].find("NBasis") >= 0) and (olines[o].find("NBE") >= 0) ):
         nbasis = int(olines[o].split()[1])
         nocc = int(olines[o].split()[3])
         nelec = nocc*2
         norb = int(olines[o+1].split()[1])
         break;


#  Find the line from where MO data starts
   for MODataBegin in range (0, len(olines)):
      if ( olines[MODataBegin].find("Molecular Orbital Coefficients") >= 0 ): break;

#  Store Basis Functions Map with Index as:
   Idx_Bas = [] # Index for Basis Function <--> Corresponding string for basis function
   icnt_atmicno = [] # Index for Basis Function <--> The Atom Symbol to which it belongs
   atm_idx = 0
   for i in range (MODataBegin+4, MODataBegin+4+nbasis):
     line = olines[i][0:22]
     words = line.split()

     if (len(words) == 4):
        Idx_Bas.append(words[3])
        atmic_no = AtmSym[words[2].upper()]
        icnt_atmicno.append(atmic_no)
     else: 
        Idx_Bas.append(words[1])
        icnt_atmicno.append(atmic_no)


#  Read and Store MOs and orbital energies
   orb = {} # Dictionary of Molecular Orbitals:: MO Index <--> MO vector
   orb_ene = {} # Dictionary of Orbital Energies:: MO Index <--> Energy

   for i in range (MODataBegin+1, len(olines), nbasis+3):
      if (olines[i].upper().find("DENSITY MATRIX") >= 0 ): break

      words = list(map(int, olines[i].split()))
      ene_words = list(map(float, olines[i+2].upper().split("EIGENVALUES --")[1].split()))
      for j in range (0, len(words)):
         orb[words[j]] = [] # Define the key (index) for the MO vector (value)
         orb_ene[words[j]] = ene_words[j]

      for k in range (i+3, i+3+nbasis):
         line = olines[k][21:71].strip("\n")

# Life was so simple, if Gaussian did not have formatted write!
#         coeff_words = list(map(float, line.split()))
         strFormat = ""
         for t in range (0, int(len(line)/10)): strFormat += "10s"

         # print("life line", line, strFormat)
         # print(struct.unpack(strFormat, bytes(line, 'utf8')))

         coeff_words = list(map(float, list(struct.unpack(strFormat, bytes(line, 'utf8')))))

         for j in range (0, len(words)):
            orb[words[j]].append(coeff_words[j])


#  Indentifying the pi-MOs
   piMOs = []
   for i in range (1, nocc+1):
      FLAG = 1
      for j in range (0, nbasis):
         eligible_for_ring = icnt_atmicno[j] in Max_Conn
         if (eligible_for_ring):
            if (abs(orb[i][j]) > 5e-4):
                if (any ((Idx_Bas[j].upper().find(S) >= 0) for S in ["1S","2S","3S"])):
                   FLAG = 0
                   break
            if (abs(orb[i][j]) > 5e-3):
                if (any ((Idx_Bas[j].upper().find(S) >= 0) for S in ["PX","PY"])):
                   FLAG = 0
                   break

      if (FLAG): piMOs.append(i)

   return piMOs, nocc


def grepPiCMO(nat, piMOs, nocc, nghost, BQ_Range, BQ_Step, outfl, outExt, distIdx):

#  Now plane is always 'XY'
   Plane = 'XY'

   olines = readFile(outfl + outExt)
   lable_string = "pi-MO#  "
   for i in range (0, len(piMOs)):
      lable_string += repr(piMOs[i]) + "    "
   f = open(outfl + ".picmo", "w")
   
   f.write(lable_string + "   " + "-Sum \n")
   dist = distIdx # BQ_Range[0]
   for i in range (nat+1, nat+nghost+1):
      for j in range (0, len(olines)):
         if ((olines[j].find("Full Cartesian NMR shielding tensor (ppm) for atom gh(" + repr(i).rjust(2) + ")") >= 0) and (olines[j+1].find("Canonical MO contributions") >= 0)): 
            break;
         elif ((olines[j].find("Full Cartesian NMR shielding tensor (ppm) for atom gh(" + repr(i).rjust(3) + ")") >= 0) and (olines[j+1].find("Canonical MO contributions") >= 0)): 
            break;

      data_string = '      '
      cmo_init = {}
      for c in range (0, len(piMOs)): 
          cmo_init[piMOs[c]] = 0.0
      prvsmo = 0
      sumval = 0.0
      data_string = "   " + fpformat.fix(dist,2)
      for k in range (0, nocc):
         if (olines[j+5+k+1].find("Total") >= 0):
            break
         else:
            words = olines[j+5+k].split()
            for l in range (prvsmo, len(piMOs)):
               if (int(float(words[0])) == int(piMOs[l])):
                  cmo_init[piMOs[l]] = float(words[9])
#                  data_string += "  " + words[9]
#                  sumval += float(words[9])
                  prvsmo = l+1
                  break
            continue
      for c in range (0, len(piMOs)):
          data_string += "  " + fpformat.fix(cmo_init[piMOs[c]],2)
          sumval += cmo_init[piMOs[c]] 
      dist += BQ_Step

      prvs = j+5+k  
      f.write(data_string + "  " + repr(-round(sumval,2)) + "\n")

   f.close()
   return dist





