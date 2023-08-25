#! /usr/bin/env python
# Author : Anuja
# 24.12.2013
# Last Updated : 24.12.2013

# This is an post-aroma run, external utility for filtering CMO-data for user-defined pi-MOs
# It takes the ".log" file as input and generates ".picmo" file with the same name.

import os
import sys
import json
import string
import fpformat

import aroma_user_constants
aroma_user_constants.loadUserConstants(sys.argv[0])
from aroma_constants import *
import aroma_constants

import aroma_util
from aroma_util import *

def grepPiCMO(piMOs, outfl, outExt):

#  Now plane is always 'XY'
   Plane = 'XY'

   olines = readFile(outfl + outExt)

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


def writeCollatedFiles(inputFileSet):

    jobType = "main"
    fileExt = "picmo"

    # get all files in the "main" run
    mainFiles = list(filter(lambda x: x["jobType"] == jobType and x["setIdx"] != "-1", inputFileSet))

    # get list of all centers
    centerList = list(set(list(map(lambda x: x["centerIdx"], mainFiles))))

    for centerIdx in centerList:
      idx = 0
      centerFileSet = list(filter(lambda x: x["centerIdx"] == centerIdx and x["jobType"] == jobType, mainFiles))

      baseprfx = centerFileSet[0]["baseprfx"] + "-center" + centerIdx

      final_file = open(externalProgram["outdir"] + baseprfx + "." + fileExt, "w")

#      collatedFileSet.append({"fileName": final_file, "flprfx": baseprfx, "jobType": jobType})
#      collatedFileSet.append({"fileName": externalProgram["outdir"] + baseprfx + "." + fileExt, "flprfx": baseprfx, "jobType": jobType, "centerIdx": centerIdx})

      for inpFil in centerFileSet:
        lines = readFile(externalProgram["outdir"] + inpFil["flprfx"] + "." + fileExt)
        for i in range(idx, len(lines)):
            final_file.write(lines[i])
        idx = 1

      final_file.close()

def main():

   piMOs = []
   if (len(sys.argv) > 3):
      for i in range (3, len(sys.argv)):
         piMOs.append(int(sys.argv[i]))
   else: print("Error: Usage is RUNTYPE filename MOs .. One of the inputs is missing, hence aborting .. "); sys.exit(10)

   runtype = sys.argv[1].upper()
   if (runtype == "GAUSSIAN" or runtype == "AROMA"): pass
   else: print("Error: The Runtype can either GAUSSIAN or AROMA. Aborting due to incorrect Runtype .. "); sys.exit(10)

   if (runtype == "AROMA"):
      flprfx  = sys.argv[2]
      outf = open(externalProgram["outdir"] + flprfx + "-inputFileSet.json", "r") # this file name needs to change, based on input
      inputFileSet = json.loads(outf.read())
      outf.close()

      # get all files in the "main" run
      jobType = "main"
      mainFiles = list(filter(lambda x: x["jobType"] == jobType and x["setIdx"] != "-1", inputFileSet))

      # get list of all centers
      centerList = list(set(list(map(lambda x: x["centerIdx"], mainFiles))))

      for centerIdx in centerList:
         idx = 0
         centerFileSet = list(filter(lambda x: x["centerIdx"] == centerIdx and x["jobType"] == jobType, mainFiles))

         for inpFil in centerFileSet:
           grepPiCMO(piMOs, externalProgram["outdir"] + inpFil["flprfx"], externalProgram["outExt"])

      writeCollatedFiles(inputFileSet)

   elif (runtype == "GAUSSIAN"):
      flnm = sys.argv[2]
      grepPiCMO(piMOs, externalProgram["outdir"] + flnm, "")

if __name__ == "__main__":
   main()
