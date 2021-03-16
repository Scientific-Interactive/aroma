#! /usr/bin/env python
# Author : Ganesh
# @ 04.02.2021

# Classes for Writing Different Inputs (com, log, chk)

import os
import sys
import string

import aroma_util
from aroma_util import *

import aroma_constants
from aroma_constants import *

class FileWriter:
   def __init__(self):
      self.hashLine = aroma_constants.externalProgram["defaultNicsKeyline"]
      self.title = "DEFAULT TITLE SET BY AROMA"

   def writeOptFile(self, flprfx, externalProgram, geom, hashLine, title, charge, mult):
      pass

   def genHeader(self, hashLine, title, charge, mult):
      pass

   def genGeomLine(self, geom):
      pass

   def genGhostAtomLine(self, pt):
      pass

   def genGhostAtomSetBreak(self):
      pass

   def genCheckpointLine(self, chkFile):
      pass

   def genTerminator(self):
      pass

class GaussianInputFileWriter(FileWriter):
   def __init__(self):
      FileWriter.__init__(self)

      self.coord_format = "{0:.5f}"

   def genHeader(self, hashLine, title, charge, mult):
      header = title 
      header = hashLine + header + "\n\n" + repr(charge) + " " + repr(mult) + "\n"
      return header

   def genGeomLine(self, geom):
      geomline = repr(geom[0]) + "   " + self.coord_format.format(geom[1]) + "   " + self.coord_format.format(geom[2]) + "   " + self.coord_format.format(geom[3]) + "\n"
      return geomline

   def genGhostAtomLine(self, pt):
      ghostAtm = 'bq     ' + self.coord_format.format(pt[0]) + "     " + self.coord_format.format(pt[1]) + "     " + self.coord_format.format(pt[2]) + "\n"
      return ghostAtm

   def genGhostAtomSetBreak(self):
      return "break"

   def genCheckpointLine(self, chkFile):
      return "%chk=" + chkFile + "\n"

   def genTerminator(self):
      return ""

   def writeOptFile(self, flprfx, externalProgram, geom, hashLine, title, charge, mult):
      # The chk file name for optimization can be same as given by the user.
      # Remove the above lines after testing this.
      hashLine_rev = hashLine

      # Generate input file for optimization
      optfl = flprfx + "-opt"
      f_opt = open(externalProgram["inpdir"] + optfl + externalProgram["inpExt"], "w")

      header = self.genHeader(hashLine_rev, title + " Optimization By Aroma ", charge, mult)
      f_opt.write(header)

      coord_format = "{0:.5f}"
      for i in range (1, len(geom)+1):
         geomline = self.genGeomLine(geom[i]) 
         f_opt.write(geomline)

      f_opt.write("\n")
      f_opt.close()

      return optfl

class OrcaInputFileWriter(FileWriter):
   def __init__(self):
      FileWriter.__init__(self)

      self.coord_format = "{0:.5f}"

   def genHeader(self, hashLine, title, charge, mult):
      header = ""
      header = hashLine + header + "\n\n*xyz " + repr(charge) + " " + repr(mult) + "\n"
      return header

   def genGeomLine(self, geom):
      geomline = repr(geom[0]) + "   " + self.coord_format.format(geom[1]) + "   " + self.coord_format.format(geom[2]) + "   " + self.coord_format.format(geom[3]) + "\n"
      return geomline

   def genGhostAtomLine(self, pt):
      ghostAtm = 'H:    ' + self.coord_format.format(pt[0]) + "     " + self.coord_format.format(pt[1]) + "     " + self.coord_format.format(pt[2]) + " newgto S 1 1 100000 1 end newauxJKgto S 1 1 200000 1 end\n"
      return ghostAtm

   def genGhostAtomSetBreak(self):
      return "break"

   def genCheckpointLine(self, chkFile):
      return ""

   def genTerminator(self):
      return "*"

   def writeOptFile(self, flprfx, externalProgram, geom, hashLine, title, charge, mult):
      # The chk file name for optimization can be same as given by the user.
      # Remove the above lines after testing this.
      hashLine_rev = hashLine

      # Generate input file for optimization
      optfl = flprfx + "-opt"
      f_opt = open(externalProgram["inpdir"] + optfl + externalProgram["inpExt"], "w")

      header = self.genHeader(hashLine_rev, title + " Optimization By Aroma ", charge, mult)
      f_opt.write(header)

      coord_format = "{0:.5f}"
      for i in range (1, len(geom)+1):
         geomline = self.genGeomLine(geom[i])
         f_opt.write(geomline)

      f_opt.write("\n" + self.genTerminator() + "\n")
      f_opt.close()

      return optfl


# Reader Function to Be Called for each Type of Format
# (Ganesh: Moved this here to remove the cyclic dependency)
GaussianSettings["writerFunctCall"] = {'geomInput':GaussianInputFileWriter()}
OrcaSettings["writerFunctCall"] = {'geomInput':OrcaInputFileWriter()}

