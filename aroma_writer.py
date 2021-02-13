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
   def __init__(self, geomfl):
      self.geomfl = geomfl
      self.geom = {}
      self.hashLine = aroma_constants.externalProgram["defaultNicsKeyline"]
      self.title = "DEFAULT TITLE SET BY AROMA"
      self.charge = 0
      self.mult = 0

   def write(self):
      pass

class GaussianInputFileWriter(FileParser):
   def __init__(self, geomfl):
      FileParser.__init__(self, geomfl)

   def write(self, flprfx, externalProgram):
      # The chk file name for optimization can be same as given by the user.
      # Remove the above lines after testing this.
      hashLine_rev = hashLine

      # Generate input file for optimization
      optfl = flprfx + "-opt"
      f_opt = open(externalProgram["inpdir"] + optfl + externalProgram["inpExt"], "w")

      title += " Optimization By Aroma "
      f_opt.write(hashLine_rev + "\n" + title + "\n\n" + repr(charge) + " " + repr(mult) + "\n")

      coord_format = "{0:.5f}"
      for i in range (1, len(geom)+1):
         geomline = repr(geom[i][0]) + "   " + coord_format.format(geom[i][1]) + "   " + coord_format.format(geom[i][2]) + "   " + coord_format.format(geom[i][3]) + "\n"
         f_opt.write(geomline)

      f_opt.write("\n")
      f_opt.close()

      print("Input file for Opitmization named as " + optfl + " is generated.")
      return optfl

# Reader Function to Be Called for each Type of Format
# (Ganesh: Moved this here to remove the cyclic dependency)
GaussianSettings["writerFunctCall"] = {'geomInput':GaussianInputFileWriter}

