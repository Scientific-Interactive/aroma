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

   def write(self):
      pass

# Reader Function to Be Called for each Type of Format
# (Ganesh: Moved this here to remove the cyclic dependency)
GaussianSettings["writerFunctCall"] = {'geomInput':GaussianInputFileWriter}

