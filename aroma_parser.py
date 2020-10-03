#! /usr/bin/env python
# Author : Anuja
# @ 07.02.2013
# Last Updated : 23.02.2014

# Classes for Parsing Different Inputs (com, log, chk)

import os
import sys
import string

import aroma_util
from aroma_util import *

import aroma_constants
from aroma_constants import *

class FileParser:
   def __init__(self, geomfl):
      self.geomfl = geomfl
      self.geom = {}
      self.hashLine = aroma_constants.DEFAULT_NICS_KEYLINE 
      self.title = "DEFAULT TITLE SET BY AROMA"
      self.charge = 0
      self.mult = 0

   def getInpData(self):
      pass

class InputFileParser(FileParser):
   def __init__(self, geomfl):
      FileParser.__init__(self, geomfl)

   def getInpData(self):
      glines = readFile(self.geomfl)

      for i in range (0, len(glines)):
        if (glines[i].find("#") >= 0):
           if (len(glines[i+4].split()) == 2):
              break;

      self.hashLine = glines[i].strip()
      self.title = glines[i+2].strip()
      self.charge, self.mult = map(string.atoi, glines[i+4].strip().split())

      # Check if its in z-matrix format or cartesian coordinates
      zmat_flag = 0
      if ( (len(glines[i+5].strip().split()) == 1) and (len(glines[i+6].strip().split()) == 3) and (len(glines[i+7].strip().split()) == 5) ): zmat_flag = 1

      if (zmat_flag):
         # Run a guess=only run to convert Z-matrx to Cartesian coordinates 
         # so that I dont have to write a separate function to do that ;)
         self.fl = self.geomfl[self.geomfl.rindex("/")+1:len(self.geomfl)+1] # chops off the path string
         self.newflprfx =  self.fl[0:self.fl.rindex(".")] + "-guessonly" # chops off the extension       
         f = open(self.geomfl[0:self.geomfl.rindex("/")+1] + self.newflprfx + aroma_constants.GaussInpExt, "w") 
         f.write("# guess=only\n")
         for j in range (i+1, len(glines)):
            f.write(glines[j] )
         f.close()
         status = os.system(aroma_constants.constructGaussCMD(self.newflprfx))
         if (status):
            print "WARNING: An error occured while running a Guess=Only job for converting Z-matrix to cartesian coordinates.\nConsequently, Aroma may terminate abnormally."
#            print "Aborting Aroma .."
#            sys.exit(10)
         
         theParser = OutputFileParser(aroma_constants.outdir + self.newflprfx + aroma_constants.GaussOutExt)
         self.geom, self.hashLine, self.title, self.charge, self.mult = theParser.getInpData()

      else:
         dig_flag = 1
         if (glines[i+5].strip().split()[0].isdigit() != 1): dig_flag = 0

         nat = 0

         if (glines[i+5].count(",") == 0):
            for j in range (i+5, len(glines)):

               words = glines[j].split()
               if ( (words[0].upper() != "BQ") and (words[0] != 0) ):
                  nat += 1
                  self.geom[nat] = []
                  if dig_flag : self.geom[nat].append(string.atoi(words[0]))
                  else: self.geom[nat].append(aroma_constants.AtmSym[words[0].upper()])

                  for k in range (1,4):
                     self.geom[nat].append(round(string.atof(words[k]),6))

               if (glines[j+1] == "\n"):  break;
         else:
            for j in range (i+5, len(glines)):

               words = glines[j].split(",")
               if ( (words[0].upper() != "BQ") and (words[0] != 0) ):
                  nat += 1
                  self.geom[nat] = []
                  if dig_flag : self.geom[nat].append(string.atoi(words[0].strip()))
                  else: self.geom[nat].append(aroma_constants.AtmSym[words[0].strip().upper()])

                  for k in range (1,4):
                     self.geom[nat].append(round(string.atof(words[k].strip()),6))

               if (glines[j+1] == "\n"):  break;

      return self.geom, self.hashLine, self.title, self.charge, self.mult

class OutputFileParser(FileParser):
   def __init__(self, geomfl):
      FileParser.__init__(self, geomfl)

   def getInpData(self):
      glines = readFile(self.geomfl)

      for i in range (0, len(glines)):
        if (glines[i].find("#") >= 0):
           if (glines[i+4].find("-----")):
              break;

      self.hashLine = glines[i].strip()

      for j in range (i+1, len(glines)):
         if ( (glines[j].find("Charge") >= 0) and (glines[j].find("Multiplicity") >= 0) ):
            break;

      self.title = glines[j-2].strip()
      words = glines[j].strip().split()
      self.charge = string.atoi(words[2])
      self.mult = string.atoi(words[5])

      # If its an output of an optimization run, then the last goemetry should be read
      # Therefore, here a reverse loop is necessary
      for i in range (len(glines)-1, -1, -1):
         if (glines[i].upper().find("ORIENTATION") >= 0 ): break;

      nat = 0
      for j in range (i+5, len(glines)):
         if ( glines[j].find("-----") >= 0 ): break;
         
         words = glines[j].split()
         nat += 1
         self.geom[nat] = []
         self.geom[nat].append(string.atoi(words[1]))

         for k in range (3,6):
            self.geom[nat].append(round(string.atof(words[k]),6))

      return self.geom, self.hashLine, self.title, self.charge, self.mult

class ChkFileParser(FileParser):
   def __init__(self, geomfl):
      FileParser.__init__(self, geomfl)

   def getInpData(self):
      os.system(FormChk_CMD + self.geomfl + " " + self.geomfl + ".fchk")
      glines = readFile(self.geomfl + ".fchk")

      self.title = glines[0].strip()
      type = glines[1][0:10].strip()
      method = glines[1][10:40].strip()
      basis = glines[1][40:70].strip()
      self.hashLine = "# " + method + "/" + basis + " " + type

      for i in range (0, len(glines)):
         if (glines[i].find("Charge") >= 0 ): self.charge = string.atoi(glines[i].split()[2])
         if (glines[i+1].find("Multiplicity") >= 0 ):
            self.mult = string.atoi(glines[i].split()[2])
            break;

      for i in range (0, len(glines)):
         if (glines[i].find("Nuclear charges") >= 0 ): break;

      nat = 0
      for j in range (i+1, len(glines)):
         if (glines[j].find("Current cartesian coordinates") >= 0): break;
         words = map(string.atof, glines[j].split())
         for k in range (0, len(words)):
            nat += 1
            self.geom[nat] = []
            self.geom[nat].append(int(words[k]))

      all_coords = []
      str = ""
      for l in range (j+1, len(glines)):
         if (glines[j].find("Force Field") >= 0): break;
         str += glines[j].strip() + "   "

      nat = 1
      words = map(string.atof, str.split())
      for k in range (0, len(words)):
         if (len(self.geom[nat]) <= 4):
            nat += 1
         # Convert to angstrom from atomic unit
         self.geom[nat].append(round(words[k]*0.529177249,6))


      return self.geom, self.hashLine, self.title, self.charge, self.mult
