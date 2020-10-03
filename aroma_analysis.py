#! /usr/bin/env python
# Author : Anuja
# 05.01.2014
# Last Updated : 27.08.2014

# Utility for Analysis of AROMA-Run

import os
import sys
import string
import fpformat
import numpy

import aroma_util
from aroma_util import *
import aroma_constants
from aroma_constants import *

def readAndCheck(fl):
   lines = readFile(outdir + fl)
   if (len(lines) == 1):
      print "AROMA output " + fl + " is empty. Hence, aborting .."
      sys.exit(10)
   else:
      return lines

def processData(lines):
#  List of list of dist, oup, inp, iso, zz
   data_dict = []
   for i in range (1, len(lines)):
       words = map(string.atof, lines[i].split())
       data_dict.append([words[0], words[1], words[4], words[5], words[8]])
       
   return data_dict 

def analyse(mfile, sfile, dist_start = DEFAULT_DISTANCE_FOR_ANALYSIS, outfl = sys.stdout):

   numpy_flag = checkNumPy()
   if (not numpy_flag): sys.exit(10)

   if (len(sys.argv) > 3): dist_start = string.atof(sys.argv[3]) 
   
   mlines = readAndCheck(mfile)
   slines = readAndCheck(sfile)

   m_dict = processData(mlines) 
   s_dict = processData(slines) 

   for i in range (0, len(m_dict)):
      if (string.atof(m_dict[i][0]) >= string.atof(dist_start)):
          dist_start = i
          break;

   if (dist_start == DEFAULT_DISTANCE_FOR_ANALYSIS):
      outfl.write("\n\nWarning: Analysis for sigma-only model is performed for BQs beyond " + repr(DEFAULT_DISTANCE_FOR_ANALYSIS) + " angstrom.\n") 
      outfl.write("There are no BQs beyond that for this job. Therefore, Aroma can not perform analysis.\n\n")
      return

   dist = []; del_oup = []; del_inp = []; del_3iso = []; del_zz = []
   for i in range (dist_start, len(m_dict)):
       dist.append(m_dict[i][0])
       del_oup.append(m_dict[i][1] - s_dict[i][1])
       del_inp.append(m_dict[i][2] - s_dict[i][2])
       del_3iso.append(3*(m_dict[i][3] - s_dict[i][3]))
       del_zz.append(m_dict[i][4] - s_dict[i][4])
  
   for i in range (0, len(del_inp)):
       if (del_inp[i] > 5.0):
           print "Warning: For some points chosen for fitting the Doop and 3Diso data, the del-inp values exceeds 5.0"
           break;

   p_oup = numpy.poly1d(numpy.polyfit(dist, del_oup, 3))
   p_3iso = numpy.poly1d(numpy.polyfit(dist, del_3iso, 3))
   p_zz = numpy.poly1d(numpy.polyfit(dist, del_zz, 3))

   nics = (p_oup(1) + p_3iso(1))/2
   err = abs(nics - p_oup(1))

   outfl.write("\n--------------------------------------------------------------------\n")
   outfl.write("Polynomials for Doop and 3Diso are :\n")
   outfl.write("\n" + str(p_oup) + "\n" + str(p_3iso) + "\n")
   outfl.write("\nThe mean NICS value is " + repr(round(nics,3)) + " with error " + repr(round(err,3)))
   outfl.write("\n--------------------------------------------------------------------\n")

if __name__ == "__main__":
   analyse(sys.argv[1], sys.argv[2])

