#! /usr/bin/env python
# Author : Anuja
# 05.01.2014
# Last Updated : 27.08.2014

# Utility for Analysis of AROMA-Run

import os
import sys
import math
import string
import fpformat
import numpy

import aroma_util
from aroma_util import *
import aroma_constants
from aroma_constants import *

def readAndCheck(fl):
   lines = readFile(externalProgram["outdir"] + fl)
   if (len(lines) == 1):
      print("AROMA output " + fl + " is empty. Hence, aborting ..")
      sys.exit(10)
   else:
      return lines

def processData(lines):
#  List of list of dist, oup, inp, iso, zz
   data_dict = []
   for i in range (1, len(lines)):
       words = list(map(float, lines[i].split()))
       data_dict.append([words[0], words[1], words[4], words[5], words[8]])
       
   return data_dict 

def processPData(lines):
   data_dict = []
   for i in range (1, len(lines)):
       words = list(map(float, lines[i].split()))
       data_dict.append(words)
       
   return data_dict 

def analyse(mfile, sfile, dist_start = DEFAULT_DISTANCE_FOR_ANALYSIS, outfl = sys.stdout, ncs_flag = 0, pfile = ""):

   numpy_flag = checkNumPy()
   if (not numpy_flag): sys.exit(10)

   if (len(sys.argv) > 3): dist_start = float(sys.argv[3]) 
   
   mlines = readAndCheck(mfile)
   slines = readAndCheck(sfile)

   m_dict = processData(mlines) 
   s_dict = processData(slines)

   plines = []
   p_dict = []
   if (ncs_flag):
      plines = readAndCheck(pfile)
      p_dict = processPData(plines)

   for i in range (0, len(m_dict)):
      if (float(m_dict[i][0]) >= float(dist_start)):
          dist_start = i
          break

   if (dist_start == DEFAULT_DISTANCE_FOR_ANALYSIS):
      outfl.write("\n\nWarning: Analysis for sigma-only model is performed for BQs beyond " + repr(DEFAULT_DISTANCE_FOR_ANALYSIS) + " angstrom.\n") 
      outfl.write("There are no BQs beyond that for this job. Therefore, Aroma can not perform analysis.\n\n")
      return

   dist = []; del_oup = []; del_inp = []; del_3iso = []; del_zz = []; p_inp = []
   for i in range (dist_start, len(m_dict)):
       dist.append(m_dict[i][0])
       del_oup.append(m_dict[i][1] - s_dict[i][1])
       del_inp.append(m_dict[i][2] - s_dict[i][2])
       del_3iso.append(3*(m_dict[i][3] - s_dict[i][3]))
       del_zz.append(m_dict[i][4] - s_dict[i][4])

       if (ncs_flag):
           p_inp.append(p_dict[i][len(p_dict[i])-1])  # take the sum column
  
   for i in range (0, len(del_inp)):
       if (del_inp[i] > 5.0):
           print("Warning: For some points chosen for fitting the Doop and 3Diso data, the del-inp values exceeds 5.0")
           break

   p_oup = numpy.poly1d(numpy.polyfit(dist, del_oup, 3))
   p_3iso = numpy.poly1d(numpy.polyfit(dist, del_3iso, 3))
   p_zz = numpy.poly1d(numpy.polyfit(dist, del_zz, 3))

   if (ncs_flag):
      p_pinp = numpy.poly1d(numpy.polyfit(dist, p_inp, 3))

   nics = (p_oup(1) + p_3iso(1))/2
   err = abs(nics - p_oup(1))

   outfl.write("\n--------------------------------------------------------------------\n")
   outfl.write("Polynomials for Doop and 3Diso are :\n")
   outfl.write("\n" + str(p_oup) + "\n" + str(p_3iso) + "\n")
   outfl.write("\nThe mean NICS value is " + repr(round(nics,3)) + " with error " + repr(round(err,3)))
   
   if (ncs_flag):
      outfl.write("\n\nPicmo polynomial fits are:\n")
      outfl.write("\n" + str(p_pinp) + "\n")
      outfl.write("\nThe NICS value using CMO analysis is " + repr(round(p_pinp(1),3)) )

   outfl.write("\n--------------------------------------------------------------------\n")

def integralnics_analyse(mfile, sfile, pfile, dist_start = 2.0, outfl = sys.stdout):


   def fitbcr(dist, dat, name):

#      p_3iso = nicsIntegralFit(dist, dat)
      p_zz = nicsIntegralFit(dist, dat)

#      nicsint_val_3iso = ((p_3iso[0]*(p_3iso[1]**100))/math.log(p_3iso[1])) - (p_3iso[0]/math.log(p_3iso[1]))
      nicsint_val_zz = ((p_zz[0]*(p_zz[1]**100))/math.log(p_zz[1])) - (p_zz[0]/math.log(p_zz[1]))

#      nics = (nicsint_val_3iso + nicsint_val_zz)/2
#      err = abs(nicsint_val_3iso - nics)

#      outfl.write("\n--------------------------------------------------------------------\n")
      outfl.write("Parameters of the AB^r curve fitting for " + name + " are :\n")
      outfl.write("\n A = " + str(p_zz[0]) + "     B = " + str(p_zz[1]) + "\n")
      outfl.write("Integral NICS " + name + " = " + str(nicsint_val_zz) + "\n\n\n")
#      outfl.write("\n A = " + str(p_3iso[0]) + "     b = " + str(p_3iso[1]) + "\n")
#      outfl.write("Integral NICS 3Diso = " + str(nicsint_val_3iso) + "\n")
#      outfl.write("\nThe mean Integral-NICS value is " + repr(round(nics,3)) + " with error " + repr(round(err,3)))
#      outfl.write("\n--------------------------------------------------------------------\n")
      return p_zz, nicsint_val_zz

   numpy_flag = checkNumPy()
   if (not numpy_flag): sys.exit(10)

   if (len(sys.argv) > 3): dist_start = float(sys.argv[3]) 
   
   mlines = readAndCheck(mfile)
   m_dict = processData(mlines) 
   if (sfile != ""): 
       slines = readAndCheck(sfile)
       s_dict = processData(slines) 

   if (pfile != ""): 
       plines = readAndCheck(pfile)
       p_dict = processPData(plines) 

   for i in range (0, len(m_dict)):
      if (float(m_dict[i][0]) >= float(dist_start)):
          dist_start = i
          break

   dist = []; m_zz = []
   for i in range (dist_start, len(m_dict)):
       dist.append(m_dict[i][0])
       m_zz.append(m_dict[i][4])

   mp_zz, nicsint_val_zz = fitbcr(dist, m_zz, "NICS-ZZ")

   if (pfile != ""):
      c_zz = []
      last = len(p_dict[1]) - 1
      for i in range (dist_start, len(p_dict)):
          c_zz.append(p_dict[i][last])
      p_cmo, nicsint_val_cmo = fitbcr(dist, c_zz, "CMO-pi")


   if (sfile != ""):
      del_oup = []; del_inp = []; del_3iso = []; del_zz = []
      for i in range (dist_start, len(m_dict)):
          del_oup.append(m_dict[i][1] - s_dict[i][1])
          del_inp.append(m_dict[i][2] - s_dict[i][2])
          del_3iso.append(3*(m_dict[i][3] - s_dict[i][3]))
          del_zz.append(m_dict[i][4] - s_dict[i][4])
  
      for i in range (0, len(del_inp)):
          if (del_inp[i] > 5.0):
              print("Warning: For some points chosen for fitting the Doop and 3Diso data, the del-inp values exceeds 5.0")
              break

      p_3iso, nicsint_val_3iso = fitbcr(dist, del_3iso, "3Diso")
      p_dzz, nicsint_val_dzz = fitbcr(dist, del_zz, "DZZ")

      nics = (nicsint_val_3iso + nicsint_val_zz)/2
      err = abs(nicsint_val_3iso - nics)
      outfl.write("\nThe mean Integral-NICS value is " + repr(round(nics,3)) + " with error " + repr(round(err,3)))
      outfl.write("\n--------------------------------------------------------------------\n")

def analyse_ncs(pfile, outfl=sys.stdout, dist_start = DEFAULT_DISTANCE_FOR_ANALYSIS):

   numpy_flag = checkNumPy()
   if (not numpy_flag): sys.exit(10)

   if (len(sys.argv) > 3): dist_start = float(sys.argv[3]) 
   
   plines = []
   p_dict = []
   plines = readAndCheck(pfile)
   p_dict = processData(plines) 

   for i in range (0, len(m_dict)):
      if (float(m_dict[i][0]) >= float(dist_start)):
          dist_start = i
          break

   if (dist_start == DEFAULT_DISTANCE_FOR_ANALYSIS):
      outfl.write("\n\nWarning: Analysis for sigma-only model is performed for BQs beyond " + repr(DEFAULT_DISTANCE_FOR_ANALYSIS) + " angstrom.\n") 
      outfl.write("There are no BQs beyond that for this job. Therefore, Aroma can not perform analysis.\n\n")
      return

   dist = []; p_inp = []
   for i in range (dist_start, len(p_dict)):
       dist.append(p_dict[i][0])
       p_inp.append(p_dict[i][len(p_dict[i])-1])  # take the sum column
  
   p_pinp = numpy.poly1d(numpy.polyfit(dist, p_inp, 3))

   outfl.write("\n\nPicmo polynomial fits are:\n")
   outfl.write("\n" + str(p_pinp) + "\n")
   outfl.write("\nThe NICS value using CMO analysis is " + repr(round(p_pinp(1),3)) )

   outfl.write("\n--------------------------------------------------------------------\n")

if __name__ == "__main__":
   analyse(sys.argv[1], sys.argv[2])

