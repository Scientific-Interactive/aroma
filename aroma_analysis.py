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
import scipy

import aroma_util
from aroma_util import *
import aroma_constants
from aroma_constants import *
from scipy.optimize import OptimizeWarning

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
   outfl.write("\nThe mean NICS value is " + str(round(nics,3)) + " with error " + str(round(err,3)))
   
   if (ncs_flag):
      outfl.write("\n\nPicmo polynomial fits are:\n")
      outfl.write("\n" + str(p_pinp) + "\n")
      outfl.write("\nThe NICS value using CMO analysis is " + repr(round(p_pinp(1),3)) )

   outfl.write("\n--------------------------------------------------------------------\n")

def integralnics_analyse(mfile, sfile, pfile, dist_start = 2.0, outfl = sys.stdout):

   outfl.write("\n\n")

   def fitabcr(dist, dat, name):

      try: 
          p_zz_abc = nicsIntegralFit_abc(dist, dat)
          nicsint_val_zz_abc1 = (p_zz_abc[0]*p_zz_abc[1]) + p_zz_abc[2]
          nicsint_val_zz_abc17 = (p_zz_abc[0]*(p_zz_abc[1]**(1.7))) + p_zz_abc[2]
          return p_zz_abc, nicsint_val_zz_abc1, nicsint_val_zz_abc17, "Success", ""
      except OptimizeWarning as o:
#          outfl.write(repr(o) + "\nHence, following values might not be accurate")
          return p_zz_abc, nicsint_val_zz_abc1, nicsint_val_zz_abc17, "Success", "Parameter values might not be accurate as: " + repr(o)
      except RuntimeError as r:
#          outfl.write(repr(r) + "\nHence, skipping this fit")
          return 0.0, 0.0, 0.0, "Failure", "The parameters are set to zeros as: " + repr(r)
      except ValueError as v:
#          outfl.write(repr(v) + "\nHence, skipping this fit")
          return 0.0, 0.0, 0.0, "Failure", "The parameters are set to zeros as: " + repr(v)


   def fitbcr(dist, dat, name):

      try:
          p_zz_ab = nicsIntegralFit(dist, dat)
          nicsint_val_zz_ab = ((p_zz_ab[0]*(p_zz_ab[1]**100))/math.log(p_zz_ab[1])) - (p_zz_ab[0]/math.log(p_zz_ab[1]))
          return p_zz_ab, nicsint_val_zz_ab, "Success", "" 
      except OptimizeWarning as o:
#          outfl.write(repr(o) + "\nHence, following values might not be accurate")
          return p_zz_ab, nicsint_val_zz_ab, "Success", "The parameter values might not be accurate as: " + repr(o)
      except RuntimeError as r:
#          outfl.write(repr(r) + "\nHence, skipping this fit")
          return 0.0, 0.0, "Failure", "The parameters are set to zero as: " + repr(r)
      except ValueError as v:
#          outfl.write(repr(v) + "\nHence, skipping this fit")
          return 0.0, 0.0, "Failure", "The parameters are set to zero as: " + repr(v)


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

   mp_zz_abc, nicsint_zz_abc1, nicsint_zz_abc17, nicsint_zz_abc_status, nicsint_zz_abc_err = fitabcr(dist, m_zz, "NICS-ZZ")
   mp_zz_ab, nicsint_zz_ab, nicsint_zz_ab_status, nicsint_zz_ab_err = fitbcr(dist, m_zz, "NICS-ZZ")

   if (pfile != ""):
      c_zz = []
      last = len(p_dict[1]) - 1
      for i in range (dist_start, len(p_dict)):
          c_zz.append(p_dict[i][last])
      p_cmo_abc, nicsint_cmo_abc1, nicsint_cmo_abc17, nicsint_cmo_abc_status, nicsint_cmo_abc_err = fitabcr(dist, c_zz, "CMO-pi")
      p_cmo_ab, nicsint_cmo_ab, nicsint_cmo_ab_status, nicsint_cmo_ab_err = fitbcr(dist, c_zz, "CMO-pi")


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

      p_3iso_abc, nicsint_3iso_abc1, nicsint_3iso_abc17, nicsint_3iso_abc_status, nicsint_3iso_abc_err = fitabcr(dist, del_3iso, "3Diso")
      p_3iso_ab, nicsint_3iso_ab, nicsint_3iso_ab_status, nicsint_3iso_ab_err = fitbcr(dist, del_3iso, "3Diso")

      p_dzz_abc, nicsint_dzz_abc1, nicsint_dzz_abc17, nicsint_dzz_abc_status, nicsint_dzz_abc_err = fitabcr(dist, del_zz, "DZZ")
      p_dzz_ab, nicsint_dzz_ab, nicsint_dzz_ab_status, nicsint_dzz_ab_err = fitbcr(dist, del_zz, "DZZ")

      outfl.write("Integral NICS: NICS(r)=A*B^r	Integral-NICS= -A/lan(B)\n\n")

      if (nicsint_zz_ab_status == "Success"):
         if (nicsint_zz_ab_err != ""): outfl.write("\n" + nicsint_zz_ab_err + "\n")
         outfl.write(("A(zz) = " + str(round(mp_zz_ab[0],4)) + "\tB(zz) = " + str(round(mp_zz_ab[1],4)) + "\tIntegral-NICS(zz) = " + str(round(nicsint_zz_ab,4)) + "\n").expandtabs(24))
      else:
         outfl.write("For ZZ: " + nicsint_zz_ab_err + "\n")

      if (pfile != ""):
         if (nicsint_cmo_ab_status == "Success"):
             if (nicsint_cmo_ab_err != ""): outfl.write("\n" + nicsint_cmo_ab_err + "\n")
             outfl.write(("A(CMO) = " + str(round(p_cmo_ab[0],4)) + "\tB(CMO) = " + str(round(p_cmo_ab[1],4)) + "\tIntegral-NICS(CMO) = " + str(round(nicsint_cmo_ab,4)) + "\n").expandtabs(24))
         else:
             outfl.write("For CMO: " + nicsint_cmo_ab_err + "\n")

      if (sfile != ""):
         if (nicsint_dzz_ab_status == "Success"):
             if (nicsint_dzz_ab_err != ""): outfl.write("\n" + nicsint_dzz_ab_err + "\n")
             outfl.write(("A(delta-ZZ) = " + str(round(p_dzz_ab[0],4)) + "\tB(delta-ZZ) = " + str(round(p_dzz_ab[1],4)) + "\tIntegral-NICS(delta-ZZ) = " + str(round(nicsint_dzz_ab,4)) + "\n").expandtabs(24))
         else:
             outfl.write("For delta-ZZ: " + nicsint_dzz_ab_err + "\n")

         if (nicsint_3iso_ab_status == "Success"):
             if (nicsint_3iso_ab_err != ""): outfl.write("\n" + nicsint_3iso_ab_err + "\n")
             outfl.write(("A(3Diso) = " + str(round(p_3iso_ab[0],4)) + "\tB(3Diso) = " + str(round(p_3iso_ab[1],4)) + "\tIntegral-NICS(3Diso) = " + str(round(nicsint_3iso_ab,4)) + "\n").expandtabs(24))
         else:
             outfl.write("For 3Diso: " + nicsint_3iso_ab_err + "\n")

         if ((nicsint_dzz_ab_status == "Success") and (nicsint_3iso_ab_status == "Success")):
             if ((nicsint_dzz_ab_err != "") or (nicsint_3iso_ab_err != "")): outfl.write("\n" + "Following value may not be accurate as the covariance of parameters could not be estimated." + "\n")
             sonly_mean = meanOf([nicsint_dzz_ab, nicsint_3iso_ab])
             outfl.write("\n                Sigma-only Integral-NICS(pi, zz) = " + str(round(sonly_mean,4)) + " \u00B1 " + str(round((sonly_mean - nicsint_dzz_ab),4)) + "\n")

      outfl.write("\nNICS(r)=A*B^r +C; NICS(1)=A*B+C	NICS(1.7)=A*B^1.7 +C\n\n")

      if (nicsint_zz_abc_status == "Success"):
         if (nicsint_zz_abc_err != ""): outfl.write("\n" + nicsint_zz_abc_err + "\n")
         outfl.write(("A(zz) = " + str(round(mp_zz_abc[0],4)) + "\tB(zz) = " + str(round(mp_zz_abc[1],4)) + "\tC(zz) = " + str(round(mp_zz_abc[2],4)) + "\tNICS(1) = " + str(round(nicsint_zz_abc1,4)) + "\tNICS(1.7) = " + str(round(nicsint_zz_abc17,4)) + "\n").expandtabs(24))
      else:
         outfl.write("For ZZ: " + nicsint_zz_abc_err + "\n")

      if (pfile != ""):
         if (nicsint_cmo_abc_status == "Success"):
             if (nicsint_cmo_abc_err != ""): outfl.write("\n" + nicsint_cmo_abc_err + "\n")
             outfl.write(("A(CMO) = " + str(round(p_cmo_abc[0],4)) + "\tB(CMO) = " + str(round(p_cmo_abc[1],4)) + "\tC(CMO) = " + str(round(p_cmo_abc[2],4)) + "\tNICS(1) = " + str(round(nicsint_cmo_abc1,4)) + "\tNICS(1.7) = " + str(round(nicsint_cmo_abc17,4)) + "\n").expandtabs(24))
         else:
             outfl.write("For CMO: " + nicsint_cmo_abc_err + "\n")

      if (sfile != ""):
         if (nicsint_dzz_abc_status == "Success"):
             if (nicsint_dzz_abc_err != ""): outfl.write("\n" + nicsint_dzz_abc_err + "\n")
             outfl.write(("A(delta-zz) = " + str(round(p_dzz_abc[0],4)) + "\tB(delta-zz) = " + str(round(p_dzz_abc[1],4)) + "\tC(delta-zz) = " + str(round(p_dzz_abc[2],4)) + "\tNICS(1) = " + str(round(nicsint_dzz_abc1,4)) + "\tNICS(1.7) = " + str(round(nicsint_dzz_abc17,4)) + "\n").expandtabs(24))
         else:
             outfl.write("For delta-ZZ: " + nicsint_dzz_abc_err + "\n")

         if (nicsint_3iso_abc_status == "Success"):
             if (nicsint_3iso_abc_err != ""): outfl.write("\n" + nicsint_3iso_abc_err + "\n")
             outfl.write(("A(3Diso) = " + str(round(p_3iso_abc[0],4)) + "\tB(3Diso) = " + str(round(p_3iso_abc[1],4)) + "\tC(3Diso) = " + str(round(p_3iso_abc[2],4)) + "\tNICS(1) = " + str(round(nicsint_3iso_abc1,4)) + "\tNICS(1.7) = " + str(round(nicsint_3iso_abc17,4)) + "\n").expandtabs(24))
         else:
             outfl.write("For 3Diso: " + nicsint_3iso_abc_err + "\n")

         if ((nicsint_dzz_ab_status == "Success") and (nicsint_3iso_ab_status == "Success")):
             if ((nicsint_dzz_ab_err != "") or (nicsint_3iso_ab_err != "")): outfl.write("\n" + "Following value may not be accurate as the covariance of parameters could not be estimated." + "\n")
             sonly_mean1 = meanOf([nicsint_dzz_abc1, nicsint_3iso_abc1])
             sonly_mean17 = meanOf([nicsint_dzz_abc17, nicsint_3iso_abc17])
             outfl.write("\n                Sigma-only Sigma-only NICS(1) = " + str(round(sonly_mean1,4)) + " \u00B1 " + str(round((sonly_mean1 - nicsint_dzz_abc1),4)) + "\n")
             outfl.write("\n                Sigma-only Sigma-only NICS(1.7) = " + str(round(sonly_mean17,4)) + " \u00B1 " + str(round((sonly_mean17 - nicsint_dzz_abc17),4)) + "\n")
      outfl.write("\n--------------------------------------------------------------------\n")

def analyse_ncs(pfile, outfl=sys.stdout, dist_start = DEFAULT_DISTANCE_FOR_ANALYSIS):

   numpy_flag = checkNumPy()
   if (not numpy_flag): sys.exit(10)

   if (len(sys.argv) > 3): dist_start = float(sys.argv[3]) 
   
   plines = []
   p_dict = []
   plines = readAndCheck(pfile)
   p_dict = processPData(plines) 

   for i in range (0, len(p_dict)):
      if (float(p_dict[i][0]) >= float(dist_start)):
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

