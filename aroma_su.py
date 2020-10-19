#! /usr/bin/env python

# Author : Anuja
# 29.04.2014
# Last Updated : 04.05.2014

# Super Script as in Wrapper to Call Aroma to run Multiple Jobs Through Single Command

import os
import sys
import time
import string

import aroma
from aroma import *
def main(suarmfile):

#   flprfx = suarmfile[0:suarmfile.rindex(".")]
   flprfx = suarmfile
   f = open(suarmfile + ".suarm", "r")
   suarmlines = f.readlines()
   f.close()

   common_str = ""
   for i in range (0, len(suarmlines)):
      if (suarmlines[i].upper().find("COMMON") >= 0 ): break;
   for j in range (i+1, len(suarmlines)):
      if (suarmlines[j].upper().find("COMMON END") >= 0 ): break;
      common_str += suarmlines[j].strip() + "\n"

   run_list = []
   for i in range (0, len(suarmlines)):
      if (suarmlines[i].upper().find("RUNNAME") >= 0 ):
         run_name = suarmlines[i].strip().split()[1]
         run_str = ""
         j = i+1
         while ( (j < len(suarmlines)) and (suarmlines[j].upper().find("RUNNAME") < 0)):
            run_str += suarmlines[j].strip() + "\n"
            j += 1
         run_list.append((run_name, common_str + run_str))
   
   for r in run_list:
      key = r[0]
      f = open (flprfx + "-" + key + ".arm", "w")
      f.write(r[1])
      f.close()

   for r in run_list:
          key = r[0]
          print(flprfx + "-" + key)
          aroma(flprfx + "-" + key)
#      try:
#          print(flprfx + "-" + key)
#          aroma(flprfx + "-" + key)
#      except:
#          print("Some Error Occured While Performing The Run Named " + key)
#          print("Moving to The Next Run ..")
#          time.sleep(10)

if __name__ == "__main__":
   main(sys.argv[1])
