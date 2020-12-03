#! /usr/bin/env python
# Author : Ganesh
# 03.12.2020
# Read external user modified constants file

import os
import sys
import importlib

from pathlib import Path

def loadUserConstants(loadDir):
  loadDir = os.path.dirname(loadDir)
  userConstantsModule = "user_aroma_constants"
  userConstantsFile = os.path.join(loadDir, userConstantsModule + ".py")
  print(userConstantsFile)

  if not Path(userConstantsFile).exists():
    print("No user specific config files available")
    return 

  try:
    del sys.modules['aroma_constants']
  except: 
    pass

  sys.path.append(loadDir)
  sys.modules['aroma_constants'] = importlib.import_module(userConstantsModule)
       

if __name__=="__main__":
  loadUserConstants(sys.argv[0])
