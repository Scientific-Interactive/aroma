#! /usr/bin/env python
# Author : Anuja
# 30.01.2013
# Last Updated : 14.10.2014

# Main Driver scipt for automation of NICS-Scan in Z and XY Directions, Sigma-Only Model and CMO-NICS 

import os
import re
import sys
import string
import fpformat

import aroma_constants
from aroma_constants import *

import aroma_util
from aroma_util import *

import aroma_parser
from aroma_parser import *

import aroma_molecule
from aroma_molecule import *

import aroma_pinbo
from aroma_pinbo import *

import aroma_analysis
from aroma_analysis import *

import aroma_ringarea
from aroma_ringarea import *

def init():
   # global flags
   global opt_flag, ncs_flag, sigma_flag, analyse_flag, area_flag, s_charge_flag, opt_external, optfl_external, xy_flag
   # global molecule-related
   global CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, xy_ref_ring_info, points, normals, exocyclic
   # global technical 
   global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, sigma_charge, analyse_dist, clear_flag, BQGuide, xy_BQ_dist

   # Initializing some global variables
   opt_flag = 0; ncs_flag = 0; sigma_flag = 0; opt_external = 0; xy_flag = 0; analyse_flag = 1; area_flag = 0
   CenterOf = {}; all_aromatic_rings = {}; exocyclic = {}; points = {}; normals = {}
   geomflext = ""; geomfl = ""; flprfx = ""; outfilename = ""
   sigma_direction = 'POSITIVE'
   n_xy_center = {1:1}; xy_ref_ring_info =[]; BQGuide = {}; xy_BQ_dist = []
   runtype = "NICSSCAN"; hashLine_opt = DEFAULT_OPTIMIZATION_KEYLINE; hashLine_nics = DEFAULT_NICS_KEYLINE; hashLine_ncs = DEFAULT_NCS_KEYLINE; hashLine_nbo = DEFAULT_NBO_KEYLINE
   BQ_Step = DEFAULT_BQ_STEP; BQ_Range = DEFAULT_BQ_RANGE
   sigma_charge = 0; s_charge_flag = 0
   analyse_dist = DEFAULT_DISTANCE_FOR_ANALYSIS 
   clear_flag = 0

def check(armfile):

   # global flags
   global opt_flag, ncs_flag, sigma_flag, xy_flag, analyse_flag, area_flag, s_charge_flag, opt_external, optfl_external
   # global molecule-related
   global CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, BQGuide, points, normals
   # global technical 
   global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, sigma_charge, analyse_dist, clear_flag

   armpath = armfile[0:armfile.rindex("/")+1] 
   flprfx = armfile[armfile.rindex("/")+1:len(armfile)]
   armlines = readFile(armpath + flprfx + ".arm")

   outfilename = outdir + flprfx + ".armlog"

   # Check for Validity of the RunType
   for i in range (0, len(armlines)):
      if (armlines[i].upper().find("RUN") >= 0): 
         runseq = re.split("[ |,]", armlines[i].upper().strip().split("=")[1])
         if (runseq.count("OPT") > 0): opt_flag = 1
         if (runseq.count("NCS") > 0): ncs_flag = 1
         if (runseq.count("SIGMA") > 0): sigma_flag = 1
         if (runseq.count("XY") > 0): xy_flag = 1

   if (opt_flag):
      for i in range (0, len(armlines)):
         if (armlines[i].upper().find("OPT_EXTERNAL") >= 0):
            opt_external = 1
            optfl_external = armlines[i].strip().split("=")[1]

   # Check for input file for Geomtry
   if (not opt_external):
      for i in range (0, len(armlines)):
         if (armlines[i].upper().find("GEOMFILE") >= 0):
             geomfl = armlines[i].strip().split("=")[1]
             break;

      geomfl_flag = os.path.exists(geomfl)
      if (not geomfl_flag):
         print(geomfl + " Could Not be Found.\nTherefore, Aborting the Run ..")
         sys.exit(10)

      geomflext = geomfl[geomfl.rindex(".")+1:len(geomfl)+1]

      valid_ext_flag = 1
      for extension in EXTENSIONS_FOR_GAUSSIAN_FILES:
         if (EXTENSIONS_FOR_GAUSSIAN_FILES[extension].count(geomflext) > 0): valid_ext_flag = 1
      if (not valid_ext_flag): print("Gaussian File with \"" + geomflext + "\" Extension Can Not be Read.\nTherefore, Aborting the Run .."); sys.exit(10)

   if (xy_flag): print("\nWARNING: You Have Requested XY-Scan. Make Sure That the Centers Are Defined in Proper Order.\n")

   for i in range (0, len(armlines)):
      if (armlines[i].upper().find("OUTFILE") >= 0): outfilename = armlines[i].strip().split("=")[1]

   # Get The Ring/Bond Info
   r_count = 0
   n_count = 0
   for i in range (0, len(armlines)):
      if (armlines[i].upper().find("CENTER") >= 0):
         r_count += 1
         CenterOf[r_count] = list(map(int, re.split("[:|=]", armlines[i].strip())[1].split(",")))
      if (armlines[i].upper().find("NORMAL") >= 0):
         n_count += 1
         normals[n_count] = list(map(float, re.split("[:|=]", armlines[i].strip())[1].split(",")))
      if (xy_flag):
         if (armlines[i].upper().find("POINT") >= 0):
            points[(r_count, r_count+1)] = list(map(float, re.split("[:|=]", armlines[i].strip())[1].split(",")))
            points[(r_count, r_count+1)].insert(0,-1)
            if (len(points[r_count, r_count+1]) != 4): print("The keyword POINT should have X, Y, Z coordinates.\n Check and Submit Again. Aborting This Run ..\n"); sys.exit(10)

   for i in range (0, len(armlines)):
      if (armlines[i].upper().find("BQSTEP") >= 0): BQ_Step = float(armlines[i].split("=")[1]); break;
   for i in range (0, len(armlines)):
      if (armlines[i].upper().find("BQRANGE") >= 0): BQ_Range = list(map(float, armlines[i].split("=")[1].split(","))); break;

   BQ_No = int((BQ_Range[1] - BQ_Range[0])/BQ_Step)
   if ( ncs_flag and (BQ_No > 100) ):
      print("The package NBO can not handle more than 100 Bqs. Aborting thr Run.\nPlease change the BQRANGE or BQSTEP and Resubmit.")
      sys.exit(10)

   if ((len(CenterOf) < 1) and (len(normals) < 1)):
      print("Rings/Bonds Are Not Defined.\nTherefore, Aborting the Run ..")
      sys.exit(10)


   # Check keywords for Gaussian for optimization and nics runs
   for i in range (0, len(armlines)):
      if (armlines[i].upper().find("KEYLINES") >= 0): break;
   for j in range (i+1, len(armlines)):
         if (armlines[j].upper().find("TEMPLATE") >= 0): 
            if (armlines[j].upper().split("_")[0] == 'OPT'):
                hashLine_opt = ""
                for k in range (j+1, len(armlines)):
                   if ( (armlines[k].upper().find("TEMPLATE") >= 0) or (armlines[k].upper().find("END KEYLINE") >= 0) ): break;
                   hashLine_opt += armlines[k]
               
            if (armlines[j].upper().split("_")[0] == 'NICSSCAN'):
                hashLine_nics = ""
                for k in range (j+1, len(armlines)):
                   if ( (armlines[k].upper().find("TEMPLATE") >= 0) or (armlines[k].upper().find("END KEYLINE") >= 0) ): break;
                   hashLine_nics += armlines[k].upper()

            if (armlines[j].upper().split("_")[0] == 'NCS'):
                hashLine_ncs = ""
                for k in range (j+1, len(armlines)):
                   if ( (armlines[k].upper().find("TEMPLATE") >= 0) or (armlines[k].upper().find("END KEYLINE") >= 0) ): break;
                   hashLine_ncs += armlines[k].upper()

            if (armlines[j].upper().split("_")[0] == 'NBO'):
                hashLine_nbo = ""
                for k in range (j+1, len(armlines)):
                   if ( (armlines[k].upper().find("TEMPLATE") >= 0) or (armlines[k].upper().find("END KEYLINE") >= 0) ): break;
                   hashLine_nbo += armlines[k].upper()

   for i in range (0, len(armlines)):
      if (armlines[i].upper().find("ANALYSE") >= 0): 
         if (armlines[i].upper().find("NO") >= 0): analyse_flag = 0; break;
         elif (armlines[i].find("=") >= 0): analyse_dist = float(armlines[i].split("=")[1])
         elif (armlines[i].upper().find("AREA") >= 0): area_flag = 1

   for i in range (0, len(armlines)):
      if (armlines[i].upper().find("CLEAR") >= 0): clear_flag = 1

   if (sigma_flag):
      for i in range (0, len(armlines)):
         if (armlines[i].upper().find("AROMATIC RING") >= 0): break;

      rings_count = 0
      for j in range (i+1, len(armlines)):
         if (armlines[j].upper().find("END") >= 0): break;
         rings_count += 1
         all_aromatic_rings[rings_count] = list(map(int, armlines[j].strip().split(",")))

      if (all_aromatic_rings == {}):
         print("\nWARNING: Aromatic Rings Are Not Provided, Therefore, Sigma-Only Model Calculations Will Not Be Performed.\n")
         sigma_flag = 0

      for i in range (0, len(armlines)):
         if (armlines[i].upper().find("EXOCYCLIC") >= 0): break;
      ec_count = 0
      for j in range (i+1, len(armlines)):
         if (armlines[j].upper().find("END") >= 0): break;
         ec_count += 1
         exocyclic[ec_count] = list(map(int, armlines[j].strip().split(",")))

      for i in range (0, len(armlines)):
         if (armlines[i].upper().find("DIRECTION") >= 0):
            words = armlines[i].upper().split()
            if ( (words[2] == 'POSITIVE' ) or (words[2] == 'NEGATIVE') ):
               sigma_direction = words[2]

      for i in range (0, len(armlines)):
         if (armlines[i].upper().find("SONLY CHARGE") >= 0):
            s_charge_flag = 1
            sigma_charge = int(armlines[i].split("=")[1])
   

def generate_Opt_Input(geom, hashLine, title, charge, mult):
    global flprfx, inpdir

    # The chk file name for optimization can be same as given by the user.
    # Remove the above lines after testing this.
    hashLine_rev = hashLine

    # Generate input file for optimization
    optfl = flprfx + "-opt" 
    f_opt = open(inpdir + optfl + GaussInpExt, "w")

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

def run_Optimization(optfl):
  
    # Run Gaussian for optimization
    # Check status and print approprate msg before proceeding 
    print("Status: Optimization Running .. ")
    status = execCmd(constructGaussCMD(optfl))
    status = 0
    print("Status: Optimization Over.")
    if (status) : 
       print("\nWARNING: Abnormal Termination of Optimization Run.")
       print("Aroma Will Continue with the Last Geometry for NICS Calculation.\n")
    


def genNicsInputs(geom, Conn, hashLine, title, charge, mult):

   # global flags
   global opt_flag, ncs_flag, sigma_flag, xy_flag, analyse_flag, area_flag, s_charge_flag, opt_external, optfl_external
   # global molecule-related
   global CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, points, normals
   # global technical 
   global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, sigma_charge, analyse_dist, clear_flag

   hashLine_rev = ''
   flag_chk = 0
   hlines = hashLine.upper().split("\n")
   for i in range (0, len(hlines)):
      if (hlines[i].find("CHK") >= 0 ):
         flag_chk = 1
      else: hashLine_rev += hlines[i] + "\n"

   if (not xy_flag):
      for ring in CenterOf:
         ring_atoms = CenterOf.get(ring)

         # The Plane is always fixed to be XY and the molecule is oriented in such a way.
         new_geom, new_points, new_normal = reorient(geom, Conn, ring_atoms)

         if (new_geom != []):
            BQs_string = generateBQs(new_geom, Conn, ring_atoms)
    
            coord_format = "{0:.5f}"
            ringf = open(inpdir + flprfx + "-center" + repr(ring) + GaussInpExt, "w")
            if (flag_chk): ringf.write("%chk=" + chkdir + flprfx + "-center" + repr(ring) + ".chk\n")
            ringf.write(hashLine_rev + title + " # Center " + repr(ring) + "\n\n" + repr(charge) + " " + repr(mult) + "\n")
            for i in range (1, len(geom)+1):
               # The dummy atoms are considered as BQs, therefore, remove them
               if (new_geom[i][0] != 0):
                  geomline = repr(new_geom[i][0]) + "   " + coord_format.format(new_geom[i][1]) + "   " + coord_format.format(new_geom[i][2]) + "   " + coord_format.format(new_geom[i][3]) + "\n"
                  ringf.write(geomline)

            ringf.write(BQs_string + "\n")

            # If NCS run is requested, then add NBO keywords at the end of the Gaussian input
            if (ncs_flag):
               ringf.write(hashLine_nbo + "\n")
   
            ringf.close()
         elif (new_geom == []):
            print("Ring no. " + repr(ring) + " is not Planar within the tolerence of " + repr(TORSION_ANGLE_TOLERANCE)  +" degrees. Therefore, ring could not be reoriented in XY plane and BQs could not be generated.")

      # For "Normal"s defined by user
      n_count = len(CenterOf)
      for norm in normals:
         n_count += 1
         ring_atoms = []
         new_geom, new_points, new_normal = reorient(geom, Conn, ring_atoms, points, normals[norm])

         if ((new_geom != []) and (new_normal != [])):
            ring_atoms = []
            BQs_string = generateBQs(new_geom, Conn, ring_atoms, new_normal)
    
            coord_format = "{0:.5f}"
            ringf = open(inpdir + flprfx + "-center" + repr(n_count) + GaussInpExt, "w")
            if (flag_chk): ringf.write("%chk=" + chkdir + flprfx + "-center" + repr(n_count) + ".chk\n")
            ringf.write(hashLine_rev + title + " # Center " + repr(n_count) + "\n\n" + repr(charge) + " " + repr(mult) + "\n")
            for i in range (1, len(geom)+1):
               # The dummy atoms are considered as BQs, therefore, remove them
               if (new_geom[i][0] != 0):
                  geomline = repr(new_geom[i][0]) + "   " + coord_format.format(new_geom[i][1]) + "   " + coord_format.format(new_geom[i][2]) + "   " + coord_format.format(new_geom[i][3]) + "\n"
                  ringf.write(geomline)

            ringf.write(BQs_string + "\n")

            # If NCS run is requested, then add NBO keywords at the end of the Gaussian input
            if (ncs_flag):
               ringf.write(hashLine_nbo + "\n")
   
            ringf.close()
         elif (new_geom == []):
            print("Ring no. " + repr(n_count) + " is not Planar within the tolerence of " + repr(TORSION_ANGLE_TOLERANCE)  +" degrees. Therefore, ring could not be reoriented in XY plane and BQs could not be generated.")

   elif (xy_flag):
      BQs_string, new_geom = generateBQs_XY(geom, Conn)
      n_fl = BQs_string.count('break') + 1
      BQs_strings = BQs_string.split("break")

      for c in range (1, n_fl+1):
         if (c > 1): n_xy_center[c] = c
         coord_format = "{0:.5f}"
         ringf = open(inpdir + flprfx + "-center" + repr(c) + GaussInpExt, "w")
         if (flag_chk): ringf.write("%chk=" + chkdir + flprfx + "-center" + repr(c) + ".chk\n")
         ringf.write(hashLine_rev + title + " # Center " + repr(c) + "\n\n" + repr(charge) + " " + repr(mult) + "\n")

         for i in range (1, len(new_geom)+1):
            # The dummy atoms are considered as BQs, therefore, remove them
            if (new_geom[i][0] != 0):
               geomline = repr(new_geom[i][0]) + "   " + coord_format.format(new_geom[i][1]) + "   " + coord_format.format(new_geom[i][2]) + "   " + coord_format.format(new_geom[i][3]) + "\n"
               ringf.write(geomline)

         ringf.write(BQs_strings[c-1] + "\n")

         # If NCS run is requested, then add NBO keywords at the end of the Gaussian input
         if (ncs_flag):
            ringf.write(hashLine_nbo + "\n")
   
         ringf.close()


# Here, the geom refers to the re-oriented geometry
# Now, the plane is always XY plane
def generateBQs(geom, Conn, ring_atoms, normal = []):

   sigma_model = 0
   direction_bq = 'POSITIVE'
   direction = 0
   H_count = 0

   for i in range (0, len(ring_atoms)):
      atm_idx = ring_atoms[i]

      if (not sigma_model):
         for j in range (0, len(Conn[atm_idx])):
            if (geom[Conn[atm_idx][j]][0] != 1):
               continue
            else: 
               if ( (geom[Conn[atm_idx][j]][3]) > 0.5):
                  H_count += 1
                  direction += 1
               elif ( (geom[Conn[atm_idx][j]][3]) < -0.5):
                  H_count += 1
                  direction -= 1

# For Normal, the detectin of Sigma-only Model is under Trial
   if ((len(ring_atoms) < 1) and (normal != [])):
       for i in range (0, len(geom)):
          atm_idx = i+1
          if (not sigma_model):
             for j in range (0, len(Conn[atm_idx])):
                if (geom[Conn[atm_idx][j]][0] != 1):
                   continue
                else: 
                   if ( (geom[Conn[atm_idx][j]][3]) > 0.5):
                      H_count += 1
                      direction += 1
                   elif ( (geom[Conn[atm_idx][j]][3]) < -0.5):
                      H_count += 1
                      direction -= 1

   if (H_count > 2): sigma_model = 1
   if (direction > 0): direction_bq = 'NEGATIVE'

   if (len(ring_atoms) != 0): cmx, cmy, cmz = getGMOfRing(geom, ring_atoms)
   else: cmx, cmy, cmz = 0.0, 0.0, 0.0

   if (sigma_model):
     print("\nSigma-Only Model detected. \nTherefore, the BQs will be generated on the opposite side of the s-only H-atoms.")

   zinc = BQ_Step
   if (direction_bq == 'NEGATIVE'): zinc = -zinc

   coord_format = "{0:.5f}"
   BQs_string = ""
   zcoord = cmz + BQ_Range[0] 
   for i in range (0, BQ_No):
      BQs_string += 'bq     ' + coord_format.format(cmx) + "     " + coord_format.format(cmy) + "     " + coord_format.format(zcoord) + "\n"
      zcoord += zinc

   return BQs_string

def generateBQs_XY(geom, Conn):

   # global flags
   global opt_flag, ncs_flag, sigma_flag, xy_flag, analyse_flag, area_flag, s_charge_flag, opt_external, optfl_external
   # global molecule-related
   global CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, xy_ref_ring_info, BQGuide, points
   # global technical 
   global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, xy_BQ_dist, sigma_charge, analyse_dist, clear_flag

   direct = ''
   sigma_model = 0
   direction_bq = 'POSITIVE'
   direction = 0
   H_count = 0
   norm_vec = {}
   ring_count = 0

   for ring in CenterOf:
       ring_atoms = getOrderedRing(Conn, CenterOf.get(ring))
       if (len(ring_atoms) > 3): 
           new_geom, new_points, new_normal = reorient(geom, Conn, ring_atoms, points)

           if (new_geom != []):
              for i in range (0, len(ring_atoms)):
                 atm_idx = ring_atoms[i]

                 if (not sigma_model):
                    for j in range (0, len(Conn[atm_idx])):
                       if (new_geom[Conn[atm_idx][j]][0] != 1):
                          continue
                       else: 
                          if ( (new_geom[Conn[atm_idx][j]][3]) > 0.5):
                             H_count += 1
                             direction += 1
                          elif ( (new_geom[Conn[atm_idx][j]][3]) < -0.5):
                             H_count += 1
                             direction -= 1
              if (H_count > 2): sigma_model = 1; print("\nSigma-Only Model detected. \nTherefore, the BQs will be generated on the opposite side of the s-only H-atoms.")
              if (direction > 0): direction_bq = 'NEGATIVE'

              cmx, cmy, cmz = getGMOfRing(new_geom, ring_atoms)
              ref_ring = ring
              xy_ref_ring_info.append(ring)
              xy_ref_ring_info.append(sorted(ring_atoms))
              normal_to_Ring = getUnitVector(getAverageNormaltoTheRing(new_geom, ring_atoms, [cmx, cmy, cmz]))
              normal_to_Ring[:] = [vc*DEFAULT_XY_DISTANCE for vc in normal_to_Ring]
              if (normal_to_Ring == [0.0, 0.0, 0.0]):
                 if (direction_bq == 'POSITIVE'):
                    direct = ('set', [0.0, 0.0, DEFAULT_XY_DISTANCE] + [cmx, cmy, cmz])
                    BQGuide[ref_ring] = ([cmx, cmy, cmz],[0.0, 0.0, DEFAULT_XY_DISTANCE])
                    norm_vec[ref_ring] = [cmx, cmy, cmz, 0.0, 0.0, DEFAULT_XY_DISTANCE]
                 elif (direction_bq == 'NEGATIVE'):
                    direct = ('set', [0.0, 0.0, -DEFAULT_XY_DISTANCE] + [cmx, cmy, cmz])
                    BQGuide[ref_ring] = ([cmx, cmy, cmz],[0.0, 0.0, -DEFAULT_XY_DISTANCE])
                    norm_vec[ref_ring] = [cmx, cmy, cmz, 0.0, 0.0, -DEFAULT_XY_DISTANCE]
              else:
                 if (direction_bq == 'POSITIVE'):
                    direct = ('set', [normal_to_Ring[0]+cmx, normal_to_Ring[1]+cmy, normal_to_Ring[2]+cmz, cmx, cmy, cmz])
                    BQGuide[ref_ring] = ([cmx, cmy, cmz], [+normal_to_Ring[0]+cmx, +normal_to_Ring[1]+cmy, +normal_to_Ring[2]+cmz])
                    norm_vec[ref_ring] = [cmx, cmy, cmz, +normal_to_Ring[0]+cmx, +normal_to_Ring[1]+cmy, +normal_to_Ring[2]+cmz]
                 elif (direction_bq == 'NEGATIVE'):
                    direct = ('set', [-normal_to_Ring[0]+cmx, -normal_to_Ring[1]+cmy, -normal_to_Ring[2]+cmz, cmx, cmy, cmz])
                    BQGuide[ref_ring] = ([cmx, cmy, cmz], [-normal_to_Ring[0]+cmx, -normal_to_Ring[1]+cmy, -normal_to_Ring[2]+cmz])
                    norm_vec[ref_ring] = [cmx, cmy, cmz, -normal_to_Ring[0]+cmx, -normal_to_Ring[1]+cmy, -normal_to_Ring[2]+cmz]
              break;

           elif (new_geom == []):
              print("None of the rings is Planar within the tolerence of " + repr(TORSION_ANGLE_TOLERANCE) + " degrees. Therefore, aboring the run .. ")
              sys.exit(10)

   xy_ref_ring_info.append(direct[1])

   for ring in CenterOf:
      if (ring != ref_ring):
         ring_atoms = getOrderedRing(Conn, CenterOf.get(ring))
         if (len(ring_atoms) > 2):
             cmx, cmy, cmz = getGMOfRing(new_geom, ring_atoms)
             normal_to_Ring = getUnitVector(getAverageNormaltoTheRing(new_geom, ring_atoms, [cmx, cmy, cmz]))
             normal_to_Ring[:] = [vc*DEFAULT_XY_DISTANCE for vc in normal_to_Ring]
             check_direction = getDihedralAngleBetween(direct[1][0:3], direct[1][3:6], [cmx, cmy, cmz], [normal_to_Ring[0]+cmx, normal_to_Ring[1]+cmy, normal_to_Ring[2]+cmz]) 
             # Tolerence of 15 degrees seem to be reasonable
             # But it is now updated to 90 which is maximum, for the case of twisted rings
             if (abs(check_direction) < 90):
                BQGuide[ring] = ([cmx, cmy, cmz], [+normal_to_Ring[0]+cmx, +normal_to_Ring[1]+cmy, +normal_to_Ring[2]+cmz])
                norm_vec[ring] = [cmx, cmy, cmz, +normal_to_Ring[0]+cmx, +normal_to_Ring[1]+cmy, +normal_to_Ring[2]+cmz]
             else:
                BQGuide[ring] = ([cmx, cmy, cmz], [-normal_to_Ring[0]+cmx, -normal_to_Ring[1]+cmy, -normal_to_Ring[2]+cmz])
                norm_vec[ring] = [cmx, cmy, cmz, -normal_to_Ring[0]+cmx, -normal_to_Ring[1]+cmy, -normal_to_Ring[2]+cmz]

   ring = 1
   if (ref_ring != 1):
      ring_atoms = getOrderedRing(Conn, CenterOf.get(ring))
      cmx, cmy, cmz = getGMOfRing(new_geom, ring_atoms)
      vec1 = norm_vec[ref_ring]
      BQGuide[1] = ([cmx, cmy, cmz], [vec1[3]-vec1[0]+cmx, vec1[4]-vec1[1]+cmy, vec1[5]-vec1[2]+cmz])
      norm_vec[ring] = [cmx, cmy, cmz, vec1[3]-vec1[0]+cmx, vec1[4]-vec1[1]+cmy, vec1[5]-vec1[2]+cmz]


   for ring in CenterOf:
      if ((ring != 1) and (ring != len(CenterOf))): 
         ring_atoms = getOrderedRing(Conn, CenterOf.get(ring))
         if (len(ring_atoms) < 3):
             cmx, cmy, cmz = getGMOfRing(new_geom, ring_atoms)

             prv = ref_ring; nxt = ref_ring
             for p in range (ring-1, 1, -1):
                if (len(CenterOf.get(p)) > 3): prv = p; break;
             for p in range (ring+1, len(CenterOf)-1):
                if (len(CenterOf.get(p)) > 3): nxt = p; break;
             if (nxt == ref_ring): nxt == prv

             vec1 = [norm_vec[prv][3] - norm_vec[prv][0] + cmx, norm_vec[prv][4] - norm_vec[prv][1] + cmy, norm_vec[prv][5] - norm_vec[prv][2] + cmz]
             vec2 = [norm_vec[nxt][3] - norm_vec[nxt][0] + cmx, norm_vec[nxt][4] - norm_vec[nxt][1] + cmy, norm_vec[nxt][5] - norm_vec[nxt][2] + cmz]
             mid_vec = [(vec1[0]+vec2[0])/2., (vec1[1]+vec2[1])/2., (vec1[2]+vec2[2])/2.]
             BQGuide[ring] = ([cmx, cmy, cmz], mid_vec)
             norm_vec[ring] = [cmx, cmy, cmz, mid_vec[0], mid_vec[1], mid_vec[2]]

   ring = len(CenterOf)
   ring_atoms = getOrderedRing(Conn, CenterOf.get(ring))
   if (len(ring_atoms) <3):
   # Did you see a heart somewhere nearby ?
      cmx, cmy, cmz = getGMOfRing(new_geom, ring_atoms)
      vec2 = norm_vec[ring-1]
      BQGuide[ring] = ([cmx, cmy, cmz], [vec2[3]-vec2[0]+cmx, vec2[4]-vec2[1]+cmy, vec2[5]-vec2[2]+cmz])
      norm_vec[ring] = [cmx, cmy, cmz, vec2[3]-vec2[0]+cmx, vec2[4]-vec2[1]+cmy, vec2[5]-vec2[2]+cmz]

   for pt in new_points:
      cmx, cmy, cmz = new_points[pt][1:4]
      prv = pt[0]
      nxt = pt[1]
      vec1 = [norm_vec[prv][3] - norm_vec[prv][0] + cmx, norm_vec[prv][4] - norm_vec[prv][1] + cmy, norm_vec[prv][5] - norm_vec[prv][2] + cmz]
      vec2 = [norm_vec[nxt][3] - norm_vec[nxt][0] + cmx, norm_vec[nxt][4] - norm_vec[nxt][1] + cmy, norm_vec[nxt][5] - norm_vec[nxt][2] + cmz]
      mid_vec = [(vec1[0]+vec2[0])/2., (vec1[1]+vec2[1])/2., (vec1[2]+vec2[2])/2.]
      BQGuide[(prv+nxt)/2.0] = (new_points[pt][1:4], mid_vec)
   
   coord_format = "{0:.5f}"
   BQs_string = ""
   bq_count = 0
   seq_BQGuide = sorted(BQGuide)
   
   bq_coord_prvs = BQGuide[seq_BQGuide[0]][1]
   for i in range (0, len(seq_BQGuide)-1):
      
      a = BQGuide[seq_BQGuide[i]][1]
      b = BQGuide[seq_BQGuide[i+1]][1]
      vec_ab = getVector(a,b)
      norm_vec_ab = vectorMagnitude(vec_ab)
      n_vec_ab = getUnitVector(vec_ab)
      n_BQ = int(round((norm_vec_ab/BQ_Step),0))
     
      for j in range (0, n_BQ):
         bq_coord = [BQ_Step*j*v for v in n_vec_ab]
         bq_coord[:] = [a[v]+bq_coord[v] for v in range (0,3)]
         xy_BQ_dist.append(round(getDistance(bq_coord, bq_coord_prvs),3))
         BQs_string += 'bq     ' + coord_format.format(bq_coord[0]) + "     " + coord_format.format(bq_coord[1]) + "     " + coord_format.format(bq_coord[2]) + "\n"
         bq_count += 1
         bq_coord_prvs = bq_coord
         if (bq_count%50 == 0): BQs_string += "break"

   # For the last BQ
   if (vectorMagnitude(getVector(bq_coord, b)) < BQ_Step):
      bq_coord = [BQ_Step*(j+1)*v for v in n_vec_ab]
      bq_coord[:] = [a[v]+bq_coord[v] for v in range (0,3)]
      xy_BQ_dist.append(round(getDistance(bq_coord, bq_coord_prvs),3))
      BQs_string += 'bq     ' + coord_format.format(bq_coord[0]) + "     " + coord_format.format(bq_coord[1]) + "     " + coord_format.format(bq_coord[2]) + "\n"
      bq_count += 1

   BQ_No = bq_count
   return BQs_string, new_geom

   
def run_Nics():

   # global flags
   global opt_flag, ncs_flag, sigma_flag, xy_flag, analyse_flag, area_flag, s_charge_flag, opt_external, optfl_external
   # global molecule-related
   global CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, normals
   # global technical 
   global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, xy_BQ_dist, sigma_charge, analyse_dist, clear_flag

   if (not xy_flag): 
      dict_cen = CenterOf.copy()
      if (len(normals) > 0):
         n_count = len(CenterOf)
         for i in normals:
            dict_cen[n_count+i] = normals[i]

   elif (xy_flag): dict_cen = n_xy_center
   for ring in dict_cen:
      flname = flprfx + "-center" + repr(ring)
      print("Job " + GaussCmd + flname + " " + flname + " running ..")
      status = execCmd(constructGaussCMD(flname))
      status = 0
      print("Job Over.")
      if (not status) : continue 
      else : 
          print("It seems that the NICS SCAN run for ring/bond number " + repr(ring) + " did not terminated normally.\n")
#          sys.exit(10)


def genSigmaModel(flprfx, geom, Conn, title, charge, mult):
   global sigma_charge, xy_flag, xy_ref_ring_info, normals, points, exocyclic

   count = len(geom)+1
   sigma_geom = {}

   # This is a complicated dictionary which consists of keys which are tuples and values which is charge on H atom to be added
   # The tuple in key denotes (Atomic number, number of connections)
   # All this is to determine the charge on the Sigma model.
   for_sigma_charge = {(5,3):-1, (6,3):0, (7,3):+1, (7,2):0, (8,2):+1, (14,3):0, (16,2):+1}
   aroma_sigma_charge = 0
   if (s_charge_flag): user_sigma_charge = sigma_charge

   for i in range (0, len(geom)):
     sigma_geom[i+1] = geom[i+1] 

   # Identify the fused bonds i.e. all pairs of fused rings
   # Check whether the pair of connected atoms is present in more than one rings.
   fused_bonds = []
   fused_rings = []
   for i in range (1, len(all_aromatic_rings)+1):
     for j in range (i+1, len(all_aromatic_rings)+1):
        intersecting_bond = list(set(all_aromatic_rings[i]) & set(all_aromatic_rings[j]))
        if ((intersecting_bond != []) and (len(intersecting_bond) == 2)): 
           fused_bonds.append(intersecting_bond)
           fused_rings.append([i,j])
           
        intersecting_bond = []

   # Identify intersections of fused bonds
   # i.e. 3 rings fused together
   fused_3 = []
   fused_3_r = []
   to_be_removed = []
   for i in range (0, len(fused_bonds)):
      for j in range (i+1, len(fused_bonds)):
         for k in range (j+1, len(fused_bonds)):
             c = set(fused_bonds[i]) & set(fused_bonds[j]) & set(fused_bonds[k])
             if ( c != set([]) ):
                fused_3.append(list(set(fused_bonds[i]) | set(fused_bonds[j]) | set(fused_bonds[k])))
                fused_3_r.append(list(set(fused_rings[i]) | set(fused_rings[j]) | set(fused_rings[k])))
                if (to_be_removed.count(fused_bonds[i]) < 1): to_be_removed.append(fused_bonds[i])
                if (to_be_removed.count(fused_bonds[j]) < 1): to_be_removed.append(fused_bonds[j])
                if (to_be_removed.count(fused_bonds[k]) < 1): to_be_removed.append(fused_bonds[k])
                
   for i,j in reversed(list(enumerate(to_be_removed))):
       fused_bonds.remove(j)
       fused_rings.remove(fused_rings[i])

   # Taking back up of All Aromatic Rings as per .arm
   all_aromatic_rings_local = all_aromatic_rings.copy() 

   indicator_for_fused_2 = len(all_aromatic_rings_local)
   r_idx = len(all_aromatic_rings_local)
   for i in range (r_idx+1, r_idx+len(fused_bonds)+1):
      all_aromatic_rings_local[i] = (fused_bonds[i-r_idx-1])

   indicator_for_fused_3 = len(all_aromatic_rings_local)
   r_idx = len(all_aromatic_rings_local)
   for i in range (r_idx+1, r_idx+len(fused_3)+1):
      all_aromatic_rings_local[i] = (fused_3[i-r_idx-1])

   ring_dummy_tuples = {}
   direct = ''
   ring_count = 0
   norm_vec = {}
   for ring in all_aromatic_rings_local:
      ring_count += 1
      if (ring_count < indicator_for_fused_3):
         ring_atoms = getOrderedRing(Conn, all_aromatic_rings_local.get(ring))
      else: ring_atoms = all_aromatic_rings_local.get(ring)
      cmx, cmy, cmz = getGMOfRing(geom, ring_atoms)

      # A dummy atom at the center of the Ring is added
      sigma_geom[count] = [0, cmx, cmy, cmz]
      count += 1
      # Another dummy atom is added from 1 angstrom distance from the CM of the ring in the user-specified direction perpendicular to the ring
      unit_normal_to_Ring = getUnitVector(getAverageNormaltoTheRing(geom, ring_atoms, [cmx, cmy, cmz]))

      if (ring_count <= indicator_for_fused_2):
         if (unit_normal_to_Ring == [0.0, 0.0, 0.0]):
            if (direct == ''):
               direct = ('set', [0.0, 0.0, 1.0] + [cmx, cmy, cmz])
               sigma_geom[count] = [0, 0.0, 0.0, 0.1]
               norm_vec[ring_count] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.1]
               if (xy_flag and (ring_atoms.sort() == xy_ref_ring_info[1])):
                  check_direction = getDihedralAngleBetween(direct[1][0:3], direct[1][3:6], xy_ref_ring_info[2][3:6], xy_ref_ring_info[2][0:3])
                  if (abs(check_direction) < 90):
                     direct = ('set', [0.0, 0.0, -1.0] + [cmx, cmy, cmz])
                     sigma_geom[count] = [0, 0.0, 0.0, -0.1]
                     norm_vec[ring_count] = [0.0, 0.0, 0.0, 0.0, 0.0, -0.1]
            else:
               sigma_geom[count] = [0, direct[1][0]+cmx, direct[1][1]+cmy, direct[1][2]+cmz]
               norm_vec[ring_count] = [cmx, cmy, cmz, direct[1][0]+cmx, direct[1][1]+cmy, direct[1][2]+cmz]

         elif (direct == ''):
            direct = ('set', [unit_normal_to_Ring[0]+cmx, unit_normal_to_Ring[1]+cmy, unit_normal_to_Ring[2]+cmz, cmx, cmy, cmz])
            sigma_geom[count] = [0, +unit_normal_to_Ring[0]+cmx, +unit_normal_to_Ring[1]+cmy, +unit_normal_to_Ring[2]+cmz]
            norm_vec[ring_count] = [cmx, cmy, cmz, +unit_normal_to_Ring[0]+cmx, +unit_normal_to_Ring[1]+cmy, +unit_normal_to_Ring[2]+cmz] 
            if (xy_flag and (ring_atoms.sort() == xy_ref_ring_info[1])):
               check_direction = getDihedralAngleBetween(direct[1][0:3], direct[1][3:6], xy_ref_ring_info[2][3:6], xy_ref_ring_info[2][0:3])
               if (abs(check_direction) < 90):
                  direct = ('set', [-unit_normal_to_Ring[0]+cmx, -unit_normal_to_Ring[1]+cmy, -unit_normal_to_Ring[2]+cmz, cmx, cmy, cmz])
                  sigma_geom[count] = [0, -unit_normal_to_Ring[0]+cmx, -unit_normal_to_Ring[1]+cmy, -unit_normal_to_Ring[2]+cmz]
                  norm_vec[ring_count] = [cmx, cmy, cmz, -unit_normal_to_Ring[0]+cmx, -unit_normal_to_Ring[1]+cmy, -unit_normal_to_Ring[2]+cmz] 

         else:
            check_direction = getDihedralAngleBetween(direct[1][0:3], direct[1][3:6], [cmx, cmy, cmz], [unit_normal_to_Ring[0]+cmx, unit_normal_to_Ring[1]+cmy, unit_normal_to_Ring[2]+cmz]) 
            # Tolerence of 15 degrees seem to be reasonable
            # But it is now updated to 90 which is maximum, for the case of twisted rings
            if (abs(check_direction) < 90):
               sigma_geom[count] = [0, +unit_normal_to_Ring[0]+cmx, +unit_normal_to_Ring[1]+cmy, +unit_normal_to_Ring[2]+cmz]
               norm_vec[ring_count] = [cmx, cmy, cmz, +unit_normal_to_Ring[0]+cmx, +unit_normal_to_Ring[1]+cmy, +unit_normal_to_Ring[2]+cmz]
            else:
               sigma_geom[count] = [0, -unit_normal_to_Ring[0]+cmx, -unit_normal_to_Ring[1]+cmy, -unit_normal_to_Ring[2]+cmz]
               norm_vec[ring_count] = [cmx, cmy, cmz, -unit_normal_to_Ring[0]+cmx, -unit_normal_to_Ring[1]+cmy, -unit_normal_to_Ring[2]+cmz]

      elif ((ring_count > indicator_for_fused_2) and (ring_count <= indicator_for_fused_3)):
         # Now all these are bonds, so unit vector is not well-defined. It has to be alightned.
         vec1 = norm_vec[fused_rings[ring_count - indicator_for_fused_2 - 1][0]]
         vec2 = norm_vec[fused_rings[ring_count - indicator_for_fused_2 - 1][1]]
         mid_vec = [(vec1[3] - vec1[0] + vec2[3] - vec2[0])/2., (vec1[4] - vec1[1] + vec2[4] - vec2[1])/2., (vec1[5] - vec1[2] + vec2[5] - vec2[2])/2.]
         sigma_geom[count] = [0, mid_vec[0] + cmx, mid_vec[1] + cmy, mid_vec[2] + cmz]

      elif (ring_count > indicator_for_fused_3):
         vec1 = norm_vec[fused_3_r[ring_count - indicator_for_fused_3 - 1][0]]
         vec2 = norm_vec[fused_3_r[ring_count - indicator_for_fused_3 - 1][1]]
         vec3 = norm_vec[fused_3_r[ring_count - indicator_for_fused_3 - 1][2]]
         mid_vec = [(vec1[3] - vec1[0] + vec2[3] - vec2[0] + vec3[3] - vec3[0])/3., (vec1[4] - vec1[1] + vec2[4] - vec2[1] + vec3[4] - vec3[1])/3., (vec1[5] - vec1[2] + vec2[5] - vec2[2] + vec3[5] - vec3[2])/3.]
         sigma_geom[count] = [0, mid_vec[0] + cmx, mid_vec[1] + cmy, mid_vec[2] + cmz]

      ring_dummy_tuples[ring] = (count-1, count)
      count += 1

#   n_dummy_couple = len(ring_dummy_tuples)
   mapec = {}
   for ec in exocyclic:
      # Check which ring the exocyclic atom is connected
      for r in all_aromatic_rings:
         ec_atom = list(set(exocyclic[ec]) - (set(exocyclic[ec])&set(all_aromatic_rings[r])))
         if ( len(ec_atom) == 1 ): 
           mapec[r] = ec_atom[0]
         # So here the mapec maps the aromatic ring with its exocyclic atom
#
#      n_dummy_couple += 1
#      ring_dummy_tuples[n_dummy_couple] = 

   # Flags to make sure that H-atoms are not added to a same atom twice !!
   flag_H = {}
   for i in range (1, len(geom)+1): flag_H[i] = 0

   sigma_conn_mat, sigma_Conn = genConnectivityMatrix(sigma_geom)
   zmat, zmat_str, zmat_idx = generateZMatrix(sigma_geom, sigma_Conn)

   no_electrons = 0
   flag_for_warning = 0
 
   ring = len(all_aromatic_rings_local) - len(fused_bonds) - len(fused_3_r) 
   for ring_atoms in fused_bonds:
      ring += 1
      for i in range (0, len(ring_atoms)):
         j = ring_atoms[i]
         if ((flag_H[j] != 1) and (len(Conn[j]) != Max_Conn[geom[j][0]])):
            zmat_str += "1  " +  repr(zmat_idx[j]) + "  " + ATM_H_BL[geom[j][0]] + "  " + repr(zmat_idx[ring_dummy_tuples[ring][0]]) + "  " + FIXED_SIGMA_ANGLE + "  " + repr(zmat_idx[ring_dummy_tuples[ring][1]]) + "  " + FIXED_SIGMA_DIHEDRAL_ANGLE + "\n" 
            flag_H[j] = 1
            no_electrons += 1 
            sigma_key = (geom[j][0], len(Conn[j]))
            if (sigma_key in for_sigma_charge):
               aroma_sigma_charge += for_sigma_charge[sigma_key]
            else:
               if (not flag_for_warning): print("\nWARNING: Aroma does not have enough data for counting charge correctly for Atom with atomic no. ", geom[j][0], "\nIt is advised to give correct charge on Sigma model externally\n"); flag_for_warning = 1
   

   for ring_atoms in fused_3:

      for i in range (0, len(ring_atoms)):
         j = ring_atoms[i]
         if ( (flag_H[j] != 1) and (len(Conn[j]) != Max_Conn[geom[j][0]])):
            # Here the angle for sigma H atom will be always 90.0 as it has to be symmetric with respect to all the fused bonds
            zmat_str += "1  " +  repr(zmat_idx[j]) + "  " + ATM_H_BL[geom[j][0]] + "  " + repr(zmat_idx[ring_dummy_tuples[ring][0]]) + "  " + "90.0" + "  " + repr(zmat_idx[ring_dummy_tuples[ring][1]]) + "  " + FIXED_SIGMA_DIHEDRAL_ANGLE + "\n"
            flag_H[j] = 1
            no_electrons += 1 
            sigma_key = (geom[j][0], len(Conn[j]))
            if (sigma_key in for_sigma_charge):
               aroma_sigma_charge += for_sigma_charge[sigma_key]
            else:
               if (not flag_for_warning): print("\nWARNING: Aroma does not have enough data for counting charge correctly for Atom with atomic no. ", geom[j][0], "\nIt is advised to give correct charge on Sigma model externally\n"); flag_for_warning = 1

   for ring in all_aromatic_rings_local:
      if (ring_count < indicator_for_fused_3):
         ring_atoms = getOrderedRing(Conn, all_aromatic_rings_local.get(ring))
      else: ring_atoms = all_aromatic_rings_local.get(ring)

      for i in range (0, len(ring_atoms)):
         j = ring_atoms[i]
         if ( (flag_H[j] != 1) and (len(Conn[j]) != Max_Conn[geom[j][0]])):
            zmat_str += "1  " +  repr(zmat_idx[j]) + "  " + ATM_H_BL[geom[j][0]] + "  " + repr(zmat_idx[ring_dummy_tuples[ring][0]]) + "  " + FIXED_SIGMA_ANGLE + "  " + repr(zmat_idx[ring_dummy_tuples[ring][1]]) + "  " + FIXED_SIGMA_DIHEDRAL_ANGLE + "\n"
            flag_H[j] = 1
            no_electrons += 1 
            sigma_key = (geom[j][0], len(Conn[j]))
            if (sigma_key in for_sigma_charge):
               aroma_sigma_charge += for_sigma_charge[sigma_key]
            else:
               if (not flag_for_warning): print("\nWARNING: Aroma does not have enough data for counting charge correctly for Atom with atomic no. ", geom[j][0], "\nIt is advised to give correct charge on Sigma model externally\n"); flag_for_warning = 1

   for r in mapec:
       zmat_str += "1  " +  repr(zmat_idx[mapec[r]]) + "  " + ATM_H_BL[geom[mapec[r]][0]] + "  " + repr(zmat_idx[ring_dummy_tuples[r][0]]) + "  " + FIXED_SIGMA_ANGLE + "  " + repr(zmat_idx[ring_dummy_tuples[r][1]]) + "  " + FIXED_SIGMA_DIHEDRAL_ANGLE + "\n"

# ADD: Z-MATRIX ELEMENTS FOR NORMALS
   for norm in normals:
      zmat_str += "0  " + repr(zmat_idx[1]) + "   " + repr(getDistance(geom[1][1:4],normals[norm][0:3])) + "   " + repr(zmat_idx[2]) + "   " + repr(getAngleBetweenPoints(normals[norm][0:3], geom[1][1:4], geom[2][1:4])) + "   " + repr(zmat_idx[3]) + "   " + repr(getDihedralAngleBetween(normals[norm][0:3], geom[1][1:4], geom[2][1:4], geom[3][1:4])) + "\n"
      zmat_str += "0  " + repr(zmat_idx[1]) + "   " + repr(getDistance(geom[1][1:4],normals[norm][3:6])) + "   " + repr(zmat_idx[2]) + "   " + repr(getAngleBetweenPoints(normals[norm][3:6], geom[1][1:4], geom[2][1:4])) + "   " + repr(zmat_idx[3]) + "   " + repr(getDihedralAngleBetween(normals[norm][3:6], geom[1][1:4], geom[2][1:4], geom[3][1:4])) + "\n"

   for p in points:
      zmat_str += "0 " + repr(zmat_idx[1]) + "   " + repr(getDistance(geom[1][1:4],points[p][0:3])) + "   " + repr(zmat_idx[2]) + "   " + repr(getAngleBetweenPoints(points[p][0:3], geom[1][1:4], geom[2][1:4])) + "   " + repr(zmat_idx[3]) + "   " + repr(getDihedralAngleBetween(points[p][0:3], geom[1][1:4], geom[2][1:4], geom[3][1:4])) + "\n"

   aroma_sigma_charge += charge
   no_electrons += -aroma_sigma_charge 
   if (no_electrons%2 != 0):
      if (charge == 0):
        aroma_sigma_charge += +1 
      elif (charge > 0):
        aroma_sigma_charge += -1
      elif (charge < 0):
        aroma_sigma_charge += +1 
   
   if (s_charge_flag):
      if (user_sigma_charge == aroma_sigma_charge):
         sigma_charge = aroma_sigma_charge
      else: 
         print("\nWARNING: User Specifed Charge on Sigma Model is " + repr(user_sigma_charge) + " while that determined by Aroma is " + repr(aroma_sigma_charge) + "\nHowever, Aroma will proceed with Charge Provided by the User.\n")
         sigma_charge = user_sigma_charge
   else: sigma_charge = aroma_sigma_charge

   sigma_flprfx = inpdir + flprfx + "-sigma" +  GaussInpExt
   f = open (sigma_flprfx, "w")
   f.write("# \n\n" + title + " sigma only model " + "\n\n" + repr(sigma_charge) + " " + repr(mult) + "\n")
   f.write(zmat_str + "\n") 
   f.close()

   return sigma_flprfx, zmat_idx

def getNewRings(geom, sigma_geom, CenterOf, zmat_idx):
   new_CenterOf = {}
   for i in CenterOf:
      org_ring = CenterOf.get(i)
      new_ring = []
      for j in range (0, len(org_ring)):
         k = zmat_idx[org_ring[j]]
         new_ring.append(k)
      new_CenterOf[i] = new_ring
   return new_CenterOf


def grepData():

      # global flags
      global opt_flag, ncs_flag, sigma_flag, xy_flag, analyse_flag, area_flag, s_charge_flag, opt_external, optfl_external
      # global molecule-related
      global CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, BQGuide
      # global technical 
      global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, xy_BQ_dist, sigma_charge, analyse_dist, clear_flag

      nBQs = []
      if (not xy_flag): 
         dict_cen = CenterOf.copy()
         if (len(normals) > 0):
            n_count = len(CenterOf)
            for i in normals:
                dict_cen[n_count+i] = normals[i]

         for i in range (0, len(dict_cen)):
            nBQs.append(BQ_No)

      elif (xy_flag): 
         dict_cen = n_xy_center
         idx = 0
         dist = 0.0
         for i in range (0, len(dict_cen)-1):
            nBQs.append(50)
         nBQs.append(BQ_No - 50*(len(nBQs)))

      for ring in dict_cen:
         Plane = 'XY'

         outlines = readFile(outdir + flprfx + "-center" + repr(ring) + GaussOutExt)

         for i in range (0, len(outlines)):
           if (outlines[i].find("NAT") >= 0 ): break;
         nat = int(outlines[i].split()[2])

         for i in range (0, len(outlines)):
           if (outlines[i].find("Magnetic shielding tensor") >= 0 ): break;

         f_out = open(outdir + flprfx + "-center" +repr(ring) + ".armdat", "w")
         f_out.write("#       oop       in1        in2       inp       iso        x         y         z\n")

         # i + 5*nat + 1, i + 5*nat + BQ_No*nat + 1 are start and end line numbers in the output file for the data for BQs
         if (not xy_flag): 
            dist = BQ_Range[0]

         for j in range (i + 5*nat + 1, i + 5*nat + nBQs[ring-1]*5 + 1, 5 ):

            BQ_data_string = ""
            if xy_flag: dist += xy_BQ_dist[idx]

            words = outlines[j].strip().split()
            if (words[1] == 'Bq'):
               # This is commented as BQ no is not important, instead distance of that BQ from GM is added (further)
               # BQ_data_string += words[0] + "   "

               # The isotropic value is in the first line
               iso = -float(words[4])

               # Then get the diagonal values for the tensor
               xx = -float(outlines[j+1].strip().split()[1])
               yy = -float(outlines[j+2].strip().split()[3])
               zz = -float(outlines[j+3].strip().split()[5])

               # Then get the three eigenvalues
               # Out-of-Plane eigenvalue is the one closest to ZZ
               e1, e2, e3 = list(map(float, outlines[j+4].split(":")[1].split()))
               sorted_e = []
               sorted_e.append(e1)
               close = abs(e1 + zz)
               if ( abs(e2 + zz) <= close): 
                  close = abs(e2 + zz)
                  sorted_e.insert(0, e2)
               elif ( abs(e2 + zz) > close):
                  sorted_e.append(e2)
               if ( abs(e3 + zz) <= close): 
                  sorted_e.insert(0, e3)
               elif ( abs(e3 + zz) > close):
                  sorted_e.append(e3)
               oup = -sorted_e[0]; inp1 = -sorted_e[1]; inp2 = -sorted_e[2]
               # In-plane chemical shift is average of the two in-plane shifts
               inp = 0.5*(inp1+inp2)

               BQ_data_string += fpformat.fix(dist,2) + "   " + fpformat.fix(oup, 4) + "   " + fpformat.fix(inp1, 4) + "   " + fpformat.fix(inp2, 4) + "   " + fpformat.fix(inp, 4) + "   " + fpformat.fix(iso, 4) + "   " + fpformat.fix(xx, 4) + "   " + fpformat.fix(yy, 4) + "   " + fpformat.fix(zz, 4) +  "\n"
               f_out.write(BQ_data_string)
               
               if (not xy_flag): dist += BQ_Step
               else : idx += 1

            else:
               print(ring)
               print("Number of Atoms and BQ number not matching .. \nTherefore, Skipping the Step Of Filtering and Storing the Data ..")
               f_out.close()
               return
            
         f_out.close()
         
      # In case of XY-Scan, Merge all armdats to one
      if (xy_flag and len(dict_cen) > 1):
         final_armdat = open(outdir + flprfx + "-center1" + ".armdat","a")
         for ring in range (1, len(dict_cen)):
            lines = readFile(outdir + flprfx + "-center" +repr(ring+1) + ".armdat")
            for i in range (1, len(lines)):
               final_armdat.write(lines[i])
         final_armdat.close()
      if (xy_flag):
         execCmd(" mv " + outdir + flprfx + "-center1" + ".armdat " + outdir + flprfx + "-allcenter" + ".armdat")


def Execute(geom, title, charge, mult, Conn):

   # global flags
   global opt_flag, ncs_flag, sigma_flag, xy_flag, analyse_flag, area_flag, s_charge_flag, opt_external, optfl_external
   # global molecule-related
   global CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, BQGuide
   # global technical 
   global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, sigma_charge, analyse_dist, clear_flag

   # Run optimization, if required
   if (hashLine_opt == "DEFAULT\n"): hashLine_opt = DEFAULT_OPTIMIZATION_KEYLINE
   if (opt_flag):
      if (not opt_external):
         optfl = generate_Opt_Input(geom, hashLine_opt, title, charge, mult)
      else:
         optfl = flprfx + "-opt"
         print("\nCopying external optimization input file to " + inpdir + optfl + GaussInpExt)
         execCmd("cp " + optfl_external + " " + inpdir + optfl + GaussInpExt)

      run_Optimization(optfl)
      theParser = ReaderFunctCall["output"](outdir + optfl + GaussOutExt)
      geom, hashLine, title, charge, mult = theParser.getInpData()
      conn_mat, Conn = genConnectivityMatrix(geom)

   if (hashLine_nics == "DEFAULT\n"): hashLine_nics = DEFAULT_NICS_KEYLINE 
   if (hashLine_ncs == "DEFAULT\n"): hashLine_ncs = DEFAULT_NCS_KEYLINE 
   if (hashLine_nbo == "DEFAULT\n"): hashLine_nbo = DEFAULT_NBO_KEYLINE 

   if (ncs_flag):
      genNicsInputs(geom, Conn, hashLine_ncs, title, charge, mult)
   else : genNicsInputs(geom, Conn, hashLine_nics, title, charge, mult)

   print("\nStatus : NICS input for all the centers generated.")

   run_Nics()
   print("\nAll the jobs are over")
   print("\nStatus : Filtering Appropriate Data .. ")

   grepData()

   if (ncs_flag):
      # pi-MOs are same in all the output files, so just identify them from the first file.
      piMOs, nocc = identifyPiMOs(outdir + flprfx + "-center1" + GaussOutExt)
      for ring in CenterOf:
         outfl = outdir + flprfx + "-center" + repr(ring)
         grepPiCMO(len(geom), piMOs, nocc, BQ_No, BQ_Range, BQ_Step, outfl, GaussOutExt)

   if (xy_flag):
      outfl = open(outfilename, "a")
      outfl.write("\n\nThe Centers for Rings/Bonds from the .arm file corresponds to following BQs\n")
      outfl.write("Ring  Distance    X       Y       Z\n")
      for i in range (1, len(BQGuide)+1):
         if i == 1: dist = 0.0
         else: dist += round(getDistance(BQGuide[i][1], BQGuide[i-1][1]),3)
         outfl.write("  " + repr(i) + "     " + fpformat.fix(dist,1) + "   " + fpformat.fix(BQGuide[i][1][0],3) + "   " + fpformat.fix(BQGuide[i][1][1],3) + "   " + fpformat.fix(BQGuide[i][1][2],3) + "\n")
      outfl.close()

def callAnalyse(flprfx, geom, CenterOf, all_aromatic_rings, analyse_dist, outfl):

   conn_mat, Conn = genConnectivityMatrix(geom)

   for ring in CenterOf:
       m_fl = flprfx + "-center" + repr(ring) + ".armdat"
       s_fl = flprfx + "-sigma-center" + repr(ring) + ".armdat"
       outfl.write("\n\nFor Center " + repr(ring))
       analyse(m_fl, s_fl, analyse_dist, outfl)

   n_count = len(CenterOf)
   if (len(normals) > 0):
       for i in normals:
          n_count += 1
          m_fl = flprfx + "-center" + repr(n_count) + ".armdat"
          s_fl = flprfx + "-sigma-center" + repr(n_count) + ".armdat"
          outfl.write("\n\nFor Center " + repr(n_count))
          analyse(m_fl, s_fl, analyse_dist, outfl)

   if (area_flag):
      area = {}
      tot_area = 0.0
      outfl.write("\nThe area in sq. ang. for each of the aromatic ring defined are: ")
      for ring in all_aromatic_rings:
          ring_atoms = all_aromatic_rings.get(ring)
          area[ring] = ring_area(inpdir + flprfx + "-center1.in", ring_atoms)
          outfl.write("\n" + repr(ring) + " : " + repr(area[ring]))
          tot_area += area[ring]
       
      outfl.write("\n\nThe total area is " + repr(tot_area) + " sq. ang.\n\n")


def aroma(armfile):

   # global flags
   global opt_flag, ncs_flag, sigma_flag, analyse_flag, area_flag, s_charge_flag, opt_external, optfl_external, xy_flag
   # global molecule-related
   global CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, normals
   # global technical 
   global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, sigma_charge, analyse_dist, clear_flag

   print("\n--------------------------------------------------------------------")
   print("                        ** Aroma Run Begins. **")
   print("--------------------------------------------------------------------\n")

   init()

   flprfx = armfile[armfile.rindex("/")+1:len(armfile)]
   check(armfile)
      
   if (not opt_external):
      for extension in EXTENSIONS_FOR_GAUSSIAN_FILES:
         if (EXTENSIONS_FOR_GAUSSIAN_FILES[extension].count(geomflext) == 1): exttype = extension

      #  Read the geometry and other data
      theParser = ReaderFunctCall[exttype](geomfl)
      geom, hashLine, title, charge, mult = theParser.getInpData() 
      conn_mat, Conn = genConnectivityMatrix(geom)
   else:
      geom = {}; title = ""; charge=""; mult=""; Conn= []

   if (xy_flag):
      print("The final output will be stored in " + outfilename)
      outfl = open(outfilename, "w")
      outfl.write("\n--------------------------------------------------------------------")
      outfl.write("\n                        ** Aroma Run **")
      outfl.write("\n--------------------------------------------------------------------\n")
      outfl.write("\nFor the Original Molecule:\n")
      outfl.close()

   Execute(geom, title, charge, mult, Conn)


   if (clear_flag):
      print("\nClearing up unnecessary files .. \n")
      execCmd("rm " + inpdir + flprfx + "-center*")
      execCmd("rm " + inpdir + flprfx + "-guessonly* " + outdir + flprfx + "-guessonly*")
      if (opt_flag):
         execCmd("rm " + inpdir + flprfx + "-opt* ")

   # Read For: For XY-Scan with Sigma-Model, just keep the final output as .armlog with r, ZZ and del-ZZ
   if (xy_flag and sigma_flag):
       xarmlogfile = outdir + flprfx + "-alldiff.armlog"
       armdatlines = readFile(outdir + flprfx + "-allcenter" + ".armdat")

   if (sigma_flag):
      if (opt_flag or opt_external):
         theParser = ReaderFunctCall["output"](outdir + flprfx + "-opt" + GaussOutExt)
         geom, hashLine, title, charge, mult = theParser.getInpData()
         conn_mat, Conn = genConnectivityMatrix(geom)

      opt_flag = 0; ncs_flag = 0
      exttype = "input" 
      
      geomfl, zmat_idx = genSigmaModel(flprfx, geom, Conn, title, charge, mult)

      # Read the geometry and other data
      theParser = ReaderFunctCall[exttype](geomfl)
      sigma_geom, hashLine, title, charge, mult = theParser.getInpData()

      # Take out the normals from the geometry, if defined
      sigma_normals = {}
      n_count = 1
      for i in range (len(sigma_geom)-(2*len(normals))+1, len(sigma_geom)+1, 2):
         sigma_normals[n_count] = sigma_geom[i][1:4] + sigma_geom[i+1][1:4]; n_count += 1
         del sigma_geom[i]
         del sigma_geom[i+1]
      org_normals = normals
      normals = sigma_normals
         
      conn_mat, Conn = genConnectivityMatrix(sigma_geom)

      new_CenterOf = getNewRings(geom, sigma_geom, CenterOf, zmat_idx)
      org_flprfx = flprfx; org_CenterOf = CenterOf
      flprfx = flprfx + "-sigma"
      CenterOf = new_CenterOf

      if (xy_flag):
         outfl = open(outfilename, "a")
         outfl.write("\n--------------------------------------------------------------------\n")
         outfl.write("\nFor the Sigma Model:\n")
         outfl.close()

      Execute(sigma_geom, title, charge, mult, Conn)

      if (clear_flag):
         print("\nClearing up unnecessary files .. \n")
         execCmd("rm " + inpdir + flprfx + "-center*")
         execCmd("rm " + inpdir + flprfx + "-guessonly* " + outdir + flprfx + "-guessonly*")

   numpy_flag = checkNumPy()
   if (sigma_flag and analyse_flag and numpy_flag and not xy_flag):
       print("The final output will be stored in " + outfilename)
       outfl = open(outfilename, "w")
       outfl.write("\n--------------------------------------------------------------------")
       outfl.write("\n                        ** Aroma Run **")
       outfl.write("\n--------------------------------------------------------------------\n")
       callAnalyse(org_flprfx, geom, org_CenterOf, all_aromatic_rings, analyse_dist, outfl)

   # For XY-Scan with Sigma-Model, just keep the final output as .armlog with r, ZZ and del-ZZ
   if (xy_flag and sigma_flag):
      armlog = open(xarmlogfile, "w")
      armlog.write("r       ZZ       Sigma-ZZ        Del-ZZ\n")
      if (xy_flag, sigma_flag):
          sarmdatlines = readFile(outdir + flprfx + "-allcenter" + ".armdat")
          lineno = min(len(armdatlines), len(sarmdatlines))
          for i in range (1, lineno):
             awords = list(map(float, armdatlines[i].split()))
             sawords = list(map(float, sarmdatlines[i].split()))
             armlog.write(fpformat.fix(awords[0],2) + "   " + fpformat.fix(awords[8],4) + "   " + fpformat.fix(sawords[8],4) + "   " + fpformat.fix(awords[8]-sawords[8], 4) + "\n" )
      armlog.close()
       

   print("\n--------------------------------------------------------------------")
   print("                        ** Aroma Run Over. **")
   print("--------------------------------------------------------------------\n")

if __name__ == "__main__":
   aroma(sys.argv[1])
