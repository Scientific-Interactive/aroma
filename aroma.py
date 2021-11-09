#! /usr/bin/env python
# Author : Anuja
# 30.01.2013
# Last Updated : 02.01.2021

# Main Driver scipt for automation of NICS-Scan in Z and XY Directions, Sigma-Only Model and CMO-NICS

import os
import re
import sys
import glob
import json
import shutil
import string
import fpformat

from aroma_constants import *
import aroma_constants

import aroma_user_constants
aroma_user_constants.loadUserConstants(sys.argv[0])
from aroma_constants import *
import aroma_constants
print(os.getcwd())

from aroma_ringarea import *
import aroma_ringarea
from aroma_analysis import *
import aroma_analysis
from aroma_pinbo import *
import aroma_pinbo
from aroma_molecule import *
import aroma_molecule
from aroma_writer import *
import aroma_writer
from aroma_parser import *
import aroma_parser
from aroma_util import *
import aroma_util


def init():
    # global flags
    global opt_flag, ncs_flag, sigma_flag, xy_flag, pointonly_flag, integralnics_flag, analyse_flag, area_flag, s_charge_flag, s_mult_flag, opt_external, optfl_external, inponly_flag, outonly_flag
    # global molecule-related
    global armpath, CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, xy_ref_ring_info, BQGuide, points, normals, exocyclic
    # global technical
    global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, xy_BQ_dist, sigma_charge, sigma_mult, analyse_dist, clear_flag, xy_extend
    # global file set - used for storing all input files generated
    global inputFileSet, collatedFileSet

    # Initializing some global variables
    opt_flag = 0
    ncs_flag = 0
    sigma_flag = 0
    opt_external = 0
    xy_flag = 0
    pointonly_flag = 0
    integralnics_flag = 0
    analyse_flag = 1
    area_flag = 0
    inponly_flag = 0
    outonly_flag = 0
    CenterOf = {}
    all_aromatic_rings = {}
    exocyclic = {}
    points = {}
    normals = {}
    geomflext = ""
    geomfl = ""
    flprfx = ""
    outfilename = ""
    sigma_direction = 'POSITIVE'
    n_xy_center = {1: 1}
    xy_ref_ring_info = []
    BQGuide = {}
    xy_BQ_dist = []
    runtype = "NICSSCAN"
    hashLine_opt = externalProgram["defaultOptimizationKeyline"]
    hashLine_nics = externalProgram["defaultNicsKeyline"]
    hashLine_ncs = externalProgram["defaultNcsKeyline"]
    hashLine_nbo = externalProgram["defaultNboKeyline"]
    BQ_Step = DEFAULT_BQ_STEP
    BQ_Range = DEFAULT_BQ_RANGE
    sigma_charge = 0
    s_charge_flag = 0
    sigma_mult = 1
    s_mult_flag = 0
    analyse_dist = DEFAULT_DISTANCE_FOR_ANALYSIS
    xy_extend = 0.0
    clear_flag = 0

    inputFileSet = []
    collatedFileSet = []


def check(armfile):

    # global flags
    global opt_flag, ncs_flag, sigma_flag, xy_flag, pointonly_flag, integralnics_flag, analyse_flag, area_flag, s_charge_flag, s_mult_flag, opt_external, optfl_external, inponly_flag, outonly_flag
    # global molecule-related
    global armpath, CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, xy_ref_ring_info, BQGuide, points, normals
    # global technical
    global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, xy_BQ_dist, sigma_charge, sigma_mult, analyse_dist, clear_flag, xy_extend

    armpath = armfile[0:armfile.rindex("/")+1]
    flprfx = armfile[armfile.rindex("/")+1:len(armfile)]
    armlines = readFile(armpath + flprfx + ".arm")

    outfilename = externalProgram["outdir"] + flprfx + ".armlog"

    # Check for Validity of the RunType
    for i in range(0, len(armlines)):
        if (armlines[i].upper().find("RUN") >= 0):
            runseq = re.split(
                "[ |,]", armlines[i].upper().strip().split("=")[1])
            if (runseq.count("OPT") > 0):
                opt_flag = 1
            if (runseq.count("NCS") > 0):
                ncs_flag = 1
            if (runseq.count("SIGMA") > 0):
                sigma_flag = 1
            if (runseq.count("XY") > 0):
                xy_flag = 1
            if (runseq.count("PTONLY") > 0):
                pointonly_flag = 1
            if (runseq.count("INTEGRALNICS") > 0):
                integralnics_flag = 1
                sigma_flag = 1
            if (runseq.count("INPONLY") > 0):
                inponly_flag = 1
            if (runseq.count("OUTONLY") > 0):
                outonly_flag = 1
            break

    if (opt_flag):
        for i in range(0, len(armlines)):
            if (armlines[i].upper().find("OPT_EXTERNAL") >= 0):
                opt_external = 1
                optfl_external = armlines[i].strip().split("=")[1]
                break

    # Check for input file for Geomtry
    if (not opt_external):
        for i in range(0, len(armlines)):
            if (armlines[i].upper().find("GEOMFILE") >= 0):
                geomfl = armlines[i].strip().split("=")[1]
                break

        geomfl_flag = os.path.exists(geomfl)
        if (not geomfl_flag):
            print(geomfl + " Could Not be Found.\nTherefore, Aborting the Run ..")
            sys.exit(10)

        geomflext = geomfl[geomfl.rindex(".")+1:len(geomfl)+1]

        valid_ext_flag = 1
        for extension in externalProgram["extensions"]:
            if (externalProgram["extensions"][extension].count(geomflext) > 0):
                valid_ext_flag = 1
        if (not valid_ext_flag):
            print("File with \"" + geomflext +
                  "\" Extension Can Not be Read.\nTherefore, Aborting the Run ..")
            sys.exit(10)

    if (xy_flag):
        print("\nWARNING: You Have Requested XY-Scan. Make Sure That the Centers Are Defined in Proper Order.\n")
        for i in range(0, len(armlines)):
            if (armlines[i].upper().find("XYEXTEND") >= 0):
                xy_extend = float(armlines[i].split("=")[1])
                break

    for i in range(0, len(armlines)):
        if (armlines[i].upper().find("OUTFILE") >= 0):
            outfilename = armlines[i].strip().split("=")[1]
            break

    # Get The Ring/Bond Info
    r_count = 0
    n_count = 0
    p_count = 0
    for i in range(0, len(armlines)):
        if (armlines[i].upper().find("CENTER") >= 0):
            r_count += 1
            CenterOf[r_count] = list(
                map(int, re.split("[:|=]", armlines[i].strip())[1].split(",")))
        if (armlines[i].upper().find("NORMAL") >= 0):
            n_count += 1
            normals[n_count] = list(
                map(float, re.split("[:|=]", armlines[i].strip())[1].split(",")))
        if (xy_flag):
            if (armlines[i].upper().find("POINT") >= 0):
                points[(r_count, r_count+1)] = list(map(float,
                                                        re.split("[:|=]", armlines[i].strip())[1].split(",")))
                points[(r_count, r_count+1)].insert(0, -1)
                if (len(points[r_count, r_count+1]) != 4):
                    print(
                        "The keyword POINT should have X, Y, Z coordinates.\n Check and Submit Again. Aborting This Run ..\n")
                    sys.exit(10)
        if (not xy_flag):
            if (armlines[i].upper().find("POINT") >= 0):
                # The 0 in tuple has no meaning. It is kept only to match the "Point" format which was defined in Aroma 1.0 version
                points[(p_count, 0)] = list(
                    map(float, re.split("[:|=]", armlines[i].strip())[1].split(",")))
                p_count += 1

    for i in range(0, len(armlines)):
        if (armlines[i].upper().find("BQSTEP") >= 0):
            BQ_Step = float(armlines[i].split("=")[1])
            break
    for i in range(0, len(armlines)):
        if (armlines[i].upper().find("BQRANGE") >= 0):
            BQ_Range = list(map(float, armlines[i].split("=")[1].split(",")))
            break

    if (not pointonly_flag):
        BQ_No = int((BQ_Range[1] - BQ_Range[0])/BQ_Step)
    else:
        BQ_No = len(points)
    if (ncs_flag and (BQ_No > 100)):
        print("The package NBO can not handle more than 100 Bqs. Aborting thr Run.\nPlease change the BQRANGE or BQSTEP and Resubmit.")
        sys.exit(10)

    if (not pointonly_flag and (len(CenterOf) < 1) and (len(normals) < 1)):
        print("Rings/Bonds Are Not Defined.\nTherefore, Aborting the Run ..")
        sys.exit(10)

    if (pointonly_flag and (len(points) < 1)):
        print("The Points (x,y,z) Are Not Defined for POINTONLY Run.\nTherefore, Aborting the Run ..")
        sys.exit(10)

    # Check keywords for Gaussian for optimization and nics runs
    for i in range(0, len(armlines)):
        if (armlines[i].upper().find("KEYLINES") >= 0):
            break
    for j in range(i+1, len(armlines)):
        if (armlines[j].upper().find("TEMPLATE") >= 0):
            if (armlines[j].upper().split("_")[0] == 'OPT'):
                hashLine_opt = ""
                for k in range(j+1, len(armlines)):
                    if ((armlines[k].upper().find("TEMPLATE") >= 0) or (armlines[k].upper().find("END KEYLINE") >= 0)):
                        break
                    hashLine_opt += armlines[k]

            if (armlines[j].upper().split("_")[0] == 'NICSSCAN'):
                hashLine_nics = ""
                for k in range(j+1, len(armlines)):
                    if ((armlines[k].upper().find("TEMPLATE") >= 0) or (armlines[k].upper().find("END KEYLINE") >= 0)):
                        break
                    hashLine_nics += armlines[k].upper()

            if (armlines[j].upper().split("_")[0] == 'NCS'):
                hashLine_ncs = ""
                for k in range(j+1, len(armlines)):
                    if ((armlines[k].upper().find("TEMPLATE") >= 0) or (armlines[k].upper().find("END KEYLINE") >= 0)):
                        break
                    hashLine_ncs += armlines[k].upper()

            if (armlines[j].upper().split("_")[0] == 'NBO'):
                hashLine_nbo = ""
                for k in range(j+1, len(armlines)):
                    if ((armlines[k].upper().find("TEMPLATE") >= 0) or (armlines[k].upper().find("END KEYLINE") >= 0)):
                        break
                    hashLine_nbo += armlines[k].upper()

    for i in range(0, len(armlines)):
        if (armlines[i].upper().find("ANALYSE") >= 0):
            if (armlines[i].upper().find("NO") >= 0):
                analyse_flag = 0
                break
            elif (armlines[i].find("=") >= 0):
                analyse_dist = float(armlines[i].split("=")[1])
            elif (armlines[i].upper().find("AREA") >= 0):
                area_flag = 1
    if integralnics_flag:
        analyse_flag = 0

    for i in range(0, len(armlines)):
        if (armlines[i].upper().find("CLEAR") >= 0):
            clear_flag = 1

    if (sigma_flag):
        for i in range(0, len(armlines)):
            if (armlines[i].upper().find("AROMATIC RING") >= 0):
                break

        rings_count = 0
        for j in range(i+1, len(armlines)):
            if (armlines[j].upper().find("END") >= 0):
                break
            rings_count += 1
            all_aromatic_rings[rings_count] = list(
                map(int, armlines[j].strip().split(",")))

        if (all_aromatic_rings == {}):
            print("\nWARNING: Aromatic Rings Are Not Provided, Therefore, Sigma-Only Model Calculations Will Not Be Performed.\n")
            sigma_flag = 0

        for i in range(0, len(armlines)):
            if (armlines[i].upper().find("EXOCYCLIC") >= 0):
                break
        ec_count = 0
        for j in range(i+1, len(armlines)):
            if (armlines[j].upper().find("END") >= 0):
                break
            ec_count += 1
            exocyclic[ec_count] = list(
                map(int, armlines[j].strip().split(",")))

        for i in range(0, len(armlines)):
            if (armlines[i].upper().find("DIRECTION") >= 0):
                words = armlines[i].upper().split()
                if ((words[2] == 'POSITIVE') or (words[2] == 'NEGATIVE')):
                    sigma_direction = words[2]

        for i in range(0, len(armlines)):
            if (armlines[i].upper().find("SONLY CHARGE") >= 0):
                s_charge_flag = 1
                sigma_charge = int(armlines[i].split("=")[1])

        for i in range(0, len(armlines)):
            if (armlines[i].upper().find("SONLY MULT") >= 0):
                s_mult_flag = 1
                sigma_mult = int(armlines[i].split("=")[1])


def generate_Opt_Input(geom, hashLine, title, charge, mult):
    global flprfx, externalProgram

    optfl = externalProgram["writerFunctCall"]["geomInput"].writeOptFile(
        flprfx, externalProgram, geom, hashLine, title, charge, mult)

    print("Input file for Optimization named as " + optfl + " is genereated.")

    return optfl


def run_Optimization(optfl):
    global externalProgram

    # Run Gaussian for optimization
    # Check status and print approprate msg before proceeding
    print("Status: Optimization Running .. ")
    status = execCmd(externalProgram["constructCmd"](optfl))
    print("Status: Optimization Over.")
    if (status):
        print("\nWARNING: Abnormal Termination of Optimization Run.")
        print("Aroma Will Continue with the Last Geometry for NICS Calculation.\n")


def writeNicsInputs(flprfx, centerIdx, flag_chk, hashLine_rev, title, charge, mult, geom, BQs_strings, jobType):
    global externalProgram
    global inputFileSet

    externalProgram["cleanupCmd"](flprfx)

    savedCenterIdx = "1"
    savedSetIdx = "-1"

    if (centerIdx.find("-set") >= 1):
       savedSetIdx = centerIdx.split("-set")[1]
       savedCenterIdx = centerIdx.split("-set")[0].split("center")[0]
    else:
       savedCenterIdx = centerIdx

    ringflName = externalProgram["inpdir"] + flprfx + \
                 "-center" + centerIdx + externalProgram["inpExt"]
    ringf = open(ringflName, "w")
    nBQs = len(BQs_strings.strip().split("\n"))
    nat = 0 
    for at in geom:
       if (geom[at][0] != 0): nat += 1

    inputFileSet.append({"filename": ringflName, "flprfx": flprfx + "-center" + centerIdx,
                        "baseprfx": flprfx,
                        "ext": externalProgram["inpExt"], "nat": nat, "nBq": nBQs, 
                        "setIdx": savedSetIdx, "centerIdx": savedCenterIdx, "jobType": jobType})

    if (flag_chk):
        ringf.write(externalProgram["writerFunctCall"]["geomInput"].genCheckpointLine(
            externalProgram["chkdir"] + flprfx + "-center" + centerIdx + ".chk"))
    # ringf.write(hashLine_rev + title + " # Center " + centerIdx + "\n\n" + repr(charge) + " " + repr(mult) + "\n")
    ringf.write(externalProgram["writerFunctCall"]["geomInput"].genHeader(
        hashLine_rev, title + " # Center " + centerIdx, charge, mult))
    for i in range(1, len(geom)+1):
        # The dummy atoms are considered as BQs, therefore, remove them
        if (geom[i][0] != 0):
            # geomline = repr(geom[i][0]) + "   " + coord_format.format(geom[i][1]) + "   " + coord_format.format(geom[i][2]) + "   " + coord_format.format(geom[i][3]) + "\n"
            geomline = externalProgram["writerFunctCall"]["geomInput"].genGeomLine(
                geom[i])
            ringf.write(geomline)
    ringf.write(
        BQs_strings + externalProgram["writerFunctCall"]["geomInput"].genTerminator() + "\n")

    # If NCS run is requested, then add NBO keywords at the end of the Gaussian input
    if (ncs_flag):
        ringf.write(hashLine_nbo + "\n")

    ringf.close()


def genNicsInputs(geom, Conn, hashLine, title, charge, mult, jobType):
    global externalProgram

    # global flags
    global opt_flag, ncs_flag, sigma_flag, xy_flag, pointonly_flag, integralnics_flag, analyse_flag, area_flag, s_charge_flag, s_mult_flag, opt_external, optfl_external, inponly_flag
    # global molecule-related
    global armpath, CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, xy_ref_ring_info, BQGuide, points, normals
    # global technical
    global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, xy_BQ_dist, sigma_charge, sigma_mult, analyse_dist, clear_flag, xy_extend

    hashLine_rev = ''
    flag_chk = 0
    hlines = hashLine.upper().split("\n")
    for i in range(0, len(hlines)):
        if (hlines[i].find("CHK") >= 0):
            flag_chk = 1
        else:
            hashLine_rev += hlines[i] + "\n"

    # POINT Keyword can add standaline points for the NICS calculation.
    if ((not xy_flag) and pointonly_flag):

        BQs_string = addBreakPoints(generateBQs_Points(points))

        # Mostly there will not be occasion where user gives more than 40-50 points, but such situation is covered.
        BQs_strings = list(map(lambda x: x.strip(), BQs_string.split("break")))

        for ring in range(0, len(BQs_strings)):
            if (BQs_strings[ring].strip() == ""):
                continue
            writeNicsInputs(flprfx, repr(ring+1), flag_chk, hashLine_rev,
                            title, charge, mult, geom, BQs_strings[ring], jobType)

    elif (not xy_flag and not pointonly_flag):  # Z-scan or Integral NICS
        for ring in CenterOf:
            ring_atoms = CenterOf.get(ring)

            # The Plane is always fixed to be XY and the molecule is oriented in such a way.
            new_geom, new_points, new_normal = reorient(geom, Conn, ring_atoms)

            if (new_geom != []):
                BQs_string = addBreakPoints(
                    generateBQs_Z(new_geom, Conn, ring_atoms))

                BQs_strings = list(
                    map(lambda x: x.strip(), BQs_string.split("break")))

                for setIdx in range(0, len(BQs_strings)):
                    if (BQs_strings[setIdx].strip() == ""):
                        continue
                    writeNicsInputs(flprfx, repr(ring) + "-set" + repr(setIdx+1), flag_chk,
                                    hashLine_rev, title, charge, mult, new_geom, 
                                    BQs_strings[setIdx], jobType)

            elif (new_geom == []):
                print("Ring no. " + repr(ring) + " is not Planar within the tolerence of " + repr(TORSION_ANGLE_TOLERANCE) +
                      " degrees. Therefore, ring could not be reoriented in XY plane and BQs could not be generated.")

        # For "Normal"s defined by user
        n_count = len(CenterOf)
        for norm in normals:
            n_count += 1
            ring_atoms = []
            new_geom, new_points, new_normal = reorient(
                geom, Conn, ring_atoms, points, normals[norm])

            if ((new_geom != []) and (new_normal != [])):
                ring_atoms = []
                BQs_string = addBreakPoints(generateBQs_Z(
                    new_geom, Conn, ring_atoms, new_normal))

                BQs_strings = list(
                    map(lambda x: x.strip(), BQs_string.split("break")))

                for setIdx in range(0, len(BQs_strings)):
                    if (BQs_strings[setIdx].strip() == ""):
                        continue
                    writeNicsInputs(flprfx, repr(n_count) + "-set" + repr(setIdx+1), flag_chk,
                                    hashLine_rev, title, charge, mult, new_geom, 
                                    BQs_strings[setIdx], jobType)

            elif (new_geom == []):
                print("Ring no. " + repr(n_count) + " is not Planar within the tolerence of " + repr(TORSION_ANGLE_TOLERANCE) +
                      " degrees. Therefore, ring could not be reoriented in XY plane and BQs could not be generated.")

    elif (xy_flag):  # XY scan
        BQs_string, new_geom = generateBQs_XY(geom, Conn)
        BQs_string = addBreakPoints(BQs_string)
        BQs_strings = list(map(lambda x: x.strip(), BQs_string.split("break")))

        for cnt in range(0, len(BQs_strings)):
            if (BQs_strings[cnt].strip() == ""):
                continue
            writeNicsInputs(flprfx, "1-set" + repr(cnt), flag_chk, hashLine_rev,
                            title, charge, mult, new_geom, BQs_strings[cnt], jobType)


def addBreakPoints(BQs_string):
    global externalProgram
    # global geom length (number of atoms)
    global nAtoms

    bqs = BQs_string.strip().split("\n")

    atmCount = nAtoms

    newBQs_string = ""

    for bq in bqs:
        newBQs_string += bq + "\n"
        atmCount += 1
        if (atmCount % MAX_BQS_IN_INPFL == 0):
            newBQs_string += externalProgram["writerFunctCall"]["geomInput"].genGhostAtomSetBreak(
            ) + "\n"
            atmCount = nAtoms

    return newBQs_string


def generateBQs_Points(points):
    global externalProgram

    BQs_string = ""
    p_count = 0
    bq_count = 0

    for p in points:
        BQs_string += externalProgram["writerFunctCall"]["geomInput"].genGhostAtomLine(
            points[p_count, 0])
        p_count += 1
        bq_count += 1
        # if (bq_count%MAX_BQS_IN_INPFL == 0): BQs_string += externalProgram["writerFunctCall"]["geomInput"].genGhostAtomSetBreak()

    return BQs_string

# Here, the geom refers to the re-oriented geometry
# Now, the plane is always XY plane


def generateBQs_Z(geom, Conn, ring_atoms, normal=[]):
    global externalProgram

    sigma_model = 0
    direction_bq = 'POSITIVE'
    direction = 0
    H_count = 0

    for i in range(0, len(ring_atoms)):
        atm_idx = ring_atoms[i]

        if (not sigma_model):
            for j in range(0, len(Conn[atm_idx])):
                if (geom[Conn[atm_idx][j]][0] != 1):
                    continue
                else:
                    if ((geom[Conn[atm_idx][j]][3]) > 0.5):
                        H_count += 1
                        direction += 1
                    elif ((geom[Conn[atm_idx][j]][3]) < -0.5):
                        H_count += 1
                        direction -= 1

    # For Normal, the detection of Sigma-only Model is under Trial
    if ((len(ring_atoms) < 1) and (normal != [])):
        for i in range(0, len(geom)):
            atm_idx = i+1
            if (not sigma_model):
                for j in range(0, len(Conn[atm_idx])):
                    if (geom[Conn[atm_idx][j]][0] != 1):
                        continue
                    else:
                        if ((geom[Conn[atm_idx][j]][3]) > 0.5):
                            H_count += 1
                            direction += 1
                        elif ((geom[Conn[atm_idx][j]][3]) < -0.5):
                            H_count += 1
                            direction -= 1

    if (H_count > 2):
        sigma_model = 1
    if (direction > 0):
        direction_bq = 'NEGATIVE'

    if (len(ring_atoms) != 0):
        cmx, cmy, cmz = getGMOfRing(geom, ring_atoms)
    else:
        cmx, cmy, cmz = 0.0, 0.0, 0.0

    if (sigma_model):
        print("\nSigma-Only Model detected. \nTherefore, the BQs will be generated on the opposite side of the s-only H-atoms.")

    zinc = BQ_Step
    if (direction_bq == 'NEGATIVE'):
        zinc = -zinc

    coord_format = "{0:.5f}"
    BQs_string = ""
    zcoord = cmz + BQ_Range[0]
    for i in range(0, BQ_No):
        bqpt = [cmx, cmy, zcoord]
        BQs_string += externalProgram["writerFunctCall"]["geomInput"].genGhostAtomLine(
            bqpt)
        zcoord += zinc

    return BQs_string


def generateBQs_XY(geom, Conn):
    global externalProgram

    # global flags
    global opt_flag, ncs_flag, sigma_flag, xy_flag, pointonly_flag, integralnics_flag, analyse_flag, area_flag, s_charge_flag, s_mult_flag, opt_external, optfl_external, inponly_flag
    # global molecule-related
    global armpath, CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, xy_ref_ring_info, BQGuide, points, normals
    # global technical
    global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, xy_BQ_dist, sigma_charge, sigma_mult, analyse_dist, clear_flag, xy_extend

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
            new_geom, new_points, new_normal = reorient(
                geom, Conn, ring_atoms, points)

            if (new_geom != []):
                for i in range(0, len(ring_atoms)):
                    atm_idx = ring_atoms[i]

                    if (not sigma_model):
                        for j in range(0, len(Conn[atm_idx])):
                            if (new_geom[Conn[atm_idx][j]][0] != 1):
                                continue
                            else:
                                if ((new_geom[Conn[atm_idx][j]][3]) > 0.5):
                                    H_count += 1
                                    direction += 1
                                elif ((new_geom[Conn[atm_idx][j]][3]) < -0.5):
                                    H_count += 1
                                    direction -= 1
                if (H_count > 2):
                    sigma_model = 1
                    print(
                        "\nSigma-Only Model detected. \nTherefore, the BQs will be generated on the opposite side of the s-only H-atoms.")
                if (direction > 0):
                    direction_bq = 'NEGATIVE'

                cmx, cmy, cmz = getGMOfRing(new_geom, ring_atoms)
                ref_ring = ring
                xy_ref_ring_info.append(ring)
                xy_ref_ring_info.append(sorted(ring_atoms))
                normal_to_Ring = getUnitVector(getAverageNormaltoTheRing(
                    new_geom, ring_atoms, [cmx, cmy, cmz]))
                normal_to_Ring[:] = [
                    vc*DEFAULT_XY_DISTANCE for vc in normal_to_Ring]
                if (normal_to_Ring == [0.0, 0.0, 0.0]):
                    if (direction_bq == 'POSITIVE'):
                        direct = (
                            'set', [0.0, 0.0, DEFAULT_XY_DISTANCE] + [cmx, cmy, cmz])
                        BQGuide[ref_ring] = (
                            [cmx, cmy, cmz], [0.0, 0.0, DEFAULT_XY_DISTANCE])
                        norm_vec[ref_ring] = [cmx, cmy, cmz,
                                              0.0, 0.0, DEFAULT_XY_DISTANCE]
                    elif (direction_bq == 'NEGATIVE'):
                        direct = (
                            'set', [0.0, 0.0, -DEFAULT_XY_DISTANCE] + [cmx, cmy, cmz])
                        BQGuide[ref_ring] = (
                            [cmx, cmy, cmz], [0.0, 0.0, -DEFAULT_XY_DISTANCE])
                        norm_vec[ref_ring] = [cmx, cmy, cmz,
                                              0.0, 0.0, -DEFAULT_XY_DISTANCE]
                else:
                    if (direction_bq == 'POSITIVE'):
                        direct = ('set', [
                                  normal_to_Ring[0]+cmx, normal_to_Ring[1]+cmy, normal_to_Ring[2]+cmz, cmx, cmy, cmz])
                        BQGuide[ref_ring] = (
                            [cmx, cmy, cmz], [+normal_to_Ring[0]+cmx, +normal_to_Ring[1]+cmy, +normal_to_Ring[2]+cmz])
                        norm_vec[ref_ring] = [cmx, cmy, cmz, +normal_to_Ring[0] +
                                              cmx, +normal_to_Ring[1]+cmy, +normal_to_Ring[2]+cmz]
                    elif (direction_bq == 'NEGATIVE'):
                        direct = (
                            'set', [-normal_to_Ring[0]+cmx, -normal_to_Ring[1]+cmy, -normal_to_Ring[2]+cmz, cmx, cmy, cmz])
                        BQGuide[ref_ring] = (
                            [cmx, cmy, cmz], [-normal_to_Ring[0]+cmx, -normal_to_Ring[1]+cmy, -normal_to_Ring[2]+cmz])
                        norm_vec[ref_ring] = [cmx, cmy, cmz, -normal_to_Ring[0] +
                                              cmx, -normal_to_Ring[1]+cmy, -normal_to_Ring[2]+cmz]
                break

            elif (new_geom == []):
                print("None of the rings is Planar within the tolerence of " +
                      repr(TORSION_ANGLE_TOLERANCE) + " degrees. Therefore, aboring the run .. ")
                sys.exit(10)

    xy_ref_ring_info.append(direct[1])

    for ring in CenterOf:
        if (ring != ref_ring):
            ring_atoms = getOrderedRing(Conn, CenterOf.get(ring))
            if (len(ring_atoms) > 2):
                cmx, cmy, cmz = getGMOfRing(new_geom, ring_atoms)
                normal_to_Ring = getUnitVector(getAverageNormaltoTheRing(
                    new_geom, ring_atoms, [cmx, cmy, cmz]))
                normal_to_Ring[:] = [
                    vc*DEFAULT_XY_DISTANCE for vc in normal_to_Ring]
                check_direction = getDihedralAngleBetween(direct[1][0:3], direct[1][3:6], [cmx, cmy, cmz], [
                                                          normal_to_Ring[0]+cmx, normal_to_Ring[1]+cmy, normal_to_Ring[2]+cmz])
                # Tolerence of 15 degrees seem to be reasonable
                # But it is now updated to 90 which is maximum, for the case of twisted rings
                if (abs(check_direction) < 90):
                    BQGuide[ring] = ([cmx, cmy, cmz], [+normal_to_Ring[0] +
                                     cmx, +normal_to_Ring[1]+cmy, +normal_to_Ring[2]+cmz])
                    norm_vec[ring] = [cmx, cmy, cmz, +normal_to_Ring[0] +
                                      cmx, +normal_to_Ring[1]+cmy, +normal_to_Ring[2]+cmz]
                else:
                    BQGuide[ring] = ([cmx, cmy, cmz], [-normal_to_Ring[0] +
                                     cmx, -normal_to_Ring[1]+cmy, -normal_to_Ring[2]+cmz])
                    norm_vec[ring] = [cmx, cmy, cmz, -normal_to_Ring[0] +
                                      cmx, -normal_to_Ring[1]+cmy, -normal_to_Ring[2]+cmz]

    ring = 1
    if (ref_ring != 1):
        ring_atoms = getOrderedRing(Conn, CenterOf.get(ring))
        cmx, cmy, cmz = getGMOfRing(new_geom, ring_atoms)
        vec1 = norm_vec[ref_ring]
        BQGuide[1] = ([cmx, cmy, cmz], [vec1[3]-vec1[0]+cmx,
                      vec1[4]-vec1[1]+cmy, vec1[5]-vec1[2]+cmz])
        norm_vec[ring] = [cmx, cmy, cmz, vec1[3]-vec1[0] +
                          cmx, vec1[4]-vec1[1]+cmy, vec1[5]-vec1[2]+cmz]

    for ring in CenterOf:
        if ((ring != 1) and (ring != len(CenterOf))):
            ring_atoms = getOrderedRing(Conn, CenterOf.get(ring))
            if (len(ring_atoms) < 3):
                cmx, cmy, cmz = getGMOfRing(new_geom, ring_atoms)

                prv = ref_ring
                nxt = ref_ring
                for p in range(ring-1, 1, -1):
                    if (len(CenterOf.get(p)) > 3):
                        prv = p
                        break
                for p in range(ring+1, len(CenterOf)-1):
                    if (len(CenterOf.get(p)) > 3):
                        nxt = p
                        break
                if (nxt == ref_ring):
                    nxt == prv

                vec1 = [norm_vec[prv][3] - norm_vec[prv][0] + cmx, norm_vec[prv][4] -
                        norm_vec[prv][1] + cmy, norm_vec[prv][5] - norm_vec[prv][2] + cmz]
                vec2 = [norm_vec[nxt][3] - norm_vec[nxt][0] + cmx, norm_vec[nxt][4] -
                        norm_vec[nxt][1] + cmy, norm_vec[nxt][5] - norm_vec[nxt][2] + cmz]
                mid_vec = [(vec1[0]+vec2[0])/2.,
                           (vec1[1]+vec2[1])/2., (vec1[2]+vec2[2])/2.]
                BQGuide[ring] = ([cmx, cmy, cmz], mid_vec)
                norm_vec[ring] = [cmx, cmy, cmz,
                                  mid_vec[0], mid_vec[1], mid_vec[2]]

    ring = len(CenterOf)
    ring_atoms = getOrderedRing(Conn, CenterOf.get(ring))
    if (len(ring_atoms) < 3):
        # Did you see a heart somewhere nearby ?
        cmx, cmy, cmz = getGMOfRing(new_geom, ring_atoms)
        vec2 = norm_vec[ring-1]
        BQGuide[ring] = ([cmx, cmy, cmz], [vec2[3]-vec2[0]+cmx,
                         vec2[4]-vec2[1]+cmy, vec2[5]-vec2[2]+cmz])
        norm_vec[ring] = [cmx, cmy, cmz, vec2[3]-vec2[0] +
                          cmx, vec2[4]-vec2[1]+cmy, vec2[5]-vec2[2]+cmz]

    for pt in new_points:
        cmx, cmy, cmz = new_points[pt][1:4]
        prv = pt[0]
        nxt = pt[1]
        vec1 = [norm_vec[prv][3] - norm_vec[prv][0] + cmx, norm_vec[prv][4] -
                norm_vec[prv][1] + cmy, norm_vec[prv][5] - norm_vec[prv][2] + cmz]
        vec2 = [norm_vec[nxt][3] - norm_vec[nxt][0] + cmx, norm_vec[nxt][4] -
                norm_vec[nxt][1] + cmy, norm_vec[nxt][5] - norm_vec[nxt][2] + cmz]
        mid_vec = [(vec1[0]+vec2[0])/2., (vec1[1]+vec2[1]) /
                   2., (vec1[2]+vec2[2])/2.]
        BQGuide[(prv+nxt)/2.0] = (new_points[pt][1:4], mid_vec)

    coord_format = "{0:.5f}"
    BQs_string = ""
    bq_count = 0
    seq_BQGuide = sorted(BQGuide)

    bq_coord_prvs = BQGuide[seq_BQGuide[0]][1]
    for i in range(0, len(seq_BQGuide)-1):

        # If the XY trajectory needs to be extended out, then xy_extend is non-zero and BQs need to be added outside of the molecule
        if ((i == 0.0) and (xy_extend != 0.0)):
            a = BQGuide[seq_BQGuide[0]][1]
            b = BQGuide[seq_BQGuide[1]][1]
            vec_ba = getVector(b, a)
            n_vec_ba = getUnitVector(vec_ba)
            # c is the point outside the molecule from where the trajectory begins
            c = [xy_extend*v for v in n_vec_ba]
            c[:] = [a[v]+c[v] for v in range(0, 3)]
            vec_ca = getVector(c, a)
            n_vec_ca = getUnitVector(vec_ca)
            n_BQ = int(round((xy_extend/BQ_Step), 0))
            bq_coord_prvs = c

            for j in range(0, n_BQ):
                bq_coord = [BQ_Step*j*v for v in n_vec_ca]
                bq_coord[:] = [c[v]+bq_coord[v] for v in range(0, 3)]
                xy_BQ_dist.append(
                    round(getDistance(bq_coord, bq_coord_prvs), 3))
                BQs_string += externalProgram["writerFunctCall"]["geomInput"].genGhostAtomLine(
                    bq_coord)
                bq_count += 1
                bq_coord_prvs = bq_coord
                # if (bq_count%MAX_BQS_IN_INPFL == 0): BQs_string += "break"

        a = BQGuide[seq_BQGuide[i]][1]
        b = BQGuide[seq_BQGuide[i+1]][1]
        vec_ab = getVector(a, b)
        norm_vec_ab = vectorMagnitude(vec_ab)
        n_vec_ab = getUnitVector(vec_ab)
        n_BQ = int(round((norm_vec_ab/BQ_Step), 0))

        for j in range(0, n_BQ):
            bq_coord = [BQ_Step*j*v for v in n_vec_ab]
            bq_coord[:] = [a[v]+bq_coord[v] for v in range(0, 3)]
            xy_BQ_dist.append(round(getDistance(bq_coord, bq_coord_prvs), 3))
            BQs_string += externalProgram["writerFunctCall"]["geomInput"].genGhostAtomLine(
                bq_coord)
            bq_count += 1
            bq_coord_prvs = bq_coord
            # if (bq_count%MAX_BQS_IN_INPFL == 0): BQs_string += "break"

    # For the last BQ
    if (vectorMagnitude(getVector(bq_coord, b)) < BQ_Step):
        bq_coord = [BQ_Step*(j+1)*v for v in n_vec_ab]
        bq_coord[:] = [a[v]+bq_coord[v] for v in range(0, 3)]
        xy_BQ_dist.append(round(getDistance(bq_coord, bq_coord_prvs), 3))
        BQs_string += externalProgram["writerFunctCall"]["geomInput"].genGhostAtomLine(
            bq_coord)
        bq_count += 1
        bq_coord_prvs = bq_coord

    # If xy_extend is non-zero then more BQs need to be added at the end of the trajectory too
    n_vec_ab = getUnitVector(vec_ab)
    n_BQ = int(round((xy_extend/BQ_Step), 0))

    for j in range(0, n_BQ):
        bq_coord = [BQ_Step*j*v for v in n_vec_ab]
        bq_coord[:] = [b[v]+bq_coord[v] for v in range(0, 3)]
        xy_BQ_dist.append(round(getDistance(bq_coord, bq_coord_prvs), 3))
        BQs_string += externalProgram["writerFunctCall"]["geomInput"].genGhostAtomLine(
            bq_coord)
        bq_count += 1
        bq_coord_prvs = bq_coord
        # if (bq_count%MAX_BQS_IN_INPFL == 0): BQs_string += "break"

    BQ_No = bq_count
    return BQs_string, new_geom


def run_Nics():
    global externalProgram

    # global flags
    global opt_flag, ncs_flag, sigma_flag, xy_flag, pointonly_flag, integralnics_flag, analyse_flag, area_flag, s_charge_flag, s_mult_flag, opt_external, optfl_external, inponly_flag
    # global molecule-related
    global armpath, CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, xy_ref_ring_info, BQGuide, points, normals
    # global technical
    global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, xy_BQ_dist, sigma_charge, sigma_mult, analyse_dist, clear_flag, xy_extend
    # global inputFileSet
    global inputFileSet

    if (not xy_flag):
        dict_cen = CenterOf.copy()
        if (len(normals) > 0):
            n_count = len(CenterOf)
            for i in normals:
                dict_cen[n_count+i] = normals[i]
    elif (xy_flag):
        dict_cen = n_xy_center

    if not pointonly_flag:
        for inpfl in inputFileSet:
            flname = inpfl["flprfx"]
            print(flname)
            print("Job " + externalProgram["extCmd"] +
                  flname + " " + flname + " running ..")
            status = execCmd(externalProgram["constructCmd"](flname))
            if (not status):
                print("Job Over.")
            else:
                print("It seems that the NICS SCAN run for ring/bond [" +
                      flname + "] did not terminated normally.\n")

    if (pointonly_flag):
        rcount = int((len(points)/MAX_BQS_IN_INPFL))+1
        for inpfl in inputFileSet:
            flname = inpfl["filename"]
            print("Job " + externalProgram["extCmd"] +
                  flname + " " + flname + " running ..")
            status = execCmd(externalProgram["constructCmd"](flname))
            if (not status):
                print("Job Over.")
            else:
                print("It seems that the NICS SCAN run for ring/bond number " +
                      repr(ring) + " did not terminated normally.\n")


def genSigmaModel(flprfx, geom, Conn, title, charge, mult):
    global sigma_charge, sigma_mult, xy_flag, xy_ref_ring_info, normals, points, exocyclic

    count = len(geom)+1
    sigma_geom = {}

    # This is a complicated dictionary which consists of keys which are tuples and values which is charge on H atom to be added
    # The tuple in key denotes (Atomic number, number of connections)
    # All this is to determine the charge on the Sigma model.
    for_sigma_charge = {(5, 3): -1, (6, 3): 0, (7, 3): +1,
                        (7, 2): 0, (8, 2): +1, (14, 3): 0, (16, 2): +1}
    aroma_sigma_charge = 0
    if (s_charge_flag):
        user_sigma_charge = sigma_charge

    for i in range(0, len(geom)):
        sigma_geom[i+1] = geom[i+1]

    # Identify the fused bonds i.e. all pairs of fused rings
    # Check whether the pair of connected atoms is present in more than one rings.
    fused_bonds = []
    fused_rings = []
    for i in range(1, len(all_aromatic_rings)+1):
        for j in range(i+1, len(all_aromatic_rings)+1):
            intersecting_bond = list(
                set(all_aromatic_rings[i]) & set(all_aromatic_rings[j]))
            if ((intersecting_bond != []) and (len(intersecting_bond) == 2)):
                fused_bonds.append(intersecting_bond)
                fused_rings.append([i, j])

            intersecting_bond = []

    # Identify intersections of fused bonds
    # i.e. 3 rings fused together
    fused_3 = []
    fused_3_r = []
    to_be_removed = []
    for i in range(0, len(fused_bonds)):
        for j in range(i+1, len(fused_bonds)):
            for k in range(j+1, len(fused_bonds)):
                c = set(fused_bonds[i]) & set(
                    fused_bonds[j]) & set(fused_bonds[k])
                if (c != set([])):
                    fused_3.append(list(set(fused_bonds[i]) | set(
                        fused_bonds[j]) | set(fused_bonds[k])))
                    fused_3_r.append(list(set(fused_rings[i]) | set(
                        fused_rings[j]) | set(fused_rings[k])))
                    if (to_be_removed.count(fused_bonds[i]) < 1):
                        to_be_removed.append(fused_bonds[i])
                    if (to_be_removed.count(fused_bonds[j]) < 1):
                        to_be_removed.append(fused_bonds[j])
                    if (to_be_removed.count(fused_bonds[k]) < 1):
                        to_be_removed.append(fused_bonds[k])

    for i, j in reversed(list(enumerate(to_be_removed))):
        fused_bonds.remove(j)
        fused_rings.remove(fused_rings[i])

    # Taking back up of All Aromatic Rings as per .arm
    all_aromatic_rings_local = all_aromatic_rings.copy()

    indicator_for_fused_2 = len(all_aromatic_rings_local)
    r_idx = len(all_aromatic_rings_local)
    for i in range(r_idx+1, r_idx+len(fused_bonds)+1):
        all_aromatic_rings_local[i] = (fused_bonds[i-r_idx-1])

    indicator_for_fused_3 = len(all_aromatic_rings_local)
    r_idx = len(all_aromatic_rings_local)
    for i in range(r_idx+1, r_idx+len(fused_3)+1):
        all_aromatic_rings_local[i] = (fused_3[i-r_idx-1])

    ring_dummy_tuples = {}
    direct = ''
    ring_count = 0
    norm_vec = {}
    for ring in all_aromatic_rings_local:
        ring_count += 1
        if (ring_count < indicator_for_fused_3):
            ring_atoms = getOrderedRing(
                Conn, all_aromatic_rings_local.get(ring))
        else:
            ring_atoms = all_aromatic_rings_local.get(ring)
        cmx, cmy, cmz = getGMOfRing(geom, ring_atoms)

        # A dummy atom at the center of the Ring is added
        sigma_geom[count] = [0, cmx, cmy, cmz]
        count += 1
        # Another dummy atom is added from 1 angstrom distance from the CM of the ring in the user-specified direction perpendicular to the ring
        unit_normal_to_Ring = getUnitVector(
            getAverageNormaltoTheRing(geom, ring_atoms, [cmx, cmy, cmz]))

        if (ring_count <= indicator_for_fused_2):
            if (unit_normal_to_Ring == [0.0, 0.0, 0.0]):
                if (direct == ''):
                    direct = ('set', [0.0, 0.0, 1.0] + [cmx, cmy, cmz])
                    sigma_geom[count] = [0, 0.0, 0.0, 0.1]
                    norm_vec[ring_count] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.1]
                    if (xy_flag and (ring_atoms.sort() == xy_ref_ring_info[1])):
                        check_direction = getDihedralAngleBetween(
                            direct[1][0:3], direct[1][3:6], xy_ref_ring_info[2][3:6], xy_ref_ring_info[2][0:3])
                        if (abs(check_direction) < 90):
                            direct = (
                                'set', [0.0, 0.0, -1.0] + [cmx, cmy, cmz])
                            sigma_geom[count] = [0, 0.0, 0.0, -0.1]
                            norm_vec[ring_count] = [
                                0.0, 0.0, 0.0, 0.0, 0.0, -0.1]
                else:
                    sigma_geom[count] = [0, direct[1][0]+cmx,
                                         direct[1][1]+cmy, direct[1][2]+cmz]
                    norm_vec[ring_count] = [cmx, cmy, cmz, direct[1]
                                            [0]+cmx, direct[1][1]+cmy, direct[1][2]+cmz]

            elif (direct == ''):
                direct = ('set', [unit_normal_to_Ring[0]+cmx, unit_normal_to_Ring[1] +
                          cmy, unit_normal_to_Ring[2]+cmz, cmx, cmy, cmz])
                sigma_geom[count] = [0, +unit_normal_to_Ring[0]+cmx, +
                                     unit_normal_to_Ring[1]+cmy, +unit_normal_to_Ring[2]+cmz]
                norm_vec[ring_count] = [cmx, cmy, cmz, +unit_normal_to_Ring[0] +
                                        cmx, +unit_normal_to_Ring[1]+cmy, +unit_normal_to_Ring[2]+cmz]
                if (xy_flag and (ring_atoms.sort() == xy_ref_ring_info[1])):
                    check_direction = getDihedralAngleBetween(
                        direct[1][0:3], direct[1][3:6], xy_ref_ring_info[2][3:6], xy_ref_ring_info[2][0:3])
                    if (abs(check_direction) < 90):
                        direct = ('set', [-unit_normal_to_Ring[0]+cmx, -unit_normal_to_Ring[1] +
                                  cmy, -unit_normal_to_Ring[2]+cmz, cmx, cmy, cmz])
                        sigma_geom[count] = [0, -unit_normal_to_Ring[0]+cmx, -
                                             unit_normal_to_Ring[1]+cmy, -unit_normal_to_Ring[2]+cmz]
                        norm_vec[ring_count] = [cmx, cmy, cmz, -unit_normal_to_Ring[0] +
                                                cmx, -unit_normal_to_Ring[1]+cmy, -unit_normal_to_Ring[2]+cmz]

            else:
                check_direction = getDihedralAngleBetween(direct[1][0:3], direct[1][3:6], [cmx, cmy, cmz], [
                                                          unit_normal_to_Ring[0]+cmx, unit_normal_to_Ring[1]+cmy, unit_normal_to_Ring[2]+cmz])
                # Tolerence of 15 degrees seem to be reasonable
                # But it is now updated to 90 which is maximum, for the case of twisted rings
                if (abs(check_direction) < 90):
                    sigma_geom[count] = [0, +unit_normal_to_Ring[0]+cmx, +
                                         unit_normal_to_Ring[1]+cmy, +unit_normal_to_Ring[2]+cmz]
                    norm_vec[ring_count] = [cmx, cmy, cmz, +unit_normal_to_Ring[0] +
                                            cmx, +unit_normal_to_Ring[1]+cmy, +unit_normal_to_Ring[2]+cmz]
                else:
                    sigma_geom[count] = [0, -unit_normal_to_Ring[0]+cmx, -
                                         unit_normal_to_Ring[1]+cmy, -unit_normal_to_Ring[2]+cmz]
                    norm_vec[ring_count] = [cmx, cmy, cmz, -unit_normal_to_Ring[0] +
                                            cmx, -unit_normal_to_Ring[1]+cmy, -unit_normal_to_Ring[2]+cmz]

        elif ((ring_count > indicator_for_fused_2) and (ring_count <= indicator_for_fused_3)):
            # Now all these are bonds, so unit vector is not well-defined. It has to be alightned.
            vec1 = norm_vec[fused_rings[ring_count -
                                        indicator_for_fused_2 - 1][0]]
            vec2 = norm_vec[fused_rings[ring_count -
                                        indicator_for_fused_2 - 1][1]]
            mid_vec = [(vec1[3] - vec1[0] + vec2[3] - vec2[0])/2., (vec1[4] - vec1[1] +
                                                                    vec2[4] - vec2[1])/2., (vec1[5] - vec1[2] + vec2[5] - vec2[2])/2.]
            sigma_geom[count] = [0, mid_vec[0] + cmx,
                                 mid_vec[1] + cmy, mid_vec[2] + cmz]

        elif (ring_count > indicator_for_fused_3):
            vec1 = norm_vec[fused_3_r[ring_count -
                                      indicator_for_fused_3 - 1][0]]
            vec2 = norm_vec[fused_3_r[ring_count -
                                      indicator_for_fused_3 - 1][1]]
            vec3 = norm_vec[fused_3_r[ring_count -
                                      indicator_for_fused_3 - 1][2]]
            mid_vec = [(vec1[3] - vec1[0] + vec2[3] - vec2[0] + vec3[3] - vec3[0])/3., (vec1[4] - vec1[1] + vec2[4] -
                                                                                        vec2[1] + vec3[4] - vec3[1])/3., (vec1[5] - vec1[2] + vec2[5] - vec2[2] + vec3[5] - vec3[2])/3.]
            sigma_geom[count] = [0, mid_vec[0] + cmx,
                                 mid_vec[1] + cmy, mid_vec[2] + cmz]

        ring_dummy_tuples[ring] = (count-1, count)
        count += 1

#   n_dummy_couple = len(ring_dummy_tuples)
    mapec = {}
    for ec in exocyclic:
        # Check which ring the exocyclic atom is connected
        for r in all_aromatic_rings:
            ec_atom = list(
                set(exocyclic[ec]) - (set(exocyclic[ec]) & set(all_aromatic_rings[r])))
            if (len(ec_atom) == 1):
                mapec[r] = ec_atom[0]
            # So here the mapec maps the aromatic ring with its exocyclic atom
#
#      n_dummy_couple += 1
#      ring_dummy_tuples[n_dummy_couple] =

    # Flags to make sure that H-atoms are not added to a same atom twice !!
    flag_H = {}
    for i in range(1, len(geom)+1):
        flag_H[i] = 0

    sigma_conn_mat, sigma_Conn = genConnectivityMatrix(sigma_geom)
    zmat, zmat_str, zmat_idx = generateZMatrix(sigma_geom, sigma_Conn)

    no_electrons = 0
    flag_for_warning = 0

    ring = len(all_aromatic_rings_local) - len(fused_bonds) - len(fused_3_r)
    for ring_atoms in fused_bonds:
        ring += 1
        for i in range(0, len(ring_atoms)):
            j = ring_atoms[i]
            if ((flag_H[j] != 1) and (len(Conn[j]) != Max_Conn[geom[j][0]])):
                zmat_str += "1  " + repr(zmat_idx[j]) + "  " + ATM_H_BL[geom[j][0]] + "  " + repr(zmat_idx[ring_dummy_tuples[ring][0]]) + \
                    "  " + FIXED_SIGMA_ANGLE + "  " + \
                    repr(zmat_idx[ring_dummy_tuples[ring][1]]) + \
                    "  " + FIXED_SIGMA_DIHEDRAL_ANGLE + "\n"
                flag_H[j] = 1
                no_electrons += 1
                sigma_key = (geom[j][0], len(Conn[j]))
                if (sigma_key in for_sigma_charge):
                    aroma_sigma_charge += for_sigma_charge[sigma_key]
                else:
                    if (not flag_for_warning):
                        print("\nWARNING: Aroma does not have enough data for counting charge correctly for Atom with atomic no. ",
                              geom[j][0], "\nIt is advised to give correct charge on Sigma model externally\n")
                        flag_for_warning = 1

    for ring_atoms in fused_3:

        for i in range(0, len(ring_atoms)):
            j = ring_atoms[i]
            if ((flag_H[j] != 1) and (len(Conn[j]) != Max_Conn[geom[j][0]])):
                # Here the angle for sigma H atom will be always 90.0 as it has to be symmetric with respect to all the fused bonds
                zmat_str += "1  " + repr(zmat_idx[j]) + "  " + ATM_H_BL[geom[j][0]] + "  " + repr(zmat_idx[ring_dummy_tuples[ring][0]]) + \
                    "  " + "90.0" + "  " + \
                    repr(zmat_idx[ring_dummy_tuples[ring][1]]) + \
                    "  " + FIXED_SIGMA_DIHEDRAL_ANGLE + "\n"
                flag_H[j] = 1
                no_electrons += 1
                sigma_key = (geom[j][0], len(Conn[j]))
                if (sigma_key in for_sigma_charge):
                    aroma_sigma_charge += for_sigma_charge[sigma_key]
                else:
                    if (not flag_for_warning):
                        print("\nWARNING: Aroma does not have enough data for counting charge correctly for Atom with atomic no. ",
                              geom[j][0], "\nIt is advised to give correct charge on Sigma model externally\n")
                        flag_for_warning = 1

    for ring in all_aromatic_rings_local:
        if (ring_count < indicator_for_fused_3):
            ring_atoms = getOrderedRing(
                Conn, all_aromatic_rings_local.get(ring))
        else:
            ring_atoms = all_aromatic_rings_local.get(ring)

        for i in range(0, len(ring_atoms)):
            j = ring_atoms[i]
            if ((flag_H[j] != 1) and (len(Conn[j]) != Max_Conn[geom[j][0]])):
                zmat_str += "1  " + repr(zmat_idx[j]) + "  " + ATM_H_BL[geom[j][0]] + "  " + repr(zmat_idx[ring_dummy_tuples[ring][0]]) + \
                    "  " + FIXED_SIGMA_ANGLE + "  " + \
                    repr(zmat_idx[ring_dummy_tuples[ring][1]]) + \
                    "  " + FIXED_SIGMA_DIHEDRAL_ANGLE + "\n"
                flag_H[j] = 1
                no_electrons += 1
                sigma_key = (geom[j][0], len(Conn[j]))
                if (sigma_key in for_sigma_charge):
                    aroma_sigma_charge += for_sigma_charge[sigma_key]
                else:
                    if (not flag_for_warning):
                        print("\nWARNING: Aroma does not have enough data for counting charge correctly for Atom with atomic no. ",
                              geom[j][0], "\nIt is advised to give correct charge on Sigma model externally\n")
                        flag_for_warning = 1

    for r in mapec:
        zmat_str += "1  " + repr(zmat_idx[mapec[r]]) + "  " + ATM_H_BL[geom[mapec[r]][0]] + "  " + repr(zmat_idx[ring_dummy_tuples[r][0]]) + \
            "  " + FIXED_SIGMA_ANGLE + "  " + \
            repr(zmat_idx[ring_dummy_tuples[r][1]]) + \
            "  " + FIXED_SIGMA_DIHEDRAL_ANGLE + "\n"

# ADD: Z-MATRIX ELEMENTS FOR NORMALS
    for norm in normals:
        zmat_str += "0  " + repr(zmat_idx[1]) + "   " + repr(getDistance(geom[1][1:4], normals[norm][0:3])) + "   " + repr(zmat_idx[2]) + "   " + repr(getAngleBetweenPoints(
            normals[norm][0:3], geom[1][1:4], geom[2][1:4])) + "   " + repr(zmat_idx[3]) + "   " + repr(getDihedralAngleBetween(normals[norm][0:3], geom[1][1:4], geom[2][1:4], geom[3][1:4])) + "\n"
        zmat_str += "0  " + repr(zmat_idx[1]) + "   " + repr(getDistance(geom[1][1:4], normals[norm][3:6])) + "   " + repr(zmat_idx[2]) + "   " + repr(getAngleBetweenPoints(
            normals[norm][3:6], geom[1][1:4], geom[2][1:4])) + "   " + repr(zmat_idx[3]) + "   " + repr(getDihedralAngleBetween(normals[norm][3:6], geom[1][1:4], geom[2][1:4], geom[3][1:4])) + "\n"

    for p in points:
        zmat_str += "0 " + repr(zmat_idx[1]) + "   " + repr(getDistance(geom[1][1:4], points[p][0:3])) + "   " + repr(zmat_idx[2]) + "   " + repr(getAngleBetweenPoints(
            points[p][0:3], geom[1][1:4], geom[2][1:4])) + "   " + repr(zmat_idx[3]) + "   " + repr(getDihedralAngleBetween(points[p][0:3], geom[1][1:4], geom[2][1:4], geom[3][1:4])) + "\n"

    aroma_sigma_charge += charge
    no_electrons += -aroma_sigma_charge
    if (no_electrons % 2 != 0):
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
            print("\nWARNING: User Specifed Charge on Sigma Model is " + repr(user_sigma_charge) + " while that determined by Aroma is " +
                  repr(aroma_sigma_charge) + "\nHowever, Aroma will proceed with Charge Provided by the User.\n")
            sigma_charge = user_sigma_charge
    else:
        sigma_charge = aroma_sigma_charge

    """
   sigma_flprfx = externalProgram["inpdir"] + flprfx + "-sigma" + externalProgram["inpExt"]
   f = open (sigma_flprfx, "w")
   f.write("# \n\n" + title + " sigma only model " + "\n\n" + repr(sigma_charge) + " " + repr(sigma_mult) + "\n")
   f.write(zmat_str + "\n") 
   f.close()
   """

    newZmat = zmat_str.split("\n")
    for idx in range(len(newZmat)):
        words = newZmat[idx].split()

        if (len(words) == 1):
            newZmat[idx] = [int(words[0])]
        elif (len(words) == 3):
            newZmat[idx] = [int(words[0]), int(words[1]), float(words[2])]
        elif (len(words) == 5):
            newZmat[idx] = [int(words[0]), int(words[1]), float(
                words[2]), int(words[3]), float(words[4])]
        elif (len(words) == 7):
            newZmat[idx] = [int(words[0]), int(words[1]), float(words[2]), int(
                words[3]), float(words[4]), int(words[5]), float(words[6])]

    sigma_geom = generateCartesianFromZmat(newZmat)

    return sigma_geom, sigma_charge, sigma_mult, zmat_idx


def getNewRings(geom, sigma_geom, CenterOf, zmat_idx):
    new_CenterOf = {}
    for i in CenterOf:
        org_ring = CenterOf.get(i)
        new_ring = []
        for j in range(0, len(org_ring)):
            k = zmat_idx[org_ring[j]]
            new_ring.append(k)
        new_CenterOf[i] = new_ring
    return new_CenterOf


def grepData():

    # global flags
    global opt_flag, ncs_flag, sigma_flag, xy_flag, pointonly_flag, integralnics_flag, analyse_flag, area_flag, s_charge_flag, s_mult_flag, opt_external, optfl_external, inponly_flag
    # global molecule-related
    global armpath, CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, xy_ref_ring_info, BQGuide, points, normals
    # global technical
    global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, xy_BQ_dist, sigma_charge, sigma_mult, analyse_dist, clear_flag, xy_extend
    # global inputFileSet
    global inputFileSet

    # this needs to be inited outside the input set loop to avoid re-initialisation
    if (not xy_flag):
       dist = BQ_Range[0]

    idx = 0

    # This loop is over each output file to filter the required tensor data
    for inpFil in inputFileSet:
        outfl = externalProgram["outdir"] + \
            inpFil["flprfx"] + externalProgram["outExt"]

        nat = inpFil["nat"]
        nBQ = inpFil["nBq"]

        theParser = externalProgram["readerFunctCall"]["output"](outfl)
        bqTensors = theParser.getMagneticTensorData(nat, nBQ, outfl)

        f_out = open(externalProgram["outdir"] +
                     inpFil["flprfx"] + ".armdat", "w")
        f_out.write(
            "#       oop       in1        in2       inp       iso        x         y         z\n")

        # i + 5*nat + 1, i + 5*nat + BQ_No*nat + 1 are start and end line numbers in the output file for the data for BQs

        for j in range(0, len(bqTensors)):  # TODO: check correctness, nBQ is replaced with len(bqTensors)

            BQ_data_string = ""
            if xy_flag:
                dist += xy_BQ_dist[idx]

            # This is commented as BQ no is not important, instead distance of that BQ from GM is added (further)
            # BQ_data_string += words[0] + "   "


            iso, xx, yy, zz, e1, e2, e3 = bqTensors[j]
            sorted_e = []
            sorted_e.append(e1)
            close = abs(e1 + zz)
            if (abs(e2 + zz) <= close):
                close = abs(e2 + zz)
                sorted_e.insert(0, e2)
            elif (abs(e2 + zz) > close):
                sorted_e.append(e2)
            if (abs(e3 + zz) <= close):
                sorted_e.insert(0, e3)
            elif (abs(e3 + zz) > close):
                sorted_e.append(e3)
            oup = -sorted_e[0]
            inp1 = -sorted_e[1]
            inp2 = -sorted_e[2]
            # In-plane chemical shift is average of the two in-plane shifts
            inp = 0.5*(inp1+inp2)

            BQ_data_string += fpformat.fix(dist, 2) + "   " + fpformat.fix(oup, 4) + "   " + fpformat.fix(inp1, 4) + "   " + fpformat.fix(inp2, 4) + "   " + fpformat.fix(
                inp, 4) + "   " + fpformat.fix(iso, 4) + "   " + fpformat.fix(xx, 4) + "   " + fpformat.fix(yy, 4) + "   " + fpformat.fix(zz, 4) + "\n"
            f_out.write(BQ_data_string)

            if (not xy_flag):
                dist += BQ_Step
            else:
                idx += 1

        f_out.close()

    # write collated set of files
    if ((len(inputFileSet) > 1)):
        writeCollatedFiles("main")
        writeCollatedFiles("sigma")

    if (xy_flag):
        shutil.move(externalProgram["outdir"] + flprfx + "-center1" + ".armdat",
                    externalProgram["outdir"] + flprfx + "-allcenter" + ".armdat")

def writeCollatedFiles(jobType):
    global inputFileSet, collatedFileSet

    # get all files in the "main" run
    mainFiles = list(filter(lambda x: x["jobType"] == jobType and x["setIdx"] != "-1", inputFileSet))

    # get list of all centers
    centerList = list(set(list(map(lambda x: x["centerIdx"], mainFiles))))

    for centerIdx in centerList:
      idx = 0
      centerFileSet = list(filter(lambda x: x["centerIdx"] == centerIdx and x["jobType"] == jobType, mainFiles))

      baseprfx = centerFileSet[0]["baseprfx"] + "-center" + centerIdx

      final_armdat = open(externalProgram["outdir"] + baseprfx + ".armdat", "w")

      collatedFileSet.append({"fileName": final_armdat, "flprfx": baseprfx, "jobType": jobType})

      for inpFil in centerFileSet:
        lines = readFile(externalProgram["outdir"] + inpFil["flprfx"] + ".armdat")
        for i in range(idx, len(lines)):
            final_armdat.write(lines[i])
        idx = 1

      final_armdat.close()


def generateAllInputs(geom, title, charge, mult, Conn, jobType):
    # global flags
    global opt_flag, ncs_flag, sigma_flag, xy_flag, pointonly_flag, integralnics_flag, analyse_flag, area_flag, s_charge_flag, s_mult_flag, opt_external, optfl_external, inponly_flag
    # global molecule-related
    global armpath, CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, xy_ref_ring_info, BQGuide, points, normals
    # global technical
    global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, xy_BQ_dist, sigma_charge, sigma_mult, analyse_dist, clear_flag, xy_extend
    # global geom length (number of atoms)
    global nAtoms

    # Run optimization, if required
    if (hashLine_opt == "DEFAULT\n"):
        hashLine_opt = externalProgram.defaultOptimizationKeyline
    if (opt_flag):
        if (not opt_external):
            optfl = generate_Opt_Input(geom, hashLine_opt, title, charge, mult)
        else:
            optfl = flprfx + "-opt"
            print("\nCopying external optimization input file to " +
                  externalProgram["inpdir"] + externalProgram["optfl"] + externalProgram["inpExt"])
            shutil.copyfile(
                optfl_external, externalProgram["inpdir"] + externalProgram["optfl"] + externalProgram["inpExt"])

        run_Optimization(optfl)
        theParser = externalProgram["readerFunctCall"]["output"](
            externalProgram["outdir"] + optfl + externalProgram["outExt"])
        geom, hashLine, title, charge, mult = theParser.getInpData()
        conn_mat, Conn = genConnectivityMatrix(geom)

    # get the number of atoms
    nAtoms = 0
    for g in geom:
       if (geom[g][0] != 0): nAtoms += 1


    if (hashLine_nics == "DEFAULT\n"):
        hashLine_nics = externalProgram.defaultNicsKeyline
    if (hashLine_ncs == "DEFAULT\n"):
        hashLine_ncs = externalProgram.defaultNcsKeyline
    if (hashLine_nbo == "DEFAULT\n"):
        hashLine_nbo = externalProgram.defaultNboKeyline

    # Generate NICS input
    if (ncs_flag):
        genNicsInputs(geom, Conn, hashLine_ncs, title, charge, mult, jobType)
    else:
        genNicsInputs(geom, Conn, hashLine_nics, title, charge, mult, jobType)

    print("\nStatus : NICS input for all the centers generated.")


def Execute():

    # global flags
    global opt_flag, ncs_flag, sigma_flag, xy_flag, pointonly_flag, integralnics_flag, analyse_flag, area_flag, s_charge_flag, s_mult_flag, opt_external, optfl_external, inponly_flag, outonly_flag
    # global molecule-related
    global armpath, CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, xy_ref_ring_info, BQGuide, points, normals
    # global technical
    global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, xy_BQ_dist, sigma_charge, sigma_mult, analyse_dist, clear_flag, xy_extend
    # global geom length (number of atoms)
    global nAtoms

    if not inponly_flag:
        if (not outonly_flag): 
            run_Nics()
            print("\nAll the jobs are over")
            
        print("\nStatus : Filtering Appropriate Data .. ")

        grepData()

        if (ncs_flag):
            # pi-MOs are same in all the output files, so just identify them from the first file.
            piMOs, nocc = identifyPiMOs(
                externalProgram["outdir"] + flprfx + "-center1-set1" + externalProgram["outExt"])
            for inpFil in inputFileSet:
                outfl = externalProgram["outdir"] + \
                        inpFil["flprfx"] # + externalProgram["outExt"]  (TODO: probable fix, need to check other flows)

                nat = inpFil["nat"]

                grepPiCMO(nat, piMOs, nocc, BQ_No, BQ_Range,
                          BQ_Step, outfl, externalProgram["outExt"])

        if (xy_flag):
            outfl = open(outfilename, "a")
            outfl.write(
                "\n\nThe Centers for Rings/Bonds from the .arm file corresponds to following BQs\n")
            outfl.write("Ring  Distance    X       Y       Z\n")
            for i in range(1, len(BQGuide)+1):
                if i == 1:
                    dist = 0.0
                else:
                    dist += round(getDistance(BQGuide[i]
                                  [1], BQGuide[i-1][1]), 3)
                outfl.write("  " + repr(i) + "     " + fpformat.fix(dist, 1) + "   " + fpformat.fix(
                    BQGuide[i][1][0], 3) + "   " + fpformat.fix(BQGuide[i][1][1], 3) + "   " + fpformat.fix(BQGuide[i][1][2], 3) + "\n")
            outfl.close()


def callAnalyse(flprfx, CenterOf, all_aromatic_rings, analyse_dist, outfl):

    #   conn_mat, Conn = genConnectivityMatrix(geom)
    global analyse_flag, integralnics_flag

    for ring in CenterOf:
        m_fl = flprfx + "-center" + repr(ring) + ".armdat"
        s_fl = flprfx + "-sigma-center" + repr(ring) + ".armdat"
        outfl.write("\n\nFor Center " + repr(ring))
        if analyse_flag:
            analyse(m_fl, s_fl, analyse_dist, outfl)
        elif integralnics_flag:
            integralnics_analyse(m_fl, s_fl, BQ_Range[0], outfl)

    n_count = len(CenterOf)
    if (len(normals) > 0):
        for i in normals:
            n_count += 1
            m_fl = flprfx + "-center" + repr(n_count) + ".armdat"
            s_fl = flprfx + "-sigma-center" + repr(n_count) + ".armdat"
            outfl.write("\n\nFor Center " + repr(n_count))
            if analyse_flag:
                analyse(m_fl, s_fl, analyse_dist, outfl)
            elif integralnics_flag:
                integralnics_analyse(m_fl, s_fl, BQ_Range[0], outfl)

    if (area_flag):
        area = {}
        tot_area = 0.0
        outfl.write(
            "\nThe area in sq. ang. for each of the aromatic ring defined are: ")
        for ring in all_aromatic_rings:
            ring_atoms = all_aromatic_rings.get(ring)
            area[ring] = ring_area(
                externalProgram["inpdir"] + flprfx + "-center1-set1" + externalProgram["inpExt"], ring_atoms)
            outfl.write("\n" + repr(ring) + " : " + repr(area[ring]))
            tot_area += area[ring]

        outfl.write("\n\nThe total area is " +
                    repr(tot_area) + " sq. ang.\n\n")


def runJobs():
    # global flags
    global opt_flag, ncs_flag, sigma_flag, xy_flag, pointonly_flag, integralnics_flag, analyse_flag, area_flag, s_charge_flag, s_mult_flag, opt_external, optfl_external, inponly_flag, outonly_flag
    # global molecule-related
    global armpath, CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, xy_ref_ring_info, BQGuide, points, normals
    # global technical
    global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, xy_BQ_dist, sigma_charge, sigma_mult, analyse_dist, clear_flag, xy_extend
    # input file set
    global inputFileSet

    if (not opt_external):  # if asked not asked for optimization, read the geometry from the input geom file
        for extension in externalProgram["extensions"]:
            if (externalProgram["extensions"][extension].count(geomflext) == 1):
                exttype = extension

        #  Read the geometry and other data
        theParser = externalProgram["readerFunctCall"][exttype](geomfl)
        geom, hashLine, title, charge, mult = theParser.getInpData()
        conn_mat, Conn = genConnectivityMatrix(geom)
    else:
        geom = {}
        title = ""
        charge = ""
        mult = ""
        Conn = []  # if asked for optimization, Execute() below will optimize and read the geom

    if (not outonly_flag):

        # generate inputs for primary run
        generateAllInputs(geom, title, charge, mult, Conn, "main")

        if (sigma_flag):
            if (opt_flag or opt_external):
                theParser = externalProgram["readerFunctCall"]["output"](
                    externalProgram["outdir"] + flprfx + "-opt" + externalProgram["outExt"])
                geom, hashLine, title, charge, mult = theParser.getInpData()
                conn_mat, Conn = genConnectivityMatrix(geom)

            opt_flag = 0
            ncs_flag = 0
            exttype = "input"

            sigma_geom, charge, mult, zmat_idx = genSigmaModel(
                flprfx, geom, Conn, title, charge, sigma_mult)

            # Read the geometry and other data
            # theParser = externalProgram["readerFunctCall"][exttype](geomfl)
            # sigma_geom, hashLine, title, charge, mult = theParser.getInpData()

            # Take out the normals from the geometry, if defined
            sigma_normals = {}
            n_count = 1
            for i in range(len(sigma_geom)-(2*len(normals))+1, len(sigma_geom)+1, 2):
                sigma_normals[n_count] = sigma_geom[i][1:4] + sigma_geom[i+1][1:4]
                n_count += 1
                del sigma_geom[i]
                del sigma_geom[i+1]
            org_normals = normals
            normals = sigma_normals

            conn_mat, Conn = genConnectivityMatrix(sigma_geom)

            new_CenterOf = getNewRings(geom, sigma_geom, CenterOf, zmat_idx)
            org_flprfx = flprfx
            org_CenterOf = CenterOf
            flprfx = flprfx + "-sigma"
            CenterOf = new_CenterOf

            if (xy_flag):
                outfl = open(outfilename, "a")
                outfl.write(
                    "\n--------------------------------------------------------------------\n")
                outfl.write("\nFor the Sigma Model:\n")
                outfl.close()

            if (s_charge_flag):
                scharge = sigma_charge
            else:
                scharge = charge

            if (s_mult_flag):
                smult = sigma_mult
            else:
                smult = mult

            # generate inputs for sigma run
            generateAllInputs(sigma_geom, title, scharge, smult, Conn, "sigma")

        # all inputs are generated at this point, dump the json
        inpf = open(externalProgram["outdir"] + flprfx + "-inputFileSet.json", "w")  # this file name needs to change, based on input
        inpf.write(json.dumps(inputFileSet))
        inpf.close()
    else:
        outf = open(externalProgram["outdir"] + flprfx + "-inputFileSet.json", "r") # this file name needs to change, based on input
        inputFileSet = json.loads(outf.read())
        outf.close()

    # execute all jobs
    Execute()

    # Read For: For XY-Scan with Sigma-Model, just keep the final output as .armlog with r, ZZ and del-ZZ
    if (xy_flag and sigma_flag):
       xarmlogfile = externalProgram["outdir"] + flprfx + "-alldiff.armlog"
       armdatlines = readFile(externalProgram["outdir"] + flprfx + "-allcenter" + ".armdat")

    # clean up
    if (clear_flag):
        print("\nClearing up unnecessary files .. \n")

        removeFiles(externalProgram["inpdir"] + flprfx + "-center*")
        removeFiles(externalProgram["inpdir"] + flprfx + "-guessonly*")
        removeFiles(externalProgram["outdir"] + flprfx + "-guessonly*")

        if (opt_flag):
            removeFiles(externalProgram["inpdir"] + flprfx + "-opt* ")

    numpy_flag = checkNumPy()
    if not inponly_flag:
        if (integralnics_flag or analyse_flag):
            if (sigma_flag and numpy_flag and not xy_flag):
                outfl = open(outfilename, "a")
                callAnalyse(org_flprfx, org_CenterOf,
                            all_aromatic_rings, analyse_dist, outfl)
                outfl.close()

        # For XY-Scan with Sigma-Model, just keep the final output as .armlog with r, ZZ and del-ZZ
        if (xy_flag and sigma_flag):
            armlog = open(xarmlogfile, "w")
            armlog.write("r       ZZ       Sigma-ZZ        Del-ZZ\n")
            if (xy_flag, sigma_flag):
                sarmdatlines = readFile(
                    externalProgram["outdir"] + flprfx + "-allcenter" + ".armdat")
                lineno = min(len(armdatlines), len(sarmdatlines))
                for i in range(1, lineno):
                    awords = list(map(float, armdatlines[i].split()))
                    sawords = list(map(float, sarmdatlines[i].split()))
                    armlog.write(fpformat.fix(awords[0], 2) + "   " + fpformat.fix(awords[8], 4) + "   " + fpformat.fix(
                        sawords[8], 4) + "   " + fpformat.fix(awords[8]-sawords[8], 4) + "\n")
            armlog.close()


def writeOutputHeader():
    # global flags
    global opt_flag, ncs_flag, sigma_flag, xy_flag, pointonly_flag, integralnics_flag, analyse_flag, area_flag, s_charge_flag, s_mult_flag, opt_external, optfl_external, inponly_flag
    # global molecule-related
    global armpath, CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, xy_ref_ring_info, BQGuide, points, normals
    # global technical
    global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, xy_BQ_dist, sigma_charge, sigma_mult, analyse_dist, clear_flag, xy_extend

    outfl = open(outfilename, "w")
    outfl.write(
        "\n--------------------------------------------------------------------")
    outfl.write("\n                        ** Aroma Run **")
    outfl.write("\n                        ** " + aromaVersion() + " **")
    outfl.write(
        "\n--------------------------------------------------------------------\n")
    outfl.write("\n --------- Input Dump ---------- \n\n")
    armlines = readFile(armpath + flprfx + ".arm")
    for aline in armlines:
        outfl.write(aline)
    outfl.write("\n --------- Input End ---------- \n")
    outfl.close()

    if (xy_flag):
        outfl = open(outfilename, "a")
        outfl.write("\nFor the Original Molecule:\n")
        outfl.close()


def aroma(armfile):
    # global flags
    global opt_flag, ncs_flag, sigma_flag, xy_flag, pointonly_flag, integralnics_flag, analyse_flag, area_flag, s_charge_flag, s_mult_flag, opt_external, optfl_external, inponly_flag, outonly_flag
    # global molecule-related
    global armpath, CenterOf, geomflext, geomfl, flprfx, outfilename, sigma_direction, all_aromatic_rings, n_xy_center, xy_ref_ring_info, BQGuide, points, normals
    # global technical
    global runtype, hashLine_nics, hashLine_opt, hashLine_ncs, hashLine_nbo, BQ_Step, BQ_Range, BQ_No, xy_BQ_dist, sigma_charge, sigma_mult, analyse_dist, clear_flag, xy_extend

    print("\n--------------------------------------------------------------------")
    print("                        ** Aroma Run Begins. **")
    print("--------------------------------------------------------------------\n")

    # init global flags
    init()

    # set the file prefix
    flprfx = armfile[armfile.rindex("/")+1:len(armfile)]

    # read aroma input file and set necessary flags and params
    check(armfile)

    # write the output header, the file where all aroma output will be spit out
    print("The final output will be stored in " + outfilename)
    writeOutputHeader()

    # run all the jobs - parent molecule and sigma model
    runJobs()

    # all over
    print("\n--------------------------------------------------------------------")
    print("                        ** Aroma Run Over. **")
    print("--------------------------------------------------------------------\n")


if __name__ == "__main__":
    aroma(sys.argv[1])
