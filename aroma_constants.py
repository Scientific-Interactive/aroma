#! /usr/bin/env python
# Author : Anuja
# 30.01.2013
# Last Updated : 05.02.2014

# A Python Script for setting the constants and paths for Package Aroma 

print("system config import")

###########################################################################################################################################
# Atom information
AtmSym = {'A':-1,'X':0, 'H':1, 'HE':2, 'LI':3, 'BE':4, 'B':5, 'C':6, 'N':7, 'O':8, 'F':9, 'AL':13,'SI':14, 'P':15, 'S':16, 'CL':17, 'TI':22}
AtmMass = {-1:0.0, 0:0.0, 1:1.007, 2:4.002, 3:6.941, 4:9.012, 5:10.810, 6:12.010, 7:14.000, 8:15.990, 9:18.990, 13:26.982, 14:28.086, 15:30.974, 16:32.065, 17:35.453, 22:47.867}
# Covalent Radii are in angstrom
AtmCovalentRadii = {-1:0.7, 0: 0.7, 1:0.23, 2:0.32, 3:0.32, 4:0.90, 5:0.82, 6:0.77, 7:0.75, 8:0.70, 9:0.71, 13:1.18, 14:1.11, 15:1.06, 16:1.05, 17:0.99, 22:1.32}

# Dictionary for Maximum number of bonds an atom can have
Max_Conn = {5:4, 6:4, 7:4, 8:3, 13:4, 14:4, 15:5, 16:3, 17:1, 22:6}

# Dictionary of key:value :: atom number : 

# Geometry Related Constants
COORDINATE_EQUALITY_TOLERENCE = 0.001 # in angstrom
COVALENT_BOND_TOLERENCE = 0.4 # in angstrom
TORSION_ANGLE_TOLERANCE = 15 # in degrees, ideal value 5
###########################################################################################################################################




###########################################################################################################################################
# Constants, Paths for setting Gaussian Runs

inpdir = "tests/"
outdir = "tests/"
chkdir = "chk/"
FormChkCmd = "/usr/local/g09/formchk "
#GaussCmd = "/usr/local/g09/g09 "
GaussCmd = "C:/Users/user/Documents/2.0Aroma/G09W/G09W/g09  "
GaussInpExt = ".in"
GaussOutExt = ".log"

# Define construction of command to run Gaussian
def constructGaussCMD(flprfx):
   return GaussCmd + inpdir + flprfx + GaussInpExt + " " + outdir + flprfx + GaussOutExt


# List of Possible extensions for a Gaussian Files
# This Dictionary can be updated as per your convience.
# Note: The Extension can be anything, however, the Format of the File should be maintained as standard
EXTENSIONS_FOR_GAUSSIAN_FILES = {'input':['com','gjf','inp','in'], 'output':['log', 'out'], 'checkpoint':['chk']}

###########################################################################################################################################




###########################################################################################################################################
# Defaults for Aroma

DEFAULT_OPTIMIZATION_KEYLINE = "%nproc=1\n%mem=1024MB\n# HF/STO-3G OPT \n"
DEFAULT_NICS_KEYLINE = "%nproc=1\n%mem=1024MB\n# HF/STO-3G NMR=GIAO INTEGRAL=(GRID=ULTRAFINE) CPHF=(GRID=FINE)\n"
DEFAULT_NCS_KEYLINE = "%nproc=1\n%mem=1024MB\n# HF/STO-3G NMR=GIAO IOP(10/46=1) POP(NBOREAD, FULL) INTEGRAL=(GRID=ULTRAFINE) CPHF=(GRID=FINE)\n"
DEFAULT_NBO_KEYLINE = "$NBO NCS=0.1 <I MO XYZ> $END\n"

# Defaults for NICS
DEFAULT_BQ_STEP = 0.1 # in angstrom
DEFAULT_BQ_RANGE = [0, 4]
# Default for distance from molecular plane in case of XY-Scan
DEFAULT_XY_DISTANCE = 1.7
# For fitting polynomials, the BQs from distance defined by following parameter onwards are considered
DEFAULT_DISTANCE_FOR_ANALYSIS = 1.1

MAX_BQS_IN_INPFL = 50
###########################################################################################################################################


###########################################################################################################################################
# Defaults for Aroma Generated Sigma-Only Model
# All the angles are in degrees and lengths in angstrom
FIXED_SIGMA_ANGLE = '95.0' 
FIXED_SIGMA_DIHEDRAL_ANGLE = '0.0' 
# Dictionary of Atom-H bond length
# Key: Value :: Atomic Number : Bond Length
ATM_H_BL = {5:'1.19', 6:'1.00', 7:'1.00', 8:'0.96', 13:'1.55', 14:'1.47',15:'1.35', 16:'1.31', 22:'1.60'}

###########################################################################################################################################


