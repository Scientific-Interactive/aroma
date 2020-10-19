#! /usr/bin/env python
# Author : Anuja
# @ 10.02.2013
# Last Updated : 20.02.2014

# Utility Functions for Package Aroma 

import os
import sys
import math
import string

import aroma_constants
from aroma_constants import *

def readFile(filename):
   f = open(filename, "r")
   lines = f.readlines()
   f.close()

   return lines

def all_same(listname):
   return all(val == listname[0] for val in listname)

def all_same_within_threshold(listname, threshold):
   return all( (abs(val-listname[0]) < threshold) for val in listname)

def detectPlane(geom):

   X = []; Y = []; Z = []
   for i in geom:
      if (geom[i][0] != 1):
         X.append(geom[i][1])
         Y.append(geom[i][2])
         Z.append(geom[i][3])

   if all_same_within_threshold(X, aroma_constants.COORDINATE_EQUALITY_TOLERENCE): plane = 'YZ'
   elif all_same_within_threshold(Y, aroma_constants.COORDINATE_EQUALITY_TOLERENCE): plane = 'XZ'
   elif all_same_within_threshold(Z, aroma_constants.COORDINATE_EQUALITY_TOLERENCE): plane = 'XY'
   else : 
      plane = ''
      print("The Molecule is not Planar. Please check the input and submit again.")

   return plane

def getCM(points, weights):

   sumx = 0.0; sumy = 0.0; sumz = 0.0
   sum_weight = 0.0
   for i in range (0, len(points)):
      atm_weight = weights[i]
      sumx += points[i][0]*atm_weight
      sumy += points[i][1]*atm_weight
      sumz += points[i][2]*atm_weight
      sum_weight += atm_weight

   cmx = sumx/sum_weight
   cmy = sumy/sum_weight
   cmz = sumz/sum_weight

   return cmx, cmy, cmz

# This is function for determining Geometrical Center
def getGM(points):

   sumx = 0.0; sumy = 0.0; sumz = 0.0
   count = len(points)
   for i in range (0, count):
      sumx += points[i][0]
      sumy += points[i][1]
      sumz += points[i][2]

   cmx = sumx/count
   cmy = sumy/count
   cmz = sumz/count

   return cmx, cmy, cmz

def meanOf(data):
   sum = 0.0
   l = len(data)
   for i in range (0, l):
     sum += data[i]

   return sum/l

# Everything about Geometry
def getDistance(a,b):
   return math.sqrt( (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2 )

def getAngleBetweenPoints(a,b,c):
   vec_ba = getVector(b,a)
   vec_bc = getVector(b,c) 
   dotprod =  dot(vec_ba, vec_bc)
   if (abs(dotprod) < 1e-5 ): 
      return 90.0
   else:
      return math.degrees(math.acos(round(dotprod/(vectorMagnitude(vec_ba)*vectorMagnitude(vec_bc)),10)))

def getAngleBetweenVec(vec_ab,vec_bc):
   norm_ab = vectorMagnitude(vec_ab)
   norm_bc = vectorMagnitude(vec_bc)
   dotprod =  dot(vec_ab, vec_bc)
   if (abs(dotprod) < 1e-5 ): return 90.0
   else:
      return math.degrees(math.acos(round(dotprod/(norm_ab*norm_bc),10)))

def getDihedralAngleBetween(a,b,c,d):
   vec_ba = getVector(b,a)
   vec_bc = getVector(b,c)
   vec_cb = getVector(c,b)
   vec_cd = getVector(c,d)
   n1 = cross(vec_ba, vec_bc)
   n2 = cross(vec_cb, vec_cd)

   sign = dot(vec_ba, cross(vec_bc, vec_cd))

   if (sign < 0):
      return (-1)*getAngleBetweenVec(n1,n2) 
   else:
      return getAngleBetweenVec(n1,n2) 

# A function to calculate a normal to the plane defined by 3 points, a, b, c
def getNormalToThePlane(a,b,c):
   ba = [a[0]-b[0], a[1]-b[1], a[0]-b[1]]
   bc = [c[0]-b[0], c[1]-b[1], c[0]-b[1]]
   normal_to_plane = cross(ba, bc)
   return normal_to_plane

def getVector(a,b):
   vec_ab = []
   for i in range (0,3): vec_ab.append(b[i]-a[i])
   return vec_ab

def vectorMagnitude(ab):
   return math.sqrt(ab[0]**2 + ab[1]**2 + ab[2]**2)

def getUnitVector(vec):
   norm_vec = vectorMagnitude(vec)
   if (norm_vec == 0): return [0.0, 0.0, 0.0]
   else: return [vec[0]/norm_vec, vec[1]/norm_vec, vec[2]/norm_vec]

def dot(ab,bc):
   return ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2]

def cross(ab,bc):
   abcrossbc = []
   abcrossbc.append(ab[1]*bc[2] - ab[2]*bc[1])
   abcrossbc.append(ab[2]*bc[0] - ab[0]*bc[2])
   abcrossbc.append(ab[0]*bc[1] - ab[1]*bc[0])
   return abcrossbc

def checkNumPy():
   try:
       import numpy
       return True
   except:
       print("\nWarning: Analysis of data can not be performed as Python library \"numpy\" is not available.")
       return False

# A function for getting area of polygon (here, ring)
# Vertices is a list of (x, y, z) for coordinates
def areaOfPolygon(vertices, ignore='z'):
   area = 0.0
   n = len(vertices)
   if (n < 3): return 0
   else:
      if (ignore.upper() == 'Z'): n1 = 0; n2 = 1
      elif (ignore.upper() == 'Y'): n1 = 0; n2 = 2
      elif (ignore.upper() == 'X'): n1 = 1; n2 = 2
      for i in range (0, n):
         if (i < n-1):
            area += 0.5*((vertices[i][n1]*vertices[i+1][n2]) - (vertices[i+1][n1]*vertices[i][n2]))
         elif (i == n-1): 
            area += 0.5*((vertices[i][n1]*vertices[0][n2]) - (vertices[0][n1]*vertices[i][n2]))

   return area

