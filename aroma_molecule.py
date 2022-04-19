#! /usr/bin/env python
# Author : Anuja
# Started: 20.02.2013
# Updated: 19.10.2014

# Processing Molecular Geometry and Identifying the Aromatic Rings for Sigma Model
# Generating Z-Matrix
# Rotating Molecule

import os
import sys
import math
import string
import itertools

from numpy import sin, cos, pi
from numpy.linalg import norm
import numpy as np

import aroma_constants
from aroma_constants import *

import aroma_util
from aroma_util import *


# ********************************************** Connections ****************************************************** #

def genConnectivityMatrix(geom):

   global Conn

   nat = len(geom)

   # Initialization of matrix
   conn_matrix = []
   Conn = {}
   for i in range (1, nat+1):
     conn_matrix.append([])
     Conn[i] = []
     for j in range (1, nat+1):
        if (i == j): conn_matrix[i-1].append(2)
        else: conn_matrix[i-1].append(0)

   for i in range (0, nat):
     for j in range (i+1, nat):
        if ((geom[i+1][0] != 0) and (geom[j+1][0] != 0)):
           dist = math.sqrt( (geom[i+1][1]-geom[j+1][1])**2 + (geom[i+1][2]-geom[j+1][2])**2 + (geom[i+1][3]-geom[j+1][3])**2 )
           if ( (dist <= AtmCovalentRadii[geom[i+1][0]] + AtmCovalentRadii[geom[j+1][0]] + COVALENT_BOND_TOLERENCE) and \
                (dist >=  AtmCovalentRadii[geom[i+1][0]] + AtmCovalentRadii[geom[j+1][0]] - COVALENT_BOND_TOLERENCE) ): 
              conn_matrix[i][j] = conn_matrix[j][i] = 1
              Conn[i+1].append(j+1); Conn[j+1].append(i+1)

   j = 1
   for i in range (0, nat):
      if (geom[i+1][0] == 0):
         conn_matrix[i][j-1] = conn_matrix[j-1][i] = 1
         Conn[j].append(i+1); Conn[i+1].append(j)

   return conn_matrix, Conn

def checkConnections(i):

   for j in range (0, len(Conn[i])):
      if ( (parentOf[Conn[i][j]] == "") and (parentOf[i] != Conn[i][j])): parentOf[Conn[i][j]] = i

      if (Conn[i][j] != parentOf[i]):
        if (Conn[i][j] == atom_to_be_found):

           child = Conn[i][j]
           flag = 1
           while (flag):
              probable_ring.append(parentOf[child])
              child = parentOf[child]
              if (child == atom_to_be_found): 
                flag = 0
           Rings.append(probable_ring)

        else:
           if (parentOf[Conn[i][j]] == i): checkConnections(Conn[i][j])

# ******************************************** Connections Ends **************************************************** #

# ********************************************** RINGS ****************************************************** #

def getCMOfRing(geom, ring_atoms):
      # Get the Center of Mass of the Ring
      points = []
      weights = []
      for i in range (0, len(ring_atoms)):
         atm_idx = ring_atoms[i]

         points.append(geom[atm_idx][1:4])
         weights.append(AtmMass[geom[atm_idx][0]])

      cmx, cmy, cmz = getCM(points, weights)
      return cmx, cmy, cmz

def getGMOfRing(geom, ring_atoms):
      # Get the Center of Mass of the Ring
      points = []
      for i in range (0, len(ring_atoms)):
         atm_idx = ring_atoms[i]

         points.append(geom[atm_idx][1:4])

      cmx, cmy, cmz = getGM(points)
      return cmx, cmy, cmz

def checkConnAndAppend(Conn, ring_atoms, ring_atoms_ordered, idx):
   for i in range (0, len(Conn[idx])):
      atmx = Conn[idx][i]
      if ( (ring_atoms.count(atmx) == 1) and (ring_atoms_ordered.count(atmx) < 1) ):
         ring_atoms_ordered.append(atmx)
         break;
   return atmx

def getOrderedRing(Conn, ring_atoms):
   # Since the rings may be specified in any random order, first they need to be ordered based on the connectivity
   nrings = len(ring_atoms)
   if (nrings > 1):
      ring_atoms_ordered = []
      idx = ring_atoms[0]
      ring_atoms_ordered.append(idx)
      while (len(ring_atoms_ordered) < nrings):
         new_idx = checkConnAndAppend(Conn, ring_atoms, ring_atoms_ordered, idx)
         idx = new_idx

      return ring_atoms_ordered
   else: return ring_atoms

def checkPlanarityOfRing(geom, ring_atoms):
   flag_planarity = 1
   nrings = len(ring_atoms)
   if (nrings < 3): return flag_planarity
   for i in range (0, nrings):
      a = ring_atoms[i]
      if (i+1 < nrings-1):
         b = ring_atoms[i+1]
      else: 
         b = ring_atoms[i+1-nrings]
      if (i+2 < nrings-1):
         c = ring_atoms[i+2]
      else: 
         c = ring_atoms[i+2-nrings]
      if (i+3 < nrings-1):
         d = ring_atoms[i+3]
      else: 
         d = ring_atoms[i+3-nrings]

      dihed = getDihedralAngleBetween(geom[a][1:4],geom[b][1:4],geom[c][1:4],geom[d][1:4])

      if (-TORSION_ANGLE_TOLERANCE < abs(dihed) < TORSION_ANGLE_TOLERANCE): 
         continue
      elif (180-TORSION_ANGLE_TOLERANCE < abs(dihed) < 180+TORSION_ANGLE_TOLERANCE): 
         continue
      else:
         flag_planarity = 0
         break;

   return flag_planarity

def angleWithZ(geom, ring_atoms_ordered, CM):

   cmx = CM[0]; cmy = CM[1]; cmz = CM[2]
   nrings = len(ring_atoms_ordered)
   # theta is a list of angles between Normal to the ring and the Z-axis
   # Normal to the Ring is determined by considering two position vectos of a pair of bonded atoms, denoted here by a and b
   # In case, of strained rings, an average angle is to be taken
   theta = []
   for i in range (0, nrings):
      a = ring_atoms_ordered[i]
      if (i+1 < nrings):
         b = ring_atoms_ordered[i+1]
      else:
         b = ring_atoms_ordered[i+1-nrings]

      # Find a normal to the plane of the ring which passes through Origin
      # Therefore, the normal is found for the Plane formed by the CM (i.e. Origin, in case of translated geometry) and atoms a, b of the Ring
      vec1 = [geom[a][1] - cmx, geom[a][2] - cmy, geom[a][3] - cmz]
      vec2 = [geom[b][1] - cmx, geom[b][2] - cmy, geom[b][3] - cmz]
      normal_to_Ring = cross(vec1, vec2)
      angle = getAngleBetweenPoints(normal_to_Ring, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0])

      theta.append(angle)

   mean_theta = math.radians(meanOf(theta))
   return mean_theta

def getAverageNormaltoTheRing(geom, ring_atoms_ordered, CM):
   cmx = CM[0]; cmy = CM[1]; cmz = CM[2]
   nrings = len(ring_atoms_ordered)

   normals = []
   for i in range (0, nrings):
      a = ring_atoms_ordered[i]
      if (i+1 < nrings):
         b = ring_atoms_ordered[i+1]
      else:
         b = ring_atoms_ordered[i+1-nrings]

      # Find a normal to the plane of the ring which passes through Origin
      # Therefore, the normal is found for the Plane formed by the CM (i.e. Origin, in case of translated geometry) and atoms a, b of the Ring
      vec1 = [geom[a][1] - cmx, geom[a][2] - cmy, geom[a][3] - cmz]
      vec2 = [geom[b][1] - cmx, geom[b][2] - cmy, geom[b][3] - cmz]
      normals.append(cross(vec1, vec2))

   xcomp = []; ycomp = []; zcomp = []
   for j in range (0, len(normals)):
      xcomp.append(normals[j][0])
      ycomp.append(normals[j][1])
      zcomp.append(normals[j][2])

   return [meanOf(xcomp), meanOf(ycomp), meanOf(zcomp)]

# Reorient Molecule so that the Ring under consideration is in XY plane
def reorient(geom, Conn, ring_atoms, points=[], normal=[]):

   nrings = len(ring_atoms)
   ring_atoms_ordered = getOrderedRing(Conn, ring_atoms)
   flag_planarity = checkPlanarityOfRing(geom, ring_atoms_ordered)

   if ((flag_planarity) and (nrings != 0)):

      # First translate the CM to origin
      cmx, cmy, cmz = getGMOfRing(geom, ring_atoms)
      translated_geom = {}
      dx = cmx; dy = cmy; dz = cmz
      for i in range(0, len(geom)):
         translated_geom[i+1] = []
         translated_geom[i+1].append(geom[i+1][0])
         translated_geom[i+1].append(geom[i+1][1] - dx)
         translated_geom[i+1].append(geom[i+1][2] - dy)
         translated_geom[i+1].append(geom[i+1][3] - dz)

      if (len(points) > 0):
         translated_points = {}
         for i in points:
           translated_points[i] = []
           translated_points[i].append(points[i][0])
           translated_points[i].append(points[i][1] - dx)
           translated_points[i].append(points[i][2] - dy)
           translated_points[i].append(points[i][3] - dz)

      mean_theta = angleWithZ(translated_geom, ring_atoms_ordered, [0.0, 0.0, 0.0])

      # Find a normal to the plane of the ring which passes through Origin
      # Therefore, the normal is found for the Plane formed by the Origin (i.e. CM) and the first two points of the Ring
      normal_to_Ring = getAverageNormaltoTheRing(translated_geom, ring_atoms_ordered, [0.0, 0.0, 0.0])

      # Average angle between the Normal to the Ring and the Z-axis, which is going to be the angle of rotation
      theta = mean_theta
      if (theta != 0.0):
         # Determine the Line of Nodes as referred in the definition of Euler Angles
         # This line is perpendicular to the Z-axis and the Normal
         # This line is going to be the axis of rotation
         line_of_nodes = cross(normal_to_Ring, [0.0, 0.0, 1.0])
         [ux, uy, uz] = getUnitVector(line_of_nodes)

         # Now rotate the whole molecule by angle theta around the line of nodes
         # This involves a complicated (!!) Rotational Matrix (R)
         costheta = math.cos(theta)
         oneminuscost = 1 - costheta
         sintheta = math.sin(theta)
         R11 = costheta + (ux**2)*oneminuscost
         R12 = ux*uy*oneminuscost - uz*sintheta
         R13 = ux*uz*oneminuscost + uy*sintheta
         R21 = ux*uy*oneminuscost + uz*sintheta
         R22 = costheta + (uy**2)*oneminuscost
         R23 = uy*uz*oneminuscost - ux*sintheta
         R31 = ux*uz*oneminuscost - uy*sintheta
         R32 = uy*uz*oneminuscost + ux*sintheta
         R33 = costheta + (uz**2)*oneminuscost

         new_geom = {}
         for i in range(0, len(geom)):
            new_geom[i+1] = []
            x1 = translated_geom[i+1][1]
            y1 = translated_geom[i+1][2]
            z1 = translated_geom[i+1][3]
            new_geom[i+1].append(geom[i+1][0])

            new_geom[i+1].append(R11*x1 + R12*y1 + R13*z1)
            new_geom[i+1].append(R21*x1 + R22*y1 + R23*z1)
            new_geom[i+1].append(R31*x1 + R32*y1 + R33*z1)

         new_points = {}
         if (len(points) > 0):
            for i in points:
               new_points[i] = []
               x1 = translated_points[i][1]
               y1 = translated_points[i][2]
               z1 = translated_points[i][3]
               new_points[i].append(points[i][0])

               new_points[i].append(R11*x1 + R12*y1 + R13*z1)
               new_points[i].append(R21*x1 + R22*y1 + R23*z1)
               new_points[i].append(R31*x1 + R32*y1 + R33*z1)

      else:
         new_geom = translated_geom
         new_points = {}
         if (len(points) > 0): new_points = translated_points

      return new_geom, new_points, normal

   elif ((nrings == 0) and (len(normal) > 0)):

      # First translate the "point" to origin
      cmx, cmy, cmz = normal[0:3]
      translated_geom = {}
      dx = cmx; dy = cmy; dz = cmz
      for i in range(0, len(geom)):
         translated_geom[i+1] = []
         translated_geom[i+1].append(geom[i+1][0])
         translated_geom[i+1].append(geom[i+1][1] - dx)
         translated_geom[i+1].append(geom[i+1][2] - dy)
         translated_geom[i+1].append(geom[i+1][3] - dz)

      translated_normal = [normal[0]-dx, normal[1]-dy, normal[2]-dz, normal[3]-dx, normal[4]-dy, normal[5]-dz]

      theta = math.radians(getAngleBetweenVec(translated_normal[3:6], [0.0, 0.0, 1.0]))
#      if ((theta != 0.0) or (theta != 180.0)):
      if (theta != 0.0):
#      if (0.0 < theta < 180.0 ):
         # Determine the Line of Nodes as referred in the definition of Euler Angles
         # This line is perpendicular to the Z-axis and the Normal
         # This line is going to be the axis of rotation
         line_of_nodes = cross(translated_normal[3:6], [0.0, 0.0, 1.0])
         [ux, uy, uz] = getUnitVector(line_of_nodes)

         # Now rotate the whole molecule by angle theta around the line of nodes
         # This involves a complicated (!!) Rotational Matrix (R)
         costheta = math.cos(theta)
         oneminuscost = 1 - costheta
         sintheta = math.sin(theta)
         R11 = costheta + (ux**2)*oneminuscost
         R12 = ux*uy*oneminuscost - uz*sintheta
         R13 = ux*uz*oneminuscost + uy*sintheta
         R21 = ux*uy*oneminuscost + uz*sintheta
         R22 = costheta + (uy**2)*oneminuscost
         R23 = uy*uz*oneminuscost - ux*sintheta
         R31 = ux*uz*oneminuscost - uy*sintheta
         R32 = uy*uz*oneminuscost + ux*sintheta
         R33 = costheta + (uz**2)*oneminuscost

         new_geom = {}
         new_normal = []
         for i in range(0, len(geom)):
            new_geom[i+1] = []
            x1 = translated_geom[i+1][1]
            y1 = translated_geom[i+1][2]
            z1 = translated_geom[i+1][3]
            new_geom[i+1].append(geom[i+1][0])

            new_geom[i+1].append(R11*x1 + R12*y1 + R13*z1)
            new_geom[i+1].append(R21*x1 + R22*y1 + R23*z1)
            new_geom[i+1].append(R31*x1 + R32*y1 + R33*z1)

         x1 = translated_normal[0]; y1 = translated_normal[1]; z1 = translated_normal[2]
         new_normal.append(R11*x1 + R12*y1 + R13*z1); new_normal.append(R21*x1 + R22*y1 + R23*z1); new_normal.append(R31*x1 + R32*y1 + R33*z1)
         x1 = translated_normal[3]; y1 = translated_normal[4]; z1 = translated_normal[5]
         new_normal.append(R11*x1 + R12*y1 + R13*z1); new_normal.append(R21*x1 + R22*y1 + R23*z1); new_normal.append(R31*x1 + R32*y1 + R33*z1)

      else:
         new_geom = translated_geom
         if (len(normal) > 0): new_normal = translated_normal

      return new_geom, points, new_normal

   else:
      return [], [], []


def identifyAromaticRings():

   global atom_to_be_found, probable_ring, Rings, parentOf

   Rings = []
   nat = len(Conn)

   parentOf = {}
   for i in range (1, nat+1): parentOf[i] = ""

   for i in range (1, nat+1):
      for j in range (1, nat+1): parentOf[j] = ""
      atom_to_be_found = i
      probable_ring = []
      checkConnections(i)


#   for r in range (0, len(Rings)): 
#      R = list(set(Rings[r]))
#      R.sort()
#      Rings[r] = R


   for r in range (0, len(Rings)): 
      Rings[r].sort()
      R =  list(Rings[r] for Rings[r],_ in itertools.groupby(Rings[r]))
      Rings[r] = R


   Rings.sort() 
   final_rings =  list(Rings for Rings,_ in itertools.groupby(Rings))

   super_rings = []
   for r1 in range (0, len(final_rings)):
      for r2 in range (r1+1, len(final_rings)): 
          superset = list(set(final_rings[r1] + final_rings[r2]))
          if ((final_rings.count(superset) > 0) and (super_rings.count(superset) < 1) ): super_rings.append(superset)

   for r in super_rings:
      final_rings.remove(r)

   print("The Final Rings:\n" + repr(final_rings))
      

# ******************************************** RINGS END **************************************************** #


def generateZMatrix(geom, Conn):
   nat = len(Conn)

   zmat = {}
   zmat_str = ""

   parent = {}
   leftouts = [] # List of dangling atoms to be added in the end
   for i in range (1, nat+1): parent[i] = ""

   # zmat_idx is a map of orginal indices to those in z-matrix
   # {key:value :: original:z-matrix index}
   zmat_idx = {}

   a = 1
   idx = 1
   zmat[a] = [a]
   zmat_idx[a] = idx
   zmat_str += repr(geom[a][0]) + "\n"
   flag = 1

   while (flag):
      flag_connections_over = 1
      for i in range (0, len(Conn[a])):
         if (not (Conn[a][i] in zmat)): 
            if (len(Conn[Conn[a][i]]) < 2): 
                leftouts.append(Conn[a][i])
                continue
            flag_connections_over = 0
            b = Conn[a][i]
            idx += 1
            zmat_idx[b] = idx
            parent[b] = a
            zmat[b] = []
            zmat[b].append(b)

            if (len(zmat[a]) == 1): 
               zmat[b].append(zmat[a][0])
            elif (len(zmat[a]) == 2): 
               zmat[b].append(zmat[a][0])
               zmat[b].append(zmat[a][1])
            else:
               zmat[b].append(zmat[a][0])
               zmat[b].append(zmat[a][1])
               zmat[b].append(zmat[a][2])

            angle = ""
            dihedral = ""
            c = ""
            d = ""
            if (len(zmat) > 3): 
               if (len(zmat[b]) == 4):
                 c = zmat[b][2]
                 d = zmat[b][3]
               if (len(zmat[b]) == 3):
                 c = zmat[b][2]
                 for j in range (0, len(Conn[a])):
                    e = Conn[a][j]
                    if ( (e != b) and (e != c) and (e in zmat_idx)): 
                       d = e
                       break;
                 zmat[b].append(d)
               if (len(zmat[b]) == 2):
                 for j in range (0, len(Conn[a])):
                    e = Conn[a][j]
                    if ((e != b) and (e in zmat_idx)): 
                       c = e 
                       break;
                 for j in range (0, len(Conn[c])):
                    e = Conn[c][j]
                    if ( (e != a) and (e != b) and (e in zmat_idx)): 
                       d = e 
                       break;
                 if (d == ''):
                    for k in range (1, idx-1):
                       e = idx - k
                       if ( (e != b) and (e != a) and (e != c) ): d = e
                     
                 zmat[b].append(c)
                 zmat[b].append(d)

               distance = getDistance(geom[b][1:4], geom[a][1:4]) 
               angle = getAngleBetweenPoints(geom[b][1:4], geom[a][1:4], geom[c][1:4])
               dihedral = -getDihedralAngleBetween(geom[b][1:4], geom[a][1:4], geom[c][1:4], geom[d][1:4])

               zmat_str += repr(geom[b][0]) + "  " + repr(zmat_idx[a]) + "   " + repr(distance) + "   " + repr(zmat_idx[c]) + "   " + repr(angle) + "   " + repr(zmat_idx[d]) + "   " + repr(dihedral) + "\n"
               break;

            elif (len(zmat) == 3):
               if (len(zmat[b]) < 3):
                  for ii in range (0, len(Conn[a])): 
                     if (Conn[a][ii] != b): break;
                  zmat[b].append(Conn[a][ii])
               c = zmat[b][2]
               distance = getDistance(geom[b][1:4], geom[a][1:4]) 
               angle = getAngleBetweenPoints(geom[b][1:4], geom[a][1:4], geom[c][1:4])

               zmat_str += repr(geom[b][0]) + "  " + repr(zmat_idx[a]) + "   " + repr(distance) + "   " + repr(zmat_idx[c]) + "   " + repr(angle) + "\n" 
               break;

            else: 
               distance = getDistance(geom[b][1:4], geom[a][1:4]) 
               zmat_str += repr(geom[b][0]) + "  " + repr(zmat_idx[a]) + "   " + repr(distance) + "\n" 
               break;

#      a = ""
      if (not flag_connections_over): a = b 
      elif (len(leftouts) > 0): a = ""
      else : a = parent[a]

      if (a == ""):
         for k in range (1, nat+1):
            if (k in zmat_idx):
               continue
            else:
               a = k
               idx += 1
               zmat_idx[a] = idx
               zmat[a] = [b, c, d]
               distance = getDistance(geom[b][1:4], geom[a][1:4]) 
               angle = getAngleBetweenPoints(geom[a][1:4], geom[b][1:4], geom[c][1:4])
               dihedral = -getDihedralAngleBetween(geom[a][1:4], geom[b][1:4], geom[c][1:4], geom[d][1:4])
               zmat_str += repr(geom[a][0]) + "  " + repr(zmat_idx[b]) + "   " + repr(distance) + "   " + repr(zmat_idx[c]) + "   " + repr(angle) + "   " + repr(zmat_idx[d]) + "   " + repr(dihedral) + "\n"
               break;

      if (len(zmat) == nat): flag = 0

   return zmat, zmat_str, zmat_idx

#
# parts of the code inspired from: https://github.com/ben-albrecht/qcl/blob/master/qcl/ccdata_xyz.py
def generateCartesianFromZmat(zmat):

   coord = []

   def toRadian(deg):
     return deg * (pi / 180.0)

   def placeLength(x): 
     atmtyp, connIdx, dist = x[:3]
     connIdx = int(connIdx)-1
     connCord = coord[connIdx]
     return [x[0], connCord[1]+float(dist), connCord[2], connCord[3]] 

   def placeAngle(x):
     atmtyp, connIdx, dist = x[:3]
     dist = float(dist)
     connIdx = int(connIdx)-1

     angleAtm, angle = x[3:5] 
     angleAtm = int(angleAtm)-1
     ang = toRadian(float(angle))

     if (len(x) > 5):
       torAtm, tor = x[5:7]
       torAtm = int(torAtm)-1
       tor = toRadian(float(tor))
     else:
       torAtm = -1
       tor = toRadian(90)
     
     avec = coord[connIdx]
     bvec = coord[angleAtm]
     
     if (torAtm == -1):
        cvec = ["X", 0.0, 1.0, 0.0]
     else:
        cvec = coord[torAtm]

     v1 = [avec[1]-bvec[1], avec[2]-bvec[2], avec[3]-bvec[3]]
     v2 = [avec[1]-cvec[1], avec[2]-cvec[2], avec[3]-cvec[3]]

     n  = np.cross(v1, v2)
     nn = np.cross(v1, n)

     n  /= norm(n)
     nn /= norm(nn)

     n  *= -sin(tor)
     nn *= cos(tor)

     v3 = n + nn
     v3 /= norm(v3)
     v3 *= dist * sin(ang) 

     v1 /= norm(v1)
     v1 *= dist * cos(ang)

     v2 = [avec[1], avec[2], avec[3]] + v3 - v1

     return [x[0], v2[0], v2[1], v2[2]]

   for zItm in zmat:
     if (len(zItm) == 0): continue

     # print(zItm)

     if (len(zItm) == 1): coord.append([zItm[0], 0, 0, 0])

     if (len(zItm) == 3):
        coord.append(placeLength(zItm))
     
     if (len(zItm) > 3):
        coord.append(placeAngle(zItm))

   # convert the coords in to 1 indexed dict
   geom = {}
   idx = 1
   for crd in coord:
     # print(' '.join(map(lambda x: repr(x) if type(x) != str else x, crd)))
     geom[idx] = crd
     idx += 1

   return geom


