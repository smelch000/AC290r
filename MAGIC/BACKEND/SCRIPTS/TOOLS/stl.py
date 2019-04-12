#!/usr/bin/env python

import sys,time,copy
import numpy
from math import *
from TOOLS.mesh import *

class Stl():

    def __init__(self):
        
        self.points = dict()
        self.points_rev = dict()
        self.tris = set()

    def readneline(self,f):
        l = f.readline()
        while not l.strip():
            l = f.readline()
        return l

    def readstl(self):

        stlf = open(INSTL,'r')

        count = 0

        line = readneline(stlf) #stlf.readline()

        while True:
            line = readneline(stlf) #stlf.readline() # facet normal
            fields = tuple(line.split())
            if (fields[0] == "endsolid"): break
	
            line = readneline(stlf) #stlf.readline() # reads "outer loop ..."
	
            line = readneline(stlf) #stlf.readline() # vertex 1
            fields = tuple(line.split())
            x = float(fields[1])
            y = float(fields[2])
            z = float(fields[3])
            p1 = (x,y,z)
            if (not p1 in self.points) :
                self.points[p1] = count
                self.points_rev[count] = p1
                count += 1

            line = readneline(stlf) #stlf.readline() # vertex 2
            fields = tuple(line.split())
            x = float(fields[1])
            y = float(fields[2])
            z = float(fields[3])
            p2 = (x,y,z)
            if (not p2 in self.points) :
                self.points[p2] = count
                self.points_rev[count] = p2
                count += 1
	
            line = readneline(stlf) #stlf.readline() # vertex 3
            fields = tuple(line.split())
            x = float(fields[1])
            y = float(fields[2])
            z = float(fields[3])
            p3 = (x,y,z)
            if (not p3 in self.points) :
                self.points[p3] = count
                self.points_rev[count] = p3
                count += 1

            self.tris.add((self.points[p1], self.points[p2], self.points[p3]))
            line = readneline(stlf) #stlf.readline() # endloop
            line = readneline(stlf) #stlf.readline() # endfacet

        stlf.close()

        return

        """
        fa=open('area.dat','w')
        stl_points2 = dict()
        stl_points_rev2 = dict()
        stl_tris2 = set()
        count2 = 0
        for vtx in self.tris:

            i1,i2,i3 = vtx

            p1 = self.points_rev[i1]
            p2 = self.points_rev[i2]
            p3 = self.points_rev[i3]

            aa = [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]
            bb = [p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2]]
            cc = [p3[0]-p2[0], p3[1]-p2[1], p3[2]-p2[2]]

            a = sqrt(aa[0]**2 + aa[1]**2 + aa[2]**2)
            b = sqrt(bb[0]**2 + bb[1]**2 + bb[2]**2)
            c = sqrt(cc[0]**2 + cc[1]**2 + cc[2]**2)

            perh = 0.5*(a+b+c) # semiperimeter
            area = sqrt(perh * (perh-a) * (perh-b) * (perh-c)) # Heron formula

            fa.write('%f \n' % area)

            if area<CUT_TRIAREA: continue

            if (not p1 in stl_points2) :
                stl_points2[p1] = count2
                stl_points_rev2[count2] = p1
                count2 += 1

            if (not p2 in stl_points2) :
                stl_points2[p2] = count2
                stl_points_rev2[count2] = p2
                count2 += 1
	
            if (not p3 in stl_points2) :
                stl_points2[p3] = count2
                stl_points_rev2[count2] = p3
                count2 += 1

            stl_tris2.add((stl_points2[p1], stl_points2[p2], stl_points2[p3]))

        fa.close()

        return stl_points_rev2, stl_tris2
        """

