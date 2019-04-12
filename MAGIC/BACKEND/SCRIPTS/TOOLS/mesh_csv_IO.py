#!/usr/bin/env python

import csv
import sys
import numpy as np
# import argparse

def readmesh_csv(fname):

    with open(fname,'r') as csvfile:
        reader = csv.reader(csvfile,delimiter=',')
        ii = 0
        for r in reader:
            if ii==0: 
                NFIELDS = len(r)

    nsz = reader.line_num-1
    csvfile.close()

    print 'reading csv...',fname,'...nsz:',nsz

    ijk = np.zeros([nsz,NFIELDS],dtype='int64')

    with open(fname,'r') as csvfile:
        reader = csv.reader(csvfile,delimiter=',')
        ii = 0
        for r in reader:
            if ii>0: 
                #
                # note how content is organized : the first 3 elements are always ijk
                # optionally the 4th is the usued for multiblock to know 
                # the original cell index, for subtraction purposes
                #
                if NFIELDS == 3:
                    ijk[ii-1,:] = int(r[0]), int(r[1]), int(r[2])

                elif NFIELDS == 4:
                    ijk[ii-1,:] = int(r[1]), int(r[2]), int(r[3]), int(r[0]) 
            ii+=1

    return ijk

########################
if __name__ == '__main__':

    ijk = readmesh_csv(f)

    imin,imax = np.amin(ijk1[:,0]), np.amax(ijk1[:,0])
    jmin,jmax = np.amin(ijk1[:,1]), np.amax(ijk1[:,1])
    kmin,kmax = np.amin(ijk1[:,2]), np.amax(ijk1[:,2])

    dims = ijk.shape
    nsz += dims[0]

    NX = max(NX,imax)
    NY = max(NY,jmax)
    NZ = max(NZ,kmax)
