#!/usr/bin/env python
#
# This script couples the geometry in an STL file with the shear 
# stress data produced by the MOEBIUS simulator (for that geometry).
# The output is gives as a VTK file.
#
import sys,os
import time
from myvtk import *
import numpy as np
import argparse
import ctypes
from TOOLS.mesh import *

###################
def main(args):

    MESHFNAME = ""
    OUTFNAME = ""
    WALLFNAME = ""

    MESHFNAME = args.infile
    OUTFNAME = args.outfile

    if args.outwall != None:
        WALLFNAME = args.outwall

    # open the .dat file
    try:
        meshf = open(MESHFNAME+'.dat','r')
    except IOError:
        print >> sys.stderr, "Cannot find file", MESHFNAME+'.dat'
        sys.exit(1)

    if (OUTFNAME == ""):
        outf = sys.stdout
    else:
        try:
            outf = open(OUTFNAME, 'w')
        except IOError:
            print >> sys.stderr, "Cannot write file", OUTFNAME
            sys.exit(1)

    if WALLFNAME:
        try:
            wallf = open(WALLFNAME, 'w')
        except IOError:
            print >> sys.stderr, "Cannot write wall file", WALLFNAME
            sys.exit(1)

    msh = Mesh()

    msh.hashcols()

    meshdata = list()
    meshdataWIO = list()

    print >> sys.stderr, 'Reading mesh file...',
    start_time = time.time()

    for line in meshf:

        fields = tuple(line.split())

        meshdata.append( (int(fields[0]), int(fields[1]), int(fields[2]), int(fields[3])) )

        if int(fields[3]) == 1: continue # exclude fluid nodes
        # if int(fields[3]) == 3: continue # exclude inlet nodes
        # if int(fields[3]) == 4: continue # exclude outlet nodes

        meshdataWIO.append( (int(fields[0]), int(fields[1]), int(fields[2]), int(fields[3])) )

    meshf.close()
    print >> sys.stderr, "done", time.time()-start_time, "secs"

    write_time = time.time()

    ijkf = np.zeros([len(meshdata),4],dtype='int')
    n = 0
    for i,j,k,f in meshdata:
        ijkf[n,:] = i,j,k,f
        n += 1

    print 'write complete vtk (unstructured points + cells file: fluid/wall/inlet/outlet)...'
    # write the complete vtk (fluid, wall, inlet, outlet)
    msh.writeVTK(OUTFNAME,ijkf)

    ijkf = np.zeros([len(meshdataWIO),4],dtype='int')
    n = 0
    for i,j,k,f in meshdataWIO:
        ijkf[n,:] = i,j,k,f
        n += 1

    if WALLFNAME:

        print 'write short vtk (legacy vtk: wall/inlet/outlet)...'

        ijkf = np.zeros([len(meshdataWIO),4],dtype='int')

        n = 0
        for i,j,k,f in meshdataWIO:
            ijkf[n,:] = i,j,k,f
            n+=1

        # write the short vtk (wall, inlet, outlet)
        msh.writeVTK(WALLFNAME,ijkf)

    print >> sys.stderr, "done", time.time()-write_time, "secs"
    print >> sys.stderr, "Total time:", time.time()-start_time, "secs"

    sys.exit(0)

###############################
def iniparser(parser):

    parser.add_argument('-i', '--infile', required=False,  help='mesh file (.dat)')
    parser.add_argument('-o', '--outfile',  required=True,  help='output complete vtk file')
    parser.add_argument('-w', '--outwall', required=False, default=None, help='output wall vtk file')

    return parser

###################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = iniparser(parser)
    args = parser.parse_args()

    main(args)
