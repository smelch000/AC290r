#!/usr/bin/env python
#
# This script couples the geometry in an STL file with the shear 
# stress data produced by the MOEBIUS simulator (for that geometry).
# The output is gives as a VTK file.
#
import sys,os
import time
import argparse
from myvtk import *
import numpy as np
import ctypes
from TOOLS.mesh import *

###################
def bgkflag2vtu(msh, MESHFNAME, OUTFNAME, WALLFNAME=None, WDIR='./'):
    
    try:
        meshf = open( os.path.join(WDIR, MESHFNAME+'.dat'), 'r')
    except IOError:
        print >> sys.stderr, "Cannot find mesh dat file", os.path.join(WDIR, MESHFNAME+'.dat')
        sys.exit(1)

    try:
        meshh = open( os.path.join(WDIR, MESHFNAME+'.hdr'), 'r')
    except IOError:
        print >> sys.stderr, "Cannot find mesh hdr file", os.path.join(WDIR, MESHFNAME+'.hdr')
        sys.exit(1)
    line = meshh.readline()
    line = meshh.readline()
    try:
        spacing = int( meshh.readline() )
    except:
        spacing = 1
    print 'Spacing:',spacing

    meshdata = list()
    meshdataWIO = list()
    # print >> sys.stderr, 'Reading mesh file...',

    msh.hashcols()

    # read the .dat file line-by-line
    for line in meshf:

        fields = tuple(line.split())

        msh.nx = max(msh.nx,int(fields[0]))
        msh.ny = max(msh.ny,int(fields[1]))
        msh.nz = max(msh.nz,int(fields[2]))

        # F,W,I,O nodes
        meshdata.append( (int(fields[0]), int(fields[1]), int(fields[2]), int(fields[3])) )

        if int(fields[3]) == msh.ID_FLUIDNODE: continue # exclude fluid nodes

        # W,I,O nodes only
        meshdataWIO.append( (int(fields[0]), int(fields[1]), int(fields[2]), int(fields[3])) )

    msh.nx += 4; msh.ny += 4; msh.nz += 4 # arbitrary offsetting

    meshf.close()

    print >> sys.stderr, '\nwrite complete vtu (unstructured mesh: points + cells)...'

    # construct the ijk array
    ijk = np.zeros([len(meshdata),4],dtype='int')
    n = 0
    for ijkf in meshdata:
        i,j,k,f = ijkf
        ijk[n,:] = i,j,k,f
        n += 1

    # write the complete vtk (fluid, wall, inlet, outlet)
    splt = OUTFNAME.split('.')
    if splt[-1] == 'vtu':
        msh.writeVTU( os.path.join(WDIR, OUTFNAME) ,ijk, spacing=spacing)
    else:
        msh.writeVTU( os.path.join(WDIR, OUTFNAME+'.vtu') ,ijk, spacing=spacing)

    del ijk

    if WALLFNAME != None:

        print >> sys.stderr, '\nwrite reduced vtk...'

        ijkf = np.zeros([len(meshdataWIO),4],dtype='int')

        n = 0
        for i,j,k,f in meshdataWIO:
            ijkf[n,:] = i,j,k,f
            n+=1

        # write the reduced vtk (wall, inlet, outlet)
        msh.writeVTK( os.path.join(WDIR, WALLFNAME) ,ijkf)

###################
def main(args):

    start_time = time.time()

    msh = Mesh()

    if args.outwall != None:
        bgkflag2vtu(msh, args.infile, args.outfile, args.outwall)
    else:
        bgkflag2vtu(msh, args.infile, args.outfile)

    print >> sys.stderr, "Total time:", time.time()-start_time, "secs"

    sys.exit(0)

###############################
def iniparser(parser):

    parser.add_argument('-i', '--infile', required=True,  help='input mesh file (.dat)')
    parser.add_argument('-o', '--outfile',  required=True,  help='output complete vtu file')
    parser.add_argument('-w', '--outwall',                 help='output wall vtk file')

    return parser

###################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = iniparser(parser)
    args = parser.parse_args()

    main(args)
