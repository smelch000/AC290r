#!/usr/bin/env python

import sys
import numpy
from scipy import ndimage
import argparse

from TOOLS.mesh import *
from TOOLS.d3q19 import *
from TOOLS.inletoutlet import *

SLOWMETHOD__ = False

def wrap_by_wall(msh):

    lat = D3q19(msh)

    # striping off previous walls
    print >> sys.stderr, 'stripping off previous walls...'
    msh.nwall = 0
    msh.itppp_w = np.zeros(msh.nwall, dtype='i8')

    itppp = np.copy(msh.itppp_f)

    print >> sys.stderr, '# of fluids:',len(msh.itppp_f)
    print >> sys.stderr, '# of inlets:',len(msh.itppp_i)
    print >> sys.stderr, '# of outlets:',len(msh.itppp_o)

    if len(msh.itppp_i)>0: itppp = np.concatenate( (itppp, msh.itppp_i), axis=0)
    if len(msh.itppp_o)>0: itppp = np.concatenate( (itppp, msh.itppp_o), axis=0)

    print >> sys.stderr, 'finding neighbors...'

    # create the neighbors layer around all nodes
    newcomers = []

    if SLOWMETHOD__ :

        for i4 in itppp:

            pnt = msh.ijk(i4)

            # find all neighbors of point (fluid, inlet or outlet)

            newneigh = lat.getneighbors_fast(pnt)

            # if pnt[2]==3: print pnt[:2] # ,':',msh.ijklist(newneigh)
            for ii4 in newneigh:

                # ii,jj,kk = msh.ijk(ii4)
                # if pnt[2]==3 and kk==3: print '    Neigh:', ii,jj

                newcomers.append(ii4)
        print >> sys.stderr, '\nwrap newcomers...',len(newcomers)

        # remove duplicates
        newcomers = list(set(newcomers))

        print >> sys.stderr, ' ====> removed duplicates, wall size:',len(newcomers)

        msh.itppp_w = numpy.zeros(len(newcomers), dtype='i8')

        msh.nwall = 0
        for i4 in newcomers:
            msh.itppp_w[msh.nwall] = i4
            msh.nwall += 1
 
        msh.itppp_w.sort()
 
    else:

        print '...dilation method valid for d3q19 lattice...'

        DEAD_NODE=0
        FLUID_NODE=1
        INLET_NODE=3
        OUTLET_NODE=4

        # numpy used in the [0:N-1] range

        # nodes = np.full([msh.nx/msh.gridspacing+1, msh.ny/msh.gridspacing+1, msh.nz/msh.gridspacing+1], \

        nodes = np.full([msh.nx/msh.gridspacing, msh.ny/msh.gridspacing, msh.nz/msh.gridspacing], \
                        DEAD_NODE, \
                        np.uint)

        imn=1e6; imx=-1e6; jmn=1e6; jmx=-1e6
        for i4 in msh.itppp_f:
            i,j,k = msh.ijk(i4)
            i-=1
            j-=1
            k-=1
            imn = min(imn,i); imx = max(imx,i); jmn = min(jmn,j); jmx = max(jmx,j)
            nodes[i/msh.gridspacing,j/msh.gridspacing,k/msh.gridspacing] = FLUID_NODE

        for i4 in msh.itppp_i:
            i,j,k = msh.ijk(i4)
            i-=1
            j-=1
            k-=1
            imn = min(imn,i); imx = max(imx,i); jmn = min(jmn,j); jmx = max(jmx,j)
            nodes[i/msh.gridspacing,j/msh.gridspacing,k/msh.gridspacing] = FLUID_NODE # done on purpose...

        for i4 in msh.itppp_o:
            i,j,k = msh.ijk(i4)
            i-=1
            j-=1
            k-=1
            imn = min(imn,i); imx = max(imx,i); jmn = min(jmn,j); jmx = max(jmx,j)
            nodes[i/msh.gridspacing,j/msh.gridspacing,k/msh.gridspacing] = FLUID_NODE # done on purpose...

        dead_nodes=(nodes == DEAD_NODE)

        kernel_mat = ndimage.generate_binary_structure(3, 2)
        assert kernel_mat.sum() == 19
        fluid_nodes_dilated = ndimage.morphology.binary_dilation(nodes==FLUID_NODE, 
                                                                 structure=kernel_mat,
                                                                 border_value=-99)

        wall_nodes = fluid_nodes_dilated & dead_nodes

        q = 0
        msh.nwall = wall_nodes.sum()
        msh.itppp_w = np.zeros(msh.nwall,np.int64)
        imn=1e6; imx=-1e6; jmn=1e6; jmx=-1e6
        for (i,j,k),w in np.ndenumerate(wall_nodes):

          if w:
            msh.itppp_w[q] =msh.i4back(i*msh.gridspacing + 1, 
                                       j*msh.gridspacing + 1, 
                                       k*msh.gridspacing + 1)
            q += 1
        msh.itppp_w.sort()

###################
"""
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input',   required=True,  help='input moebius mesh file')
    parser.add_argument('-o', '--output',  required=True,  help='output moebius mesh file')
    parser.add_argument('-v', '--vtmoutput',               help='output multiblock VTM file')
    parser.add_argument('-c', '--conf',                    help='output .xyz file')
    args = parser.parse_args()

"""

###################
def wrapit(infile,outfile,vtmoutput=None,conf=None):

    # input files
    INHDR = infile+'.hdr'
    INDAT = infile+'.dat'
    INIOS = infile+'.ios'

    #output files
    OUTPUT  = outfile # root of .hdr/.dat filename
    OUTVTM  = vtmoutput
    OUTCONF = conf

    print >> sys.stderr  
    print >> sys.stderr, 'Running wrap_by_wall...'
    print >> sys.stderr, 'Inputs :',INHDR,' / ',INDAT,' / ',INIOS
    print >> sys.stderr, 'Outputs:',OUTPUT+'.hdr',' / ',OUTPUT+'.dat',' / ',OUTPUT+'.ios',
    if OUTVTM: print >> sys.stderr, ' / ',OUTVTM,
    if OUTCONF: print >> sys.stderr, ' / ',OUTCONF,
    print >> sys.stderr, '\n'

    msh = Mesh()

    msh.loadMOEBIUSinput(INHDR,INDAT)

    lat = D3q19(msh)

    # striping off previous walls
    print >> sys.stderr, 'stripping off previous walls...'
    msh.nwall = 0
    msh.itppp_w = np.zeros(msh.nwall, dtype='i8')

    itppp = np.copy(msh.fluid)

    print >> sys.stderr, '# of inlets:',len(msh.inlet)
    print >> sys.stderr, '# of outlets:',len(msh.outlet)

    if len(msh.inlet)>0:  itppp = np.concatenate( (itppp, msh.inlet),  axis=0)
    if len(msh.outlet)>0: itppp = np.concatenate( (itppp, msh.outlet), axis=0)

    print >> sys.stderr, 'finding neighbors...please be patient...it can take 30 minutes or longer'

    # create the neighbors layer around all nodes
    newcomers = []

    for pnt in itppp:
        # find all neighbors of point (fluid, inlet or outlet)
        newneigh = lat.getneighbors(pnt)

        for i4 in newneigh:
            newcomers.append(i4)

    print >> sys.stderr, '\nwrap newcomers...',len(newcomers)

    # remove duplicates
    newcomers = list(set(newcomers))

    print >> sys.stderr, ' ====> removed duplicates, len:',len(newcomers)

    msh.itppp_w = numpy.zeros(len(newcomers), dtype='i8')

    msh.nwall = 0
    for i4 in newcomers:
        msh.itppp_w[msh.nwall] = i4
        msh.nwall += 1
 
    msh.itppp_w.sort()
 

    #############################

    # dump .hdr/.dat files
    io_missing = msh.writeMOEBIUSinput(OUTPUT)

    # dump inletoutlet.inp file only if inlet / outlet nodes are present
    if not io_missing:
        io = Inletoutlet(msh)

        isfile = io.loadMOEBIUSinput(INIOS)

        if isfile:

            iohd = []
            iohd.append('1  inlet TBA 0 0 1   0. 0. 1.   0.\n')
            iohd.append('2 outlet TBA 0 0 1   0. 0. 1.   0.\n')

            io.writeMOEBIUSinput(OUTPUT,iohead=iohd)

    # dump .xyz
    if OUTCONF:
        msh.writexyz(OUTCONF)

    # dump .vtm
    if OUTVTM:
        msh.writeVTM(OUTVTM)

###################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--infile',   required=True,  help='input moebius mesh file root')
    parser.add_argument('-o', '--outfile',  required=True,  help='output moebius mesh file root')
    parser.add_argument('-v', '--vtmoutput',               help='output multiblock VTM file')
    parser.add_argument('-c', '--conf',                    help='output .xyz file')
    args = parser.parse_args()

    wrapit(args.infile,args.outfile,args.vtmoutput,args.conf)

