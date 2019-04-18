#!/usr/bin/env python

import sys
import numpy
from TOOLS.mesh import *
from TOOLS.d3q19 import *
from TOOLS.inletoutlet import *
import argparse

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

