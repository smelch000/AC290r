#!/usr/bin/env python

import csv
import sys
import math
import numpy as np
import argparse
from mesh_csv_IO import *

def i4back(i,j,k):
    global LNX,LNY,LNXY
    return long(k)*LNY*LNX + long(j)*LNX + long(i)

def ijk(i4):
    global LNX,LNY,LNXY

    k = i4 / LNY / LNX
    j = (i4 - long(k)*LNY*LNX) / LNX
    i = i4 - long(k)*LNY*LNX - long(j) * LNX

    return i,j,k

def check_sorting(itppp):

    dims = itppp.shape

    for ii in xrange(dims[0]-1):
        # if itppp[ii,0] >= itppp[ii+1,0] : 
        if itppp[ii] >= itppp[ii+1] : 
            print 'ordering not correct:',itppp[ii],itppp[ii+1]
            sys.exit(1)

########################
def main(args):

    infiles = args.infiles
    outfile = args.outfile

    for f in infiles:
        print '<---- infile:',f
    print '----> outfile:',outfile
    print 

    NX = 0; NY = 0; NZ = 0
    ijks = []
    nsz = 0
    for f in infiles:

        ijk1 = readmesh_csv(f)
        ijks.append( ijk1 )

        imin,imax = np.amin(ijk1[:,0]), np.amax(ijk1[:,0])
        jmin,jmax = np.amin(ijk1[:,1]), np.amax(ijk1[:,1])
        kmin,kmax = np.amin(ijk1[:,2]), np.amax(ijk1[:,2])

        dims = ijk1.shape
        nsz += dims[0]

        # print 'nsz:',nsz,'     imax,jmax,kmax:', imax,jmax,kmax
        print '   imin,jmin,kmin:', imin,jmin,kmin,'   imax,jmax,kmax:', imax,jmax,kmax,'   nsz:', nsz

        NX = max(NX,imax)
        NY = max(NY,jmax)
        NZ = max(NZ,kmax)

    NX += 4; NY += 4; NZ += 4 # some offsetting

    LNX = long(NX)
    LNY = long(NY)
    LNXY = LNX * LNY

    print '\nNX,NY,NZ:',NX,NY,NZ

    itppp = np.zeros(nsz,dtype='int64')

    jj = 0
    for ijk1 in ijks:

        ii = 0
        for i,j,k in ijk1:
            if not k*NY*NX == long(k)*long(NY)*long(NX): print 'AHAH'

            itppp[jj] = k*NY*NX + j*NX + i
            ii += 1
            jj += 1

        # itppp = np.concatenate((itp,itppp),axis=0)

    # sanity check
    for i4 in itppp:
        i,j,k = ijk(i4)
        if not i4back(i,j,k) == i4:
            print 'fail1'
            sys.exit(1)
        if not (i,j,k) == ijk(i4):
            print 'fail2', (i,j,k), ijk(i4)
            sys.exit(1)

    itppp = np.sort(itppp,0)

    # check sorting
    check_sorting(itppp)

    print '....Writing files:', outfile+'.hdr', 'and', outfile+'.dat'

    fh = open(outfile+'.hdr','w')
    fd = open(outfile+'.dat','w')

    fh.write('%d %d %d \n' % (NX,NY,NZ))
    fh.write('%d %d %d %d %d \n' % (0,0,0,0,0))
    for i4 in itppp:

        i,j,k = ijk(i4)
        fd.write('%d %d %d %d\n' % (i,j,k,i4))

    fh.close()
    fd.close()

    print '....done!'

###############################
def iniparser(parser):

    parser.add_argument('-i', '--infile', required=True, nargs='+', help='input csv files')
    parser.add_argument('-o', '--outfile', required=True, default='bgkflag', help='output mesh prefix')

    return parser

########################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = iniparser(parser)
    args = parser.parse_args()

    main(args)

