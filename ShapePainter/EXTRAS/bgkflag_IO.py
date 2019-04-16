#!/usr/bin/env python

import sys
import numpy as np

def bgkflag_hdr_read(fname):

    f = open(fname,'r')
   
    line = f.readline()

    line = line.split()

    nx,ny,nz = int(line[0]), int(line[1]), int(line[2])

    line = f.readline()

    line = line.split()

    iskip = int(line[0])
    nfluid = int(line[1])
    nwall = int(line[2])
    ninlet = int(line[3])
    noutlet = int(line[4])

    return nx,ny,nz,iskip,nfluid,nwall,ninlet,noutlet

def bgkflag_dat_read(fname):

    f = open(fname,'r')
   
    itmp = []
    ifla = []
    nfluid = 0; nwall = 0; ninlet = 0; noutlet = 0
    for line in f.readlines():
        line = line.split()

        i = int(line[0]); j = int(line[1]); k = int(line[2]); fl = int(line[3])

        itmp.append([i-1,j-1,k-1])
        ifla.append(fl)

        if fl==1:
            nfluid += 1
        elif fl==2:
            nwall += 1
        elif fl==3:
            ninlet += 1
        elif fl==4:
            noutlet += 1

    ijkf = np.zeros((len(itmp),4),dtype='int')

    for n in xrange(len(itmp)):
        i,j,k = itmp[n]
        fl = ifla[n]
        ijkf[n,:] = i,j,k,fl

    return ijkf,nfluid,nwall,ninlet,noutlet


def bgkflag_dat_memwaste_read(fname,nx,ny,nz):

    f = open(fname,'r')
   
    print 'NX,NY,NZ:',nx,ny,nz
    iflg = np.zeros((nx,ny,nz),dtype='int')

    nfluid = 0; nwall = 0; ninlet = 0; noutlet = 0
    for line in f.readlines():
        line = line.split()

        i = int(line[0]); j = int(line[1]); k = int(line[2]); fl = int(line[3])
        iflg[i-1,j-1,k-1] = fl

        if i<1 or i>nx: print 'Error in i range:',i,nx
        if j<1 or j>ny: print 'Error in j range:',j,ny
        if k<1 or k>nz: print 'Error in k range:',k,nz

        if fl==1:
            nfluid += 1
        elif fl==2:
            nwall += 1
        elif fl==3:
            ninlet += 1
        elif fl==4:
            noutlet += 1
        else:
            print __name__,'flag not found in bgkflag_dat_memwaste_read:',fl
            sys.exit(1)

    return iflg,nfluid,nwall,ninlet,noutlet



