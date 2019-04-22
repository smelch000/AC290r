#!/usr/bin/env python

import os, sys, math
import random
from ctypes import *
from numpy import empty
import argparse

import TOOLS.mesh as msh
from TOOLS.wrap_by_wall import *
import numpy as np

# RBC model parameters
__VRBC = 2.1 # volume of RBC
__DMIN_RBC = 1.0    # minimum diameter RBC's
__DMAX_RBC = 2.0    # maximum diameter RBC's

# __CUTOFF_INSERTION = 1.0 # 2.0 # cutoff for random insertion of particles
__CUTOFF_INSERTION = 1.2
__CUTOFF_WALL = .5
__GAMMAT = 0.0001
__GAMMAR = 0.0001

def generate_cylinder(RADIUS,NX,NY,NZ):

    print 'CYLINDER RADIUS:',RADIUS,'(l.u.) = ',RADIUS * 4.,'(micron)'
    print 'CYLINDER LENGTH:',NX,'(l.u.) = ',NX * 4.,'(micron)'
    print 'NX:',NX,'    NY,NZ:',NY,NZ

    MSH = msh.Mesh(nx=NX,ny=NY,nz=NZ)

    RADIUS2 = RADIUS * RADIUS

    itp_f, itp_w, itp_i, itp_o = [], [], [], []
    ii_f, jj_f, kk_f = [], [], []
    ii_w, jj_w, kk_w = [], [], []

    #
    # paint the mesh
    for k in xrange(1,MSH.nz+1):
        for j in xrange(1,MSH.ny+1):
            for i in xrange(1,MSH.nx+1):

                d = (j-MSH.ny/2)**2 + (k-MSH.nz/2)**2

                if d < RADIUS2:
                    itp_f.append( MSH.i4back(i,j,k) )
                    ii_f.append(i)
                    jj_f.append(j)
                    kk_f.append(k)
                else:
                    ii_w.append(i)
                    jj_w.append(j)
                    kk_w.append(k)

    MSH.specifyNodes(itp_f,itp_w,itp_i,itp_o)

    # complete by wrapping
    wrap_by_wall(MSH)

    print 'Nfluid nodes:',MSH.nfluid

    # mesh is complete
    # print 'writing mesh.xyz'
    # MSH.writexyz('mesh.xyz')

    MSH.writeMOEBIUSinput('bgkflag')
    return MSH, ii_f, jj_f, kk_f, ii_w, jj_w, kk_w, itp_f, itp_w, itp_i, itp_o


def insert_RBC(mesh, RADIUS, HCT, alignment):

    if mesh != None:

        MSH = Mesh()
        MSH.loadMOEBIUSinput(mesh+'.hdr', mesh+'.dat')

    else:
        NX = int(RADIUS * 6)
        NY = int(2*RADIUS) + 10
        NZ = int(2*RADIUS) + 10

        # generate mesh files for a test cylinder
        MSH, ii_f, jj_f, kk_f, ii_w, jj_w, kk_w, itp_f, itp_w, itp_i, itp_o = generate_cylinder(RADIUS,NX,NY,NZ)

    __NRBC = int(HCT  * MSH.nfluid  / __VRBC)

    # Spheres model parameter
    __NSPH = 0
    __D_SPH    = 1.0    # diameter spheres

    print 'Hematocrit:',HCT
    print 'NRBC:',__NRBC, 'Volume of RBC:', __VRBC,'l.u.'
    print 'NSPH = ',__NSPH

    # write atom.inp
    a = open('atom.inp','w')
    print >> a, \
"""
MODEL
# title
molecular types    3
_Ellipsoids
nummols    %d
atoms    1
E  1. 0.  1  0  .5 1. 1. 0   ! name, mass, chge, rept, frzt, inertia x/y/z, frzr
hydro
.false. 0.04   1.6  2.0   2.0  0.0      ! passive scalar, tumbling coeff, visc_enhancer, csi, oblate, smooth
%f    %f  ! gammas
finish
_Spheroids
nummols    %d
atoms    1
S  %f 0.  1  0  1.0 1.0 1.0 0   ! name, mass, chge, rept, frzt, inertia x/y/z, frzr
hydro
.false. 0.04   1.6  2.0   2.0  0.0      ! passive scalar, tumbling coeff, visc_enhancer, csi, oblate, smooth
%f    %f  ! gammas
finish
Wall
nummols    1
atoms    1
W
finish
vdw    6
E E gb    1.e-4   1.0  2.0  2.0    0.0
E S gb    1.e-4   %f %f   %f   0.0
S S gb    1.e-4   %f %f %f    0.0
E W lj    1.e-4   1.0  0.0  0.0    0.0
S W lj    1.e-4   %f  0.0  0.0    0.0
W W lj    1.e-4   1.0  0.0  0.0    0.0
CLOSEMODEL
""" % (__NRBC,__GAMMAT,__GAMMAR,__NSPH,(8*(__D_SPH/2)**3),__GAMMAR,__GAMMAT, \
      (__D_SPH + __DMIN_RBC)/2., (__D_SPH + __DMAX_RBC)/2., (__D_SPH + __DMAX_RBC)/2., \
      __D_SPH, __D_SPH, __D_SPH, \
      __D_SPH)

    print >> a, "CONFIG"
    print >> a, __NRBC + __NSPH

    """
    i4rand = random.sample(itp_f, __NRBC + __NSPH)
    print 'atom position uniqueness:',len(i4rand) == len(set(i4rand))
    atompos = []
    for i4 in i4rand:
        atompos.append( MSH.ijk(i4) )
    print >> a, "CONFIG"
    print >> a, len(atompos)
    V = 0.001
    for x,y,z in atompos:
        # print >> a, x,y,z,0.,0.,0.
        print >> a, x,y,z,random.uniform(-V,V),random.uniform(-V,V),random.uniform(-V,V)
    """

    # now start insersion of RBC (+ SPH)

    def dist2(a,b):

        def NINT(f):
            return int(f + 0.5)

        dx = a[0] - b[0]; dy = a[1] - b[1]; dz = a[2] - b[2]
        dx = dx - NINT(dx/NX)*NX
        # dy = dy - NINT(dy/BOX)*BOX
        # dz = dz - NINT(dz/BOX)*BOX
        return dx**2 + dy**2 + dz**2

    cut2 = __CUTOFF_INSERTION**2

    p = open('RBC.xyz','w')
    print >> p, __NRBC
    print >> p

    USE_SO_LIBRARY = True

    if USE_SO_LIBRARY:


        fnr = os.path.join(os.getenv("MOEBIUS_ROOT"), 'BACKEND/SCRIPTS/TOOLS/insert_RBC_part')
        # os.system('gfortran -cpp -DLINKCELL -g -shared -fPIC %s.f90 -o %s.so'%(fnr,fnr))
        os.system('gfortran -cpp -g -shared -fPIC %s.f90 -o %s.so'%(fnr,fnr))

        # LIB = CDLL('./preproc_part.so', RTLD_GLOBAL)
        LIB = CDLL('%s.so'%fnr, RTLD_GLOBAL)

        # NSZ = MSH.nfluid + MSH.nwall

        ii_f_l = empty(MSH.nfluid, dtype="int64")
        jj_f_l = empty(MSH.nfluid, dtype="int64")
        kk_f_l = empty(MSH.nfluid, dtype="int64")

        ii_w_l = empty(MSH.nwall, dtype="int64")
        jj_w_l = empty(MSH.nwall, dtype="int64")
        kk_w_l = empty(MSH.nwall, dtype="int64")

        if mesh != None:
            for ifl in range(MSH.nfluid):
                i4 = MSH.itppp_f[ifl]
                i,j,k = MSH.ijk(i4)
                ii_f_l[ifl] = i
                jj_f_l[ifl] = j
                kk_f_l[ifl] = k

            for ifl in range(MSH.nwall):
                i4 = MSH.itppp_w[ifl]
                i,j,k = MSH.ijk(i4)
                ii_w_l[ifl] = i
                jj_w_l[ifl] = j
                kk_w_l[ifl] = k

        else:
            ii_f_l = np.copy(ii_f)
            jj_f_l = np.copy(jj_f)
            kk_f_l = np.copy(kk_f)
            ii_w_l = np.copy(ii_w)
            jj_w_l = np.copy(jj_w)
            kk_w_l = np.copy(kk_w)


        xo = empty(MSH.nfluid, dtype="double")
        yo = empty(MSH.nfluid, dtype="double")
        zo = empty(MSH.nfluid, dtype="double")

        if alignment == 'x': ialign = 1
        if alignment == 'y': ialign = 2
        if alignment == 'z': ialign = 3

        LIB.preproc_part(c_int(ialign), \
                         c_int(MSH.nx), c_int(MSH.ny), c_int(MSH.nz), \
                         c_int(MSH.nfluid), c_int(MSH.nwall), \
                         c_int(__NRBC), c_int(__NSPH), \
                         c_float(__CUTOFF_INSERTION), c_float(__CUTOFF_WALL), \
                         ii_f_l.ctypes.data_as(POINTER(c_int64)), \
                         jj_f_l.ctypes.data_as(POINTER(c_int64)), \
                         kk_f_l.ctypes.data_as(POINTER(c_int64)), \
                         ii_w_l.ctypes.data_as(POINTER(c_int64)), \
                         jj_w_l.ctypes.data_as(POINTER(c_int64)), \
                         kk_w_l.ctypes.data_as(POINTER(c_int64)), \
                         xo.ctypes.data_as(POINTER(c_double)), \
                         yo.ctypes.data_as(POINTER(c_double)), \
                         zo.ctypes.data_as(POINTER(c_double)))

        for i in xrange(__NRBC + __NSPH):
            print >> a, xo[i],yo[i],zo[i],0.,0.,0.,i+1
            print >> p, 'C', xo[i],yo[i],zo[i]

    else:

        rro = []
        while True:
            i4 = random.sample(itp_f, 1)
            rr = MSH.ijk(i4[0])

            tooclose = False

            d = math.sqrt( (rr[1]-MSH.ny/2)**2 + (rr[2]-MSH.nz/2)**2 )
            # if d > RADIUS - __CUTOFF_INSERTION:
            if d > RADIUS - 0.50:
                tooclose = True

            if not tooclose:
                for ro in rro:
                    if dist2(rr,ro) < cut2:
                        tooclose = True
                        break

            if not tooclose:
                # dd = 1.e5
                # for ro in rro:
                #     dd =  min(dd,math.sqrt(dist2(rr,ro)))
                # print 'accepted ....mindist : ', dd

                rro.append(rr)
                if len(rro)%1000 == 0: print len(rro),'/',__NRBC + __NSPH
                if len(rro) >= __NRBC + __NSPH: 
                    break

        for i in xrange(__NRBC + __NSPH):
            r = rro[i]
            print >> a, r[0],r[1],r[2],0.,0.,0.,i+1
            print >> p, 'C', r[0],r[1],r[2]

    for i in xrange(__NRBC + __NSPH):
        if i%3==0:
            print >> a, '1. 0. 0. 0. 0. 0.',i+1
        elif i%3==1:
            print >> a, '0. 1. 0. 0. 0. 0.',i+1
        elif i%3==2:
            print >> a, '0. 0. 1. 0. 0. 0.',i+1

    print >> a, 'CLOSECONFIG'

    a.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-H', '--hematocrit', required=True,  help='hematocrit in percentile (0.0: ~0.50 range)')
    parser.add_argument('-m', '--mesh', required=False,  help='input mesh file')
    parser.add_argument('-R', '--radius', required=False,  help='cylinder radius (l.u.)')
    parser.add_argument('-A', '--alignment', required=True,  help='alignment of system (x, y or z)')

    args = parser.parse_args()

    HCT = float(args.hematocrit)
    if HCT < 0. or HCT > 0.60:
        print 'hematrocit level unreasonable:', HCT, '...set in the [0-0.50] range and rerun'
        sys.exit(1)

    alignment = args.alignment
    if alignment not in ['x','y','z']:
        print 'alignment can only be x,y or z'
        sys.exit(1)

    if args.mesh == None:
        if args.radius == None:
            print 'For the test case, insert cylinder radius'
            sys.exit(1)

        RADIUS = float(radius)

    insert_RBC(args.mesh, RADIUS, HCT, alignment)

