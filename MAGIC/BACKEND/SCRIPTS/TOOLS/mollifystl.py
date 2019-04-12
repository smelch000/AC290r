#!/usr/bin/env python

import sys,time,copy
import numpy
from math import *

from TOOLS.mesh import *
from TOOLS.stl import *

# input STL file
INSTL = 'mdc_2_0.stl'
# distance from inlet/outlet wehre IB nodes are frozen
CUT_IOFREEZE = 12.0 
# force constant for tethering to initial position
KTETH = 1.e-6 
# force constant for bond forces
KBOND = 0.2 

# mesh input files
INHDR = 'bgkflag.hdr'
INDAT = 'bgkflag.dat'

# mesh input files
OUTATOM = 'atom.inp'

# CUT_TRIAREA = 0.12

#################################
msh = Mesh()

msh.loadMOEBIUSinput(INHDR,INDAT)

scale,transl = msh.loadtransformfile('mesh_transform.inp')

stl = Stl()

stl.points_rev,stl.tris = stl.readfile()

g = open(OUTATOM,'w')

g.write(
'CONTROL\n\
   1   1   0.2        ! MD sub-sub-iterations,sub-iterations,smallest timestep\n\
.false.     ! rotational motion\n\
5          ! 1:NGP 2:IB 4:BOUZIDI 5:DELTA\n\
0.   0.   0.   ! force\n\
.false.     ! linsert\n\
CLOSECONTROL\n\
\n\
MODEL\n\
# test case\n\
molecular types    2\n\
Test_Particles\n\
nummols    1 \n')

npoints = len(stl.points_rev)

tag = numpy.zeros([msh.nx,msh.ny,msh.nz], dtype='i')

NR=3
FR=10
for m in xrange(msh.ninlet):
    i,j,k = msh.ijk( msh.itppp_i[m] )
    tag[i-FR:i+FR,j-FR:j+FR,k-FR:k+FR] = 1 # zone around inlet

for m in xrange(msh.ninlet):
    i,j,k = msh.ijk( msh.itppp_i[m] )
    tag[i-NR:i+NR,j-NR:j+NR,k-NR:k+NR] = 2 # zone inside inlet

for m in xrange(msh.noutlet):
    i,j,k = msh.ijk( msh.itppp_o[m] )
    tag[i-FR:i+FR,j-FR:j+FR,k-FR:k+FR] = 1 # zone around outlet

for m in xrange(msh.noutlet):
    i,j,k = msh.ijk( msh.itppp_o[m] )
    tag[i-NR:i+NR,j-NR:j+NR,k-NR:k+NR] = 2 # zone inside outlet

atoms=[]
present=[]
n2=0; n1=0
stl_tris3 = set()
# stl_tris3 = copy.deepcopy(stl.tris)
count = 0
for n in xrange(npoints):

    p = stl.points_rev[n]

    x = p[0]*scale + transl[0]
    y = p[1]*scale + transl[1]
    z = p[2]*scale + transl[2]

    i,j,k = int(x),int(y),int(z)

    present.append(-1)

    if tag[i,j,k]==2: 
        n2+=1
        present[n] = -1
        continue # exclude atoms in inlet and outlet
    elif tag[i,j,k]==1: 
        n1+=1

    atoms.append([x,y,z])
    present[n] = count
    count += 1

for vtx in stl.tris:
    i1,i2,i3=vtx
    # print i1,i2,i3,present[i1],present[i2],present[i3]

    if present[i1]<0 or present[i2]<0 or present[i3]<0: continue

    j1,j2,j3 = present[i1],present[i2],present[i3]
    stl_tris3.add((j1,j2,j3))

natms = len(atoms)

print 'NVERTEX:',natms,'(near i/o:',n1,'(far i/o:',n2


g.write('atoms %d\n' % (natms))

xyz = open('STL.xyz','w')
xyz.write('%d\n\n' % natms)

FAST = True

nfreeze = 0
freeze = []
c2 = CUT_IOFREEZE**2
for n in xrange(natms):

    if n%100==0:
        print n,'.',
        sys.stdout.flush()

    x,y,z = atoms[n]

    frozen = False
    if FAST:
        i,j,k = int(x),int(y),int(z)
        frozen = (tag[i,j,k]>0)

    else:
        for m in xrange(msh.ninlet):
            i,j,k = msh.ijk( msh.itppp_i[m] )

            d = (x-i)**2 + (y-j)**2 + (z-k)**2
            if d < c2:
                frozen = True
                break

        if not frozen:
            for m in xrange(msh.noutlet):
                i,j,k = msh.ijk( msh.itppp_o[m] )

                d = (x-i)**2 + (y-j)**2 + (z-k)**2
                if d < c2:
                    frozen = True
                    break

    if frozen:
        freeze.append(True)
        nfreeze += 1
        g.write('B  1.  0.  1  1  1. 1. 1. 0  ! name, mass, chge, rept, frzt, inertia x/y/z, frzr\n')
        xyz.write('O %f %f %f\n' % (x,y,z))

    else:
        freeze.append(False)
        g.write('A  1.  0.  1  0  1. 1. 1. 0  ! name, mass, chge, rept, frzt, inertia x/y/z, frzr\n')
        xyz.write('C %f %f %f\n' % (x,y,z))

xyz.close()

print 'frozen ',nfreeze,'/', natms

l = open('links.tcl','w')
for vtx in stl_tris3:
    i1,i2,i3 = vtx
    l.write('topo addbond %d %d\n' % (i1,i2))
    l.write('topo addbond %d %d\n' % (i1,i3))
    l.write('topo addbond %d %d\n' % (i2,i3))
l.close()

######## BONDS
# bonds=[]
bonds=set()
for vtx in stl_tris3:

    i1,i2,i3 = vtx

    p1 = atoms[i1]
    p2 = atoms[i2]
    p3 = atoms[i3]

    aa = [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]
    bb = [p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2]]
    cc = [p3[0]-p2[0], p3[1]-p2[1], p3[2]-p2[2]]

    a = sqrt(aa[0]**2 + aa[1]**2 + aa[2]**2)
    b = sqrt(bb[0]**2 + bb[1]**2 + bb[2]**2)
    c = sqrt(cc[0]**2 + cc[1]**2 + cc[2]**2)

    # if not freeze[i1] and not freeze[i2]: bonds.append([i1,i2,a])
    # if not freeze[i1] and not freeze[i3]: bonds.append([i1,i3,b])
    # if not freeze[i2] and not freeze[i3]: bonds.append([i2,i3,c])
    if not (freeze[i1] and freeze[i2]): bonds.add((i1,i2,a))
    if not (freeze[i1] and freeze[i3]): bonds.add((i1,i3,b))
    if not (freeze[i2] and freeze[i3]): bonds.add((i2,i3,c))

nbonds = len(bonds)

g.write('bonds %d\n' % nbonds)
for b in bonds:
    i1,i2,dist = b
    g.write('harm %d %d %f %f\n' % (i1+1,i2+1,KBOND,dist))

######## TETHERS
g.write('tether %d\n' % (natms-nfreeze))
for n in range(natms):
    x,y,z = atoms[n]
    if not freeze[n]:
       g.write('harm %d %f %f %f %f\n' % (n+1,KTETH,x,y,z))


######## HYDRO and VDW

g.write('\
hydro\n\
.false.   0.   1.  1.   1.  0.  ! passive scalar,tumbling coeff,visc_enhancer,csi,oblate,smooth\n\
0.5    0.0   ! gammas: translational and rotational coupling coeffs\n\
finish\n\
Wall\n\
nummols    1\n\
atoms    1\n\
W\n\
finish\n\
\n\
vdw      6\n\
A   A    lj 0.     0.  0.   0.   0.  ! epsilon/sigma/null/null/cutoff\n\
A   B    lj 0.     0.  0.   0.   0. \n\
B   B    lj 0.     0.  0.   0.   0. \n\
A   W    lj 0.     .5  0.   0.   .99 \n\
B   W    lj 0.     .5  0.   0.   .99 \n\
W   W    lj 0.     0.  0.   0.   0.\n\
CLOSEMODEL\n\
\n\
')

####### CONFIG #########
g.write('\n\n')
g.write('CONFIG  \n')
g.write(str(natms)+'\n')

for n in xrange(natms):

    x,y,z = atoms[n]

    g.write("%10.5f %10.5f %10.5f 0. 0. 0. %d \n" % (x,y,z,n))

g.write('ENDCONFIG  \n')
g.close()
