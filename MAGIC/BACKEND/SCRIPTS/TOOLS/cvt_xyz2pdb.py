#!/usr/bin/env python

import numpy as np
import sys

ftem = open('template.pdb','r')
if len(sys.argv) > 1:
    fxyz = open(sys.argv[1],'r')
    fpdb = open(sys.argv[1]+'.pdb','w')
else:
    fxyz = open('CONF.xyz','r')
    fpdb = open('CONF.pdb','w')

# nrept = int(raw_input('number of clones -> \n'))

"""
box0 = float(raw_input('box edge-> \n'))
box=[box0,box0,box0]
def refold(x,y,z,xr,yr,zr):
    global box

    dx = x - xr
    dx = dx - round(dx/box[0]) * box[0]

    dy = y - yr
    dy = dy - round(dy/box[1]) * box[1]

    dz = z - zr
    dz = dz - round(dz/box[2]) * box[2]

    return xr+dx, yr+dy, zr+dz
"""

#x0 = np.zeros(nrept)
#y0 = np.zeros(nrept)
#z0 = np.zeros(nrept)

frame = 0
while True:

    natms = int(fxyz.readline())
    fxyz.readline()

    # natprot = natms / nrept
    print 'new frame',frame+1,'natms:',natms # ,'natprot:',natprot

    fpdb.write('MODEL\n')
    fpdb.write('HEADER\n')

    if frame==0:
        lab = []
        ftem.readline()
        for l in ftem.readlines():
            if l.split()[0] == 'TER':
                continue
            lab.append(l)

    xx = np.zeros(natms,np.float32)
    yy = np.zeros(natms,np.float32)
    zz = np.zeros(natms,np.float32)
    id = np.zeros(natms,np.int32)

    for i in xrange(natms):
       l = fxyz.readline() 
       l = l.split()

       id[i] = int(l[4]) - 1

       xx[id[i]] = float(l[1])
       yy[id[i]] = float(l[2])
       zz[id[i]] = float(l[3])

    """
    for n in xrange(natms):
        for ip in xrange(natprot):
            i = natprot*n + ip

            if ip==0:
                if frame>0: 
                    xx[i],yy[i],zz[i] = refold(xx[i],yy[i],zz[i],x0[n],y0[n],z0[n])
                x0[n] = xx[i]; y0[n] = yy[i]; z0[n] = zz[i]
            else:
                xx[i],yy[i],zz[i] = refold(xx[i],yy[i],zz[i],xx[i-1],yy[i-1],zz[i-1])
    """

    print 'natms:',natms
    for i in xrange(natms):
        fpdb.write('%.30s%8.3f%8.3f%8.3f\n' % (lab[i], xx[i], yy[i], zz[i]))
        # fpdb.write('%.30s%8.2f%8.2f%8.2f\n' % (lab[i%natprot], xx[i], yy[i], zz[i]))

    fpdb.write('ENDMDL\n')

    frame+=1
        




