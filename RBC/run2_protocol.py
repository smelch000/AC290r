#!/usr/bin/env python

import sys
from MagicUniverse import *

MagicBegins()

C0 = ...
C1 = ...
GROWTIME = ...

RESTART = False

u.decorate()

# initial thin band of drug
nx,ny,nz = m.getBox()
nx = int(nx); ny = int(ny); nz = int(nz)

myid = get_myproc()

if not RESTART:
    profile2 = c.getArray(nx*ny*nz)
    for k in range(1,nz+1):
        for j in range(1,ny+1):
            for i in range(1,nx+1):

                ifl = m.getLocator(i,j,k)

                if i > nx/8 - 3 and i < nx/8. + 3: # create BOLUS
                    profile2[ifl] = C1
                else:
                    profile2[ifl] = C0

    c.setDensityProfile(profile2)

    # freeze the fluid and drug until released
    f.setFreeze(True)
    c.setFreeze(True)

    # cap RBC forces to a large roof to avoid instabilities
    a.setCapForces(True, forcecap=1.e4, torquecap=1.e4, velcap=0.4, angvelcap=0.6)
    # set a robust friction coefficient for initial equilibration of RBC
    a.setGamma(gammaT=0.1, gammaR=0.1)
    a.setZeroVelocity()

for itime in u.cycle():

  if not RESTART:

    # gently increase excluded volume
    if itime <= GROWTIME:
        tscale = 0.5 * (1. + itime/GROWTIME) 
        if myid==0 and itime%10==0: 
            print 'Rescaling interaction to', tscale
            sys.stdout.flush()
        a.scaleVdwParameters(tscale)

    #  increase plasma-RBC coupling
    elif itime == int(1.5*GROWTIME):
        f.setFreeze(False)
        a.setGamma(gammaT=0.001, gammaR=0.001)
        a.setZeroVelocity()

    # RBC free to move with right coupling
    elif itime == int(2.*GROWTIME):
        a.setGamma(gammaT=0.01, gammaR=0.01)
        a.setZeroVelocity()

    # allow the drug to move 
    elif itime == int(2.5*GROWTIME):
        c.setFreeze(False)

  # BC: stripe of given drug density
  c.setDensityStripe('x', nx, C0) 

  u.animate()

MagicEnds()
