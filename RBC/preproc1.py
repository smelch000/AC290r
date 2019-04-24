#!/usr/bin/env python

from ShapePainter import *

Re = 2
visc = 0.1
u = 0.01

RADIUS = int(Re * visc / (2.*u))

# RADIUS = 10
LENGTH = 20 * RADIUS

print ('RADIUS:',RADIUS,'LENGTH:',LENGTH)


SP( NX = 2*RADIUS + 4, 
    NY = 2*RADIUS + 4, 
    NZ = LENGTH, 
    RADIUS = RADIUS,
    INBC = 'flow', OUTBC = 'pressure',
    INVAL = '0.01', OUTVAL = '0.0',
    OUTSTL = 'SP.stl',
    VISUAL = False,
    PERIODIC = 'z')

