#!/usr/bin/env python

Re = 2
visc = 0.1
u = 0.01

RADIUS = int(Re * visc / (2.*u))

# RADIUS = 10
LENGTH = 20 * RADIUS

print ('RADIUS:',RADIUS,'LENGTH:',LENGTH)

from TOOLS.insert_RBC import *

insert_RBC('bgkflag', RADIUS, 0.30, 'z')

