#!/usr/bin/env python

from MagicUniverse import *
import sys

# preprocess_parallel_mesh(4); sys.exit(1)

NU = 0.16666
RHO = 0.22
NSTEP = 10000
NDIAG = 100
NVTKFREQ = 500

MagicBegins()

# define universe

u = Universe()
s = Scale()
m = Mesh()
f = Fluid()
t = Tracker()

u.addItems([s,m,f,t])

u.setTitle('Pisa/Massa 000_patient_m')
u.setNumberOfSteps(NSTEP)
u.setTemperature(0.0)
u.setStateRestart(False)
u.setStateDumpFrequency(-1)

u.create()

# set params

s.set(name='MonoScale', mesh=m,  actors=[f,t])

m.setRegularMesh(False)

m.setPeriodicity('000')
m.setDomainDecomposition(8)

t.setDiagnosticFrequency(NDIAG)
# t.setDataShow(velocity=True)
# t.setMapDirections('zx')

t.setVtkDump(True, meshtype='unstructured', frequency=NVTKFREQ)

f.setName('BloodFlow')

f.setDensityUniform(RHO)
f.setViscosity(NU)

f.setInletOutletMethod('zouhe')
# f.setInletOutletMethod('equilibrium')
# f.setInletOutletMethod('closed')

f.setStabilizeLB(True)

u.decorate()

for itime in u.cycle():

	u.animate()

MagicEnds()

