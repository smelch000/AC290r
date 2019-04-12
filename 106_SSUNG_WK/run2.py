#!/usr/bin/env python

try:
    from MagicUniverse import *
except:
    from muphy2wrapper import *

try:
    MagicBegins()
except:
    muphy2parser()

# define universe
u = Universe()
s = Scale()
m = Mesh()
f = Fluid()
t = Tracker()

u.addItems([s,m,f,t])

u.setTitle('106_SSUNG_WK')
u.setNumberOfSteps(5000)
u.setStateRestart(False)
u.setStateDumpFrequency(10000)

u.create()

# set params
s.set(name='MonoScale', mesh=m,  actors=[f,t])

m.setRegularMesh(False)
m.setPeriodicity('000')
# m.setSystemExchange('inlet/outlet')
m.setDomainDecomposition(8)

t.setDiagnosticFrequency(100)
# t.setDataShow(velocity=True)
# t.setMapDirections('zx')

t.setVtkDump(True, meshtype='unstructured', frequency=100)

f.setName('BGK')
f.setCollisionType('BGK')

f.setDensityUniform(1.0)
f.setViscosity(0.16667)

# f.setInletOutletMethod('zouhe')
f.setInletOutletMethod('equilibrium')

f.setStabilizeLB(True)

u.decorate()

# prepare.check(f)

for itime in u.cycle():

    u.animate()

try:
    MagicEnds()
except:
    muphy2wrapper_finish()
