#!/usr/bin/env python

from __future__ import division

import pint
from uncertainties import ufloat
from uncertainties.umath import *
import numpy as np

class Units():
    def __init__(self):
        self.ureg = pint.UnitRegistry()
        self.Q_ = self.ureg.Quantity

        self.dx = None
        self.dt = None
        self.dm = None
        self.Uchar = None
        self.Machar = None

        self.u = self.Q_( 0.1 ) # sweet spot velocity

        self.rho = self.Q_( 1.0 ) # default density
        
    def setMeshSpacing(self, dx): 
        self.dx = self.Q_(dx).to('m')
    def getMeshSpacing(self): return self.dx

    # def getTimeStep(self): return self.dx * self.u / self.Uchar
    def getTimeStep(self): 
        return self.dx * self.u / self.Uchar

    def setCharacteristicVelocity(self, Uchar): self.Uchar = self.Q_(Uchar).to('m/s')
    def getCharacteristicVelocity(self): return self.Uchar

    def setCharacteristicMass(self, Mchar): self.Mchar = self.Q_(Mchar).to('kg')
    def getCharacteristicMass(self): return self.Mchar

    def setCharacteristicPressure(self, Pchar): self.Pchar = self.Q_(Pchar).to('Pa')
    def getCharacteristicPressure(self): return self.Pchar

    def setCharacteristicTemperature(self, Tchar): self.Tchar = self.Q_(Tchar).to('C')
    def getCharacteristicTemperature(self): return self.Tchar

    def setCharacteristicMach(self, Machar): self.Machar = self.Q_(Machar)
    def getCharacteristicMach(self): return self.Machar

    def setPhysicalViscosity(self, nuPhys): self.nuPhys = self.Q_(nuPhys).to('m^2/s')
    def getPhysicalViscosity(self): return self.nuPhys

    def setPhysicalDensity(self, rhoPhys): self.rhoPhys = self.Q_(rhoPhys).to('kg/m^3')
    def getPhysicalDensity(self): return self.rhoPhys

    def setSweetSpotVelocity(self, u): 
        self.u = self.Q_(u) # adimensional
        self.dt = self.u * self.dx / self.Uchar
    def getSweetSpotVelocity(self): return self.u
    def getTimestep(self): 
        self.dt = self.u * self.dx / self.Uchar
        return self.dt

    def getVelocity(self):
        return self.Uchar * self.dt / self.dx

    def getViscosity(self):
        return self.nuPhys * self.dt / self.dx**2

    def setMass(self, dm): 
        self.dm = self.Q_(dm).to('kg') # overrides default value

    def getMass(self): 
        self.dm = self.rhoPhys * self.dx**3
        return self.dm

    def getDensity(self):
        if self.dm != None:
            return self.rhoPhys * self.dx**3 / self.dm
        else:
            return self.rho

if __name__ == '__main__':

    u = Units()

    u.setMeshSpacing('0.1 cm')
    print 'mesh spacing:',u.getMeshSpacing()

    u.setCharacteristicVelocity('3 cm/s')
    print 'char velocity:',u.getCharacteristicVelocity()

    u.setSweetSpotVelocity('0.1')
    print 'sweet spot velocity:',u.getSweetSpotVelocity()

    print '-> timestep:', u.getTimestep()
    print

    u.setCharacteristicPressure('10. kg/m/s^2')
    print 'Pchar:',u.getCharacteristicPressure()
    u.setCharacteristicPressure('10. Pa')
    print 'Pchar:',u.getCharacteristicPressure(), '=', u.getCharacteristicPressure().to('millibar')
    print

    u.setPhysicalViscosity('1.e-5 m^2/s')
    print 'nuPhys:',u.getPhysicalViscosity()
    print '-> nu:',u.getViscosity()
    print

    u.setPhysicalDensity('1.e3 kg/m^3')
    print 'rhoPhys:',u.getPhysicalDensity()
    # u.setMass('1.e-2 kg') # optional
    print '-> density:',u.getDensity(), '-> mass:',u.getMass()
    print


