#!/usr/bin/env python3

import sys,os,stat

import ctypes
from ctypes import *
import socket
import pickle
import math
import subprocess
import shlex
import tempfile
import argparse
import time
import threading
import resource
import io
import Mlogging

try:
    import muphy2graphics as M2G
except:
    M2G = None

try: # c_bool defined from python 2.5 on
    c_bool()
except:
    c_bool=c_int

"""
module for muphy2
"""

class Empty():
    pass

Magic = None

M = Empty()

##########################################

# def MagicBegins(serialrun=False):
class MagicBegins():

    def __init__(self, serialrun=False):
        global Magic

        """
        parse and init
        """

        Magic = self

        self.XPU = 'cpu'
        self.SOCKET = False
        self.MUPHYorMOBIUS = '__None__'

        f = open('.MagicPid','w')
        # print >> f, os.getpid()
        f.write( str(os.getpid()) )
        f.close()

        # Completely replace the stdout file descriptor for this process and any 
        # subprocesses, we can even capture output from child processes.
        # os.dup2(stdout, 1)
        # print read_pipe()

        parser = argparse.ArgumentParser()

        parser.add_argument('-x', '--xpu',     required=False, default='cpu', help='device [cpu/gpu]')
        parser.add_argument('-d', '--debug',  required=False, default='no',  help='debug output [yes/no]')
        parser.add_argument('-p', '--preproc', required=False, default=None,  help='parallel preprocessing on given number of tasks [int]')
        # parser.add_argument('-s', '--socket',  required=False, default='no',  help='socket port active [yes/no]')
        args = parser.parse_args()

        Mlogging.Mlogging(args.debug)

        # done with logging defs

        self.XPU = args.xpu

        try:
            MENV = os.path.join(os.environ.get('MUPHY_ROOT'), 'EXECUTE')
            self.MUPHYorMOBIUS = 'Muphy'
        except:
            MENV = os.path.join(os.environ.get('MOEBIUS_ROOT'), 'BACKEND', 'SHOP')
            self.MUPHYorMOBIUS = 'Moebius'

        if self.XPU != 'cpu' and self.XPU != 'gpu':
            Mlogging.Mlogging.debug( 'xpu:' + self.XPU + '...not of the right kind!' )
            sys.exit(1)

        elif self.XPU == 'cpu':
            if self.MUPHYorMOBIUS == 'Muphy':
                THELIB = os.path.join( MENV, 'libmuphy2.so')
            elif self.MUPHYorMOBIUS == 'Moebius':
                THELIB = os.path.join( MENV, 'libmoebius.so')

        elif self.XPU == 'gpu':
            if self.MUPHYorMOBIUS == 'Muphy':
                THELIB = os.path.join( MENV, 'libmuphy2.gpu.so')
            elif self.MUPHYorMOBIUS == 'Moebius':
                THELIB = os.path.join( MENV, 'libmoebius.gpu.so')

        if not os.path.isfile(THELIB):
            print('Shared library does not exist:' + THELIB,' ... Halting')
            Mlogging.Mlogging.debug( 'Shared library does not exist:' + THELIB )
            sys.exit(1)

        M.MUPHY = CDLL(THELIB, RTLD_GLOBAL)

        #f = io.BytesIO()
        #with Mlogging.stdout_redirector(f):
        #    print('foobar')
        #    print(12)
        #    M.MUPHY.testprint()
        #    os.system('echo and this is from echo')

        self.SOCKET = False
        #if args.socket=='yes':
        #    self.SOCKET = True

        if args.preproc:
            preprocess_parallel_mesh( int(args.preproc) )
            sys.exit(1)

        muphy2wrapper_init(mycomm=0,serialrun=c_int(serialrun))

##########################################

def muphy2wrapper_init(mycomm=0,serialrun=False):
    """
    initialize wrapper and sockers
    -> M.MUPHY.muphy2wrapper_init
    """

    global Magic

    Mlogging.Mlogging.debug( 'Run starts on host: ' + socket.gethostname() )

    # M.fill_muphy_object(version='2', xpu=Magic.XPU)

    try:
        if Magic.SOCKET:
            M.set_socket_params()
    except:
        Magic.SOCKET=False
        pass

    # M.MUPHY.muphy2wrapper_init(mycomm,len(str),str)

    M.MUPHY.muphy2wrapper_init(c_int(mycomm),c_bool(serialrun))

    if get_numprocs() >  1: return

    if Magic.SOCKET:
        try:
            M.req_q = mxqueue.MXQueue()
            sock = MySocket(M.req_q)
        except:
            print('socket already in use....not connecting')
            Magic.SOCKET=False

def print_resource_usage(msg): 

    myid_l = (c_int*1)()
    M.MUPHY.getMyproc(byref(myid_l))
    myid = myid_l[0]
    if myid != 0: return

    usage = resource.getrusage(resource.RUSAGE_SELF)

    # for name, desc in [
    #     ('ru_utime', 'User time'),
    #     ('ru_stime', 'System time'),
    #     ('ru_maxrss', 'Max. Resident Set Size'),
    #     ('ru_ixrss', 'Shared Memory Size'),
    #     ('ru_idrss', 'Unshared Memory Size'),
    #     #('ru_isrss', 'Stack Size'),
    #     #('ru_inblock', 'Block inputs'),
    #     #('ru_oublock', 'Block outputs'),
    #     ]:
    #     print '%-25s (%-10s) = %s' % (desc, name, getattr(usage, name))

    allmsg = '%-20s %-20s: %d MB      %-13s: %d MB     %-15s: %d MB' % \
                        (msg, 'Max. Resident Memory', int(getattr(usage, 'ru_maxrss')/1048576.),
                         'Shared Memory', int(getattr(usage, 'ru_ixrss')/1048576.),
                         'Unshared Memory', int(getattr(usage, 'ru_idrss')/1048576.))

    print(100*'~')
    print(allmsg)
    print(100*'~')

    Mlogging.Mlogging.info( allmsg )

def MagicEnds(): 
    """
    finalize run
    -> M.MUPHY.muphy2wrapper_finish()
    """

    print_resource_usage('Final')

    wrapper_finish = M.MUPHY.muphy2wrapper_finish()

    del M.MUPHY

    print('Magic ends, bye.')

##########################################
def get_myproc():
    """
    get my processor
    """
    myproc = (c_int*1)()

    M.MUPHY.getMyproc(byref(myproc))

    return myproc[0]
##########################################
def get_numprocs():
    """
    get number of processors
    """
    numprocs = (c_int*1)()
    M.MUPHY.getNumprocs(byref(numprocs))
    return numprocs[0]

##########################################
def sumpara(x):
    """
    parallel sum
    """
    if type(x) == int:
        x_ = (c_int*1)()
        x_[0] = x
        M.MUPHY.sum_dd_i(byref(x_))
        return x_[0]

    elif type(x) == float:
        x_ = (c_float*1)()
        x_[0] = x
        M.MUPHY.sum_dd_f(byref(x_))
        return x_[0]

    elif type(x) == list:

        n = len(x)
        if type(x[0]) == int:
            x_ = (c_int*n)()
            for i in range(n):
                x_[i] = x[i]
            M.MUPHY.sum_dd_ai(c_int(n), byref(x_))

        elif type(x[0]) == float:
            x_ = (c_float*n)()
            for i in range(n):
                x_[i] = x[i]
            M.MUPHY.sum_dd_a(c_int(n), byref(x_))

        return x_

    else:
        print('Hard Stop in sumpara')
        sys.exit(1)

################################
class MuArray(object):
    
    def __init__(self,mesh,allocate=True):

        self.ptr = (c_void_p*1)()
        self.mesh = mesh

        self.PRC = (c_void_p*1)()
        M.MUPHY.getFloatPrecision(byref(self.PRC))

        if self.PRC[0] == 4:
            self.c_PRC = c_float

        elif self.PRC[0] == 8:
            self.c_PRC = c_double

        self.lbound, self.ubound = mesh.getFluidNodesBounds()

        if allocate:
            M.MUPHY.initMuArray_float(c_int(self.lbound), c_int(self.ubound), self.ptr)
        else:
            M.MUPHY.initMuArray_float(c_int(self.lbound), c_int(self.ubound), self.ptr)
            #M.MUPHY.XinitMuArray_float(c_int(self.lbound), c_int(self.ubound), self.ptr)
            # M.MUPHY.initMuArray_float(c_int(0), c_int(1), self.ptr)
            # M.MUPHY.initMuArray_double(c_int(lb), c_int(ub), self.ptr) # TBD !!!

    def pointMuArray(self,parr):
        M.MUPHY.pointMuArray(self.ptr, parr)

    def __add__(self,other):
        new = MuArray(self.mesh)
        M.MUPHY.addMuArray_float(self.ptr, other.ptr, new.ptr)
        return new

    def __mul__(self,f):
        new = MuArray(self.mesh)
        M.MUPHY.mulMuArray_float(self.ptr, self.c_PRC(f), new.ptr)
        return new

    def __rmul__(self, f):
        return self.__mul__(f)

    def __getitem__(self,key):
        if isinstance(key, int) :
            if key < 0: key += len(self)
            if key < 0 or key >= len(self): 
                raise IndexError() # , "Index (%d) out of range" % key
            v = (self.c_PRC*1)()
            M.MUPHY.getElementMuArray_float(self.ptr, c_int(key), byref(v))
            return v[0]

        elif isinstance(key, slice) :
            v = (self.c_PRC*1)()
            vv = []
            for index in range(*key.indices(len(self))):
                M.MUPHY.getElementMuArray_float(self.ptr, c_int(index), byref(v))
                vv.append(v[0])
            return vv

    def __setitem__(self,key,v):
        if isinstance(key, int) :
            if key < 0: key += len(self)
            if key < 0 or key >= len(self): 
                raise IndexError() # , "Index (%d) out of range" % key
            M.MUPHY.setElementMuArray_float(self.ptr, c_int(key), self.c_PRC(v))

        elif isinstance(key, slice) :
            for index in range(*key.indices(len(self))):
                M.MUPHY.setElementMuArray_float(self.ptr, c_int(index), self.c_PRC(v))

    def __len__(self):
        sz = (c_int*1)()
        M.MUPHY.getSizeMuArray(self.ptr, byref(sz))
        return sz[0]

    def __del__(self):
        M.MUPHY.delMuArray(self.ptr)

################################

class Universe(object):
    """
    class for universe
    """

    DEFAULT = dict(NAME = 'Eden',
                   NUMBEROFSTEPS = 20,
                   TEMPERATURE = 0.,
                   RESTARTFLAG = False,
                   DUMPFREQUENCY = -1,
                   RESOLUTION = 1.0,
                   PHYSICALMASSDENSITY = 1.0,
                   PHYSICALVISCOSITY = 1.e-6,
                   CHARLENGTH = 1.0,
                   CHARVELOCITY = 1.0,
                   CHARPRESSURE = 1.0,
                   CHARMACH = 1.0,
                   HEARTBEAT = 1,
                   PERIOD = 1)

    def setMergeMAP(self,frequency,inputs,output):
        self.mergeMAP.active = True
        self.mergeMAP.frequency = frequency
        self.mergeMAP.inputs = inputs
        self.mergeMAP.output = output

    def setMergeVTK(self,frequency,influids,outputprefix):
        self.mergeVTK.active = True
        self.mergeVTK.frequency = frequency
        self.mergeVTK.influids = influids
        self.mergeVTK.outputprefix = outputprefix

    def __init__(self,unitsystem='MKS'):
        """
        self-explain
        -> M.MUPHY.getReferenceUniverse()
        """

        self.mergeMAPthread = None
        self.mergeMAP = Empty()
        self.mergeMAP.active = False

        self.mergeVTKthread = None
        self.mergeVTK = Empty()
        self.mergeVTK.active = False

        self.connectAllMeshesDone = False
        
        if unitsystem != 'MKS':
            print('Only MKS system allowed')
            sys.exit(1)

        self.unitsystem = unitsystem
        self.numberofsteps = None

        self.Lchar = 1.
        self.Dx = 1.
        self.Vchar = 1.
        self.Pchar = 1.
        self.PhysicalMassDensity = 1.
        self.PhysicalViscosity = 1.
        self.Resolution = 1.
        self.Mach = 1.
        self.cs = 1./math.sqrt(3.0)

        self.scalePath = []

        self.scaleList = []
        self.meshList = []
        self.fluidList = []
        self.atomList = []
        self.ODEList = []
        self.trackList = []
        self.crosstrackList = []

        # self.ref = (py_object*1)() # address of object
        # self.ref = (c_int*1)() # address of object
        self.ref = (c_void_p*1)() # address of object

        # return the MUPHY address and store it in self.ref
        M.MUPHY.getReferenceUniverse(byref(self.ref))

        self.setPeriod(1)
        
        #PM the following string identifies the algorithm used to treat linked scales
        #PM This is specified by the user by invoking self.linkScales()
        #PM As for now, the only method available is "MultiMesh" (default)
        self.linkedScalesMethod = ""  

    def addItems(self,items=[]): # add items to the universe
        """
        self-explanatory
        """
        for item in items:
            self.addItem(item)

    def addItem(self,item): # add items to the universe
        """
        self-explanatory
        """

        if isinstance(item,Mesh):
            self.meshList.append(item)

        elif isinstance(item,Scale):
            self.scaleList.append(item)

        elif isinstance(item,Fluid):
            self.fluidList.append(item)

        elif isinstance(item,Atom):
            self.atomList.append(item)

        elif isinstance(item,ODE):
            self.ODEList.append(item)

        elif isinstance(item,Tracker):

            if isinstance(item,CrossTracker):
                self.crosstrackList.append(item) # cross-scale tracker
            else:
                self.trackList.append(item) # in-scale tracker

        else:
            print('addItem: Item not recognized')

    def setPeriod(self,val): 
        """
        set period of universe clock
        -> M.MUPHY.setPeriodUniverse()
        """
        M.MUPHY.setPeriodUniverse(self.ref,c_int(val))

    def setHeartBeat(self,val): 
        """
        set heartbead of universe
        -> M.MUPHY.setHeartBeatUniverse()
        """
        M.MUPHY.setHeartBeatUniverse(self.ref, c_int(val))

    def getPeriod(self): 
        """
        get heartbead of universe
        -> M.MUPHY.getPeriodUniverse()
        """
        val = (c_int*1)()
        M.MUPHY.getPeriodUniverse(self.ref, byref(val))
        return int(val[0])

    def setTitle(self,name): 
        """
        set title of universe
        """
        self.setName(name)
    def setName(self,name=DEFAULT['NAME']): 
        """
        set name of universe
        -> M.MUPHY.setName()
        """
        self.name = name
        M.MUPHY.setName(self.ref, c_int(len(name)), name)
    def setNumberOfSteps(self,val=DEFAULT['NUMBEROFSTEPS']): 
        """
        self-explain
        -> M.MUPHY.setNumberOfSteps()
        """
        self.numberofsteps = val
        M.MUPHY.setNumberOfSteps(self.ref, c_int(val))
    def getNumberOfSteps(self):
        """
        self-explain
        -> M.MUPHY.getNumberOfSteps()
        """
        return M.MUPHY.getNumberOfSteps(self.ref)
    def setTemperature(self,val=DEFAULT['TEMPERATURE']): 
        """
        set universe temeprature. Can be overridden for specific actors
        -> M.MUPHY.setTemperature()
        """
        M.MUPHY.setTemperature(self.ref, c_float(val))
    def getTemperature(self): 
        """
        self-explain
        -> M.MUPHY.setTemperature()
        """
        return M.MUPHY.getTemperature(self.ref)
    def setStateRestart(self,val=DEFAULT['RESTARTFLAG'],resettimer=False): 
        """
        flag for restart
        -> M.MUPHY.setStateRestart()
        """
        M.MUPHY.setStateRestart(self.ref, c_bool(val), c_bool(resettimer))
    def getStateRestart(self):
        """
        flag for restart
        -> M.MUPHY.getStateRestart()
        """
        return M.MUPHY.getStateRestart(self.ref)
    def setStateDumpFrequency(self,val=DEFAULT['DUMPFREQUENCY']): 
        """
        frequency for restart
        -> M.MUPHY.setStateDumpFrequency()
        """
        M.MUPHY.setStateDumpFrequency(self.ref, c_int(val))
    def getStateDumpFrequency(self):
        """
        frequency for restart
        -> M.MUPHY.getStateDumpFrequency()
        """
        return M.MUPHY.getStateDumpFrequency(self.ref)

    def setLchar(self,val=DEFAULT['CHARLENGTH']):
        """
        set char length in m
        """
        self.Lchar = val

    def setResolution(self,voxm=None,voxcm=None,voxmm=None,voxLchar=None):
        """
        resolution in voxel per m, cm, mm or per Lchar
        """

        num_args = sum(0 if arg is None else 1 for arg in [voxm, voxcm, voxmm, voxLchar])

        # if voxm + voxcm + voxmm + voxLchar > 1:
        if num_args > 1:
            print('multiple conflicting arguments in setResolution')
            sys.exit(1)
        # elif voxm + voxcm + voxmm + voxLchar == 0:
        elif num_args == 0:
            print('no argument in setResolution')
            sys.exit(1)

        if voxm != None:

            self.Resolution = voxm

        elif voxcm != None:

            self.Resolution = 100 * voxcm

        elif voxmm != None:

            self.Resolution = 1000 * voxmm

        elif voxLchar != None:

            if self.Lchar == None:
                print('set characteristic length first !')
                sys.exit(1)

            self.Resolution = voxLchar / self.Lchar 

        else:
            print('specify resolution in voxm,voxcm,voxmm or voxLchar !')
            sys.exit(1)

        self.Dx = 1. / self.Resolution # in MKS

    def setPhysicalMassDensity(self,val=DEFAULT['PHYSICALMASSDENSITY']):
        """
        physical density in Kg/m3
        """
        self.PhysicalMassDensity = val

    def getPhysicalMassDensity(self):

        return self.PhysicalMassDensity

    def setPhysicalViscosity(self,val=DEFAULT['PHYSICALVISCOSITY']):
        """
        set viscosity in m2/s
        """
        self.PhysicalViscosity = val

    def getSpacing(self):
        """
        get spatial resolution
        """
        return self.Dx

    def getTimeStep(self):
        """
        get time step
        """
        return self.Dt

    def getLatticeViscosity(self):
        """
        get lattice velocity
        """
        return self.LatticeViscosity

    def setCharacteristicVelocity(self,val=DEFAULT['CHARVELOCITY']):
        """
        set char velocity in m/s
        """
        self.Vchar = val

    def setCharacteristicPressure_mmHg(self,val=DEFAULT['CHARPRESSURE'] * 0.007500616827042):
        """
        set char pressure in mmHg
        """
        self.Pchar_mmHg = val
        self.Pchar = val / 0.007500616827042

    def setCharacteristicPressure(self,val=DEFAULT['CHARPRESSURE']):
        """
        set char pressure in Pa
        """
        self.Pchar = val
        self.Pchar_mmHg = val * 0.007500616827042

    def getCharacteristicPressure_mmHg(self):
        """
        get char pressure in mmHg
        """
        return self.Pchar_mmHg

    def getCharacteristicPressure(self):
        """
        get char pressure in mmHg
        """
        return self.Pchar

    def setCharacteristicMach(self,val=DEFAULT['CHARMACH']):
        """
        char Mach
        """
        self.Mach = val

    def getCharacteristicMach(self):
        return self.Mach

    def setRandomRepeatableSequence(self,val):
        M.MUPHY.setRandomRepeatableSequence(c_bool(val))

    def setRandomMarsagliaZhang(self,val):
        M.MUPHY.setRandomMarsagliaZhang(c_bool(val))

    def create(self): # init the whole universe
        """
        create the universe
        -> M.MUPHY.banner()
        -> M.MUPHY.general_output()
        -> M.MUPHY.allocateMesh()
        -> M.MUPHY.allocateFluid()
        -> M.MUPHY.allocateAtom()
        -> M.MUPHY.allocateODE()
        -> M.MUPHY.allocateTracker()
        -> M.MUPHY.prepareDump()
        """

        self.setHeartBeat(0)

        # multiscale / multimesh case
        self.nscale = len(self.scaleList)
        # allocate scalelist within MUPHY (not scalepath)
        M.MUPHY.allocateScale(c_int(self.nscale))
        
        if self.nscale == 1 and len(self.scalePath) == 0: 
          a=self.scaleList[0]
          #get MUPHY scale reference in scaleList
          a.reference(1)
          a.id = 1
          a.setName('SCALE1')

          self.linkScales([self.scaleList[0]],method="")

        else:
          
          id=0
          for a in self.scaleList:
              id+=1
              a.reference(id)
              a.id = id
              a.setName('SCALE'+str(id))
          #>>>> the user MUST call linkScales to complete scales linkage
        

        # write out banner
        M.MUPHY.banner()

        self.Dt = 1.0
        self.LatticeViscosity = 0.0

        # compute and write out units for fluids
        if len(self.fluidList) > 0:

            self.Dt = self.Dx * self.cs * self.Mach / self.Vchar
            self.LatticeViscosity = self.PhysicalViscosity * self.Dt / self.Dx**2

            myid = self.getMyproc()

            if myid == 0:

                print
                print( '============================================================================================ ')
                print
                print( '                              CFD Characteristic Quantities')
                print('Characteristic Length       %10.4g [m]    ' % self.Lchar)
                print('Voxelization Density        %10g [Vox/Lch]' % (self.Resolution * self.Lchar),)
                print(                             '%10g [Vox/m]  ' % self.Resolution,)
                print(                             '%10g [Vox/cm] ' % (self.Resolution * 0.01))
                print('Spatial Resolution (Dx)     %10.4g [m]    ' % self.Dx)
                print('Physical Density            %10.4g [Kg/m3]' % self.PhysicalMassDensity)
                print('Characteristic Velocity     %10.4g [m/s]  ' % self.Vchar)
                print('Physical Viscosity          %10.4g [m2/s] ' % self.PhysicalViscosity)
                print('Simulated Mach              %10.4g        ' % self.Mach)
                # print('Characteristic Mach       %10.4g        ' % self.Mach)
                # print('Lattice Sound Speed       %10.4g        ' % self.cs)
                print
                print( '                                  Derived Quantities')
                print('Timestep                    %10.4g [s]    ' % self.Dt)
                if self.numberofsteps != None:
                    print('Total Simulation Time       %10.4g [s]    ' % (self.Dt * self.numberofsteps))
                print('Lattice Viscosity           %10.4g [LU]   ' % self.LatticeViscosity)
                print('Lattice Char. Velocity      %10.4g [LU]  ' % ((self.Dt / self.Dx) * self.Vchar))
                print
                print( '============================================================================================ ')
                print

            self.Mchar = self.PhysicalMassDensity * self.Dx**3

            M.MUPHY.setCharacteristicValues(self.ref, \
                                        c_float(self.Lchar), \
                                        c_float(self.Vchar), \
                                        c_float(self.Pchar), \
                                        c_float(self.Mchar),  \
                                        c_float(self.Mach))

        M.MUPHY.general_output()

        # treat meshes
        self.nmesh = len(self.meshList)
        M.MUPHY.allocateMesh(c_int(self.nmesh))

        id=0
        for a in self.meshList:
            id+=1
            a.reference(id)
            a.id = id
            a.setName('MESH'+str(id)) # can be overridden by user
            a.setGridSpacing(self.Dx) 
            a.setTimestep(self.Dt) 

        # treat fluids
        self.nFluid = len(self.fluidList)
        M.MUPHY.allocateFluid(c_int(self.nFluid))

        id=0
        for a in self.fluidList:
            id+=1
            a.reference(id)
            a.id = id
            a.universe = self
            a.setName('FLUID'+str(id)) # can be overridden by user

        # treat atoms
        self.nAtom = len(self.atomList)
        M.MUPHY.allocateAtom(c_int(self.nAtom))

        id=0
        for a in self.atomList:
            id+=1
            a.reference(id)
            a.id = id
            a.universe = self
            a.setName('ATOM'+str(id)) # can be overridden by user

        # treat ODEs
        self.nODE = len(self.ODEList)
        M.MUPHY.allocateODE(c_int(self.nODE))

        id=0
        for a in self.ODEList:
            id+=1
            a.reference(id)
            a.id = id
            a.universe = self
            a.setName('ODE'+str(id)) # can be overridden by user

        # prepare to save the state of each actor
        M.MUPHY.prepareDump()

        # treat trackers and crosstrackers
        self.nTracker = len(self.trackList)
        self.nCrossTracker = len(self.crosstrackList)

        M.MUPHY.allocateTracker(c_int(self.nTracker + self.nCrossTracker))

        id=0
        for a in self.trackList:
            id+=1
            a.reference(id)
            a.id = id
            a.universe = self
            a.setName('TRK'+str(id)) # can be overridden by user

        for a in self.crosstrackList:
            id+=1
            a.reference(id)
            a.id = id
            a.universe = self
            a.setName('CrossTRK'+str(id)) # can be overridden by user

        print_resource_usage('Created_Universe')

    def decorate(self): 
        """
        decorate the whole universe
        """
        # each call will allocate the prope arrays (scale,mesh,tracks,fluidsall,atomsall,ODEall)
        # for a in self.scaleList: a.prepare()

        #find the coarsest mesh and store its ref
        maxsp=0
        for a in self.meshList: #PM

          a.prepare()

          if a.getSpacing()>maxsp:
            maxsp = a.getSpacing()
            self.CoarsestMesh = a
          
        if maxsp == 0: 
            print('>>> ERROR: undefined mesh spacing. Check Universe.decorate() call. Aborting...\n')
            sys.exit(1)

        if self.linkedScalesMethod == "MultiMesh" : #PM
            
          """
          for i in range(len(self.scalePath)-1):
            if self.scalePath[i].mesh.getSpacing()<self.scalePath[i+1].mesh.getSpacing():
              print('In MultiMesh method, scales in linkScales must be given in spacing descending order!\nAborting...')
              sys.exit(1)
          """
          #PM NOTE: the order in which scales are stored within scalePath is now IRRELEVANT
          self.connectAllMeshes()

        for a in self.fluidList: a.prepare()
        for a in self.atomList: a.prepare()
        for a in self.ODEList: a.prepare()
        for a in self.trackList: a.prepare()
        for a in self.crosstrackList: a.prepare()
        for a in self.scaleList: a.prepare()

        for a in self.meshList:
          a.finalprepare()

        M.MUPHY.gpu_init()

        M.MUPHY.resetTimer()

        M.MUPHY.startSimulationTime()

        print_resource_usage('Decorate')

    def animate(self): 
        """
        animate the whole universe for one step (i.e. the "major" step, the mesh's longest one in MG)
        """

        if self.mergeMAPthread != None:
            self.mergeMAPthread.join() # thread for graphics
        if self.mergeVTKthread != None:
            self.mergeVTKthread.join() # thread for graphics

        for scale in self.scalePath:
              scale.preupdateCommon()

        #PM questo loop ha senso se linkedScalesMethod=="" oppure linkedScalesMethod=="MultiMesh"
        #PM Per altri metodi di avanzamento bisogna implementarli opportunamente
        
        #advance all meshes until the coarsest one is fully saturated.
        #That will mark the end of the "major" step
        heartbeat=0
        while True:

            self.setHeartBeat(heartbeat)

            for scale in self.scalePath:
              scale.update(self.linkedScalesMethod)
              # scale.update_PM()
              
            if self.CoarsestMesh.isSaturated(): break
            heartbeat += 1

        for tt in self.crosstrackList:
            tt.update()

        if self.getMyproc() == 0 and M2G != None:

            if self.mergeMAP.active and self.itime % self.mergeMAP.frequency == 0:
                self.mergeMAPthread = threading.Thread(target=M2G.MergeMAPWorker, args=(self.itime, self.mergeMAP))
                self.mergeMAPthread.start()

            if self.mergeVTK.active and self.itime % self.mergeVTK.frequency == 0:

                self.mergeVTK.inputs = []
                print(self.mergeVTK.influids)
                for fl in self.mergeVTK.influids:

                    self.mergeVTK.inputs.append( 'DIRDATA_' + fl.name + '/VTK/T' + str(self.itime).zfill(10) + '.pvtu' )

                self.mergeVTK.output = self.mergeVTK.outputprefix + str(self.itime).zfill(10) + '.vtu'

                self.mergeVTKthread = threading.Thread(target=M2G.MergeVTKWorker, args=(self.itime, self.mergeVTK))
                self.mergeVTKthread.start()
        """
        for heartbeat in range(self.getPeriod()):

            self.setHeartBeat(heartbeat)

            for scale in self.scalePath:
                if scale.timechart[heartbeat] == 0: #do nothing
                    continue
                elif scale.timechart[heartbeat] == 1: #do c&s.
                    #In MG:  the first time will be done on saturated nodes only,
                    # while the 2nd time will be done on insaturated nodes (through inter-mesh saturation)
                    scale.update()
                elif scale.timechart[heartbeat] == 2 and scale.isMG == True: 
                    #MG ONLY: do complete c&s (i.e. saturated nodes c&s first, then saturation on the unsaturated ones)
                    scale.update()
                    scale.update()                    
                else: 
                    print('Timechart exception! ...aborting job')
                    sys.exit(1)

        """
                  
        """
        for heartbeat in range(self.getPeriod()):
            self.setHeartBeat(heartbeat)

            for scale in self.scalePath:

                if scale.timechart[heartbeat] == 0: 
                    scale.update2()
        """

    def getMyproc(self):
        """
        get my processor
        """
        myproc = (c_int*1)()
        M.MUPHY.getMyproc(byref(myproc))
        return myproc[0]

    def getNumprocs(self):
        """
        get my processor
        """
        numprocs = (c_int*1)()
        M.MUPHY.getNumprocs(byref(numprocs))
        return numprocs[0]

    def getItime(self):
        """
        get initial and final timestep
        """

        itime_start = (c_int*1)()
        ncycle = (c_int*1)()
        dummy = (c_int*1)()

        M.MUPHY.get_itime(byref(itime_start),byref(ncycle),byref(dummy))

        return itime_start[0], ncycle[0]
    def setItime(self,itime): 
        """
        set current time
        """
        M.MUPHY.set_itime(byref(c_int(itime)))

    def cycle(self):
        """
        time cycling
        """

        itime_start,ncycle = self.getItime()
        #WARNING in a MG context this makes that the number of complete c&s steps executed
        # on the COARSEST grid is equal to the steps set by setNumberOfSteps plus one.
        # E.g. if one sets u=Universe; u.setNumberOfSteps(0) then 1 complete step on the
        # COARSEST mesh will be done.
        for self.itime in range(itime_start, itime_start + ncycle+1):
            self.setItime(self.itime)
            yield self.itime

    # to be called after scale.set(...)
    # define "scales linkage": set and allocate (in MUPHY) ScalePath;
    # NOT ScaleList which is set and allocated within Universe.create()
    def linkScales(self,scalePath,method="MultiMesh"):
        if not hasattr(scalePath[0],'ref'):
            print('>>> ERROR: linkScales has to be called after Universe.create(). Aborting...\n')
            sys.exit(1)
        if hasattr(self,'scalepathref'):
            print('>>> WARNING! Scales already linked. linkScales call ignored...')
            return
        #self.checkMeshConnections = True
        """
        link scales
        -> M.MUPHY.allocateScale()
        -> M.MUPHY.allocateScalepath()
        -> M.MUPHY.getReferenceScalepath()
        -> M.MUPHY.linkScalesScalepath()
        """

        ##
        self.scalePath = scalePath
        #PM in MultiMesh context, sort scalePath in mesh spacing descending order
        #PM (useful in Universe.animate and connectAllMeshes)
        if method != "MultiMesh" and method != "":
            print('>>> ERROR: undefined mesh spacing. Check Universe.decorate() call. Aborting...\n')
            sys.exit(1)
            
            
        self.linkedScalesMethod = method
        """ #PM
        self.nscale = len(self.scaleList)
        M.MUPHY.allocateScale(c_int(self.nscale))
        id=0
        for a in self.scaleList:
            id+=1
            a.reference(id)
            a.id = id
            a.setName('SCALE'+str(id))
        """
        ##
        M.MUPHY.allocateScalepath()

        self.scalepathref = (c_void_p*1)() # address of object

        # return the MUPHY address and store it in self.scalepathref
        M.MUPHY.getReferenceScalepath(byref(self.scalepathref))

        # link the scales
        a = get_refs(self.scalePath)
        M.MUPHY.linkScalesScalepath(self.scalepathref, c_int(len(a)), byref(a))

    def flushout(self):
        """
        deallocate all actors
        -> M.MUPHY.deallocateMesh()
        -> M.MUPHY.deallocateFluid()
        -> M.MUPHY.deallocateAtom()
        -> M.MUPHY.deallocateODE()
        -> M.MUPHY.deallocateTracker()
        -> M.MUPHY.barrier()
        """
        M.MUPHY.deallocateMesh()
        M.MUPHY.deallocateFluid()
        M.MUPHY.deallocateAtom()
        M.MUPHY.deallocateODE()
        M.MUPHY.deallocateTracker()
        M.MUPHY.barrier()

    def connectTwoMeshes(self,s1,s2):
        # -> M.MUPHY.connectTwoMeshes(s1,s2)
        #s1.enableMG()
        #s2.enableMG()
        M.MUPHY.connectTwoMeshes(s1.ref,s2.ref)
        
        # do all mesh pairs connections.
        # Every mesh is connected to both the one having
        # the smallest spacing among those with higher resolution
        # and the one with the greatest spacing among those
        # with lower resolution.
        # Then checkConnections is called for each mesh after the connections setup
        # to check actual status of each node and communicate various infos

    def connectAllMeshes(self): #PM

        if self.connectAllMeshesDone: return
      
        orderedMeshes={}
        for s in self.scalePath:
          orderedMeshes[s.mesh.getSpacing()]=s
        spacings=orderedMeshes.keys()
        spacings.sort()
        
        #DO NOT change connection order! It was carefully arranged!
        # note that self-connections is done last for each mesh
        for i in range(0,len(spacings)-1):
          self.connectTwoMeshes(orderedMeshes[spacings[i]],orderedMeshes[spacings[i+1]])
          self.connectTwoMeshes(orderedMeshes[spacings[i+1]],orderedMeshes[spacings[i]])
          #the following does mesh self-connection if enabled in fortran code, otherwise does nothing
          self.connectTwoMeshes(orderedMeshes[spacings[i]],orderedMeshes[spacings[i]])
          
        #the following does mesh self-connection if enabled in fortran code, otherwise does nothing
        self.connectTwoMeshes(orderedMeshes[spacings[-1]],orderedMeshes[spacings[-1]])
        for i in range(len(self.scalePath)):
          M.MUPHY.checkConnections(self.scalePath[i].ref)
            
        self.connectAllMeshesDone = True

class Scale(object):

    def __init__(self):
        self.actors = []
        self.timechart = []

    def reference(self,id):
        """
        make a reference
        -> M.MUPHY.getReferenceScale()
        """
        # self.ref = (py_object*1)() # address of object
        # self.ref = (c_int*1)() # address of object
        self.ref = (c_void_p*1)() # address of object

        M.MUPHY.getReferenceScale(c_int(id), \
                                byref(self.ref))
    def prepare(self): 
        """
        prepare scale
        -> M.MUPHY.prepareScale()
        """
        M.MUPHY.prepareScale(self.ref)

    def set(self,name='GenScale',mesh=None,timechart=[1],actors=[],temperature=0.):
        """
        set the scale with mesh, timechart, actors and temperature
        """

        self.setName(name)
        self.setMesh(mesh)
        self.setTimeChart(timechart)
        self.addActors(actors)
        self.setTemperature(temperature)
        #self.isMG=False #PM

        for a in actors: # back reference from the actors
            a.scale = self


    def setName(self,name): 
        """
        set name of scale 
        -> M.MUPHY.setNameScale()
        """
        self.name = name
        M.MUPHY.setNameScale(self.ref, c_int(len(name)), name)

    def setMesh(self,mesh): 
        """
        set mesh of scale 
        -> M.MUPHY.setMeshScale()
        """
        self.mesh = mesh
        M.MUPHY.setMeshScale(self.ref, self.mesh.ref)

    def setTimeChart(self,timechart): 
        """
        set timechart of scale 
        -> M.MUPHY.setTimeChartScale()
        """
        self.timechart = timechart
        tc = (c_int*len(timechart))()
        for i in range(len(timechart)): tc[i] = timechart[i]
        M.MUPHY.setTimeChartScale(self.ref, c_int(len(timechart)), byref(tc))

    def setTemperature(self,val): 
        """
        set temperature of scale 
        -> M.MUPHY.setTemperatureScale()
        """
        M.MUPHY.setTemperatureScale(self.ref, c_float(val))

    def getName(self): 
        """
        get name of scale
        """
        return self.name

    def addActors(self,actors): # add actors to scale. scale.mesh should be already referenced
        """
        add actors to scale
        -> M.MUPHY.addFluidsScale()
        -> M.MUPHY.addAtomsScale()
        -> M.MUPHY.addODEScale()
        -> M.MUPHY.addTrackersScale()
        """

        # all actors
        self.fluids = []
        self.atoms = []
        self.ODEs = []
        trks = []

        for a in actors:

            if isinstance(a,Fluid): 

                self.fluids.append(a)

            elif isinstance(a,Atom):

                self.atoms.append(a)

            elif isinstance(a,ODE): 

                self.ODEs.append(a)

            elif isinstance(a,Tracker): 

                trks.append(a)

            else:
                print('error for actors')
                sys.exit(1)

        self.tracker = None
        if len(trks) == 1:
            self.tracker = trks[0]
            self.tracker.addActorsTracker(actors)

        elif len(trks) > 1:
            print('error: a single scale cannot have more than one tracker')
            sys.exit(1)

        f = get_refs(self.fluids)
        M.MUPHY.addFluidsScale(self.ref, c_int(len(f)), byref(f))

        a = get_refs(self.atoms)
        M.MUPHY.addAtomsScale(self.ref, c_int(len(a)), byref(a))

        o = get_refs(self.ODEs)
        M.MUPHY.addODEsScale(self.ref, c_int(len(o)), byref(o))

        t = get_refs([self.tracker])
        M.MUPHY.addTrackersScale(self.ref, c_int(len(t)), byref(t))

    # update common to all actors in scale
    def preupdateCommon(self):
        """
        update common to all actors in scale
        """

        if len(self.fluids)>0: 
            self.fluids[0].preupdateCommon(self)

        if len(self.atoms)>0: 
            self.atoms[0].preupdateCommon(self)

        if self.tracker and self.tracker != None: 
            self.tracker.preupdateCommon(self)

    # update all scale actors
    def update(self,linkedScalesMethod):
        """
        update scale
        """

        # no MG: returns always true 
        #    MG: exit if this mesh hasn't to be updated 

        updated = []
        if linkedScalesMethod == "MultiMesh":
          #PM check is this scale's fluids can be updated or not
          #   It can always be updated in a single-mesh context
          if self.isToBeUpdated():
            M.MUPHY.updateFluids(self.ref)
            updated.append(True)
          else:
            updated.append(False)

        elif linkedScalesMethod == "":
          M.MUPHY.updateFluids(self.ref)
          updated.append(True)

        """
        for a in self.fluids: 
            a.update(self) 
        """
        for a in self.atoms: 
            updated.append ( a.update(self) ) 

        for a in self.ODEs: 
            a.update(self)

        if not all(updated): return

        self.tracker.update(self)


    def printIt(self): 
        """
        self-explain
        -> M.MUPHY.printItScale()
        """
        M.MUPHY.printItScale(self.ref)
        
    def isToBeUpdated(self):
       u = (c_int*1)()
       M.MUPHY.isToBeUpdatedMesh(self.ref,u)
       return u[0]==1
       
      

class Mesh(object):
    """
    class for mesh
    """

    DEFAULT = dict(NAME = 'Cubic',
                   PARALLELEPIPEDBOX = [10,10,10],
                   FILEROOT = 'bgkflag',
                   FILES = ['bgkflag.hdr','bgkflag.dat','bgkflag.ios'],
                   MESHFREE = False,
                   REGULARMESH = False,
                   READPREPROCESSEDMESH = False,
                   WRITEPREPROCESSEDMESH = False,
                   WORKASPREPROCESSOR = False,
                   PERIODICITYSTRING = '111',
                   PERIODICITY = [True,True,True],
                   SYSTEMEXCHANGE = None,
                   GRIDSPACING = 1.0,
                   TIMESTEP = 1.0,
                   DOMAINDECOMPOSITION = ['Z_Slabs',
                                          'Y_Slabs',
                                          'X_Slabs',
                                          'ZY_Slabs',
                                          'ZX_Slabs',
                                          'XY_Slabs',
                                          'Cubes',
                                          'Read_in',
                                          'Broadcasted'],
                   PARTITIONALONGXYZ = [2,2,2],
                   SYNCHRONOUSSEND = False)

    dd2int = dict(Z_Slabs = 1, 
                  Y_Slabs = 2,
                  X_Slabs = 3,
                  ZY_Slabs = 4,
                  ZX_Slabs = 5,
                  YX_Slabs = 6,
                  Cubes = 7,
                  Read_in = 8,
                  Broadcasted = 9)

    def __init__(self):
        self.boxsetexplicit=False
        pass

    def reference(self,id):
        """
        self-explain
        -> M.MUPHY.getReferenceMesh()
        """
        # self.ref = (py_object*1)() # address of object
        # self.ref = (c_int*1)() # address of object
        self.ref = (c_void_p*1)() # address of object

        M.MUPHY.getReferenceMesh(c_int(id), \
                                byref(self.ref))

    def setName(self,name=DEFAULT['NAME']): 
        """
        self-explain
        -> M.MUPHY.setNameMesh()
        """
        self.name = name
        M.MUPHY.setNameMesh(self.ref, c_int(len(name)), name)

    def createWorkArray(self):
        arr = MuArray() # allocate ctypes pointer
        M.MUPHY.createWorkArrayMesh(self.ref, 
                                byref(arr.ptr), 
                                arr.size, arr.lbound, arr.ubound)
        arr.size = arr.size[0]
        arr.lbound = arr.lbound[0]
        arr.ubound = arr.ubound[0]
        if arr.lbound==None: arr.lbound=0
        if arr.ubound==None: arr.ubound=0
        return arr

    def destroyWorkArray(self,arr):
        M.MUPHY.destroyWorkArrayMesh(self.ref, arr.ptr)
        del arr

    def getLaplacian(self, arr, arrin):
        M.MUPHY.getLaplacianMesh(self.ref, 
                                 byref(arr.ptr), 
                                 c_int(arr.size), c_int(arr.lbound), c_int(arr.ubound),
                                 byref(arrin.ptr))
        return arr

    def setParallelepipedBox(self,nx=DEFAULT['PARALLELEPIPEDBOX'][0], \
                                  ny=DEFAULT['PARALLELEPIPEDBOX'][1], \
                                  nz=DEFAULT['PARALLELEPIPEDBOX'][2]) :
        """
        set box enclosing mesh
        -> M.MUPHY.getMyproc()
        """
        # myid = M.MUPHY.getMyproc()
        # myid = universe.getMyproc()  # universe not defined...
        myid_l = (c_int*1)()
        M.MUPHY.getMyproc(byref(myid_l))
        myid = myid_l[0]

        prefix = '_box.'+str(self.id)

        if myid==0:

            f = open(prefix+'.hdr','w')
            print >> f, int(nx),int(ny),int(nz)
            # f.write( int(nx),int(ny),int(nz) )
            print >> f, 1,int(nx)*int(ny)*int(nz),0,0,0
            # f.write( 1,int(nx)*int(ny)*int(nz),0,0,0 )
            #f.write('{0:d} {1:d} {2:d}\n'.format(int(nx),int(ny),int(nz)))
            #f.write('{0:d} {1:d} 0 0 0\n'.format(1,int(nx)*int(ny)*int(nz)))
            print >>f, 1
            f.close()

            f = open(prefix+'.dat','w')
            # print >> f,'-1 -1 -1 -1'
            f.write( '-1 -1 -1 -1' )
            #f.write('-1 -1 -1 -1\n')
            f.close()

        self.setFiles(prefix+'.hdr',prefix+'.dat')

        self.boxsetexplicit=True

    def setFileroot(self,nameroot=DEFAULT['FILEROOT']): 
        """
        set file root of mesh
        -> M.MUPHY.setFileroot()
        """
        if self.boxsetexplicit:
            print('conflicting options....box already set explicitly')
            return
        M.MUPHY.setFileroot(self.ref, \
                            c_int(len(nameroot)), nameroot)

    def getArray(self,n): 

        array = (c_float*n)()
        return array

    def setFiles(self,namehdr=DEFAULT['FILES'][0], \
                      namedat=DEFAULT['FILES'][1], \
                      nameios=DEFAULT['FILES'][2]): 
        """
        set file names for mesh
        -> M.MUPHY.setFilesMesh()
        """
        if self.boxsetexplicit:
            print('conflicting options....box already set explicitly')
            return
        M.MUPHY.setFilesMesh(self.ref, \
                            c_int(len(namehdr)), namehdr, \
                            c_int(len(namedat)), namedat, \
                            c_int(len(nameios)), nameios)


    def setMeshFree(self,val=DEFAULT['MESHFREE']): 
        """
        specify that the mesh is free: basically an empty box for MD
        -> M.MUPHY.setMeshFreeMesh()
        """
        M.MUPHY.setMeshFreeMesh(self.ref, c_bool(val))
    def setRegularMesh(self,val=DEFAULT['REGULARMESH']): 
        """
        specify that the mesh is a regular one
        -> M.MUPHY.setRegularMesh()
        """
        M.MUPHY.setRegularMesh(self.ref, c_bool(val))
    def setReadPreprocessedMesh(self,val=DEFAULT['READPREPROCESSEDMESH']): 
        """
        specify that the mesh is read in as preprocessed
        -> M.MUPHY.setReadPreprocessedMesh()
        """
        M.MUPHY.setReadPreprocessedMesh(self.ref, c_bool(val))
    def setWritePreprocessedMesh(self,val=DEFAULT['WRITEPREPROCESSEDMESH']):
        """
        specify that the mesh is written as preprocessed
        -> M.MUPHY.setWritePreprocessedMesh()
        """
        M.MUPHY.setWritePreprocessedMesh(self.ref, c_bool(val))
    def setWorkAsPreprocessor(self,val=DEFAULT['WORKASPREPROCESSOR']): 
        """
        specify that the run acts as the preprocessor
        -> M.MUPHY.setWorkAsPreprocessor()
        """
        M.MUPHY.setWorkAsPreprocessor(self.ref, c_bool(val))
    def setPeriodicity(self,name=DEFAULT['PERIODICITYSTRING']):
        """
        set periodicity
        -> M.MUPHY.setPeriodicityMesh()
        """
        M.MUPHY.setPeriodicityMesh(self.ref, c_int(len(name)), name)
    def setSystemExchange(self,name=DEFAULT['SYSTEMEXCHANGE']):
        """
        set open/close mesh
        -> M.MUPHY.setSystemExchangeMesh()
        """
        print('(py) Inactive....call instead setInletOutletMethod')
        M.MUPHY.setSystemExchangeMesh(self.ref, c_int(len(name)), name)
    def setGridSpacing(self,val=DEFAULT['GRIDSPACING']):
        """
        set grid spacing
        -> M.MUPHY.setGridSpacingMesh()
        """
        M.MUPHY.setGridSpacingMesh(self.ref, c_float(val))
    def setTimestep(self,val=DEFAULT['TIMESTEP']):
        """
        set grid spacing
        -> M.MUPHY.setTimestepMesh()
        """
        M.MUPHY.setTimestepMesh(self.ref, c_float(val))
    def setDomainDecomposition(self,val=DEFAULT['DOMAINDECOMPOSITION']):
        """
        set domain decomposition
        -> M.MUPHY.setDomainDecompositionMesh()
        """
        if isinstance(val, int):
            ival = val
        else:
            ival = dd2int[val]

        M.MUPHY.setDomainDecompositionMesh(self.ref, c_int(ival))
    def setPartitionAlongXYZ(self,
                             v1=DEFAULT['PARTITIONALONGXYZ'][0],
                             v2=DEFAULT['PARTITIONALONGXYZ'][1],
                             v3=DEFAULT['PARTITIONALONGXYZ'][2]):
        """
        set decomposition along X Y and Z
        -> M.MUPHY.setPartitionAlongXYZMesh()
        """
        M.MUPHY.setPartitionAlongXYZMesh(self.ref, c_int(v1),c_int(v2),c_int(v3))
    def setSynchronousSend(self,val=DEFAULT['SYNCHRONOUSSEND']):
        """
        set sync vs async communications
        -> M.MUPHY.setSynchronousSendMesh()
        """
        M.MUPHY.setSynchronousSendMesh(self.ref, c_bool(val))

    def setPoints(self,fluidnodes,inletnodes,outletnodes): 
        """
        set mesh points
        -> M.MUPHY.setPointsMesh()
        """
        ntot = len(fluidnodes + inletnodes + outletnodes)

        nwguess = (c_int*1)()
        M.MUPHY.setAllocatePointsMesh(self.ref, c_int(ntot), byref(nwguess))
        ntot = ntot + nwguess
        ipnts_l = (c_int*n)(*ipnts); jpnts_l = (c_int*n)(*jpnts); kpnts_l = (c_int*n)(*kpnts)
        M.MUPHY.setWallPointsMesh(self.ref, c_int(ntot), byref(ipnts_l), byref(jpnts_l), byref(kpnts_l))

        M.MUPHY.setFluidPointsMesh(self.ref, c_int(n), byref(ipnts_l), byref(jpnts_l), byref(kpnts_l))

    def setFrame(self,nx,ny,nz): 
        """
        set box frame and strides
        -> M.MUPHY.setFrameMesh()
        -> M.MUPHY.setBoxMesh()
        """
        M.MUPHY.setFrameMesh(self.ref, c_int(nx), c_int(ny), c_int(nz))
        bx = 1.*nx
        by = 1.*ny
        bz = 1.*nz
        M.MUPHY.setBoxMesh(self.ref, c_float(bx), c_float(by), c_float(bz))

    def setBox(self,bx,by,bz): 
        """
        set box
        -> M.MUPHY.setBoxMesh()
        """
        M.MUPHY.setBoxMesh(self.ref, c_float(bx), c_float(by), c_float(bz))
    def getBox(self): 
        """
        get box
        -> M.MUPHY.getBoxMesh()
        """
        bx = (c_float*1)()
        by = (c_float*1)()
        bz = (c_float*1)()
        M.MUPHY.getBoxMesh(self.ref, byref(bx), byref(by), byref(bz))
        return [bx[0],by[0],bz[0]]

    def getSpacing(self): #PM return this mesh spacing 
        c = (c_int*1)()
        M.MUPHY.getSpacingMesh(self.ref,byref(c))
        return c[0]

    def isSaturated(self): #PM return True if current mesh is fully saturated 
        c = (c_int*1)()
        M.MUPHY.IsSaturatedMesh(self.ref,byref(c))
        return c[0]==1

    def prepare(self):
        """
        prepare mesh
        -> M.MUPHY.prepareMesh(byref()
        -> M.MUPHY.prepareMeshConnector()
        """
        id = 1
        M.MUPHY.prepareMesh(byref(self.ref))

        M.MUPHY.prepareMeshConnector(byref(self.ref))

    def finalprepare(self):
        """
        prepare mesh
        -> M.MUPHY.finalprepareMesh(byref()
        """
        M.MUPHY.finalprepareMesh(byref(self.ref))

    def getFluidNodesBounds(self): 

        lb = (c_int*1)()
        ub = (c_int*1)()
        M.MUPHY.getFluidNodesBoundsMesh(self.ref, byref(lb), byref(ub))

        return lb[0],ub[0]

    def getLocator(self,i,j,k): 

        ifl = (c_int*1)()
        M.MUPHY.getLocatorMesh(self.ref, c_int(i), c_int(j), c_int(k), byref(ifl))

        return ifl[0]

    def getNodePosition(self,ifl): 

        i = (c_int*1)()
        j = (c_int*1)()
        k = (c_int*1)()
        M.MUPHY.getNodePositionMesh(self.ref, c_int(ifl), byref(i), byref(j), byref(k))

        return i[0],j[0],k[0]

    def getNodeType(self,ifl): 

        itype = (c_int*1)()
        M.MUPHY.getNodeTypeMesh(self.ref, c_int(ifl), byref(itype))

        return itype[0]

    def getLocatorNodeType(self,i,j,k): 

        itype = (c_int*1)()
        M.MUPHY.getLocatorNodeTypeMesh(self.ref, c_int(i), c_int(j), c_int(k), byref(itype))

        if itype[0] == 1:

            return 'FLUID'

        elif itype[0] == 2:

            return 'WALL'

        elif itype[0] == 3:

            return 'INLET'

        elif itype[0] == 4:

            return 'OUTLET'

        else:
            print('getLocatorNodeType: node type not recognized:',itype[0])
            sys.exit(1)

class Actor(object):
    scale = None
    mesh = None
    pass


class Fluid(Actor):
    """
    class for fluid
    """

    DEFAULT = dict(NAME = 'Air',
                   DUMPINFORMAT = 'unformatted',
                   DUMPOUTFORMAT = 'unformatted',
                   ADR = False,
                   MOMENTUMFREEZE = False,
                   FREEZE = False,
                   INERT = False,
                   STABILIZEDLB = True,
                   STABILIZEDLBVMAX = 0.9,
                   REGULARIZEDLB = False,
                   ENTROPICLB = False,
                   ENTROPICLBITERATION = 20,
                   ENTROPICLBTOLERANCE = 1.e-6,
                   COLLISIONTYPE = ['BGK','HSMIX'],
                   VISCOSITY = 1./6.,
                   DENSITY = 0.22,
                   MASS = 1.0,
                   DIFFUSIVITY = 1./6.,
                   CHARGE = 1.,
                   HSDIAMETER = 0.,
                   FLUIDONPARTICLE = True,
                   ADDNOISE = False,
                   TEMPERATURE = 0.0,
                   HOMOGENEOUSFORCE = (0.,0.,0.),
                   HOMOGENEOUSELECTRICFORCE = (0.,0.,0.),
                   INITIALVELOCITY = (0.,0.,0.),
                   SHEAR = False,
                   SHEARDIRS = (1,2),
                   SHEARFORCING = 0.0,
                   PHYSICALVISCOSITY = 1.,
                   PHYSICALDENSITY = 1.,
                   PHYSICALMASSDENSITY = 1.,
                   STOPCOM = False,
                   AUXILIARYDISTRIBUTION = False,
                   EDMINTEGRATION = False,
                   TRAPEZOIDALINTEGRATION = False,
                   HSDECORRELATE = False,
                   CAPFORCES = False,
                   INLETOUTLETMETHOD = ['closed','equilibrium','zouhe'],
                   INLETOUTLETFILE = 'bgkflag.ios',
                   WBCFILE = 'bgkflag.wbc',
                   DENSITYRANDOMFACTOR = 1.0,
                   SHANCHEN = True,
                   SHANCHENWALL = True,
                   SHANCHENMIX = True,
                   WALLLEAKAGEFACTOR = 0.0)

    def __init__(self):
        """
        self-explain
        """
        self.outletmethod = None
        self.monitorInletOutletInited = None

    def reference(self,id):
        """
        a few references, such as gradient and laplacian
        -> M.MUPHY.getReferenceFluid()
        """
        # self.ref = (py_object*1)() # address of object
        # self.ref = (c_int*1)() # address of object
        self.ref = (c_void_p*1)() # address of object

        # self.gradient = (py_object*1)() # object function
        # self.gradient = (c_int*1)() # object function
        self.gradient = (c_void_p*1)() # object function

        # self.laplacian = (py_object*1)() # object function
        # self.laplacian = (c_int*1)() # object function
        self.laplacian = (c_void_p*1)() # object function

        M.MUPHY.getReferenceFluid(c_int(id), \
                                    byref(self.ref), \
                                    byref(self.gradient), \
                                    byref(self.laplacian))

    def setName(self,name=DEFAULT['NAME']): 
        """
        set name of fluid
        -> M.MUPHY.setNameFluid()
        """
        self.name = name
        M.MUPHY.setNameFluid(self.ref, c_int(len(name)), name)

    def setDumpInFormat(self,fmt=DEFAULT['DUMPINFORMAT']): 
        """
        set name of fluid
        -> M.MUPHY.setDumpInFormatFluid()
        """
        M.MUPHY.setDumpInFormatFluid(self.ref, c_int(len(fmt)), fmt)

    def setDumpOutFormat(self,fmt=DEFAULT['DUMPOUTFORMAT']): 
        """
        set name of fluid
        -> M.MUPHY.setDumpOutFormatFluid()
        """
        M.MUPHY.setDumpOutFormatFluid(self.ref, c_int(len(fmt)), fmt)

    def setAdvector(self,f): 
        """
        set freeze flag for fluid
        -> M.MUPHY.setAdvectorFluid()
        """
        M.MUPHY.setAdvectorFluid(self.ref, f.ref)

    def setADR(self,val=DEFAULT['ADR']): 
        """
        set ADR flag for fluid
        -> M.MUPHY.setADRFluid()
        """
        M.MUPHY.setADRFluid(self.ref, c_bool(val))

    def setMomentumFreeze(self,val=DEFAULT['MOMENTUMFREEZE']): 
        """
        set momentum freeze flag for fluid
        -> M.MUPHY.setMomentumFreeze()
        """
        M.MUPHY.setMomentumFreezeFluid(self.ref, c_bool(val))

    def setFreeze(self,val=DEFAULT['FREEZE']): 
        """
        set freeze flag for fluid
        -> M.MUPHY.setFreezeFluid()
        """
        M.MUPHY.setFreezeFluid(self.ref, c_bool(val))

    def setInert(self,val=DEFAULT['INERT']): 
        """
        set inert flag for fluid
        -> M.MUPHY.setInertFluid()
        """
        M.MUPHY.setInertFluid(self.ref, c_bool(val))

    def setStabilizeLB(self,val=DEFAULT['STABILIZEDLB'],vmax=DEFAULT['STABILIZEDLBVMAX']): 
        """
        set stabilization for fluid
        -> M.MUPHY.setStabilizeLBFluid()
        """
        M.MUPHY.setStabilizeLBFluid(self.ref, c_bool(val), c_float(vmax))
    def setRegularizedLB(self,val=DEFAULT['REGULARIZEDLB']): 
        """
        set regularization for fluid : no ghosts modes
        -> M.MUPHY.setRegularizedLBFluid()
        """
        M.MUPHY.setRegularizedLBFluid(self.ref, c_bool(val))
    def setEntropicLB(self,val=DEFAULT['ENTROPICLB'],itermax=DEFAULT['ENTROPICLBITERATION'],nrtol=DEFAULT['ENTROPICLBTOLERANCE']): 
        """
        set stabilization for fluid
        -> M.MUPHY.setEntropicLBFluid()
        """
        M.MUPHY.setEntropicLBFluid(self.ref, c_bool(val), c_int(itermax), c_float(nrtol))
    def setCollisionType(self,name=DEFAULT['COLLISIONTYPE'][0]): 
        """
        set collision type for fluid
        -> M.MUPHY.setCollisionTypeFluid()
        """
        M.MUPHY.setCollisionTypeFluid(self.ref, c_int(len(name)), name)
    def setLocalViscosity(self,ifl,val): 
        """
        set lattice viscosity for fluid
        -> M.MUPHY.setViscosityFluid()
        """
        if ifl<0: return
        M.MUPHY.setLocalViscosityFluid(self.ref, c_int(ifl), c_float(val))
    def setViscosity(self,val=DEFAULT['VISCOSITY']): 
        """
        set lattice viscosity for fluid
        -> M.MUPHY.setViscosityFluid()
        """
        M.MUPHY.setViscosityFluid(self.ref, c_float(val))

    def setDensityUniform(self,val=DEFAULT['DENSITY']): 
        """
        set initial density for fluid
        -> M.MUPHY.setDensityFluid()
        """
        M.MUPHY.setDensityFluid(self.ref, c_float(val))

    def getMeanDensity(self,val): 
        """
        get density for fluid
        -> M.MUPHY.getMeanDensityFluid()
        """
        val = (c_float*1)()
        M.MUPHY.getMeanDensityFluid(self.mesh.ref, self.ref, byref(val))

        return val[0]

    def getDensity(self):

        parr = (c_void_p*1)()
        M.MUPHY.getDensityFluid(self.ref, byref(parr))

        # muarray = MuArray(self.mesh, allocate=False)
        muarray = MuArray(self.mesh, allocate=True)

        M.MUPHY.pointMuArray(muarray.ptr, parr)

        # print 'DENSITY:', len(muarray), (muarray.lbound, muarray.ubound), muarray[:3]
        return muarray

    def setDensityTMP(self,parr):

        #print '\n\n',parr.ptr,'\n\n'
        M.MUPHY.setDensityTMPFluid(self.ref, parr.ptr)

    def setEquilibriumPopulation(self,prho,pu,pv,pw): 
        """
        set initial density for fluid
        -> M.MUPHY.setEquilibriumPopulationFluid()
        """
        M.MUPHY.setEquilibriumPopulationFluid(self.ref, prho.ptr, pu.ptr, pv.ptr, pw.ptr)

    def getArray(self,n): 
        array = (c_float*n)() # allocating n elements !
        return array

    def setDensityProfile(self,rho): 
        """
        set initial density for fluid
        -> M.MUPHY.setDensityProfileFluid()
        """
        n = len(rho)

        M.MUPHY.setDensityProfileFluid(self.ref, c_int(n), rho)

    def getDensityProfile(self,rho):  # MAKES A COPY OF THE ARRAY !!!
        """
        get instantaneous density of fluid
        -> M.MUPHY.getDensityProfileFluid()
        """
        n = len(rho)

        M.MUPHY.getDensityProfileFluid(self.ref, c_int(n), rho)

        return rho

    def readInitialProfile(self,name): 
        """
        read initial profile from file for fluid
        -> M.MUPHY.readInitialProfileFluid()
        """
        M.MUPHY.readInitialProfileFluid(self.ref, c_int(len(name)), name)

    def setShanChenWallCouplingProfile(self,couplingarray): 
        """
        set initial density for fluid
        -> M.MUPHY.setDensityProfileFluid()
        """
        n = len(couplingarray)

        M.MUPHY.setShanChenWallCouplingProfileFluid(self.ref, c_int(n), couplingarray)

    def setMass(self,val=DEFAULT['MASS']): 
        """
        set mass for fluid
        -> M.MUPHY.setMassFluid()
        """
        M.MUPHY.setMassFluid(self.ref, c_float(val))
    def setDiffusivity(self,val=DEFAULT['DIFFUSIVITY']): 
        """
        set mass diffusivity for fluid
        -> M.MUPHY.setDiffusivityFluid()
        """
        M.MUPHY.setDiffusivityFluid(self.ref, c_float(val))
    def setCharge(self,val=DEFAULT['CHARGE']): 
        """
        set charge for fluid
        -> M.MUPHY.setChargeFluid()
        """
        M.MUPHY.setChargeFluid(self.ref, c_float(val))
    def setHSDiameter(self,val=DEFAULT['HSDIAMETER']): 
        """
        set hard sphere diameter for fluid
        -> M.MUPHY.setHSDiameterFluid()
        """
        M.MUPHY.setHSDiameterFluid(self.ref, c_float(val))
    def setFluidOnParticle(self,val=DEFAULT['FLUIDONPARTICLE']): 
        """
        set if fluid acts on particle
        -> M.MUPHY.setFluidOnParticle()
        """
        M.MUPHY.setFluidOnParticle(self.ref, c_bool(val))
    def setAddNoise(self,val=DEFAULT['ADDNOISE']): 
        """
        set if fluid has noise
        -> M.MUPHY.setAddNoiseFluid()
        """
        M.MUPHY.setAddNoiseFluid(self.ref, c_bool(val))
    def setTemperature(self,val=DEFAULT['TEMPERATURE']): 
        """
        set temperature for fluid
        -> M.MUPHY.setTemperatureFluid()
        """
        M.MUPHY.setTemperatureFluid(self.ref, c_float(val))
    def setHomogeneousForce(self,v1=DEFAULT['HOMOGENEOUSFORCE'][0], 
                                 v2=DEFAULT['HOMOGENEOUSFORCE'][1], 
                                 v3=DEFAULT['HOMOGENEOUSFORCE'][2]) :
        """
        set homogeneous force for fluid
        -> M.MUPHY.setHomogeneousForceFluid()
        """
        M.MUPHY.setHomogeneousForceFluid(self.ref, c_float(v1), c_float(v2), c_float(v3))
    def setHomogeneousElectricForce(self,v1=DEFAULT['HOMOGENEOUSELECTRICFORCE'][0],
                                         v2=DEFAULT['HOMOGENEOUSELECTRICFORCE'][1],
                                         v3=DEFAULT['HOMOGENEOUSELECTRICFORCE'][2]):
        """
        set homogeneous electric force for fluid
        -> M.MUPHY.setHomogeneousElectricForceFluid()
        """
        M.MUPHY.setHomogeneousElectricForceFluid(self.ref, c_float(v1), c_float(v2), c_float(v3))
    def setInitialVelocity(self,v1=DEFAULT['INITIALVELOCITY'][0],
                                v2=DEFAULT['INITIALVELOCITY'][1],
                                v3=DEFAULT['INITIALVELOCITY'][2]): 
        """
        set initial velocity for fluid
        -> M.MUPHY.setInitialVelocityFluid()
        """
        M.MUPHY.setInitialVelocityFluid(self.ref, c_float(v1), c_float(v2), c_float(v3))
    def setShear(self,flag=DEFAULT['SHEAR'],dirs=DEFAULT['SHEARDIRS'],forcing=DEFAULT['SHEARFORCING']): 
        """
        set shear force for fluid
        -> M.MUPHY.setShearFluid()
        """
        M.MUPHY.setShearFluid(self.ref, c_bool(flag), c_int(2), dirs, c_float(forcing))
    def setPhysicalViscosity(self,val=DEFAULT['PHYSICALVISCOSITY']): 
        """
        set physical viscosity for fluid
        -> M.MUPHY.setPhysicalViscosityFluid()
        """
        M.MUPHY.setPhysicalViscosityFluid(self.ref, c_float(val))

    def setPhysicalDensity(self,val=DEFAULT['PHYSICALDENSITY']): 
        """
        set physical mass density for fluid
        """
        self.setPhysicalMassDensity(val)

    def setPhysicalMassDensity(self,val=DEFAULT['PHYSICALMASSDENSITY']): 
        """
        set physical mass density for fluid
        -> M.MUPHY.setPhysicalMassDensityFluid()
        """
        M.MUPHY.setPhysicalMassDensityFluid(self.ref, c_float(val))
    def setStopCOM(self,val=DEFAULT['STOPCOM']): 
        """
        set flag to stop center of mass for fluid
        -> M.MUPHY.setStopCOMFluid()
        """
        M.MUPHY.setStopCOMFluid(self.ref, c_bool(val))
    def setAuxiliaryDistribution(self,val=DEFAULT['AUXILIARYDISTRIBUTION']): 
        """
        set flag to use auxiliary distribution method for fluid
        -> M.MUPHY.setAuxiliaryDistributionFluid()
        """
        M.MUPHY.setAuxiliaryDistributionFluid(self.ref, c_bool(val))
    def setEDMIntegration(self,val=DEFAULT['EDMINTEGRATION']): 
        """
        set flag to use exact difference method for fluid
        -> M.MUPHY.setExactDiffMethIntegrationFluid()
        """
        M.MUPHY.setExactDiffMethIntegrationFluid(self.ref, c_bool(val))
    def setTrapezoidalIntegration(self,val=DEFAULT['TRAPEZOIDALINTEGRATION']): 
        """
        set flag to use trapezoidal method for fluid
        -> M.MUPHY.setTrapezoidalIntegrationFluid()
        """
        M.MUPHY.setTrapezoidalIntegrationFluid(self.ref, c_bool(val))
    def setHSDecorrelate(self,val=DEFAULT['HSDECORRELATE']): 
        """
        set flag to decorrelate equilibrium for hard spheres (ENSKOG), equals to ORTHOGONAL equilibrium 
        -> M.MUPHY.setHSDecorrelateFluid()
        """
        M.MUPHY.setHSDecorrelateFluid(self.ref, c_bool(val))
    def setCapVelocities(self,val=False,roof=1.): 
        """
        set params for capping forces for fluid
        -> M.MUPHY.setCapForcesFluid()
        """
        M.MUPHY.setCapVelocitiesFluid(self.ref, c_bool(val), c_float(roof))
    def setCapForces(self,val=DEFAULT['CAPFORCES'],roof=1.): 
        """
        set params for capping forces for fluid
        -> M.MUPHY.setCapForcesFluid()
        """
        M.MUPHY.setCapForcesFluid(self.ref, c_bool(val), c_float(roof))
    def setInletOutletMethod(self,name=DEFAULT['INLETOUTLETMETHOD'][1]): 
        """
        set inlet/outlet boundary method for fluid
        -> M.MUPHY.setInletOutletMethod()
        """
        M.MUPHY.setInletOutletMethodFluid(self.ref, c_int(len(name)), name)

    def setInletOutletFile(self,name=DEFAULT['INLETOUTLETFILE']): 
        """
        set inlet/outlet file for fluid
        -> M.MUPHY.setInletOutletFileFluid()
        """
        M.MUPHY.setInletOutletFileFluid(self.ref, c_int(len(name)), name)

    def setWBCFile(self,name=DEFAULT['WBCFILE']): 
        """
        set inlet/outlet file for fluid
        -> M.MUPHY.setWBCFileFluid()
        """
        M.MUPHY.setWBCFileFluid(self.ref, c_int(len(name)), name)

    def setInitialDensityRandomFactor(self,rdf=DEFAULT['DENSITYRANDOMFACTOR']): 
        """
        set randomization factor for initial density
        -> M.MUPHY.setInitialDensityRandomFactorFluid()
        """
        M.MUPHY.setInitialDensityRandomFactorFluid(self.ref, c_float(rdf))

    def setShanChen(self,flag=DEFAULT['SHANCHEN'],coupling=0.0,psi_type=1,a=1.0,b=4.0,t=1.0,omega=1.0): 
        """
        set flag for shan-chen method
        -> M.MUPHY.setShanChenFluid()
        """
        M.MUPHY.setShanChenFluid(self.ref, c_bool(flag), c_float(coupling), \
         c_int(psi_type), c_float(a), c_float(b), c_float(t), c_float(omega))

    def setShanChenWall(self,flag=DEFAULT['SHANCHENWALL'],wallcoupling=0.0,walldensity=1.0): 
        """
        set flag for shan-chen method
        -> M.MUPHY.setShanChenWallFluid()
        """
        M.MUPHY.setShanChenWallFluid(self.ref, c_bool(flag), c_float(wallcoupling), c_float(walldensity))

    def setShanChenMix(self,fluidfriend,flag=DEFAULT['SHANCHENMIX'],coupling=0.0,psiAB_type=3): 
        """
        set flag for shan-chen method
        -> M.MUPHY.setShanChenMixFluid()
        """
        M.MUPHY.setShanChenMixFluid(self.ref, fluidfriend.ref, c_bool(flag), c_float(coupling), c_int(psiAB_type))

    def resetShanChenMix(self,fluidfriend,flag,coupling=0.0,psiAB_type=3): 
        """
        set flag for shan-chen method
        -> M.MUPHY.resetShanChenMixFluid()
        """
        M.MUPHY.resetShanChenMixFluid(self.ref, fluidfriend.ref, c_bool(flag), c_float(coupling), c_int(psiAB_type))

    def getSizeArray(self):
        """
        get size of array for fluid
        -> M.MUPHY.getSizeArrayFluid()
        """
        n = (c_int*1)()
        # n = (c_long*1)()
        M.MUPHY.getSizeArrayFluid(self.ref, byref(n))
        return n[0]

    def getLocalDensity(self,i,j,k):
        """
        get local density of fluid on point
        -> M.MUPHY.getLocalDensityFluid()
        """
        rho_l = (c_float*1)()
        M.MUPHY.getLocalDensityFluid(self.mesh.ref,self.ref,c_int(i),c_int(j),c_int(k),byref(rho_l))
        return rho_l[0]

    def getLocalVelocity(self,i,j,k):
        """
        get local velocity of fluid on point
        -> M.MUPHY.getLocalVelocityFluid()
        """
        vx_l = (c_float*1)(); vy_l = (c_float*1)(); vz_l = (c_float*1)()
        M.MUPHY.getLocalVelocityFluid(self.mesh.ref,self.ref,c_int(i),c_int(j),c_int(k),byref(vx_l),byref(vy_l),byref(vz_l))
        return vx_l[0],vy_l[0],vz_l[0]

    def getLocalCurrent(self,i,j,k):
        """
        get local current of fluid on point
        -> M.MUPHY.getLocalCurrentFluid()
        """
        jx_l = (c_float*1)(); jy_l = (c_float*1)(); jz_l = (c_float*1)()
        M.MUPHY.getLocalCurrentFluid(self.mesh.ref,self.ref,c_int(i),c_int(j),c_int(k),byref(jx_l),byref(jy_l),byref(jz_l))
        return jx_l[0],jy_l[0],jz_l[0]

    def getLocalForce(self,i,j,k):
        """
        get local force of fluid on point
        -> M.MUPHY.getLocalForceFluid()
        """
        fx_l = (c_float*1)(); fy_l = (c_float*1)(); fz_l = (c_float*1)()
        M.MUPHY.getLocalForceFluid(self.mesh.ref,self.ref,c_int(i),c_int(j),c_int(k),byref(fx_l),byref(fy_l),byref(fz_l))
        return fx_l[0],fy_l[0],fz_l[0]

    def setLocalDensity(self,i,j,k,val):
        """
        get local density of fluid on point
        -> M.MUPHY.setLocalDensityFluid()
        """
        M.MUPHY.setLocalDensityFluid(self.mesh.ref,self.ref,c_int(i),c_int(j),c_int(k),c_float(val))

    def setLocalForce(self,i,j,k,fx,fy,fz):
        """
        set local force of fluid on point
        -> M.MUPHY.setLocalForceFluid()
        """
        M.MUPHY.setLocalForceFluid(self.mesh.ref,self.ref,c_int(i),c_int(j),c_int(k),c_float(fx),c_float(fy),c_float(fz))

    def getIOBCType(self,ioname,ioid): 
        """
        get inlet/outlet flow value
        -> M.MUPHY.getIOBCType()
        """
        bctype = ' '*256
        M.MUPHY.getIOBCType(self.ref, c_int(len(ioname)), ioname, c_int(ioid), c_int(len(bctype)), bctype)
        return bctype

    def getIOFlow(self,ioname,ioid): 
        """
        get inlet/outlet flow value
        -> M.MUPHY.getIOFlow()
        """
        val = (c_float*1)()
        M.MUPHY.getIOFlow(self.ref, c_int(len(ioname)), ioname, c_int(ioid), byref(val))
        return val[0]

    def getIOArea(self,ioname,ioid): 
        """
        get inlet/outlet area
        -> M.MUPHY.getIOArea()
        """
        val = (c_float*1)()
        M.MUPHY.getIOArea(self.ref, c_int(len(ioname)), ioname, c_int(ioid), byref(val))
        return val[0]

    def getIOPressure(self,ioname,ioid): 
        """
        get inlet/outlet flow value
        -> M.MUPHY.getIOPressure()
        """
        val = (c_float*1)()
        M.MUPHY.getIOPressure(self.ref, c_int(len(ioname)), ioname, c_int(ioid), byref(val))
        return val[0]

    def getIOValue(self,ioname,ioid): 
        """
        get inlet/outlet value
        -> M.MUPHY.getIOValue()
        """
        val = (c_float*1)()
        M.MUPHY.getIOValue(self.ref, c_int(len(ioname)), ioname, c_int(ioid), byref(val))
        return val[0]

    def setInletOutletBC(self,ioname,ioid,bctype,iodir=None,val=0.,val_dir=None): 
        """
        set inlet/outlet value
        -> M.MUPHY.setInletOutletBC()
        """
        val_dir_ = (c_float*3)()
        if val_dir==None:
            if iodir != None:
                val_dir_[0] = iodir[0] * 1.
                val_dir_[1] = iodir[1] * 1.
                val_dir_[2] = iodir[2] * 1.
            else:
                val_dir_[0] = -9999.
                val_dir_[1] = -9999.
                val_dir_[2] = -9999.

        iodir_ = (c_int*3)()
        if iodir==None:
            iodir_[0] = -9999
            iodir_[1] = -9999
            iodir_[2] = -9999
        else:
            iodir_[0] = iodir[0]
            iodir_[1] = iodir[1]
            iodir_[2] = iodir[2]

        M.MUPHY.setInletOutletBC(self.ref, c_int(len(ioname)),ioname, 
                        c_int(ioid), 
                        c_int(len(bctype)),bctype,
                        byref(iodir_),
                        c_float(val),
                        byref(val_dir_))

    def setIOValue(self,ioname,ioid,val): 
        """
        set inlet/outlet value
        -> M.MUPHY.setIOValue()
        """
        M.MUPHY.setIOValue(self.ref, c_int(len(ioname)),ioname, c_int(ioid), c_float(val))

    def setInflowVelocity(self,val): 
        """
        set inlet/outlet inflow value
        -> M.MUPHY.setInflowVelocity()
        """
        M.MUPHY.setInflowVelocity(self.ref, c_float(val))

    def setWallLeakageFactor(self,val=DEFAULT['WALLLEAKAGEFACTOR']): 
        """
        set wall leakage factor
        -> M.MUPHY.setWallLeakageFactor()
        """
        M.MUPHY.setWallLeakageFactor(self.ref, c_float(val))

    def setEquilibrium(self,i,j,k,rho,u,v,w):
        """
        set populations of fluid
        -> M.MUPHY.setEquilibriumFluid)
        """
        M.MUPHY.setEquilibriumFluid(self.ref, \
                                   c_int(i),c_int(j),c_int(k), \
                                   c_float(rho),c_float(u),c_float(v),c_float(w))

    def getPopulation(self):
        """
        get populations of fluid
        -> M.MUPHY.getPopulationFluid()
        """
        n = self.getSizeArray()
        p0 = (c_float*n)()
        p1 = (c_float*n)()
        p2 = (c_float*n)()
        p3 = (c_float*n)()
        p4 = (c_float*n)()
        p5 = (c_float*n)()
        p6 = (c_float*n)()
        p7 = (c_float*n)()
        p8 = (c_float*n)()
        p9 = (c_float*n)()
        p10 = (c_float*n)()
        p11 = (c_float*n)()
        p12 = (c_float*n)()
        p13 = (c_float*n)()
        p14 = (c_float*n)()
        p15 = (c_float*n)()
        p16 = (c_float*n)()
        p17 = (c_float*n)()
        p18 = (c_float*n)()
        M.MUPHY.getPopulationFluid(self.ref,c_int(n), \
                                    byref(p0), \
                                    byref(p1), \
                                    byref(p2), \
                                    byref(p3), \
                                    byref(p4), \
                                    byref(p5), \
                                    byref(p6), \
                                    byref(p7), \
                                    byref(p8), \
                                    byref(p9), \
                                    byref(p10), \
                                    byref(p11), \
                                    byref(p12), \
                                    byref(p13), \
                                    byref(p14), \
                                    byref(p15), \
                                    byref(p16), \
                                    byref(p17), \
                                    byref(p18))
        return [p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18]

    def setPopulation(self,pop):
        """
        set populations of fluid
        -> M.MUPHY.setPopulationFluid)
        """

        """
        p0 = pop[0]
        p1 = pop[1]
        p2 = pop[2]
        p3 = pop[3]
        p4 = pop[4]
        p5 = pop[5]
        p6 = pop[6]
        p7 = pop[7]
        p8 = pop[8]
        p9 = pop[9]
        p10 = pop[10]
        p11 = pop[11]
        p12 = pop[12]
        p13 = pop[13]
        p14 = pop[14]
        p15 = pop[15]
        p16 = pop[16]
        p17 = pop[17]
        p18 = pop[18]
        """

        n = len(pop[0])

        p0_l = (c_float*n)(*pop[0])
        p1_l = (c_float*n)(*pop[1])
        p2_l = (c_float*n)(*pop[2])
        p3_l = (c_float*n)(*pop[3])
        p4_l = (c_float*n)(*pop[4])
        p5_l = (c_float*n)(*pop[5])
        p6_l = (c_float*n)(*pop[6])
        p7_l = (c_float*n)(*pop[7])
        p8_l = (c_float*n)(*pop[8])
        p9_l = (c_float*n)(*pop[9])
        p10_l = (c_float*n)(*pop[10])
        p11_l = (c_float*n)(*pop[11])
        p12_l = (c_float*n)(*pop[12])
        p13_l = (c_float*n)(*pop[13])
        p14_l = (c_float*n)(*pop[14])
        p15_l = (c_float*n)(*pop[15])
        p16_l = (c_float*n)(*pop[16])
        p17_l = (c_float*n)(*pop[17])
        p18_l = (c_float*n)(*pop[18])
        M.MUPHY.setPopulationFluid(self.ref,c_int(n), \
                                    byref(p0_l), \
                                    byref(p1_l), \
                                    byref(p2_l), \
                                    byref(p3_l), \
                                    byref(p4_l), \
                                    byref(p5_l), \
                                    byref(p6_l), \
                                    byref(p7_l), \
                                    byref(p8_l), \
                                    byref(p9_l), \
                                    byref(p10_l), \
                                    byref(p11_l), \
                                    byref(p12_l), \
                                    byref(p13_l), \
                                    byref(p14_l), \
                                    byref(p15_l), \
                                    byref(p16_l), \
                                    byref(p17_l), \
                                    byref(p18_l))

    # prepare Fluids, passing context of this scale
    def prepare(self):
        """
        prepare fluid
        -> M.MUPHY.prepareFluid()
        """

        self.mesh = self.scale.mesh
        M.MUPHY.linkToMeshFluid(self.ref, self.mesh.ref)

        self.track = self.scale.tracker

        M.MUPHY.prepareFluid(self.ref, self.scale.ref)

    def preupdateCommon(self,scale): 
        """
        propagate all common fluids
        -> M.MUPHY.preupdateCommonFluid()
        """
        M.MUPHY.preupdateCommonFluid(scale.ref)

    def update(self,scale): 

        """
        propagate fluid by one step

        MG: the call automatically recognize the right time of updating. 
            If the updating is done then return True (u[0]=1), otherwise False (u[0]=0):
               1st call only saturated nodes are treated
               2nd call only insaturated nodes are treated
        thus the c&s on mesh is complete after 2 updating calls

        -> M.MUPHY.updateFluid()
        """

        M.MUPHY.updateFluid(self.ref,scale.ref)

    #def update2(self,scale): 
        #"""
        #propagate fluid by one step
        #-> M.MUPHY.update2Fluid()
        #"""
        #M.MUPHY.update2Fluid(self.ref,scale.ref)

    def pressuremmHg(self):
        """
        print info on fluid
        -> M.MUPHY.setPressureMmHg()
        """
        M.MUPHY.setPressureMmHg(self.ref)

    def printIt(self): 
        """
        print info on fluid
        -> M.MUPHY.printItFluid()
        """
        M.MUPHY.printItFluid(self.ref)

    def setBodyForce(self,o,field):
        """
        set body force on fluid
        -> M.MUPHY.setBodyForce()
        """
        M.MUPHY.setBodyForce(byref(self.ref), byref(o.ref), byref(field))

    inlets = {}
    def addInlet(self,id=None,aream2=None,areamm2=None,velocity=None):
        i = Empty()
        i.Id = id

        if areamm2 != None:
            i.PhysicalArea = areamm2 * 1.e-6 # convert from mm2 to m2
        elif aream2 != None:
            i.PhysicalArea = aream2
        else:
            print('addIn/Outlet interface changed....use aream2 or areamm2 instead')
            sys.exit(1)

        i.PhysicalVelocity = None

        if velocity != None:
            i.PhysicalVelocity = velocity # m/s
            i.PhysicalFlow = i.PhysicalArea * velocity # m3/s
            self.inlets[id] = i
        else:
            self.inlets[id] = i

    outlets = {}
    def addOutlet(self,id=None,aream2=None,areamm2=None,velocity=None):
        o = Empty()
        o.Id = id

        if areamm2 != None:
            o.PhysicalArea = areamm2 * 1.e-6 # convert from mm2 to m2
        elif aream2 != None:
            o.PhysicalArea = aream2
        else:
            print('addIn/Outlet interface changed....use aream2 or areamm2 instead')
            sys.exit(1)

        o.PhysicalVelocity = None

        if velocity != None:
            o.PhysicalVelocity = velocity # m/s
            o.PhysicalFlow = o.PhysicalArea * velocity # m3/s

        self.outlets[id] = o

    def getInletPhysicalArea(self,id=None):
        return self.inlets[id].PhysicalArea

    def getInletPhysicalVelocity(self,id=None):
        return self.inlets[id].PhysicalVelocity

    def getInletLatticeVelocity(self,id=None):
        return self.inlets[id].PhysicalVelocity * self.universe.Dt / self.universe.Dx

    def getOutletPhysicalVelocity(self,id=None):
        return self.outlets[id].PhysicalVelocity

    def getOutletLatticeVelocity(self,id=None):
        return self.outlets[id].PhysicalVelocity * self.universe.Dt / self.universe.Dx

    def setOutletMethodMurray(self):
        self.outletmethod = 'murray'

    def IOworkout(self,silent=False):

        # grab the inlet
        if not len(self.inlets)==1:
            print('only a single inlet allowed!')
            sys.exit(1)
        else:
            ikey = self.inlets.keys()[0]
            i = self.inlets[ikey]
            if not silent:
                print('Inlet Velocity:')
                print('     ',ikey, '    vel:', i.PhysicalVelocity, 'm/s                 ', \
                          (self.universe.Dt / self.universe.Dx ) * i.PhysicalVelocity, 'l.u.')
                print

        if self.outletmethod == None:
            pass

        elif self.outletmethod == 'murray':

            sm = 0
            for okey in self.outlets.keys():
                diam = math.sqrt(4. * self.outlets[okey].PhysicalArea / math.pi) * 1000. # in mm !!!
                sm += diam**2.2 # beware for the type of method used here !!!

            alpha = i.PhysicalFlow / sm # (m3/s) / m2.2

            # ck = 0
            if not silent:
                print('Outlet Velocities:')

            for okey in self.outlets.keys():

                diam = math.sqrt(4. * self.outlets[okey].PhysicalArea / math.pi) * 1000. # in mm !!!
                flow = alpha * diam**2.2

                o = self.outlets[okey]
                o.PhysicalFlow = flow
                o.PhysicalVelocity = flow / o.PhysicalArea
                # ck += flow

                if not silent:
                    print('     ',okey, \
                        'vel:', self.outlets[okey].PhysicalVelocity, 'm/s    ', \
                        (self.universe.Dt / self.universe.Dx ) * o.PhysicalVelocity,'l.u.',\
                        '   flow:', i.PhysicalFlow, self.outlets[okey].PhysicalFlow)

            # if not silent: print '\nCheck for equality:',ck,'=',i.PhysicalFlow,'m3/s','\n'

        else:
            print('outlet method unknown:',self.outletmethod)
            sys.exit(1)


    def monitorInletOutlet(self,itime,iofreq,printfreq):

        myid = self.universe.getMyproc()

        if self.monitorInletOutletInited == None:

            if myid==0:
                self.monitorfp = open('outpressure.dat','w')
                self.monitorff = open('outflow.dat','w')

            self.monitorInletOutletInited = True

        if myid==0 and itime%iofreq == 0:
            self.monitorfp.write(' %d' % itime)
            self.monitorff.write(' %d' % itime)

        # Inlet
        for key in self.inlets.keys():

            id = self.inlets[key].Id

            p_in = self.getIOPressure('inlet',id)
            f_in = self.getIOFlow('inlet',id)
            a_in = self.getIOArea('inlet',id)

            if myid == 0 and itime%printfreq == 0: 
                print('Step:',itime,'id:',id, \
                      ' Inlet Press:',p_in, \
                      'Flow:',f_in, \
                      'Area:',a_in, \
                      'Mean Velocity:',f_in/max(1.e-8,a_in))

            if myid ==0 and itime%iofreq == 0:
                self.monitorfp.write(' %g' % p_in); self.monitorff.write(' %g' % f_in)

        # Outlets
        for key in self.outlets.keys():

            id = self.outlets[key].Id

            p_out = self.getIOPressure('outlet',id)
            f_out = self.getIOFlow('outlet',id)
            a_out = self.getIOArea('outlet',id)

            if myid==0 and itime%printfreq==0: 
                print('Step:',itime,'id:',id, \
                      ' Outlet Press:',p_out, \
                      'Flow:',f_out, \
                      'Area:',a_out, \
                      'Mean Velocity:',f_out/max(1.e-8,a_out))

            if myid==0 and itime%iofreq == 0:
                self.monitorfp.write(' %g' % (p_out)); self.monitorff.write(' %g' % (f_out))

        if myid == 0 and itime%iofreq == 0:
            self.monitorfp.write('\n'); self.monitorff.write('\n'); 
            self.monitorfp.flush(); self.monitorff.flush(); 

    # check that *.ios header is ok: pressure for inlet, flow for outlets
    def checkhemosetup(self):

        for key in self.inlets.keys():
            id = self.inlets[key].Id
            bctype = self.getIOBCType('inlet',id)
            if bctype[:8] != 'pressure':
                print('Inlet BC must be pressure!')
                sys.exit(1)

        for key in self.outlets.keys():
            id = self.outlets[key].Id
            bctype = self.getIOBCType('outlet',id)
            if bctype[:4] != 'flow':
                print('Outlet BC must be flow!')
                sys.exit(1)

    def snapshot(self,type,qty,fout):

        myid = Universe().getMyproc()

        if myid != 0:
            return

        qqty = qty
        if qty == 'density':
            qqty = 'dens'
        elif qty == 'velocity':
            qqty = 'vel'

        if type == 'probe':

            print('Not implemented yet...just lazy')

        elif type == 'dat' or type == 'line':

            print('Not implemented yet...just lazy')

        elif type == 'map':

            # subprocess.Popen('gmuphy gmap -i DIRDATA_%s/%s.map -o %s -A T -b T -O T'%(self.name,qqty,fout), shell=True)
            subprocess.Popen('gmuphy gmap -i DIRDATA_%s/%s.map -o %s -b T -O T'%(self.name,qqty,fout), shell=True)
        
        elif type == 'time':

            print('Not implemented yet...just lazy')

class Atom(Actor):
    """
    class for particles
    """

    DEFAULT = dict(NAME = 'Air',
                   IDENTIFIER = ['molecule','dna'],
                   TEMPERATURE = 0.0,
                   PROPAGATIONTYPE = ['MD','MD_1ST','AM','BD','BD_HYBRID'],
                   CONFIGURATIONFILE = 'atom.inp',
                   ENERGYFILE = 'emol.dat',
                   LOOP = [1,1,1.0],
                   ENVIRONMENTCOUPLINGMETHOD = ['nearest grid point','delta particle','bouzidi','immersed boundary'],
                   ROTATIONALMOTION = False,
                   RANDOMINSERTION = False,
                   FORCEALLPAIRS = False,
                   MINIMIZATION = False,
                   RDATCONRELATIVEFACTOR = 1.0,
                   RDATANGRELATIVEFACTOR = 1.0,
                   RDATDIHRELATIVEFACTOR = 1.0,
                   PATCHINLET = False,
                   INTERATOMICFORCES = True,
                   BLOCKCOORDINATE = 'No',
                   UNBLOCKCOORDINATE = 'No',
                   CAPFORCES = False,
                   CAPFORCESFLUID = False,
                   HOMOGENEOUSFORCE = (0.,0.,0.),
                   HOMOGENEOUSELECTRICFORCE = (0.,0.,0.),
                   INLETOUTLETMAGICDOORS = False,
                   GAMMAT = 0.1,
                   GAMMAR = 0.0,
                   PASSIVESCALAR = False,
                   VISCENHANCER = 0.0,
                   BOX = (0.,0.,0.),
                   DUMPINFORMAT = 'unformatted',
                   DUMPOUTFORMAT = 'unformatted',
                   STOPCOM = False)

    def __init__(self):
        """
        self-explain
        """
        pass

    def reference(self,id):
        """
        self-explain
        -> M.MUPHY.getReferenceAtom()
        """
        # self.ref = (py_object*1)() # address of object
        # self.ref = (c_int*1)() # address of object
        self.ref = (c_void_p*1)() # address of object

        M.MUPHY.getReferenceAtom(c_int(id), \
                                byref(self.ref))

    def setName(self,name=DEFAULT['NAME']): 
        """
        self-explain
        -> M.MUPHY.setNameAtom()
        """
        self.name = name
        M.MUPHY.setNameAtom(self.ref, c_int(len(name)), name)

    def setIdentifier(self,name=DEFAULT['IDENTIFIER']): 
        """
        TBF
        -> M.MUPHY.setIdentifiedAtom()
        """
        M.MUPHY.setIdentifiedAtom(self.ref, c_int(len(name)), name)

    def setTemperature(self,val=DEFAULT['TEMPERATURE']): 
        """
        self-explain
        -> M.MUPHY.setTemperatureAtom()
        """
        M.MUPHY.setTemperatureAtom(self.ref, c_float(val))
    def getTemperature(self):
        """
        self-explain
        -> M.MUPHY.getTemperatureAtom()
        """
        val = (c_float*1)()
        M.MUPHY.getTemperatureAtom(self.ref, byref(val))
        return val[0]

    def setPropagationType(self,name=DEFAULT['PROPAGATIONTYPE']): 
        """
        self-explain
        -> M.MUPHY.setPropagationTypeAtom()
        """
        M.MUPHY.setPropagationTypeAtom(self.ref, c_int(len(name)), name)

    def setConfigurationFile(self,name=DEFAULT['CONFIGURATIONFILE']): 
        """
        self-explain
        -> M.MUPHY.setConfigurationFileAtom()
        """
        M.MUPHY.setConfigurationFileAtom(self.ref, c_int(len(name)), name)

    def setConfigurationDump(self,start=0,frequency=-1,name='CONF.xyz'):
        """
        self-explain
        -> M.MUPHY.setConfigurationDumpAtom()
        """
        M.MUPHY.setConfigurationDumpAtom(self.ref, c_int(len(name)),name,c_int(start),c_int(frequency))

    def setVelocityDump(self,start=0,frequency=-1,name='VEL.xyz'):
        """
        self-explain
        -> M.MUPHY.setVelocityDumpAtom()
        """
        M.MUPHY.setVelocityDumpAtom(self.ref, c_int(len(name)),name,c_int(start),c_int(frequency))

    def setSynergyDump(self,start=0,frequency=-1,name='SYN.xyz'):
        """
        self-explain
        -> M.MUPHY.setSynergyDumpAtom()
        """
        M.MUPHY.setSynergyDumpAtom(self.ref, c_int(len(name)),name,c_int(start),c_int(frequency))

    def setEnergyFile(self,name=DEFAULT['ENERGYFILE']): 
        """
        self-explain
        -> M.MUPHY.setEnergyFileAtom()
        """
        M.MUPHY.setEnergyFileAtom(self.ref, c_int(len(name)),name)

    def setConfigurationVelocityDump(self,start=0,frequency=-1,name='VEL.xyz'):
        """
        self-explain
        -> M.MUPHY.setConfigurationVelocityDumpAtom()
        """
        M.MUPHY.setConfigurationVelocityDumpAtom(self.ref, c_int(len(name)),name,c_int(start),c_int(frequency))

    def setLoop(self,outercycle=DEFAULT['LOOP'][0],
                     innercycle=DEFAULT['LOOP'][1],
                     innertimestep=DEFAULT['LOOP'][2]): 
        """
        outer,inner cycle and inntertimestep for MD looping
        -> M.MUPHY.setLoopAtom()
        """
        M.MUPHY.setLoopAtom(self.ref, c_int(outercycle),c_int(innercycle),c_float(innertimestep))

    def setEnvironmentCouplingMethod(self,name=DEFAULT['ENVIRONMENTCOUPLINGMETHOD']): 
        """
        self-explain
        -> M.MUPHY.setEnvironmentCouplingMethodAtom()
        """
        M.MUPHY.setEnvironmentCouplingMethodAtom(self.ref, c_int(len(name)), name)

    def setEnvironmentCouplingMethodMG(self,scale,name):
        """
        self-explain
        -> M.MUPHY.setEnvironmentCouplingMethodMGAtom()
        """
        M.MUPHY.setEnvironmentCouplingMethodMGAtom(self.ref, scale.ref, c_int(len(name)), name)

    def setConfigurationAngularVelocityDump(self,start,frequency): 
        """
        self-explain
        -> M.MUPHY.setConfigurationAngularVelocityDumpAtom()
        """
        M.MUPHY.setConfigurationAngularVelocityDumpAtom(self.ref, c_int(start),c_int(frequency))

    def setRotationalMotion(self,val=DEFAULT['ROTATIONALMOTION']): 
        """
        self-explain
        -> M.MUPHY.setRotationalMotionAtom()
        """
        M.MUPHY.setRotationalMotionAtom(self.ref, c_bool(val))

    def setRandomInsertion(self,val=DEFAULT['RANDOMINSERTION']): 
        """
        self-explain
        -> M.MUPHY.setRandomInsertion()
        """
        M.MUPHY.setRandomInsertion(self.ref, c_bool(val))

    def setForceAllPairs(self,val=DEFAULT['FORCEALLPAIRS']): 
        """
        self-explain
        -> M.MUPHY.setForceAllPairs()
        """
        M.MUPHY.setForceAllPairs(self.ref, c_bool(val))

    def setBennettMinimization(self,val=DEFAULT['MINIMIZATION']): 
        """
        self-explain
        -> M.MUPHY.setBennettMinimizationAtom()
        """
        M.MUPHY.setBennettMinimizationAtom(self.ref, c_bool(val))

    def scaleRdatconRelativeFactor(self,val=DEFAULT['RDATCONRELATIVEFACTOR']): 
        """
        self-explain
        -> M.MUPHY.scaleRdatconRelativeFactorAtom()
        """
        M.MUPHY.scaleRdatconRelativeFactorAtom(self.ref, c_bool(val))

    def scaleRdatangRelativeFactor(self,val=DEFAULT['RDATANGRELATIVEFACTOR']): 
        """
        self-explain
        -> M.MUPHY.scaleRdatangRelativeFactorAtom()
        """
        M.MUPHY.scaleRdatangRelativeFactorAtom(self.ref, c_bool(val))

    def scaleRdatdihRelativeFactor(self,val=DEFAULT['RDATDIHRELATIVEFACTOR']): 
        """
        self-explain
        -> M.MUPHY.scaleRdatdihRelativeFactorAtom()
        """
        M.MUPHY.scaleRdatdihRelativeFactorAtom(self.ref, c_bool(val))

    def setPatchInlet(self,val=DEFAULT['PATCHINLET']): 
        """
        self-explain
        -> M.MUPHY.setPatchInletAtom()
        """
        M.MUPHY.setPatchInletAtom(self.ref, c_bool(val))

    def setInterAtomicForces(self,val=DEFAULT['INTERATOMICFORCES']): 
        """
        self-explain
        -> M.MUPHY.setInterAtomicForcesAtom()
        """
        M.MUPHY.setInterAtomicForcesAtom(self.ref, c_bool(val))

    def setFreeze(self,id): 
        """
        self-explain
        -> M.MUPHY.setFreezeAtom()
        """
        M.MUPHY.setFreezeAtom(self.ref, c_int(id))
    def setParticleFreeze(self,ipf):
        """
        self-explain
        -> M.MUPHY.setParticleFreezeAtom()
        """
        n = len(ipf)
        ipf_l = (c_int * len(ipf))(*ipf)
        M.MUPHY.setParticleFreezeAtom(self.ref, c_int(n), byref(ipf_l))

    def setParticleFreezeX(self,ipf):

        self.setParticleFreezeDirection(ipf,1)

    def setParticleFreezeY(self,ipf):

        self.setParticleFreezeDirection(ipf,2)

    def setParticleFreezeZ(self,ipf):

        self.setParticleFreezeDirection(ipf,3)

    def setParticleFreezeDirection(self,ipf,idir):
        """
        self-explain
        -> M.MUPHY.setParticleFreezeDirectionAtom()
        """
        n = len(ipf)
        ipf_l = (c_int * len(ipf))(*ipf)
        M.MUPHY.setParticleFreezeDirectionAtom(self.ref, c_int(n), byref(ipf_l), c_int(idir))

    def setBlockCoordinate(self,axe=DEFAULT['BLOCKCOORDINATE']): 
        """
        blocking motion for atom along given direction 
        -> M.MUPHY.setBlockCoordinateAtom()
        """
        M.MUPHY.setBlockCoordinateAtom(self.ref, axe)
   
    def setUnFreeze(self,id): 
        """
        self-explain
        -> M.MUPHY.setUnfreezeAtom()
        """
        M.MUPHY.setUnfreezeAtom(self.ref, c_int(id))
    def printFrozen(self): 
        """
        self-explain
        -> M.MUPHY.printFrozenAtom()
        """
        M.MUPHY.printFrozenAtom(self.ref)
    def setUnblockCoordinate(self,axe=DEFAULT['UNBLOCKCOORDINATE']): 
        """
        blocking motion for atom along given direction 
        -> M.MUPHY.setUnblockCoordinateAtom()
        """
        M.MUPHY.setUnblockCoordinateAtom(self.ref, axe)
   
    def setActiveMatterFrictionTemperature(self,valf,valt): 
        """
        friction for active matter
        -> M.MUPHY.setActiveMatterFrictionAtom()
        """
        M.MUPHY.setActiveMatterFrictionTemperatureAtom(self.ref, c_float(valf), c_float(valt))

    def setCapForces(self,val=DEFAULT['CAPFORCES'],forcecap=1.,torquecap=1.,velcap=1.,angvelcap=1.): 
        """
        self-explain
        -> M.MUPHY.setCapForcesAtom()
        """
        M.MUPHY.setCapForcesAtom(self.ref, c_bool(val),
                                            c_float(forcecap),c_float(torquecap),
                                            c_float(velcap),c_float(angvelcap))
    def setCapForcesFluid(self,val=DEFAULT['CAPFORCESFLUID'],momroof=1.,angmomroof=1.): 
        """
        self-explain
        -> M.MUPHY.setCapForcesFluidAtom()
        """
        M.MUPHY.setCapForcesFluidAtom(self.ref, c_bool(val),
                                                c_float(momroof),c_float(angmomroof))
    def setHomogeneousForce(self,f1=DEFAULT['HOMOGENEOUSFORCE'][0],
                                 f2=DEFAULT['HOMOGENEOUSFORCE'][1],
                                 f3=DEFAULT['HOMOGENEOUSFORCE'][2]): 
        """
        self-explain
        -> M.MUPHY.setHomogeneousForceAtom()
        """
        M.MUPHY.setHomogeneousForceAtom(self.ref,c_float(f1),c_float(f2),c_float(f3))

    def setHomogeneousElectricForce(self,f1=DEFAULT['HOMOGENEOUSELECTRICFORCE'][0],
                                         f2=DEFAULT['HOMOGENEOUSELECTRICFORCE'][1],
                                         f3=DEFAULT['HOMOGENEOUSELECTRICFORCE'][2]): 
        """
        self-explain
        -> M.MUPHY.setHomogeneousElectricForceAtom()
        """
        M.MUPHY.setHomogeneousElectricForceAtom(self.ref,c_float(f1),c_float(f2),c_float(f3))

    def setInletOutletMagicDoors(self,val=DEFAULT['INLETOUTLETMAGICDOORS']): 
        """
        self-explain
        -> M.MUPHY.setInletOutletMagicDoors()
        """
        M.MUPHY.setInletOutletMagicDoors(self.ref, c_bool(val))

    def setSolvation(self,ia,isp,fact): 
        """
        self-explain
        -> M.MUPHY.setSolvation()
        -> M.MUPHY.setGammaRAtom()
        """
        M.MUPHY.setSolvationAtom(self.ref, c_int(ia), c_int(isp), c_float(fact))

    def setGamma(self,ia,gammaT=DEFAULT['GAMMAT'],gammaR=DEFAULT['GAMMAR']): 
        """
        self-explain
        -> M.MUPHY.setGammaTAtom()
        -> M.MUPHY.setGammaRAtom()
        """
        if gammaT != None:
            M.MUPHY.setGammaTAtom(self.ref, c_int(ia), c_float(gammaT))

        if gammaR != None:
            M.MUPHY.setGammaRAtom(self.ref, c_int(ia), c_float(gammaR))

    def setPassiveScalar(self,isp,val=DEFAULT['PASSIVESCALAR']): 
        """
        self-explain
        -> M.MUPHY.setPassiveScalarAtom()
        """
        M.MUPHY.setPassiveScalarAtom(self.ref, c_int(isp), c_bool(val))
    def setViscEnhance(self,isp,val=DEFAULT['VISCENHANCER']): 
        """
        self-explain
        -> M.MUPHY.setViscEnhance()
        """
        M.MUPHY.setViscEnhance(self.ref, c_int(isp), c_float(val))
    def scaleVdwParameters(self,val): 
        """
        self-explain
        -> M.MUPHY.scaleVdwParameters()
        """
        M.MUPHY.scaleVdwParameters(self.ref, c_float(val))
    def scaleVelocity(self,val): 
        """
        self-explain
        -> M.MUPHY.scaleVelocityAtom()
        """
        M.MUPHY.scaleVelocityAtom(self.ref, c_float(val))
    def scaleVelocityToTemperature(self,val): 
        """
        self-explain
        -> M.MUPHY.scaleVelocityToTemperatureAtom()
        """
        M.MUPHY.scaleVelocityToTemperatureAtom(self.ref, c_float(val))

    def getPotentialEnergy(self): 
        """
        self-explain
        -> M.MUPHY.getPotentialEnergyAtom()
        """
        val = (c_float*1)()
        M.MUPHY.getPotentialEnergyAtom(self.ref,byref(val))
        return val[0]

    def getInternalPressure(self): 
        """
        self-explain
        -> M.MUPHY.getInternalPressureAtom()
        """
        val = (c_float*1)()
        M.MUPHY.getInternalPressureAtom(self.ref,byref(val))
        return val[0]

    def getKineticEnergy(self): 
        """
        self-explain
        -> M.MUPHY.getKineticEnergyAtom()
        """
        val = (c_float*1)()
        M.MUPHY.getKineticEnergyAtom(self.ref,byref(val))
        return val[0]

    def getNumberTotal(self):
        """
        self-explain
        -> M.MUPHY.getNumberTotalAtom()
        """
        natms_tot = (c_int*1)()
        M.MUPHY.getNumberTotalAtom(self.ref,byref(natms_tot))
        return natms_tot[0]

    def getNumber(self):
        """
        self-explain
        -> M.MUPHY.getNumberAtom()
        """
        natms = (c_int*1)()
        M.MUPHY.getNumberAtom(self.ref,byref(natms))
        return natms[0]

    def getNumberExtended(self):
        """
        self-explain
        -> M.MUPHY.getNumberExtendedAtom()
        """
        natms_ext = (c_int*1)()
        M.MUPHY.getNumberExtendedAtom(self.ref,byref(natms_ext))
        return natms_ext[0]

    def getUid(self):
        """
        self-explain
        -> M.MUPHY.getUidAtom()
        """
        # n = self.getNumber(); px_l = (c_int*n)()
        n = self.getNumberExtended(); px_l = (c_int*n)()
        M.MUPHY.getUidAtom(self.ref,c_int(n),byref(px_l))
        return px_l

    def getPosition(self):
        """
        self-explain
        -> M.MUPHY.getPositionAtom()
        """
        # n = self.getNumber(); px_l = (c_float*n)(); py_l = (c_float*n)(); pz_l = (c_float*n)()
        n = self.getNumberExtended(); px_l = (c_float*n)(); py_l = (c_float*n)(); pz_l = (c_float*n)()
        M.MUPHY.getPositionAtom(self.ref,c_int(n),byref(px_l),byref(py_l),byref(pz_l))
        return px_l,py_l,pz_l

    def getVelocity(self):
        """
        self-explain
        -> M.MUPHY.getVelocityAtom()
        """
        # no remote velocities are transferred from neighbor procs
        n = self.getNumber(); px_l = (c_float*n)(); py_l = (c_float*n)(); pz_l = (c_float*n)()
        M.MUPHY.getVelocityAtom(self.ref,c_int(n),byref(px_l),byref(py_l),byref(pz_l))
        return px_l,py_l,pz_l

    def getForce(self):
        """
        self-explain
        -> M.MUPHY.getForceAtom()
        """
        # no remote forces are transferred from neighbor procs
        n = self.getNumber(); px_l = (c_float*n)(); py_l = (c_float*n)(); pz_l = (c_float*n)()
        M.MUPHY.getForceAtom(self.ref,c_int(n),byref(px_l),byref(py_l),byref(pz_l))
        return px_l,py_l,pz_l

    def getElectricForce(self):
        """
        self-explain
        -> M.MUPHY.getElectricForceAtom()
        """
        # no remote forces are transferred from neighbor procs
        n = self.getNumber(); px_l = (c_float*n)(); py_l = (c_float*n)(); pz_l = (c_float*n)()
        M.MUPHY.getElectricForceAtom(self.ref,c_int(n),byref(px_l),byref(py_l),byref(pz_l))
        return px_l,py_l,pz_l

    def getMass(self):
        """
        self-explain
        -> M.MUPHY.getMassAtom()
        """
        # no remote forces are transferred from neighbor procs
        n = self.getNumber(); pm_l = (c_float*n)()
        M.MUPHY.getMassAtom(self.ref,c_int(n),byref(pm_l))
        return pm_l

    def getShearstress(self):
        """
        self-explain
        -> M.MUPHY.getStressAtom()
        """
        n = self.getNumber(); px_l = (c_float*n)()
        M.MUPHY.getShearstressAtom(self.ref,c_int(n),byref(px_l))
        return px_l

    def getForceFast(self):
        """
        self-explain
        -> M.MUPHY.getForceFastAtom()
        """
        n = self.getNumber(); px_l = (c_float*n)(); py_l = (c_float*n)(); pz_l = (c_float*n)()
        M.MUPHY.getForceFastAtom(self.ref,c_int(n),byref(px_l),byref(py_l),byref(pz_l))
        return px_l,py_l,pz_l

    #PM 11-Apr-18
    # It returns the array whose i-th element is the delta density 
    # at dd of the i-th atom (with support radius=2)
    # Use only cpl = COUPLE_DELTA or cpl=COUPLE_IMMERSEDBOUNDARY
    def getDeltaIsoAtom(self,dd,cpl="COUPLE_DELTA"):
        """
        self-explain
        -> M.MUPHY.getStressAtom()
        """
        if cpl=="COUPLE_DELTA":
          keycpl=2
        elif cpl=="COUPLE_IMMERSEDBOUNDARY":
          keycpl=5
        else:
          print('>>> ERROR: undefined coupling method for getDeltaIsoAtom. Aborting...\n')
          sys.exit(1)
        
        is_inside=(c_int*1)()
        ddd=(c_float*3)()
        dens=[]
        n = self.getNumberExtended(); px_l = (c_float*n)(); py_l = (c_float*n)(); pz_l = (c_float*n)()
        M.MUPHY.getPositionAtom(self.ref,c_int(n),byref(px_l),byref(py_l),byref(pz_l))
        for i in range(n):
          ddd[0]=dd[0]-px_l[i]
          ddd[1]=dd[1]-py_l[i]
          ddd[2]=dd[2]-pz_l[i]
          out=(c_float*1)()
          M.MUPHY.getDeltaIsoAtom(c_int(keycpl),byref(ddd),byref(is_inside),byref(out))
          dens.append(out[0])
          
        return dens

    ##########

    def setPosition(self,px,py,pz):
        """
        self-explain
        -> M.MUPHY.setPositionAtom()
        """
        n = len(px)
        px_l = (c_float*n)(*px); py_l = (c_float*n)(*py); pz_l = (c_float*n)(*pz)
        M.MUPHY.setPositionAtom(self.ref,c_int(n), byref(px_l), byref(py_l), byref(pz_l))
   
    def setVelocity(self,px,py,pz):
        """
        self-explain
        -> M.MUPHY.setVelocityAtom()
        """
        n = len(px)
        px_l = (c_float*n)(*px); py_l = (c_float*n)(*py); pz_l = (c_float*n)(*pz)
        M.MUPHY.setVelocityAtom(self.ref,c_int(n), byref(px_l), byref(py_l), byref(pz_l))
   
    def setForce(self,px,py,pz):
        """
        self-explain
        -> M.MUPHY.setForceAtom()
        """
        n = len(px)
        px_l = (c_float*n)(*px); py_l = (c_float*n)(*py); pz_l = (c_float*n)(*pz)
        M.MUPHY.setForceAtom(self.ref,c_int(n), byref(px_l), byref(py_l), byref(pz_l))
   
    def setBox(self,bx=DEFAULT['BOX'][0],
                    by=DEFAULT['BOX'][1],
                    bz=DEFAULT['BOX'][2]): 
        """
        self-explain
        -> M.MUPHY.setBoxAtom()
        """
        M.MUPHY.setBoxAtom(self.ref, c_float(bx), c_float(by), c_float(bz))

    def getBox(self): 
        """
        self-explain
        -> M.MUPHY.getBoxAtom()
        """
        bx = (c_float*1)()
        by = (c_float*1)()
        bz = (c_float*1)()
        M.MUPHY.getBoxAtom(self.ref, byref(bx), byref(by), byref(bz))
        return [bx[0],by[0],bz[0]]

    def setDumpInFormat(self,fmt=DEFAULT['DUMPINFORMAT']): 
        """
        set dump input format of atom
        -> M.MUPHY.setDumpInFormatAtom()
        """
        M.MUPHY.setDumpInFormatAtom(self.ref, c_int(len(fmt)), fmt)

    def setDumpOutFormat(self,fmt=DEFAULT['DUMPOUTFORMAT']): 
        """
        set dump output format of atom
        -> M.MUPHY.setDumpOutFormatAtom()
        """
        M.MUPHY.setDumpOutFormatAtom(self.ref, c_int(len(fmt)), fmt)

    def setStopCOM(self,val=DEFAULT['STOPCOM']): 
        """
        set flag to stop center of mass for atom
        -> M.MUPHY.setStopCOMAtom()
        """
        M.MUPHY.setStopCOMAtom(self.ref, c_bool(val))

    # prepare Atoms, passing context of this scale
    def prepare(self):
        """
        self-explain
        -> M.MUPHY.prepareAtom()
        """

        # self.linkToMesh(self.scale.mesh)
        self.mesh = self.scale.mesh
        M.MUPHY.linkToMeshAtom(self.ref, self.mesh.ref)

        self.track = self.scale.tracker

        M.MUPHY.prepareAtom(self.ref, self.scale.ref)

    def preupdateCommon(self,scale): 
        """
        propagate all common atoms
        -> M.MUPHY.preupdateCommonAtom()
        """
        M.MUPHY.preupdateCommonAtom(scale.ref)

    def update(self,scale): 
        """
        self-explain
        -> M.MUPHY.updateAtom()
        """
        u = (c_int*1)()
        M.MUPHY.updateAtom(self.ref,scale.ref,u)
        return u[0]==1

    def getCOM(self):
        cx = (c_float*1)()
        cy = (c_float*1)()
        cz = (c_float*1)()
        M.MUPHY.getCOMAtom(self.ref, byref(cx), byref(cy), byref(cz))

        return(cx[0],cy[0],cz[0])

    def getMolecularCOM(self,indmol=1):
        cx = (c_float*1)()
        cy = (c_float*1)()
        cz = (c_float*1)()
        M.MUPHY.getMolecularCOMAtom(self.ref, c_int(indmol), byref(cx), byref(cy), byref(cz))

        return(cx[0],cy[0],cz[0])

    def getMolecularVelocity(self,indmol=1):

        vx = (c_float*1)()
        vy = (c_float*1)()
        vz = (c_float*1)()
        M.MUPHY.getMolecularVelocityAtom(self.ref, c_int(indmol), byref(vx), byref(vy), byref(vz))

        return(vx[0],vy[0],vz[0])

    def rescaleMolecularVelocity(self,indmol=1,action='scale'):

        if action == 'scale':
            iscale = 1
        elif action == 'gaussian':
            iscale = 2
        else:
            print('rescale Molecular Velocity can only be scale/gaussian, got:',action)
            sys.exit(1)

        M.MUPHY.rescaleMolecularVelocityAtom(self.ref, c_int(indmol), c_int(iscale))

    def rescaleMolecularVelocities(self,action='scale'):

        if action == 'scale':
            iscale = 1
        elif action == 'gaussian':
            iscale = 2
        else:
            print('rescale Molecular Velocity can only be scale/gaussian, got:',action)
            sys.exit(1)

        M.MUPHY.rescaleMolecularVelocitiesAtom(self.ref, c_int(iscale))

    def rescaleMolecularAngularVelocity(self,indmol=1,action='scale'):

        if action == 'scale':
            iscale = 1
        elif action == 'gaussian':
            iscale = 2
        else:
            print('rescale Molecular Angular Velocity can only be scale/gaussian, got:',action)
            sys.exit(1)

        M.MUPHY.rescaleMolecularAngularVelocityAtom(self.ref, c_int(indmol), c_int(iscale))

    def rescaleMolecularAngularVelocities(self,action='scale'):

        if action == 'scale':
            iscale = 1
        elif action == 'gaussian':
            iscale = 2
        else:
            print('rescale Molecular Angular Velocity can only be scale/gaussian, got:',action)
            sys.exit(1)

        M.MUPHY.rescaleMolecularAngularVelocitiesAtom(self.ref, c_int(iscale))

    def rescaleVelocities(self,action='scale'):

        if action == 'scale':
            iscale = 1
        elif action == 'gaussian':
            iscale = 2
        else:
            print('rescale  Velocity can only be scale/gaussian, got:',action)
            sys.exit(1)

        M.MUPHY.rescaleVelocitiesAtom(self.ref, c_int(iscale))

    def printIt(self): 
        """
        self-explain
        -> M.MUPHY.printItAtom()
        """
        M.MUPHY.printItAtom(self.ref)

class ODE(Actor):
    """
    self-explain
    """

    def __init__(self):
        """
        self-explain
        """
        pass

    def reference(self,id):
        """
        set self-reference and solve, field and fieldgradient
        -> M.MUPHY.getReferenceODE()
        """
        # self.ref = (py_object*1)() # address of object
        # self.ref = (c_int*1)() # address of object
        self.ref = (c_void_p*1)() # address of object

        # self.solve = (py_object*1)() # object function
        # self.solve = (c_int*1)() # object function
        self.solve = (c_void_p*1)() # object function

        # self.field = (py_object*1)() # object function
        # self.field = (c_int*1)() # object function
        self.field = (c_void_p*1)() # object function

        # self.fieldgradient = (py_object*1)() # object function
        # self.fieldgradient = (c_int*1)() # object function
        self.fieldgradient = (c_void_p*1)() # object function

        M.MUPHY.getReferenceODE(c_int(id), \
                                byref(self.ref), \
                                byref(self.solve), \
                                byref(self.field), \
                                byref(self.fieldgradient))

    def setName(self,name): 
        """
        self-explain
        -> M.MUPHY.setNameODE()
        """
        self.name = name
        M.MUPHY.setNameODE(self.ref, c_int(len(name)), name)

    def setGeneric(self,bctype,action,diffusion=True,advection=True,initialvalue=1.0,diffusivity=1.0,mixedbclength=1.e9): 
        """
        self-explain
        -> M.MUPHY.setGenericODE()
        """

        M.MUPHY.setGenericODE(self.ref, c_int(len(bctype)), bctype, 
                                        c_int(len(action)), action,
                                        c_bool(diffusion), 
                                        c_bool(advection),
                                        c_float(initialvalue),
                                        c_float(diffusivity),
                                        c_float(mixedbclength))

    def setElectro(self,bjerrum=1.,debye=1.,fix_debye=True): 
        """
        self-explain
        -> M.MUPHY.setElectroODE()
        """

        M.MUPHY.setElectroODE(self.ref, 
                            c_float(bjerrum), 
                            c_float(debye),
                            c_bool(fix_debye))

    def setElectroBjerrum(self,bjerrum=1.): 
        """
        self-explain
        -> M.MUPHY.setElectroX()
        """

        M.MUPHY.setElectroBjerrum(self.ref, 
                            c_float(bjerrum))

    def setDensity(self,tolerance=0.00001,maxiter=10000): 
        """
        self-explain
        -> M.MUPHY.setDensityODE()
        """

        M.MUPHY.setDensityODE(self.ref, 
                            c_float(tolerance), 
                            c_int(maxiter))

    def setLocalMixedbclength(self,arr): 
        """
        set initial density for fluid
        -> M.MUPHY.setLocalMixedbclengthODE()
        """
        n = len(arr)

        M.MUPHY.setLocalMixedbclengthODE(self.ref, c_int(n), arr)

    def setSOR(self,iterations=100,omega=0.3,tolerance=0.001,**kwargs): 
        """
        self-explain
        -> M.MUPHY.setSORODE()
        """
        M.MUPHY.setSORODE(self.ref, c_int(iterations), c_float(omega), c_float(tolerance))

    def setSORslow(self): 
        """
        self-explain
        -> M.MUPHY.setSORslowODE()
        """
        M.MUPHY.setSORslowODE(self.ref)

    def printIt(self): 
        """
        self-explain
        -> M.MUPHY.printItODE()
        """
        M.MUPHY.printItODE(self.ref)

    # prep ODEs, passing context of this scale
    def prepare(self):
        """
        self-explain
        -> M.MUPHY.prepareODE()
        """

        self.mesh = self.scale.mesh
        M.MUPHY.linkToMeshAtom(self.ref, self.mesh.ref)

        self.track = self.scale.tracker

        M.MUPHY.prepareODE(self.ref, self.scale.ref)

    def preupdateCommon(self,scale): 
        """
        propagate all common ODEs
        -> M.MUPHY.preupdateCommonODE()
        """
        M.MUPHY.preupdateCommonODE(scale.ref)

    def update(self,scale): 
        """
        self-explain
        -> M.MUPHY.updateODE()
        """
        M.MUPHY.updateODE(self.ref,scale.ref)

    def setFieldProfile(self,field): 
        """
        set profile for field quantity
        -> M.MUPHY.setFieldProfileFluid()
        """
        n = len(field)

        M.MUPHY.setFieldProfileODE(self.ref, c_int(n), field)

    def setFreeze(self,val): 
        """
        set freeze flag for ODE
        -> M.MUPHY.setFreezeODE()
        """
        M.MUPHY.setFreezeODE(self.ref, c_bool(val))


    # def setOperator(self,owner,op):
    #     self.operatorOwner = owner
    #     self.operator = op

    # def solveODE(self,f):
    #     M.MUPHY.solveODE(byref(self.ref), byref(self.operatorOwner.ref), byref(self.operator))

class Frame(Actor):
    """
    self-explain
    """

    def __init__(self):
        """
        self-explain
        """
        pass

    def reference(self,id):
        """
        -> M.MUPHY.getReferenceFrame()
        """
        self.ref = (c_void_p*1)() # address of object

        M.MUPHY.getReferenceFrame(c_int(id), \
                                  byref(self.ref))

    def setName(self,name): 
        """
        self-explain
        -> M.MUPHY.setNameODE()
        """
        self.name = name
        M.MUPHY.setNameFrame(self.ref, c_int(len(name)), name)

class Tracker(Actor):
    """
    class for tracker
    """

    DEFAULT = dict(NAME = 'Track',
                   DIAGNOSTICFREQUENCY = 10,
                   MEASUREMENT = 'density',
                   MEASUREMENTTYPE = 'map',
                   DATASHOWDENSITY = False,
                   DATASHOWVELOCITY = False,
                   DATASHOWCURRENT = False,
                   DATASHOWFLOWRATE = False,
                   DATASHOWELECTRO = False,
                   DATASHOWTEMPERATURE = False,
                   DATASHOWPRESSURE = False,
                   MAPDIRECTIONS = 'xy',
                   MERGEDUMP = False,
                   SAMPLEANALYSIS = False,
                   SAMPLEFREQUENCYANALYSIS = False,
                   SAMPLEFILES = '',
                   SAMPLETRANSFORMFILE = '',
                   MEMORYTRACKER = True,
                   VTKDUMPFLAG = False, 
                   VTKMESHTYPE = 'unstructured', 
                   VTKDUMPSTART = 1,
                   VTKDUMPSTOP = -1,
                   VTKDUMPFREQUENCY = 10,
                   VTKPOINTFREQUENCY = 1)

    def __init__(self):
        """
        self-explain
        """
        pass

    def reference(self,id):
        """
        self-explain
        -> M.MUPHY.getReferenceTracker()
        """

        # self.ref = (py_object*1)() # reference to object
        # self.ref = (c_int*1)() # reference to object
        self.ref = (c_void_p*1)() # reference to object
        M.MUPHY.getReferenceTracker(c_int(id), \
                                    byref(self.ref))

    def prepare(self):
        """
        self-explain
        -> M.MUPHY.prepareTracker()
        """

        M.MUPHY.prepareTracker(self.ref, self.scale.ref)

    def setName(self,name=DEFAULT['NAME']): 
        """
        self-explain
        -> M.MUPHY.setNameTracker()
        """
        self.name = name
        M.MUPHY.setNameTracker(self.ref, c_int(len(name)), name)

    def setDiagnosticFrequency(self,val=DEFAULT['DIAGNOSTICFREQUENCY']): 
        """
        self-explain
        -> M.MUPHY.setDiagnosticFrequencyTracker()
        """
        M.MUPHY.setDiagnosticFrequencyTracker(self.ref, c_int(val))

    def setMeasurement(self, fluid,
                       qty=DEFAULT['MEASUREMENT'],
                       ftype=DEFAULT['MEASUREMENTTYPE'],
                       fname='unknown',
                       ijk=[-1,-1,-1],
                       action='ontime', start=0, stop=int(1e6)): 
        """
        self-explain
        -> M.MUPHY.setMeasurementTracker()
        """

        if not isinstance(fluid,Fluid):
            print('setMeasurement requires a Fluid object as argument')
            sys.exit(1)

        if action == 'ontime':
            action_ = 0
        elif action == 'average':
            action_ = 1
        elif action == 'rms':
            action_ = 2
        else:
            print('action can only be: ontime/average/rms ! got :',action)
            sys.exit(1)

        # convert list to ctypes array
        _ijk = (c_int * len(ijk))(*ijk)

        M.MUPHY.setMeasurementTracker(self.ref, fluid.ref, 
                            c_int(len(qty)), qty, 
                            c_int(len(ftype)), ftype, 
                            c_int(len(fname)), fname, 
                            _ijk, 
                            c_int(action_), c_int(start), c_int(stop))

    def setDataShow(self,density=DEFAULT['DATASHOWDENSITY'],
                         velocity=DEFAULT['DATASHOWVELOCITY'],
                         current=DEFAULT['DATASHOWCURRENT'],
                         flowrate=DEFAULT['DATASHOWFLOWRATE'],
                         electro=DEFAULT['DATASHOWELECTRO'],
                         temperature=DEFAULT['DATASHOWTEMPERATURE'],
                         pressure=DEFAULT['DATASHOWPRESSURE']): 
        """
        self-explain
        -> M.MUPHY.setDataShowTracker()
        """
        M.MUPHY.setDataShowTracker(self.ref, 
                            c_bool(density),
                            c_bool(velocity),
                            c_bool(current),
                            c_bool(flowrate),
                            c_bool(electro),
                            c_bool(temperature),
                            c_bool(pressure))

    def setMapDirections(self, name=DEFAULT['MAPDIRECTIONS'], midpoint=-999): 
        """
        self-explain
        -> M.MUPHY.setMapDirectionsTracker()
        """
        M.MUPHY.setMapDirectionsTracker(self.ref,c_int(len(name)), name, c_int(midpoint))

    def setMergeDump(self, val=DEFAULT['MERGEDUMP']): 
        """
        self-explain
        -> M.MUPHY.setMergeDumpTracker()
        """
        M.MUPHY.setMergeDumpTracker(self.ref, c_bool(val))

    def setSampleAnalysis(self,val=DEFAULT['SAMPLEANALYSIS']): 
        """
        self-explain
        -> M.MUPHY.setSampleAnalysis()
        """
        M.MUPHY.setSampleAnalysis(self.ref,c_bool(val))

    def setSampleFrequencyAnalysis(self,val=DEFAULT['SAMPLEFREQUENCYANALYSIS']): 
        """
        self-explain
        -> M.MUPHY.setSampleFrequencyAnalysis()
        """
        M.MUPHY.setSampleFrequencyAnalysis(self.ref,c_int(val))

    def setSampleFiles(self,names=DEFAULT['SAMPLEFILES']): 
        """
        self-explain
        -> M.MUPHY.setSampleFiles()
        """
        M.MUPHY.setSampleFiles(self.ref,c_int(len(names)), names)

    def setSampleTransformFile(self,name=DEFAULT['SAMPLETRANSFORMFILE']): 
        """
        self-explain
        -> M.MUPHY.setSampleTransformFile()
        """
        M.MUPHY.setSampleTransformFile(self.ref,c_int(len(name)), name)

    def getMemoryTracker(self=DEFAULT['MEMORYTRACKER']):
        """
        self-explain
        -> M.MUPHY.getMemoryTracker()
        """
        mymem = (c_float*1)()
        M.MUPHY.getMemoryTracker(self.ref,byref(mymem))
        return mymem[0]

    def setVtkDump(self,flag=DEFAULT['VTKDUMPFLAG'],
                        meshtype=DEFAULT['VTKMESHTYPE'],
                        start=DEFAULT['VTKDUMPSTART'],
                        stop=DEFAULT['VTKDUMPSTOP'],
                        frequency=DEFAULT['VTKDUMPFREQUENCY'],
                        pointfrequency=DEFAULT['VTKPOINTFREQUENCY']): 
        """
        self-explain
        -> M.MUPHY.setVtkDumpTracker()
        """
        M.MUPHY.setVtkDumpTracker(self.ref, c_bool(flag))

        if meshtype != None: self.setVtkMeshType(meshtype)
        if start != None: self.setVtkDumpStart(start)
        if stop != None: self.setVtkDumpStart(stop)
        if frequency != None: self.setVtkDumpFrequency(frequency)
        if pointfrequency != None: self.setVtkDumpPointFrequency(pointfrequency)

    def setVtkMeshType(self,name=DEFAULT['VTKMESHTYPE']): 
        """
        self-explain
        -> M.MUPHY.setVtkMeshTypeTracker()
        """
        M.MUPHY.setVtkMeshTypeTracker(self.ref,c_int(len(name)),name)
    def setVtkDumpStart(self,val=DEFAULT['VTKDUMPSTART']): 
        """
        self-explain
        -> M.MUPHY.setVtkDumpStartTracker()
        """
        M.MUPHY.setVtkDumpStartTracker(self.ref, c_int(val))
    def setVtkDumpStop(self,val=DEFAULT['VTKDUMPSTOP']): 
        """
        self-explain
        -> M.MUPHY.setVtkDumpStopTracker()
        """
        M.MUPHY.setVtkDumpStopTracker(self.ref, c_int(val))
    def setVtkDumpFrequency(self,val=DEFAULT['VTKDUMPFREQUENCY']): 
        """
        self-explain
        -> M.MUPHY.setVtkDumpFrequencyTracker()
        """
        M.MUPHY.setVtkDumpFrequencyTracker(self.ref, c_int(val))
    def setVtkDumpPointFrequency(self,val=DEFAULT['VTKPOINTFREQUENCY']): 
        """
        self-explain
        -> M.MUPHY.setVtkDumpPointFrequencyTracker()
        """
        M.MUPHY.setVtkDumpPointFrequencyTracker(self.ref, c_int(val))

    def preupdateCommon(self,scale): 
        """
        pre-update all common trackers in a scale
        -> M.MUPHY.preupdateCommonTracker()
        """
        M.MUPHY.preupdateCommonTracker(scale.ref)

    def update(self,scale): 
        """
        update all common trackers in a scale
        -> M.MUPHY.updateTracker)
        """
        M.MUPHY.updateTracker(self.ref,scale.ref)

    def printIt(self): 
        """
        self-explain
        -> M.MUPHY.printItTracker()
        """
        M.MUPHY.printItTracker(self.ref)

    def addActorsTracker(self,actors):
        """
        self-explain
        """

        self.fluids = []
        self.atoms = []
        self.ODEs = []
        for a in actors:

            if isinstance(a,Tracker): # skip itself
                continue

            if isinstance(a,Fluid):
                self.fluids.append(a)

            elif isinstance(a,Atom):
                self.atoms.append(a)

            elif isinstance(a,ODE):
                self.ODEs.append(a)

class Interact(Actor):

    def __init__(self,A1,A2):

        if isinstance(A1,Atom):
            self.a = A1
        elif isinstance(A2,Atom):
            self.a = A2

        if isinstance(A1,Fluid):
            self.f = A1
        elif isinstance(A2,Fluid):
            self.f = A2

        if isinstance(A1,ODE):
            self.o = A1
        elif isinstance(A2,ODE):
            self.o = A2

    def setCouplings(self, ids=None, drag=None, solvation=None): 

        try:
            if drag != None:
                for ia in ids:
                    self.a.setGamma(ia, drag)

            if solvation != None:
                for ia in ids:
                    self.a.setSolvation(ia, self.f.id, solvation)
        except:
            print('Incorrect Atom - Fluid coupling !!',ids,drag,solvation)

        # ... add similar snippets for Fluid - ODE, Atom - ODE, etc.

class CrossTracker(Tracker):
    """
    class for cross-scale tracker 
    """

    def setCrossTrack(self,actors):

        # check that all followees are homogeneous actors
        F = 0; A = 0; O = 0
        for a in actors:
            if   isinstance(a,Fluid): F = 1
            elif isinstance(a,Atom): A = 1
            elif isinstance(a,ODE): O = 1

        if F + A + O > 1:
            print('error: crosstrack followees are not homogeneous...',f,a,o )
            sys.exit(1)

        self.superfollowees = []
        for a in actors:
            self.superfollowees.append(a)

    def prepare(self):
        """
        self-explain
        -> M.MUPHY.prepareCrossTracker()
        """

        ff = []
        aa = []
        oo = []
        for a in self.superfollowees:
            if   isinstance(a,Fluid): ff.append(a)
            elif isinstance(a,Atom): aa.append(a)
            elif isinstance(a,ODE): oo.append(a)
        f = get_refs(ff)
        a = get_refs(aa)
        o = get_refs(oo)

        # print 'SSSS:',self.name
        # for a in self.superfollowees:
        #     if   isinstance(a,Fluid): 
        #         print 'followee :',a.name,a.track.name
        M.MUPHY.prepareCrossTracker(self.ref, 
                                    c_int(len(f)), byref(f),
                                    c_int(len(a)), byref(a),
                                    c_int(len(o)), byref(o) )

    def update(self): 
        """
        self-explain
        -> M.MUPHY.updateCrossTracker()
        """
        M.MUPHY.updateCrossTracker(self.ref)

'''
def list(count, p_items):
    """Returns a python list for the given items represented by a pointer and the number of items"""
    items = []
    for i in range(count):
        items.append(p_items[i])
    return items

def p_list(items):
    """Returns a pointer to a list of items"""
    c_items = (type(items[0])*len(items))(*items)
    p_items = cast(c_items, POINTER(type(items[0])))

    return p_items
'''

def ptr_add(ptr, offset):
    """
    returns pointer address
    """
    address = addressof(ptr.contents) + offset
    return pointer(type(ptr.contents).from_address(address))

def get_refs(items):
    """
    get reference of argument
    """

    a = (c_void_p*len(items))()
    # print list(a)
    # for i in items: print 'NAME is:',i.name

    ptr = cast(a, POINTER(c_void_p))
    # print 'CONTENTS:',ptr.contents

    for i in range(len(items)):
        # ptr = ptr_add(ptr, i*(sizeof(c_void_p)))
        ptr.contents.value = items[i].ref[0]
        ptr = ptr_add(ptr, (sizeof(c_void_p)))

    # print 'LIST(a):',list(a)

    return a

#########################################

def preprocess_parallel_mesh(ntasks, infileroot='bgkflag'):

    for i in range(ntasks):

        f = infileroot + '_'+str(i)+'.dat'
        if os.path.isfile(f): os.remove(f)

        f = infileroot + '_'+str(i)+'.hdr'
        if os.path.isfile(f): os.remove(f)

    # muphy2wrapper_init(mycomm=0) 
    # MagicBegins(mycomm=0) 

    u = Universe()
    s = Scale()
    m = Mesh()
    f = Fluid()
    t = Tracker()

    u.addItems([s,m,f,t])

    u.setTitle('Pre-processing for irregular domain decomposition')
    u.setNumberOfSteps(1)

    u.create()

    s.set(name='MonoScale', mesh=m,  actors=[f,t])

    m.setFiles(infileroot+'.hdr', infileroot+'.dat', infileroot+'.ios')
    m.setRegularMesh(False)
    m.setWorkAsPreprocessor(True)
    m.setWritePreprocessedMesh(False)

    m.setPeriodicity('000')
    m.setDomainDecomposition(8)

    s.prepare()
    m.prepare()

    tmpdir = tempfile.mkdtemp(prefix='muphy2TMP') + '/'

    # write all side scripts
    # blddeco_sh   = write_blddeco_sh(tmpdir)
    # convgraph_sh = write_convgraph_sh(tmpdir)
    # goscotch_sh  = write_goscotch_sh(tmpdir)
    # ...now run all
    # run_command( blddeco_sh )
    # run_command( convgraph_sh )
    # fi = open(os.path.join(tmpdir, 'deco'),'r')
    # fo = open(os.path.join(tmpdir, 'deco.proscotch'),'w')

    ft = open('deco', 'w')
    fi = open('nl.out', 'r') # file reports only the number of verices and edges in graph
    print >> ft, fi.readline().strip()
    fi.close()
    fi = open('com.out', 'r') # file reports the i4 of each the 18 node neighbors
    for line in fi.readlines():
        print >> ft, line.strip()
    fi.close()
    ft.close()

    METIS = True
    SCOTCH = False

    if   METIS:
        proc = 'gpmetis %s %d'%('deco', ntasks)
        print ('Running Metis ...: ', proc)
        subprocess.call(proc, shell=True)
        print('...done')

        # paste i4.out deco.part | sort -T. -n > nodeownr.inp
        fi1 = open('i4.out', 'r')
        fi2 = open('deco.part.%d'%ntasks , 'r')
        dct = {}
        while (True):
            try:
                i4 = int( fi1.readline().rstrip() )
                iproc = int( fi2.readline().rstrip() )
                dct[i4] = iproc
            except:
                break
        fi1.close()
        fi2.close()

        fo = open('nodeownr.inp' , 'w')
        for i4 in sorted (dct.keys()):
            print >> fo, i4, dct[i4]
        fo.close()


    elif SCOTCH:
        ft = open('deco', 'r')
        fo = open('deco.proscotch', 'w') # adds a header and no. of edges per node, needed by scotch
        line = ft.readline().split()
        nvert = int(line[0])
        nedges = int(line[1])
        print >> fo, '0'
        print >> fo, nvert, nedges
        print >> fo, '1 000'
        for line in ft.readlines():
            l = line.split()
            print >> fo, len(l), line.rstrip()
        ft.close()
        fo.close()

        proc = 'dgpart %d %s %s'%(ntasks,  'deco.proscotch', 'deco.%d'%ntasks )
        print('Scotch Running...:', proc)
        subprocess.Popen(proc, shell=True)

        makeownr_sh  = write_makeownr_sh(tmpdir)
        run_command( makeownr_sh + ' deco.%d i4.out nodeownr.%d'%(ntasks,ntasks) )


    h = open(infileroot + '.hdr', 'r')
    l = h.readline().split()
    nx,ny,nz = int(l[0]), int(l[1]), int(l[2])

    l = h.readline().split()
    nmaj,nfl,nwl,nin,nou = int(l[0]), int(l[1]), int(l[2]), int(l[3]), int(l[4])
    h.close()

    if nmaj != 5:
        print('Only operates for majority nodes as deadnodes (5)...halting')
        sys.exit(1)

    nodes = {}
    for i in range(ntasks):
        nodes[i] = []

    nn = nfl + nwl + nin + nou
    # print 'nfl, nwl, nin, nou:', nfl, nwl, nin, nou, 'ntotal:', nn

    fi = open('nodeownr.inp', 'r')
    fd = open(infileroot + '.dat', 'r')
    for n in range(nn) : 

        line = fi.readline().split()
        i4 = int(line[0])
        itask = int(line[1])

        k = i4 / ((nx+2)*(ny+2))

        i4 = i4 - (nx+2)*(ny+2)*k
        j = i4 / (nx+2)

        i = i4 - (nx+2)*j

        line = fd.readline().split()
        ii,jj,kk,flg = int(line[0]), int(line[1]), int(line[2]), int(line[3])

        if   ii != i: 
            print('i mismatch', ii,i); sys.exit(1)
        elif jj != j: 
            print('j mismatch', jj,j); sys.exit(1)
        elif kk != k: 
            print('k mismatch', kk,k); sys.exit(1)

        # print >> n, i4, itask
        nodes[itask].append([i,j,k,flg])
    fi.close()

    for itask in range(ntasks):

        nfl,nwl,nin,nou = 0,0,0,0
        for node in nodes[itask]:
            i,j,k,flg = node
            if flg == 1: nfl += 1
            if flg == 2: nwl += 1
            if flg == 3: nin += 1
            if flg == 4: nou += 1

        bh = open( infileroot + '_' + str(itask) + '.hdr', 'w')
        print >> bh, nx, ny, nz
        print >> bh, 5, nfl, nwl, nin, nou
        bh.close()

        bd = open( infileroot + '_' + str(itask) + '.dat', 'w')
        for node in nodes[itask]:
            i,j,k,flg = node
            print >> bd,i,j,k,flg
        bd.close()

    #mkflgown_pl  = write_mkflgown_pl(tmpdir)
    #run_command( mkflgown_pl + ' %d bgkflag.hdr bgkflag.dat nodeownr.%d'%(ntasks,ntasks) )
    #run_command( 'mv nodeownr.%d nodeownr.inp'%(ntasks) )

    run_command( 'rm nl.out com.out i4.out deco' )
    if METIS:
        run_command( 'rm deco.part.%d'%(ntasks) )
    if SCOTCH:
        run_command( 'rm deco.%d deco.proscotch'%(ntasks) )

    try:
        import TOOLS.mesh as _mesh
        # import TOOLS.bgkflag2vtk as cvt

        pds = []
        msh = _mesh.Mesh()
        for itask in range(ntasks):

            pd = msh.meshfileToPolyData('bgkflag_%d'%itask, fieldval=itask)
            pds.append( pd )

        msh.writeCombinedVTP('bgkflag_combined.vtp', pds)
    
    except:
        return

def run_command(cmd):

    time.sleep(1)
    print(cmd)
    proc = subprocess.call(shlex.split(cmd))
    time.sleep(2)

def write_blddeco_sh(tmpdir):

    text=r"""#!/bin/bash

# echo 0 > fakempi.inp
# echo 1 >> fakempi.inp
# rm flagdd_0.bin

[ ! -f com.out ] && echo "Missing file com.out...quitting" && exit 0
[ ! -f nl.out ] && echo "Missing file nl.out...quitting" && exit 0

cat nl.out com.out > deco
"""

    blddeco_sh = os.path.join(tmpdir,'blddeco.sh')
    print('...writing file:',blddeco_sh)
    tmpfl = open(blddeco_sh,'w')
    # print >> tmpfl, text
    tmpfl.write( text )
    tmpfl.close()
    try:
        os.chmod(blddeco_sh, stat.S_IEXEC)
    except:
        pass
        # os.chmod(blddeco_sh,0744)
    return blddeco_sh

def write_convgraph_sh(tmpdir):

    text=r"""#!/bin/bash
echo -n "Running convgraph.pl < deco > deco.proscotch..."
OUTPUT=$(%s/convgraph.pl < deco > deco.proscotch)
""" % (tmpdir)

    convgraph_sh = os.path.join(tmpdir,'convgraph.sh')
    tmpfl = open(convgraph_sh,'w')
    # print >> tmpfl, text
    tmpfl.write( text )
    tmpfl.close()
    try:
        os.chmod(convgraph_sh,stat.S_IEXEC)
    except:
        pass
        # os.chmod(convgraph_sh,0744)
    # os.chmod(convgraph_sh,744)

    text=r"""#!/usr/bin/perl
use strict;
my ($line, $n, $nv, $nedges, $ne);
my @fields;
print "0\n";
$line=<STDIN>;
chomp($line); #removes any trailing string corresponding to the record separator
($nv,$nedges)=split/ /,$line;
print "$nv $nedges\n";
print "1 000\n";
while(<STDIN>) {
    chomp($_);
    @fields=split/ /,$_;
    $ne=$#fields+1;
    print "$ne @fields\n";
}
"""

    convgraph_pl = os.path.join(tmpdir,'convgraph.pl')
    print('...writing file:',convgraph_pl)
    tmpfl = open(convgraph_pl,'w')
    # print >> tmpfl, text
    tmpfl.write( text )
    tmpfl.close()
    try:
        os.chmod(convgraph_pl,stat.S_IEXEC)
    except:
        pass
        # os.chmod(convgraph_pl,0744)
    # os.chmod(convgraph_pl,744)

    return convgraph_sh

def write_mkflgown_pl(tmpdir):

    text=r"""#!/usr/bin/perl
use strict;
my @nf=(); # my: declares the variables in LIST to be scoped within the block
my @nw=();
my @ni=();
my @no=();
my @fields=();
my ($bgkline,$ownrline,$i4,$p,$c,$n,$flag);
$c=0;
for($p=0; $p<$ARGV[0]; $p++) {
  $nf[$p]=$nw[$p]=$ni[$p]=$no[$p]=0;
}
# die "Usage: $0 number_of_tasks bgkflag [nodeownr]" if(@ARGV<2);
die "Usage: $0 number_of_tasks bgkflag [nodeownr]" if(@ARGV<3);

# open(BGKFLAG,"<$ARGV[1]") || die "Could not open $ARGV[1]: $!";
open(BGKHDR,"<$ARGV[1]") || die "Could not open $ARGV[1]: $!";
open(BGKDAT,"<$ARGV[2]") || die "Could not open $ARGV[2]: $!";

# $bgkline=<BGKFLAG>;
$bgkline=<BGKHDR>;
chomp($bgkline);
$bgkline=~s/^\s+//;
my ($nx, $ny, $nz) = split/\s+/,$bgkline;

# $bgkline=<BGKFLAG>;
$bgkline=<BGKHDR>;

# if(@ARGV>2) {
if(@ARGV>3) {
  # open(NODEOWNR,"<$ARGV[2]") || die "Could not open $ARGV[2]: $!";
  open(NODEOWNR,"<$ARGV[3]") || die "Could not open $ARGV[3]: $!";
  $n="bgkflag_" . $c . ".dat";
} else {
  $p=0;
  $n="bgkflag.dat";
}
open(BGKFLAGP,">>$n") || die "Could not open $n: $!";
# while(<BGKFLAG>) {
while(<BGKDAT>) {
  chomp($_);
  $bgkline=$_;
  $bgkline=~s/^\s+//;
  @fields=split/\s+/,$bgkline;
  $flag=$fields[$#fields];
  # if(@ARGV>2) {
  if(@ARGV>3) {
    $ownrline=<NODEOWNR>;
    chomp($ownrline);
    ($i4, $p)=split/\s+/,$ownrline;
    if($p<0 || $p >=$ARGV[0]) {
      die "Invalid task: $p\n";
    }
  }
  if($flag == 1) {
    $nf[$p]++;
  } elsif($flag == 2) {
    $nw[$p]++;
  } elsif($flag == 3) {
    $ni[$p]++;
  } elsif($flag == 4) {
    $no[$p]++;
  } elsif($flag == 5) {
  } else {
    die "Invalid flag value: $flag\n";
  }

  if($p != $c) {
    close(BGKFLAGP);
    $n="bgkflag_" . $p . ".dat";
    open(BGKFLAGP,">>$n") || die "Could not open $n: $!";
    $c=$p;
  }
  print BGKFLAGP "$bgkline\n";
}
for($p=0; $p<$ARGV[0]; $p++) {
  # if(@ARGV>2) {
  if(@ARGV>3) {
    $n="bgkflag_" . $p . ".hdr";
  } else {
    $n="bgkflag.hdr";
  }
  open(BGKFLAGH,">>$n") || die "Could not open $n: $!";
  print BGKFLAGH "$nx $ny $nz\n";
  print BGKFLAGH "5 $nf[$p] $nw[$p] $ni[$p] $no[$p]\n";
  close BGKFLAGH;
}
"""
    mkflgown_pl = os.path.join(tmpdir,'mkflgown.pl')
    print('...writing file:',mkflgown_pl)
    tmpfl = open(mkflgown_pl,'w')
    # print >> tmpfl, text
    tmpfl.write( text )
    tmpfl.close()
    # os.chmod(mkflgown_pl,0744)
    try:
        os.chmod(mkflgown_pl,stat.S_IEXEC)
    except:
        pass
        # os.chmod(mkflgown_pl,0744)

    return mkflgown_pl

def write_makeownr_sh(tmpdir):

    text=r"""#!/bin/bash
if [ $# -ne 3 ]
then
    echo Usage $0 partition i4 targetfile
    exit 1
fi
part=$1
i4file=$2
tf=$3
tail -n +2 $part | sort -n -T. | cut -f 2 | paste $i4file - | sort -T. -n > $tf
"""

    makeownr_sh = os.path.join(tmpdir,'makeownr.sh')
    print('...writing file:',makeownr_sh)
    tmpfl = open(makeownr_sh,'w')
    # print >> tmpfl, text
    tmpfl.write( text )
    tmpfl.close()
    try:
        os.chmod(makeownr_sh,stat.S_IEXEC)
    except:
        pass
        # os.chmod(makeownr_sh,0744)
    # os.chmod(makeownr_sh,744)

    return makeownr_sh

def write_goscotch_sh(tmpdir):

    if Magic.MUPHYorMOBIUS == 'Muphy':
        executable = os.path.join(os.environ.get('MUPHY_ROOT'), 'EXECUTE', 'dgpart')

    elif Magic.MUPHYorMOBIUS == 'Moebius':
        executable = os.path.join(os.environ.get('MOEBIUS_ROOT'), 'BACKEND', 'SHOP', 'dgpart')
    else:
        print ('var undefined', Magic.MUPHYorMOBIUS)
        sys.exit(1)

    text=r"""#!/bin/bash

if [ $# -ne 3 ]
then
    echo Usage $0 n inputfile outputfile
    exit 1
fi
n=$1
input=$2
output=$3

# mpirun -np 4 -hostfile ./hostlist ./dgpart $n $input $output
# ./dgpart $n $input $output
%s $n $input $output
""" % executable

# $MUPHY_BIN/dgpart $n $input $output

    goscotch_sh = os.path.join(tmpdir,'goscotch.sh')
    print('...writing file:',goscotch_sh)
    tmpfl = open(goscotch_sh,'w')
    # print >> tmpfl, text
    tmpfl.write( text )
    tmpfl.close()
    try:
        os.chmod(goscotch_sh,stat.S_IEXEC)
    except:
        pass
        # os.chmod(goscotch_sh,0744)
    # os.chmod(goscotch_sh,744)

    return goscotch_sh


def splitParallelDomainsEqualSlabs(nproc=1, splitdir='x', infileroot='bgkflag'):
    """
    Method to split a mesh file into x, y or z regions with almost-equal split zones.
    It does not guarantee the same number of voxels in each slab.
    Input:  *.hdr/.dat mesh file
    Output: *_PE.hdr/*_PE.dat files and
            nodeownr.inp
    """

    print ('Splitting:', nproc,'parallel domains', \
          'along:', splitdir, 'direction', \
          'for mesh files:', infileroot, '.*')

    h = open(infileroot + '.hdr', 'r')
    l = h.readline()
    l = l.split()
    nx,ny,nz = int(l[0]),int(l[1]),int(l[2])
    h.close()

    n = open('nodeownr.inp','w')

    nodes = {}
    for i in range(nproc):
        nodes[i] = []

    iproc = 0
    mesh = {}
    nn = 0
    imn,jmn,kmn = +1e6,+1e6,+1e6
    imx,jmx,kmx = -1e6,-1e6,-1e6

    d = open(infileroot + '.dat', 'r')
    for line in d.readlines():
        l = line.split()
        i = int(l[0]); j = int(l[1]); k = int(l[2])
        imn, imx = min(i,imn), max(i,imx)
        jmn, jmx = min(j,jmn), max(j,jmx)
        kmn, kmx = min(k,kmn), max(k,kmx)
    d.close()
    print ('imn,imx:',imn,imx)
    print ('jmn,jmx:',jmn,jmx)
    print ('kmn,kmx:',kmn,kmx)

    d = open(infileroot + '.dat', 'r')

    for line in d.readlines():

        l = line.split()
        i = int(l[0]); j = int(l[1]); k = int(l[2])

        flg = int(l[3])

        if splitdir == 'x':
            iproc = (i-imn) / (nx/nproc )

        elif splitdir == 'y':
            iproc = (j-jmn) / (ny/nproc )

        elif splitdir == 'z':
            iproc = (k-kmn) / (nz/nproc )

        else:
            print ('Unrecognized split direction :',splitdir)

        i4 = (nx+2)*(ny+2)*k + (nx+2)*j + i
        print >> n, i4, iproc
        nodes[iproc].append([i,j,k,flg])
        nn += 1

    n.close()
    d.close()

    label = ['C','N','O','H','X','Y','K','I']
    fxyz = open('msh.xyz','w')
    print >> fxyz, nn
    print >> fxyz
    for iproc in range(nproc):

        nfl,nwl,nin,nou = 0,0,0,0
        for node in nodes[iproc]:
            i,j,k,flg = node
            if flg == 1: nfl += 1
            if flg == 2: nwl += 1
            if flg == 3: nin += 1
            if flg == 4: nou += 1

        bh = open( infileroot + '_' + str(iproc) + '.hdr', 'w')
        print >> bh,nx,ny,nz
        print >> bh, 5, nfl,nwl,nin,nou
        bh.close()

        bd = open( infileroot + '_' + str(iproc) + '.dat', 'w')
        for node in nodes[iproc]:
            i,j,k,flg = node
            print >> bd,i,j,k,flg
            print >> fxyz,label[iproc],i,j,k

