#!/usr/bin/env python

import os,sys

# sys.path.append(os.environ["MUPHY_ROOT"]+"/GMUPHY")
# sys.path.append(os.environ["MUPHY_ROOT"]+"/ECOMUPHY")
# sys.path.append(os.environ["MUPHY_ROOT"]+"/ECOMUPHY/ICONS")

from ctypes import *
import pickle
import math
import argparse

"""
module for muphy1
"""

try:
    from muphyserver import *
    import mxqueue
    SOCKET=True    
except ImportError:
    SOCKET=False    

try: # c_bool defined from python 2.5 on
        c_bool()
except:
        c_bool=c_int

class MuphyP:
    """
    main muphy class
    """

    global SOCKET

    def fill_muphy_object(self,version='1',xpu='cpu',mycomm=0):
        """
        various starting options
        """

        self.version = version
        self.xpu = xpu

        EXEDIR = os.environ.get('MOEBIUS_ROOT')
        #EXEDIR = os.environ.get('MUPHY_BIN')

        THELIB = EXEDIR + 'BACKEND/SHOP/libmoebius'

        if(self.xpu != 'cpu' and self.xpu != 'gpu'):
            print 'xpu:',xpu,'...not of the right kind!'
            sys.exit(1)

        if xpu == 'gpu':
            THELIB = THELIB + '.gpu'

        THELIB = THELIB + '.so'

        if not os.path.isfile(THELIB):
            print 'libray does not exist:',THELIB
            sys.exit(1)

        self.library = THELIB

        # lib1 = cdll.LoadLibrary("/usr/local/cuda/lib/libcudart.dylib")
        # lib1 = CDLL("/usr/local/cuda/lib/libcudart.dylib")
        # lib2 = cdll.LoadLibrary("/usr/local/cuda/lib/libcurand.dylib")
        # lib2 = CDLL("/usr/local/cuda/lib/libcurand.dylib")
        
        self.MUPHY = CDLL(self.library, RTLD_GLOBAL)
        # self.MUPHY = cdll.LoadLibrary(self.library)

    def check_parallelism(self):
        """
        inquire library and set the var self.PARALLEL
        """

        self.PARALLEL = False

        if get_numprocs() > 1:
            self.PARALLEL = True
            return True

    def set_socket_params(self):
        """
        set parameters for socket
        """

        if not SOCKET: return

        self.remote_obj = {}

        self.remote_obj['get_fluid_sample_points']         = ['Sampler',['points'],['1']]
        self.remote_obj['get_fluid_sample_points_density'] = ['Sample' ,['scalar'],['1']]
        self.remote_obj['get_fluid_sample_points_velocity']= ['Sample' ,['vector'],['1']]

        self.remote_obj['get_fluid_map_density']           = ['Map'    ,['scalar'],['none']]
        self.remote_obj['get_fluid_map_vel_x']             = ['Map'    ,['scalar'],['none']]
        self.remote_obj['get_fluid_map_velmod']            = ['Map'    ,['scalar'],['none']]

        self.remote_obj['get_fluid_momentum_total']        = ['xy'     ,['vector'],['none']]
        self.remote_obj['get_fluid_number_total']          = ['xy'     ,['scalar'],['none']]
        self.remote_obj['get_total_charge']                = ['xy'     ,['scalar'],['none']]
        self.remote_obj['get_particles_pot_energy']        = ['xy'     ,['scalar'],['none']]
        self.remote_obj['get_particles_kin_energy']        = ['xy'     ,['scalar'],['none']]
        self.remote_obj['get_fluid_density_profile']       = ['Nxy'    ,['scalar'],['none']]
        self.remote_obj['get_fluid_vel_profile']           = ['Nxy'    ,['vector'],['none']]
        self.remote_obj['get_fluid_current_profile']       = ['Nxy'    ,['vector'],['none']]
        self.remote_obj['get_fluid_density_centerline']    = ['Nxy'    ,['scalar'],['1','2','3','4']]
        self.remote_obj['get_fluid_velocity_centerline']    = ['Nxy'    ,['scalar'],['1','2','3','4']]

# instance the muphy object
M = MuphyP()

##########################################

# def muphywrapper_init(xpu='cpu',mycomm=0):
# def muphywrapper_init(xpu='cpu',mycomm=0):
def muphywrapper_init(*args):
    """
    muphy init: run parser and do socket ops
    """

    global SOCKET

    XPU = 'cpu'
    mycomm = 0

    for a in args:
        if a == 'cpu':
            XPU = 'cpu'
        elif a == 'gpu':
            XPU = 'gpu'

    # the command args will override 
    parser = argparse.ArgumentParser()

    parser.add_argument('-x', '--xpu',  default='cpu', help='device [cpu/gpu]')
    parser.add_argument('-s', '--socket',  default='no', help='socket port active [yes/no]')
    args = parser.parse_args()

    XPU = args.xpu

    M.fill_muphy_object(version='1',xpu=XPU)

    if SOCKET:
        M.set_socket_params()

    str = ' Xpu: '+XPU+'\n'+ \
        ' Muphy shared library: '+M.library+'\n' + \
        ' Muphy version: '+M.version

    wrapper_init = M.MUPHY.muphywrapper_init(mycomm,len(str),str)

    r = M.check_parallelism()

    if r: return

    if SOCKET:
        try:
            M.req_q = mxqueue.MXQueue()
            sock = MySocket(M.req_q)
        except:
            print 'socket already in use....not connecting'
            SOCKET=False

def tracker_finish(): 
    """
    tracker ending step
    """
    M.MUPHY.tracker_finish()

def muphywrapper_finish(): 
    """
    muphy ending step
    """
    wrapper_finish = M.MUPHY.muphywrapper_finish()
    print 'Muphy ends, bye.'

##########################################
def get_itime():
    """
    get initial and final timestep
    """

    itime_start = (c_int*1)()
    ncycle = (c_int*1)()
    dummy = (c_int*1)()

    M.MUPHY.get_itime(byref(itime_start),byref(ncycle),byref(dummy))

    return itime_start[0], ncycle[0]

##########################################
def xchange_fluid_cube_sampler(line):
    """
    get fluid cube sampler number and points
    """

    line = line.split(',')

    o = [float(line[1]),float(line[2]),float(line[3])]
    u = [float(line[4]),float(line[5]),float(line[6])]
    v = [float(line[7]),float(line[8]),float(line[9])]
    w = [float(line[10]),float(line[11]),float(line[12])]

    n_l = (c_int*1)()
    o_l = (c_float*3)()
    u_l = (c_float*3)()
    v_l = (c_float*3)()
    w_l = (c_float*3)()

    o_l[0] = o[0]; o_l[1] = o[1]; o_l[2] = o[2]
    u_l[0] = u[0]; u_l[1] = u[1]; u_l[2] = u[2]
    v_l[0] = v[0]; v_l[1] = v[1]; v_l[2] = v[2]
    w_l[0] = w[0]; w_l[1] = w[1]; w_l[2] = w[2]

    M.MUPHY.get_fluid_cube_sampler_number(o_l,u_l,v_l,w_l,n_l)
    n = n_l[0]

    px_l = (c_float*n)()
    py_l = (c_float*n)()
    pz_l = (c_float*n)()

    M.MUPHY.get_fluid_cube_sampler_points(byref(n_l), \
                                        byref(o_l), byref(u_l), byref(v_l), byref(w_l), \
                                        byref(px_l), byref(py_l), byref(pz_l))

    px=[]; py=[]; pz=[]
    for i in range(n):
        px.append(px_l[i])
        py.append(py_l[i])
        pz.append(pz_l[i])

    return n,px,py,pz

##########################################
def get_itime3():
    """
    get initial,final and current timestep
    """

    itime_start = (c_int*1)()
    ncycle = (c_int*1)()
    itime_current = (c_int*1)()

    M.MUPHY.get_itime(byref(itime_start),byref(ncycle),byref(itime_current))

    return itime_start[0], ncycle[0], itime_current[0]

##########################################
def get_myproc():
    """
    get my processor
    """

    myproc = (c_int*1)()

    M.MUPHY.get_myproc(byref(myproc))

    return myproc[0]

##########################################
def get_numprocs():
    """
    get number of processors
    """

    numprocs = (c_int*1)()

    M.MUPHY.get_numprocs(byref(numprocs))

    return numprocs[0]

##########################################
def set_itime(itime): 
    """
    set current time
    """
    M.MUPHY.set_itime(byref(c_int(itime)))
    
def timestep_bgk_init(): 
    """
    init single fluid 
    """
    M.MUPHY.timestep_bgk_init()

def timestep_mixture_init(): 
    """
    init multi fluids
    """
    M.MUPHY.timestep_mixture_init()

def timestep_molecular_init(): 
    """
    init molecules
    """
    M.MUPHY.timestep_molecular_init()

def timestep_bgk(): 
    """
    propagate one step for single fluid
    """
    M.MUPHY.timestep_bgk()

def timestep_mixture(): 
    """
    propagate one step for multi fluid
    """
    M.MUPHY.timestep_mixture()

def timestep_molecular(): 
    """
    propagate one step for multi molecules
    """
    M.MUPHY.timestep_molecular()

##########################################
def tracker():
    """
    track every actor
    """
    global SOCKET

    M.MUPHY.tracker()

    if M.PARALLEL: return
    if not SOCKET: return

    # queue handling
    reqs = M.req_q.get_queue()

    if not reqs: return

    for rstring,q in reqs:

        # the request string is either a single name or a 
        # comma-separated list of names: var,type,rank component,reduction
        rlist = rstring.split(',')
        r = rlist[0]
        
        if (r == 'REMOTE_OBJS'):
            serialized = pickle.dumps(M.remote_obj)
            q.enqueue(serialized)
            continue

        '''
        if (r == 'get_all_remote_vars'):    # old style to send back a string
            response = M.all_remote_vars
            q.enqueue(response)
            continue
        '''

        if (r[0:26] == 'xchange_fluid_cube_sampler'):
            databody = xchange_fluid_cube_sampler(rstring)
            serialized = pickle.dumps(databody)
            q.enqueue(serialized)
            continue

        elif (r == 'get_time'):
            iii = get_itime3()
            response = str(iii[2])
            q.enqueue(response)
            continue

        if len(rlist)==1: continue

        databody = []
        for (ff,ttt) in M.remote_obj.iteritems():

            tt = ttt[0]

            var = type = rankcomp = reduce = None
            if len(rlist)==1:
                var = rlist[0]
            else:
                var, type, rankcomp, reduce = rlist

            if var == ff:

                # get the method from the namespace
                methodToCall = globals()[var]

                databody.append(tt) # message envelope: data type

                if tt == 'xy':

                    # convention of message: [step,...data...]
                    databody.append(get_itime3()[2])    # timestep
                    databody.append(methodToCall(type,rankcomp,reduce))        # value

                elif tt == 'Nxy':

                    ax,ay = methodToCall(type,rankcomp,reduce)

                    # convention of message: [...data...]
                    for (x,y) in zip(ax,ay): 
                        databody.append([x,y])

                elif tt == 'Map':

                    # convention of message: [n1,n2,...data...]
                    ax,ay,vals = methodToCall(type,rankcomp,reduce)
                    n1 = len(ax)
                    n2 = len(ay)

                    databody.append([n1,n2])
                    for i1 in range(n1): 
                        for i2 in range(n2): 
                            databody.append(vals[i1*n2 + i2])

                elif tt == 'Sampler':

                    # convention of message: [size,points_vector]
                    n,px,py,pz,np,i1,i2,i3 = methodToCall(type,rankcomp,reduce)

                    # points
                    databody.append(n)
                    for i in range(n): 
                        databody.append([px[i],py[i],pz[i]])

                    # triangles
                    databody.append(np)
                    for i in range(np): 
                        databody.append([i1[i],i2[i],i3[i]])

                elif tt == 'Sample':

                    # convention of message: [size,...data...]
                    n,v = methodToCall(type,rankcomp,reduce)

                    databody.append(n)
                    for i in range(n): 
                        databody.append(v[i])

                # pack serialized response
                serialized = pickle.dumps(databody)

                q.enqueue(serialized)

                break # handle only one request

##########################################
def set_restart(flag):
    """
    set flag for restarting
    """

    cflag = c_bool(flag)

    M.MUPHY.set_restart(cflag)

##########################################
def set_surface_stress(flag):
    """
    set flag for computing surface stress
    """

    cflag = c_bool(flag)

    M.MUPHY.set_surface_stress(cflag)

##########################################
def scale_prmvdw(factor):
    """
    scale nonbonding parameters by some factor
    """

    cfactor = c_float(factor)

    M.MUPHY.scale_prmvdw(cfactor)

##########################################
def scale_prmbnd(factor):
    """
    scale bonding parameters by some factor
    """

    cfactor = c_float(factor)

    M.MUPHY.scale_rdatcon_relativefactor(cfactor)

##########################################
def scale_prmang(factor):
    """
    scale angular parameters by some factor
    """

    cfactor = c_float(factor)

    M.MUPHY.scale_rdatang_relativefactor(cfactor)

##########################################
def scale_prmdih(factor):
    """
    scale dihedral parameters by some factor
    """

    cfactor = c_float(factor)

    M.MUPHY.scale_rdatdih_relativefactor(cfactor)


##########################################
def set_cap_force(flag,f1,f2,f3,f4):
    """
    cap forces by max absolute vals
    """

    cflag = c_bool(flag)
    cf1 = c_float(f1)
    cf2 = c_float(f2)
    cf3 = c_float(f3)
    cf4 = c_float(f4)

    M.MUPHY.set_cap_force(cflag,cf1,cf2,cf3,cf4)


##########################################
def set_cap_force_fluid(flag,f1,f2):
    """
    cap fluid force by max absolute vals
    """

    cflag = c_bool(flag)
    cf1 = c_float(f1)
    cf2 = c_float(f2)

    M.MUPHY.set_cap_force_fluid(cflag,cf1,cf2)

#########################################
def set_cap_force_fluid_BGK(flag,f1):
    """
    cap fluid force of single fluid by max absolute vals
    """

    cflag = c_bool(flag)
    cf1 = c_float(f1)

    M.MUPHY.set_cap_force_fluid_BGK(cflag,cf1)

##########################################
def set_bennett(flag):
    """
    flag to activate bennett minimization for particles
    """

    cflag = c_bool(flag)

    M.MUPHY.set_bennett(cflag)


##########################################
def scale_velocity(factor):
    """
    factor to scale velocities for particles
    """

    cfactor = c_float(factor)

    M.MUPHY.scale_velocity(cfactor)


##########################################
def mergedumps(numproc):
    """
    merge dumps of different processors
    """

    M.MUPHY.mergedumps(c_int(numproc))

##########################################
def tracker_minimal_serial():
    """
    minimal tracking of a serial run
    """

    M.MUPHY.tracker_minimal_serial()

##########################################
def stop_mixture_com(isp):
    """
    stop center of mass of multi fluid
    """

    cisp = c_int(isp)

    M.MUPHY.stop_com_mix_specie(cisp)

##########################################
def get_box():
    """
    get simulation box
    """

    nx = (c_int*1)()
    ny = (c_int*1)()
    nz = (c_int*1)()

    M.MUPHY.get_box(byref(nx),byref(ny),byref(nz))

    return nx[0],ny[0],nz[0]

##########################################
def get_fluid_nodesize():
    """
    get number of fluid nodes
    """

    nf = (c_int*1)()
    M.MUPHY.get_fluid_nodesize(byref(nf))
    return nf[0]

##########################################
def get_wall_nodesize():
    """
    get number of wall nodes
    """

    nf = (c_int*1)()
    M.MUPHY.get_wall_nodesize(byref(nf))
    return nf[0]

##########################################
def get_inlet_nodesize():
    """
    get number of inlet nodes
    """

    nf = (c_int*1)()
    M.MUPHY.get_inlet_nodesize(byref(nf))
    return nf[0]

##########################################
def get_outlet_nodesize():
    """
    get number of outlet nodes
    """

    nf = (c_int*1)()
    M.MUPHY.get_outlet_nodesize(byref(nf))
    return nf[0]

##########################################
def get_particle_number():
    """
    get number of particles on processor
    """

    natms = (c_int*1)()

    M.MUPHY.get_particle_number(byref(natms))

    return natms[0]

##########################################
def get_particle_number_total():
    """
    get global number of particles
    """

    natms_tot = (c_int*1)()

    M.MUPHY.get_particle_number_total(byref(natms_tot))

    return natms_tot[0]

##########################################
def get_particle_uid():
    """
    get universal ids of particles on processor
    """

    n = get_particle_number()

    uid_l = (c_int*n)()
   
    M.MUPHY.get_particle_uid(c_int(n), byref(uid_l))
   
    return uid_l

##########################################
def get_particle_mass():
    """
    get mass of particles on processor
    """

    n = get_particle_number()

    mass_l = (c_float*n)()
   
    M.MUPHY.get_particle_uid(c_int(n), byref(mass_l))
   
    return mass_l

##########################################
def get_particle_freeze():
    """
    get freeze flag of particles on processor
    """

    n = get_particle_number()

    freeze_l = (c_int*n)()
   
    M.MUPHY.get_particle_freeze(c_int(n), byref(freeze_l))
   
    return freeze_l

##########################################
def get_particle_position():
    """
    get position of particles on processor
    """

    n = get_particle_number()

    px_l = (c_float*n)()
    py_l = (c_float*n)()
    pz_l = (c_float*n)()
   
    M.MUPHY.get_particle_position(c_int(n), byref(px_l), byref(py_l), byref(pz_l))
   
    return px_l,py_l,pz_l

##########################################
def get_particle_velocity():
    """
    get velocity of particles on processor
    """

    n = get_particle_number()

    px_l = (c_float*n)()
    py_l = (c_float*n)()
    pz_l = (c_float*n)()
   
    M.MUPHY.get_particle_velocity(c_int(n), byref(px_l), byref(py_l), byref(pz_l))

    return px_l,py_l,pz_l

##########################################
def get_particle_force():
    """
    get force of particles on processor
    """

    n = get_particle_number()

    px_l = (c_float*n)()
    py_l = (c_float*n)()
    pz_l = (c_float*n)()

    M.MUPHY.get_particle_force(c_int(n), byref(px_l), byref(py_l), byref(pz_l))

    return px_l,py_l,pz_l

##########################################
def broadcast_array(iproc,p):
    """
    broadcast array to all processors
    """

    n = len(p)

    if isinstance(p[0],float):

        p_l = (c_float * len(p))(*p)
        M.MUPHY.broadcast_array_float(c_int(iproc), c_int(n), byref(p_l))
        return p_l

    elif isinstance(p[0],int):

        p_l = (c_int * len(p))(*p)
        M.MUPHY.broadcast_array_int(c_int(iproc), c_int(n), byref(p_l))
        return p_l

    else:
        return None

##########################################
def set_fluid_force(fx,fy,fz):
    """
    set force on fluid
    """

    M.MUPHY.set_fluid_force(c_float(fx), c_float(fy), c_float(fz))
   
##########################################
def get_fluid_force():
    """
    get force on fluid
    """

    fx_l = (c_float*1)()
    fy_l = (c_float*1)()
    fz_l = (c_float*1)()

    M.MUPHY.get_fluid_force(byref(fx_l), byref(fy_l), byref(fz_l))

    return fx_l[0],fy_l[0],fz_l[0]

##########################################
def set_particle_position(px,py,pz):
    """
    set position of particles
    """

    n = len(px)

    #px_l = (c_float*n)()
    #py_l = (c_float*n)()
    #pz_l = (c_float*n)()
    #for i in range(n):
    #   px_l[i] = px[i]
    #   py_l[i] = py[i]
    #   pz_l[i] = pz[i]
    px_l = (c_float * len(px))(*px)
    py_l = (c_float * len(px))(*py)
    pz_l = (c_float * len(px))(*pz)
   
    M.MUPHY.set_particle_position(c_int(n), byref(px_l), byref(py_l), byref(pz_l))
   
##########################################
def set_particle_velocity(px,py,pz):
    """
    set velocity of particles
    """

    n = len(px)

    #px_l = (c_float*n)()
    #py_l = (c_float*n)()
    #pz_l = (c_float*n)()
    #for i in range(n):
    #   px_l[i] = px[i]
    #   py_l[i] = py[i]
    #   pz_l[i] = pz[i]
    px_l = (c_float * len(px))(*px)
    py_l = (c_float * len(px))(*py)
    pz_l = (c_float * len(px))(*pz)
   
    M.MUPHY.set_particle_velocity(c_int(n), byref(px_l), byref(py_l), byref(pz_l))

##########################################
def set_particle_force(px,py,pz):
    """
    set force of particles
    """

    n = len(px)

    #px_l = (c_float*n)()
    #py_l = (c_float*n)()
    #pz_l = (c_float*n)()
    #for i in range(n):
    #   px_l[i] = px[i]
    #   py_l[i] = py[i]
    #   pz_l[i] = pz[i]
    px_l = (c_float * len(px))(*px)
    py_l = (c_float * len(px))(*py)
    pz_l = (c_float * len(px))(*pz)
   
    M.MUPHY.set_particle_force(c_int(n), byref(px_l), byref(py_l), byref(pz_l))

##########################################
def set_particle_freeze(ipf):
    """
    set freeze flag of particles
    """

    n = len(ipf)

    ipf_l = (c_int * len(ipf))(*ipf)
   
    M.MUPHY.set_particle_freeze(c_int(n), byref(ipf_l))
   
##########################################
def get_fluid_density(n,px,py,pz):
    """
    get density of fluid on given point
    """

    #   # numpy way
    #   rho = empty(n, dtype="float32")
    #   M.MUPHY.sample_fluid_density(c_int(n),
    #                              px.ctypes.data_as(POINTER(c_float)), 
    #                              py.ctypes.data_as(POINTER(c_float)), 
    #                              pz.ctypes.data_as(POINTER(c_float)),
    #                              rho.ctypes.data_as(POINTER(c_float)))

    # px,py,pz already allocated. Just bind to them.
    px_l = (c_float * len(px))(*px)
    py_l = (c_float * len(px))(*py)
    pz_l = (c_float * len(px))(*pz)
    # allocate rho
    # rho = [1]*len(px)
    # bind to it
    # rho_l = (c_float * len(px))(*rho) # this does not work on exit from function !!!
    rho_l = (c_float * len(px))()

    M.MUPHY.sample_fluid_density(c_int(n),byref(px_l), byref(py_l), byref(pz_l), byref(rho_l))

    return rho_l

##########################################
def get_fluid_velocity(n,px,py,pz):
    """
    get velocity of fluid on given point
    """

    # px,py,pz already allocate. Just bind to them.
    px_l = (c_float * len(px))(*px)
    py_l = (c_float * len(px))(*py)
    pz_l = (c_float * len(px))(*pz)

    # # allocate vels # the following does  ot work on function exit !!!!
    # u = [0]*len(px)
    # v = [0]*len(px)
    # w = [0]*len(px)
    # # bind to them
    # u_l = (c_float * len(px))(*u)
    # v_l = (c_float * len(px))(*v)
    # w_l = (c_float * len(px))(*w)
    u_l = (c_float * len(px))()
    v_l = (c_float * len(px))()
    w_l = (c_float * len(px))()

    M.MUPHY.sample_fluid_velocity(c_int(n),
                              byref(px_l), byref(py_l), byref(pz_l), 
                              byref(u_l), byref(v_l), byref(w_l))

    # return u,v,w
    return u_l,v_l,w_l

##########################################
def get_fluid_vorticity(n,px,py,pz):
    """
    get vorticity of fluid on given point
    """

    # px,py,pz already allocate. Just bind to them.
    px_l = (c_float * len(px))(*px)
    py_l = (c_float * len(px))(*py)
    pz_l = (c_float * len(px))(*pz)

    # # allocate vels # the following does  ot work on function exit !!!!
    # u = [0]*len(px)
    # v = [0]*len(px)
    # w = [0]*len(px)
    # # bind to them
    # u_l = (c_float * len(px))(*u)
    # v_l = (c_float * len(px))(*v)
    # w_l = (c_float * len(px))(*w)
    u_l = (c_float * len(px))()
    v_l = (c_float * len(px))()
    w_l = (c_float * len(px))()

    M.MUPHY.sample_fluid_vorticity(c_int(n),
                              byref(px_l), byref(py_l), byref(pz_l), 
                              byref(u_l), byref(v_l), byref(w_l))

    # return u,v,w
    return u_l,v_l,w_l

##########################################
def get_fluid_velocity_on_mesh():
    """
    get velocity of fluid on all mesh
    """

    nx,ny,nz = get_box()

    ss = (nx+2)*(ny+2)*(nz+2)
    # ss = nx*ny*nz

    # u = [0]*ss
    # v = [0]*ss
    # w = [0]*ss
    # u_l = (c_float * ss)(*u)
    # v_l = (c_float * ss)(*v)
    # w_l = (c_float * ss)(*w)

    u_l = (c_float * ss)()
    v_l = (c_float * ss)()
    w_l = (c_float * ss)()

    M.MUPHY.get_fluid_velocity_on_mesh(c_int(ss),byref(u_l), byref(v_l), byref(w_l))

    # for i in range(ss):
    #    u[i] = u_l[i]
    #    v[i] = v_l[i]
    #    w[i] = w_l[i]

    # return u,v,w # binding should be already sufficient
    return u_l,v_l,w_l

##########################################
def get_fluid_density_on_mesh():
    """
    get density of fluid on all mesh
    """

    nx,ny,nz = get_box()

    ss = (nx+2)*(ny+2)*(nz+2)
    # ss = nx*ny*nz

    rho_l = (c_float * ss)()

    M.MUPHY.get_fluid_density_on_mesh(c_int(ss),byref(rho_l))

    return rho_l

##########################################
def get_fluid_density_profile(*args):
    """
    get density of fluid as profile
    """

    nx,ny,nz = get_box()

    # this is a cut in the z-plane at midplane !!! to be generalize
    u_l = (c_float*nz)()

    z_l=[]
    for i in range(0,nz): z_l.append(i)

    if(len(args)>2): 

        M.MUPHY.get_fluid_density_profile(c_int(nz), byref(u_l))

    else:
        M.MUPHY.get_fluid_density_profile(c_int(nz), byref(u_l))

    return z_l,u_l

##########################################
def get_fluid_vel_profile(*args):
    """
    get velocity of fluid as profile
    """

    nx,ny,nz = get_box()

    # this is a cut in the z-plane at midplane !!! to be generalize
    u_l = (c_float*nz)()

    z_l=[]
    for i in range(0,nz): z_l.append(i)

    if(len(args)>2): 

        if args[1]=='v.x':
            M.MUPHY.get_fluid_vel_profile(c_int(1), c_int(nz), byref(u_l))

        elif args[1]=='v.y':
            M.MUPHY.get_fluid_vel_profile(c_int(2), c_int(nz), byref(u_l))

        elif args[1]=='v.z':
            M.MUPHY.get_fluid_vel_profile(c_int(3), c_int(nz), byref(u_l))

        elif args[1]=='|v|' or args[1]=='|v|^2':
            M.MUPHY.get_fluid_vel_profile(c_int(1), c_int(nz), byref(u_l))

            v_l = (c_float*nz)()
            M.MUPHY.get_fluid_vel_profile(c_int(2), c_int(nz), byref(v_l))

            w_l = (c_float*nz)()
            M.MUPHY.get_fluid_vel_profile(c_int(3), c_int(nz), byref(w_l))

            for i in range(len(u_l)):
                u_l[i] = u_l[i]*u_l[i] + v_l[i]*v_l[i] + w_l[i]*w_l[i]

            if args[1]=='|v|' :
                for i in range(len(u_l)):
                    u_l[i] = math.sqrt(u_l[i])

    else:
        M.MUPHY.get_fluid_vel_profile(c_int(1), c_int(nz), byref(u_l))

    return z_l,u_l


##########################################
def get_fluid_current_profile(*args):
    """
    get current of fluid as profile
    """

    nx,ny,nz = get_box()

    # this is a cut in the z-plane at midplane !!! to be generalize
    u_l = (c_float*nz)()

    z_l=[]
    for i in range(0,nz): z_l.append(i)

    if(len(args)>2): 

        if args[1]=='v.x':
            M.MUPHY.get_fluid_current_profile(c_int(1), c_int(nz), byref(u_l))

        elif args[1]=='v.y':
            M.MUPHY.get_fluid_current_profile(c_int(2), c_int(nz), byref(u_l))

        elif args[1]=='v.z':
            M.MUPHY.get_fluid_current_profile(c_int(3), c_int(nz), byref(u_l))

        elif args[1]=='|v|' or args[1]=='|v|^2':
            M.MUPHY.get_fluid_current_profile(c_int(1), c_int(nz), byref(u_l))

            v_l = (c_float*nz)()
            M.MUPHY.get_fluid_current_profile(c_int(2), c_int(nz), byref(v_l))

            w_l = (c_float*nz)()
            M.MUPHY.get_fluid_current_profile(c_int(3), c_int(nz), byref(w_l))

            for i in range(len(u_l)):
                u_l[i] = u_l[i]*u_l[i] + v_l[i]*v_l[i] + w_l[i]*w_l[i]

            if args[1]=='|v|' :
                for i in range(len(u_l)):
                    u_l[i] = math.sqrt(u_l[i])

    else:
        M.MUPHY.get_fluid_current_profile(c_int(1), c_int(nz), byref(u_l))

    return z_l,u_l


##########################################
def get_fluid_density_centerline(*args):
    """
    get density of fluid on centerline
    """

    if(len(args)<3): return 

    n_l = (c_int*1)()
    idfile = int(args[2])

    M.MUPHY.get_fluid_centerline_size(c_int(idfile),byref(n_l))

    n = n_l[0]

    # this is a cut in the z-plane at midplane !!! to be generalize
    u_l = (c_float*n)()

    z_l=[]
    for i in range(0,n): z_l.append(i)

    if(len(args)>2): 

        M.MUPHY.get_fluid_density_centerline(c_int(idfile), c_int(n), byref(u_l))

    else:
        M.MUPHY.get_fluid_density_centerline(c_int(idfile), c_int(n), byref(u_l))

    return z_l,u_l

##########################################
def get_fluid_velocity_centerline(*args):
    """
    get velocity of fluid on centerline
    """

    if(len(args)<3): return 

    n_l = (c_int*1)()
    idfile = int(args[2])

    M.MUPHY.get_fluid_centerline_size(c_int(idfile),byref(n_l))

    n = n_l[0]

    # this is a cut in the z-plane at midplane !!! to be generalize
    u_l = (c_float*n)()

    z_l=[]
    for i in range(0,n): z_l.append(i)

    if(len(args)>2): 

        M.MUPHY.get_fluid_velocity_centerline(c_int(idfile), c_int(n), byref(u_l))

    else:
        M.MUPHY.get_fluid_velocity_centerline(c_int(idfile), c_int(n), byref(u_l))

    return z_l,u_l

##########################################
def get_fluid_sample_size(args):
    """
    get size of sample for fluid
    """

    id = int(args[2])
    n = (c_int*1)()
    np = (c_int*1)()

    M.MUPHY.get_fluid_sample_size(c_int(id), byref(n), byref(np))

    return id,n[0],np[0]
##########################################
def get_fluid_sample_points(*args):
    """
    get sample for fluid
    """

    argp = args
    id,n,np = get_fluid_sample_size(argp)

    px_l = (c_float*n)()
    py_l = (c_float*n)()
    pz_l = (c_float*n)()

    M.MUPHY.get_fluid_sample_points(c_int(id), c_int(n), byref(px_l), byref(py_l), byref(pz_l))

    i1_l = (c_int*np)()
    i2_l = (c_int*np)()
    i3_l = (c_int*np)()

    M.MUPHY.get_fluid_sample_lineprimitives(c_int(id), c_int(np), byref(i1_l), byref(i2_l), byref(i3_l))

    # print 'NNN:',n
    # for i in range(len(px_l)):
    #     if i%100==0: print 'Sample Point:',i,px_l[i],py_l[i],pz_l[i]

    return n,px_l,py_l,pz_l,np,i1_l,i2_l,i3_l

##########################################
def get_fluid_sample_points_density(*args):
    """
    get density on sample for fluid
    """

    argp = args
    id,n,np = get_fluid_sample_size(argp)

    v_l = (c_float*n)()

    M.MUPHY.get_fluid_sample_points_density(c_int(id), c_int(n), byref(v_l))

    return n,v_l

##########################################
def get_fluid_sample_points_velocity(*args):
    """
    get velocity on sample for fluid
    """

    argp = args
    id,n,np = get_fluid_sample_size(argp)

    v_l = (c_float*n)()

    if(len(args)>2): 

        if args[1]=='v.x':
            M.MUPHY.get_fluid_sample_points_velocity(c_int(id), c_int(1), c_int(n), byref(v_l))

        elif args[1]=='v.y':
            M.MUPHY.get_fluid_sample_points_velocity(c_int(id), c_int(2), c_int(n), byref(v_l))

        elif args[1]=='v.z':
            M.MUPHY.get_fluid_sample_points_velocity(c_int(id), c_int(3), c_int(n), byref(v_l))

        elif args[1]=='|v|':
            M.MUPHY.get_fluid_sample_points_velocity(c_int(id), c_int(4), c_int(n), byref(v_l))

        elif args[1]=='|v|^2':
            M.MUPHY.get_fluid_sample_points_velocity(c_int(id), c_int(5), c_int(n), byref(v_l))

    return n,v_l

##########################################
def get_fluid_map_density(*args):
    """
    get density map for fluid
    """

    nx,ny,nz = get_box()

    i_l = (c_float*nx)()
    j_l = (c_float*nz)()

    for i in range(0,nx): i_l[i] = i
    for j in range(0,nz): j_l[j] = j
    # for i in range(0,nx): i_l.append(i)
    # for j in range(0,nz): j_l.append(j)

    v_l = (c_float*(nx*nz))()

    M.MUPHY.get_fluid_map_density(c_int(nx), c_int(nz), byref(v_l))

    return i_l,j_l,v_l

##########################################
def get_fluid_map_vel_x(*args):
    """
    get velocity map for fluid
    """

    nx,ny,nz = get_box()

    i_l = (c_float*nx)()
    j_l = (c_float*nz)()

    for i in range(0,nx): i_l[i] = i
    for j in range(0,nz): j_l[j] = j
    # for i in range(0,nx): i_l.append(i)
    # for j in range(0,nz): j_l.append(j)

    v_l = (c_float*(nx*nz))()

    M.MUPHY.get_fluid_map_vel_x(c_int(nx), c_int(nz), byref(v_l))

    return i_l,j_l,v_l

##########################################
def get_fluid_map_velmod(*args):
    """
    get velmod map for fluid
    """

    nx,ny,nz = get_box()

    i_l = (c_float*nx)()
    j_l = (c_float*nz)()

    for i in range(0,nx): i_l[i] = i
    for j in range(0,nz): j_l[j] = j
    # for i in range(0,nx): i_l.append(i)
    # for j in range(0,nz): j_l.append(j)

    v_l = (c_float*(nx*nz))()

    M.MUPHY.get_fluid_map_velmod(c_int(nx), c_int(nz), byref(v_l))

    return i_l,j_l,v_l

##########################################
def get_fluid_momentum_total(*args):
    """
    get total momentum for fluid
    """

    nsp = get_number_fluid_specie()
    if nsp == 0: return 0.0

    j = 0.0

    if(len(args)>2): 
        if args[1]=='v.x':
            for isp in range(nsp): j += get_fluid_momentum_total_specie(1,isp+1)

        elif args[1]=='v.y':
            for isp in range(nsp): j += get_fluid_momentum_total_specie(2,isp+1)

        elif args[1]=='v.z':
            for isp in range(nsp): j += get_fluid_momentum_total_specie(3,isp+1)

        elif args[1]=='|v|' or args[1]=='|v|^2':
            jx = jy = jz = 0
            for isp in range(nsp):
                jx += get_fluid_momentum_total_specie(1,isp+1)
                jy += get_fluid_momentum_total_specie(2,isp+1)
                jz += get_fluid_momentum_total_specie(3,isp+1)
            j = jx**2  + jy**2 + jz**2
            if args[1]=='|v|' : j = math.sqrt(j)

    else: # default is x-direction
        for isp in range(nsp): j += get_fluid_momentum_total_specie(1,isp+1)

    return j


##########################################
def get_fluid_momentum_total_specie(idir,isp):
    """
    get total momentum for specie of fluid
    """

    ju = (c_float*1)()
    jv = (c_float*1)()
    jw = (c_float*1)()

    M.MUPHY.get_fluid_momentum_total(c_int(isp), byref(ju), byref(jv), byref(jw))

    if idir==1: return ju[0]
    if idir==2: return jv[0]
    if idir==3: return jw[0]

##########################################
def get_number_fluid_specie():
    """
    get number of fluid specie
    """

    inum = (c_int*1)()
    M.MUPHY.get_number_fluid_specie(byref(inum))
    return inum[0]

##########################################
def get_fluid_number_total(*args):
    """
    get number of fluid specie...TBF
    """

    nsp = get_number_fluid_specie()
    if nsp == 0: return 0

    n = 0
    for isp in range(nsp): n += get_fluid_number_total_specie(isp+1)

    return n

##########################################
def get_fluid_number_total_specie(isp):
    """
    ...TBF
    """

    fnum = (c_float*1)()

    M.MUPHY.get_fluid_number_total(c_int(isp), byref(fnum))

    return fnum[0]

##########################################
def get_total_charge(*args):
    """
    get total charge of fluids
    """

    fnum = (c_float*1)()

    M.MUPHY.get_total_charge(byref(fnum))

    return fnum[0]

##########################################
def get_particles_pot_energy(*args):
    """
    get potential energy of particles
    """

    fnum = (c_float*1)()

    M.MUPHY.get_particles_pot_energy(byref(fnum))

    return fnum[0]

##########################################
def get_particles_kin_energy(*args):
    """
    get kinetic energy of particles
    """

    fnum = (c_float*1)()

    M.MUPHY.get_particles_kin_energy(byref(fnum))

    return fnum[0]

##########################################
def get_io_value(ioname,ioid): 
    """
    get inlet/outlet value
    """
    val = (c_float*1)()
    M.MUPHY.get_io_value(c_int(len(ioname)), ioname, c_int(ioid), byref(val))
    return val[0]

##########################################
def get_io_flow(ioname,ioid): 
    """
    get inlet/outlet flow value
    """
    val = (c_float*1)()
    M.MUPHY.get_io_flow(c_int(len(ioname)), ioname, c_int(ioid), byref(val))
    return val[0]

##########################################
def set_io_value(ioname,ioid,val): 
    """
    set inlet/outlet value
    """
    M.MUPHY.set_io_value(c_int(len(ioname)),ioname, c_int(ioid), c_float(val))

##########################################
def set_inflow_velocity(val): 
    """
    set inlet/outlet inflow value
    """
    M.MUPHY.set_inflow_velocity(c_float(val))

##########################################
def set_shearforcing(val): 
    """
    set shear forcing value
    """
    M.MUPHY.set_shearforcing(c_float(val))

##########################################
def set_global_temperature(val): 
    """
    set global temperature of fluid and particles
    """
    M.MUPHY.set_global_temperature(c_float(val))

##########################################
def get_global_temperature():
    """
    get global temperature of fluid and particles
    """

    temperature = (c_float*1)()

    M.MUPHY.get_global_temperature(byref(temperature))

    return temperature[0]

##########################################
def set_temperature(val): 
    """
    get temperature of fluid
    """
    M.MUPHY.set_temperature(c_float(val))
##########################################
def set_particle_temperature(val): 
    """
    get temperature of particles
    """
    M.MUPHY.set_particle_temperature(c_float(val))
##########################################
def set_wall_leakage_factor(val): 
    """
    set wall leakage factor
    """
    M.MUPHY.set_wall_leakage_factor(c_float(val))
##########################################
def refold(px,py,pz):
    """
    refold position within PBC
    """
    nx,ny,nz = get_box()

    for i in range(0,len(px)):
        if px[i]-px[i-1]>+nx: px[i] = px[i]-nx
        if px[i]-px[i-1]<-nx: px[i] = px[i]+nx
        if py[i]-py[i-1]>+ny: py[i] = py[i]-ny
        if py[i]-py[i-1]<-ny: py[i] = py[i]+ny
        if pz[i]-pz[i-1]>+nz: pz[i] = pz[i]-nz
        if pz[i]-pz[i-1]<-nz: pz[i] = pz[i]+nz
        
##########################################
def taylor(px,py,pz):
    """
    taylor a group of particles to minimum distance w respect to the 1st particle
    """

    nx,ny,nz = get_box()

    tx = []; ty = []; tz = []

    p0x = px[0]; p0y = py[0]; p0z = pz[0]

    tx.append(p0x); ty.append(p0y); tz.append(p0z)

    for i in range(1,len(px)):

        dx = px[i] - p0x
        dy = py[i] - p0y
        dz = pz[i] - p0z

        dx = dx - NINT(dx/nx)*nx
        dy = dy - NINT(dy/ny)*ny
        dz = dz - NINT(dz/nz)*nz

        tx.append(p0x + dx)
        ty.append(p0y + dy)
        tz.append(p0z + dz)

    return tx,ty,tz
        
##########################################
def NINT(a): 
    """
    nearest integer
    """
    return int(a+0.5)
##########################################
def inbox(px,py,pz):
    """
    reset positiong within central box
    """

    nx,ny,nz = get_box()

    for i in range(0,len(px)):
        px[i] = minimage(nx,px[i]) + nx/2.
        py[i] = minimage(ny,py[i]) + ny/2.
        pz[i] = minimage(nz,pz[i]) + nz/2.

##########################################
def minimage(L,d):
    """
    minimum image
    """
    return d - int(d/L+0.5)*L

##########################################
def write_current_vtk():
    """
    write vtk file
    """

    M.MUPHY.write_current_vtk()

