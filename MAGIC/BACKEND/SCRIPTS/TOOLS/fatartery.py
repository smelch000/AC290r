#!/usr/bin/env python

import sys
import numpy
from mesh import *
from d3q19 import *
from inletoutlet import *

# number of iterations for fattening (~ 10 mesh points)
NITER = 10 

# input files
INHDR = 'bgkflag.orig.hdr'
INDAT = 'bgkflag.orig.dat'
INIO = 'inletoutlet.orig.inp'

#output files
OUTPUT = 'bgkflag'
OUTIO = 'inletoutlet.inp'

##################
print 'reading inlet/outlet file',INIO

msh = Mesh()

msh.loadMOEBIUSinput(INHDR,INDAT)

lat = D3q19(msh)

io = Inletoutlet(msh)
io.loadMOEBIUSinput(INIO)

msh.writexyz('CONF.init.xyz')

########################
ninlet_old = msh.ninlet
itppp_i_old = numpy.copy(msh.itppp_i)
itppp_i_id_old = numpy.copy(io.itppp_i_id)

noutlet_old = msh.noutlet
itppp_o_old = numpy.copy(msh.itppp_o)
itppp_o_id_old = numpy.copy(io.itppp_o_id)

########################
print '\nGrowing Inlet / Outlet'
for iter in range(NITER):

    # msh.writexyz('CONF.'+str(iter)+'.xyz')

    print 
    print 'iter:',iter, 'ninlet+noutlet:',msh.ninlet+msh.noutlet

    for id in xrange(io.nin_id + io.nout_id):

        if io.type[id]=='inlet':
            print 'ninlet:',msh.ninlet,
        else:
            print 'noutlet:',msh.noutlet,

        newcomers = []
        for i4,ifl in io.iodict[id].iteritems():

            f = msh.ijk(i4)

            newneigh = lat.getneighbors_aniso(f,io.dir[id])
            for i4 in newneigh:
                newcomers.append(i4)

        # print ' id:',id,io.type[id],' remove duplicates, len:',len(newcomers)
        # remove duplicates
        newcomers = list(set(newcomers))

        print ' add new neighbors len:',len(newcomers),

        if io.type[id]=='inlet':
            jtppp = numpy.zeros(msh.ninlet + len(newcomers), dtype='i8')
            for n in xrange(msh.ninlet): jtppp[n] = msh.itppp_i[n]
            for i4 in newcomers:
                jtppp[msh.ninlet] = i4
                msh.ninlet += 1

        else:
            jtppp = numpy.zeros(msh.noutlet + len(newcomers), dtype='i8')
            for n in xrange(msh.noutlet): jtppp[n] = msh.itppp_o[n]
            for i4 in newcomers:
                jtppp[msh.noutlet] = i4
                msh.noutlet += 1

        jtppp.sort()

        if io.type[id]=='inlet':
            inlet=[]
            for n in xrange(msh.ninlet):

                msh.itppp_i[n] = jtppp[n]

                i,j,k = msh.ijk(msh.itppp_i[n])
                inlet.append([i,j,k])
                io.iodict[id][msh.itppp_i[n]] = 0

            print ' total len:',len(inlet)
        else:
            outlet=[]
            for n in xrange(msh.noutlet):

                msh.itppp_o[n] = jtppp[n]

                i,j,k = msh.ijk(msh.itppp_o[n])
                outlet.append([i,j,k])
                io.iodict[id][msh.itppp_o[n]] = 0

            print ' total len:',len(outlet)


########################
print '\nGrowing fluid'
for iter in range(NITER):

    # msh.writexyz('CONF.'+str(iter)+'.xyz')

    print 
    print 'iter:',iter, 'nfluid:',msh.nfluid

    newcomers = []
    for f in msh.fluid:
        newneigh = lat.getneighbors(f)
        for i4 in newneigh:
            newcomers.append(i4)
    print ' remove duplicates, len:',len(newcomers),
    # remove duplicates
    newcomers = list(set(newcomers))

    if iter<NITER-1:

        jtppp = numpy.zeros(msh.nfluid + len(newcomers), dtype='i8')
        for n in xrange(msh.nfluid):
            jtppp[n] = msh.itppp_f[n]

        print '  add new neighbors len:',len(newcomers),
        for i4 in newcomers:
            jtppp[msh.nfluid] = i4
            msh.nfluid += 1
     
        print '   sort ...new nfluid:', msh.nfluid
        jtppp.sort()
     
        msh.fluid=[]
        for n in xrange(msh.nfluid):
     
            msh.itppp_f[n] = jtppp[n]
     
            i,j,k = msh.ijk(msh.itppp_f[n])
            msh.fluid.append([i,j,k])

    else:

        jtppp = []
        for i4 in newcomers:
            jtppp.append(i4)
            msh.nwall += 1
        print '  ...new wall:', msh.nwall
     

######
print '\nRecovering original inlet / outlet nodes'
ntmp = 0
for n in xrange(msh.ninlet):
    i4 = msh.itppp_i[n]
    ir = msh.binsearch(ninlet_old, itppp_i_old, i4)

    if itppp_i_old[ir]==i4 :
        msh.itppp_i[ntmp] = i4
        io.itppp_i_id[ntmp] = itppp_i_id_old[ir]
        ntmp += 1 
    else:
        # msh.itppp_w[msh.nwall] = i4
        jtppp.append(i4)
        msh.nwall += 1 
msh.ninlet = ntmp

ntmp = 0
for n in xrange(msh.noutlet):
    i4 = msh.itppp_o[n]
    ir = msh.binsearch(noutlet_old, itppp_o_old, i4)

    if itppp_o_old[ir]==i4 :
        msh.itppp_o[ntmp] = i4
        io.itppp_o_id[ntmp] = itppp_o_id_old[ir]
        ntmp += 1 
    else:
        jtppp.append(i4)
        msh.nwall += 1 
msh.noutlet = ntmp

jtppp = list(set(jtppp)) # list of uniques
jtppp.sort()

i = 0
for i4 in jtppp:
    msh.itppp_w[i] = i4
    i+=1

#############################
msh.writexyz('CONF.final.xyz')

msh.writeMOEBIUSinput(OUTPUT)

io.writeMOEBIUSinput(OUTIO)

