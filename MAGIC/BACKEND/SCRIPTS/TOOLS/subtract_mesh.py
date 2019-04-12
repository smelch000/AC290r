#!/usr/bin/env python
#
import sys,os
import numpy as np
import argparse
from TOOLS.mesh import *
from TOOLS.inletoutlet import *

###################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--subtrahend_mesh_csvfile',  required=True,            help='subtrahend mesh csv file')
    parser.add_argument('-x', '--subtractors_mesh_csvfile', required=True, nargs='+', help='subtractors mesh csvfile ')
    parser.add_argument('-t', '--nodetypes_subtractors',    required=True, nargs='+', help='node types of subtractors : (wall,inlet,outlet) ')
    parser.add_argument('-o', '--output',                   required=True,            help='output MOEBIUS files')
    parser.add_argument('-v', '--vtmoutput',                required=True,            help='output multiblock (vtm) file')
    args = parser.parse_args()

    subtrahend_mesh_csvfile = args.subtrahend_mesh_csvfile
    subtractors_mesh_csvfile = args.subtractors_mesh_csvfile
    nodetypes = args.nodetypes_subtractors
    fileout = args.output
    vtmfileout = args.vtmoutput

    print
    print 'Outputs:', fileout+'.hdr',' / ',fileout+'.dat',' / ',fileout+'.ios',
    print ' / ',vtmfileout
    print '\n'

    if not subtrahend_mesh_csvfile:
        print 'fluid_input_csvfile missing'
        sys.exit(1)

    if not subtractors_mesh_csvfile:
        print 'subtractors_mesh_csvfile missing'
        sys.exit(1)

    msh = Mesh()

    io = Inletoutlet(msh)

    # ijkroot = readmesh_csv(subtrahend_mesh_csvfile)
    ijkroot = msh.loadfile_csv(subtrahend_mesh_csvfile)

    shaperoot = ijkroot.shape
    tags = np.zeros(shaperoot[0],dtype='int')

    #msh.ijk_w = np.zeros([1,3],dtype='int')
    #msh.ijk_i = np.zeros([1,3],dtype='int')
    #msh.ijk_o = np.zeros([1,3],dtype='int')
    msh.nwall = 0; msh.ninlet = 0; msh.noutlet = 0

    inlet_id = {}
    outlet_id = {}
    io.iohead = []

    print 

    # read and subtract from the root nodes the sequence of subtractor files
    for n in range(len(subtractors_mesh_csvfile)):

        fname = subtractors_mesh_csvfile[n]

        msh.nx = 0; msh.ny = 0; msh.nz = 0

        ijk = msh.loadfile_csv(fname)

        msh.nx = max(msh.nx,np.amax(ijk[:,0])); msh.ny = max(msh.ny,np.amax(ijk[:,1])); msh.nz = max(msh.nz,np.amax(ijk[:,2]))

        msh.nxy2 = msh.nx * msh.ny; msh.nx2 = msh.nx # bad hack

        if nodetypes[n][0:4] == 'wall':

            print '... declared as wall (shlld not be used with impunity...)'

            sys.exit(1)

            msh.nwall = ijk.shape[0]
            msh.ijk_w = np.copy(ijk[:,0:3])

        elif nodetypes[n][0:5] == 'inlet':

            io_id = int(nodetypes[n][5:])
            print '... declared as inlet ... id:',io_id
            msh.ninlet = ijk.shape[0]
            msh.ijk_i = np.copy(ijk[:,0:3])

            for i,j,k,f in ijk:
                i4 = msh.i4back(i,j,k)
                inlet_id[i4] = io_id

            io.iohead.append('%d  inlet pressure 1 0 0   0. 0. 0.   0.\n' % (io_id))

        elif nodetypes[n][0:6] == 'outlet':

            io_id = int(nodetypes[n][6:])
            print '... declared as outlet ... id:',io_id
            msh.noutlet = ijk.shape[0]
            msh.ijk_o = np.copy(ijk[:,0:3])

            for i,j,k,f in ijk:
                i4 = msh.i4back(i,j,k)
                outlet_id[i4] = io_id

            io.iohead.append('%d outlet pressure 1 0 0   0. 0. 0.   0.\n' % (io_id))

        else:
            print 'nodetype not recognized:',nodetypes[n]
            sys.exit(1)

        print 


        # tag nodes to remove later
        for i,j,k,idx in ijk:

            ii,jj,kk = ijkroot[idx]

            if i==ii and j==jj and k==kk: # paranoid test
                tags[idx] = 1
            else:
                print 'internal error in subtract_mesh'
                sys.exit(1)
        
    print 

    # remove tagged nodes
    msh.nfluid = 0
    for n in xrange(shaperoot[0]):

        if tags[n] == 1: continue

        ijkroot[msh.nfluid,0:3] = ijkroot[n,0:3]
        msh.nfluid += 1

    msh.ijk_f = np.zeros([msh.nfluid,3],dtype='int')

    for n in xrange(msh.nfluid):
        msh.nx = max(msh.nx,ijkroot[n,0]); msh.ny = max(msh.ny,ijkroot[n,1]); msh.nz = max(msh.nz,ijkroot[n,2])
        msh.ijk_f[n,:] = ijkroot[n,:]
    del ijkroot

    # add arbitrary offset
    msh.nx += 1; msh.ny += 1; msh.nz += 1
    msh.nxy2 = msh.nx * msh.ny; msh.nx2 = msh.nx # hack

    msh.writeMOEBIUSinput(fileout)

    # handle inlet/outlet

    #for h in iohead:
    #    string = 'ccc %d' % (h)
    #    print string
    #    io.iohead.append(string+'\n')

    io.itppp_i_id = np.zeros(len(inlet_id),dtype='i8')
    n = 0
    for i4 in inlet_id:
        io.itppp_i_id[n] = inlet_id[i4]
        n += 1
        
    io.itppp_o_id = np.zeros(len(outlet_id),dtype='i8')
    n = 0
    for i4 in outlet_id:
        io.itppp_o_id[n] = outlet_id[i4]
        n += 1
        
    io.writeMOEBIUSinput(fileout)

    # dump VTM (multiblock) file
    msh.writeVTM(vtmfileout)

