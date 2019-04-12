#!/usr/bin/env python

import numpy as np
import sys, argparse
import logging
import vtk


def MergeVTKWorker(itime,args):

    appendF = vtk.vtkAppendFilter()

    for inpts in args.inputs:

        # g = vtkXMLUnstructuredGridReader() # VTU file
        g = vtk.vtkXMLPUnstructuredGridReader() # PVTU file
        g.SetFileName( inpts )
        g.Update()
        appendF.AddInputData( g.GetOutput() )

    # print appendF.GetMergePoints()
    appendF.MergePointsOn()
    # print appendF.GetMergePoints()
    appendF.Update()

    W = vtk.vtkXMLUnstructuredGridWriter()
    W.SetFileName( args.output )
    W.SetInputData(appendF.GetOutput())
    W.Update()

def probemapfile(fn):

    typ = None

    f = open(fn,'r')
    Ix = +100000; Iy = +100000
    Nx = -100000; Ny = -100000
    odd = False; even = False
    for line in f.readlines():
        line = line.split()

        if line[0] == '#': continue

        # i = int(float(line[0]))
        # j = int(float(line[1]))
        i = int(float(line[1]))
        j = int(float(line[0]))

        Ix = min(Ix,i); Iy = min(Iy,j)
        Nx = max(Nx,i); Ny = max(Ny,j)

        # probe
        odd = odd or i%2 == 0
        even = even or (i+1)%2 == 0

        if len(line) == 3:
            typ = 'scalar'
        elif len(line) == 6:
            typ = 'vector'

    finesse = (odd and even)

    f.close()

    return Nx,Ny,Ix,Iy,typ,finesse

def loadmapfile(Ixt,Iyt,Nx,Ny,finesse,fn):

    if finesse:
        vinit = np.nan
    else:
        vinit = 0.

    u = np.zeros( (Nx,Ny,4) )
    u[:,:] = vinit

    f = open(fn,'r')
    for line in f.readlines():
        line = line.split()

        if not line: continue # Skip empty lines

        if line[0][0] == '#': continue

        # ii = int(float(line[0]))
        # jj = int(float(line[1]))
        ii = int(float(line[1]))
        jj = int(float(line[0]))

        iii = ii-Ixt
        jjj = jj-Iyt

        if len(line) == 3:
            u[iii,jjj] = float(line[2])
        elif len(line) == 6:
            u[iii,jjj,:] = [float(line[2]),float(line[3]),float(line[4]),float(line[5])]
            # if not finesse: print Ixt,ii,iii,line[2]
    f.close()

    return u

def interpolate(Nx,Ny,Ix,Iy,u):

    o2 = 1./2.
    o4 = 1./4.

    for j in xrange(0,Ny-Iy,2):
        for i in xrange(0,Nx-Ix,2):
    # for j in xrange(0,Ny-Iy-1,2):
    #     for i in xrange(0,Nx-Ix-1,2):

            for icart in [0,1,2,3]:
                v = u[i,j][icart]
                u[i+1, j  ][icart] += o2 * v
                u[i-1, j  ][icart] += o2 * v
                u[i,   j+1][icart] += o2 * v
                u[i,   j-1][icart] += o2 * v
                u[i+1, j+1][icart] += o4 * v
                u[i-1, j+1][icart] += o4 * v
                u[i+1, j-1][icart] += o4 * v
                u[i-1, j-1][icart] += o4 * v

    return u

def MergeMAPWorker(itime,args):

    if itime % args.frequency != 0:
        return

    Ix = +10000; Iy = +10000
    Nx = -10000; Ny = -10000
    for inpts in args.inputs:

        Nx0,Ny0,Ix0,Iy0,typ0,finesse0 = probemapfile(inpts)
        Ix = min(Ix0,Ix); Iy = min(Iy0,Iy)
        Nx = max(Nx0,Nx); Ny = max(Ny0,Ny)

    Ux = {}; Uy = {}; Uz = {}; UU = {}
    for inpts in args.inputs:

        f = open(inpts,'r')

        ux = {}; uy = {}; uz = {}; uu = {}
        for line in f.readlines():
            line = line.split()
            if not line: continue # Skip empty lines
            if line[0][0] == '#': continue
            ii = float(line[1])
            jj = float(line[0])
            if typ0 == 'scalar':
                uu[iii,jjj] = float(line[2])
            elif typ0 == 'vector':
                ux[ii,jj] = float(line[2])
                uy[ii,jj] = float(line[3])
                uz[ii,jj] = float(line[4])
                uu[ii,jj] = float(line[5])
        # for rr in uu.keys():
        #     ii,jj = rr
        #     Ux[ii,jj] = Ux[ii,jj]

        f.close()

    f = open(args.output,'w')
    for rr in uu.keys():
        if typ0 == 'scalar':
            print >> f, rr[0], rr[1], ux[rr]
        elif typ0 == 'vector':
            print >> f, rr[0], rr[1], ux[rr], uy[rr], uz[rr], uu[rr]
    f.close()
        

class MergeMap():

    def __init__(self,args,interp=True):

        Nx0,Ny0,Ix0,Iy0,typ0,finesse0 = probemapfile(args.inputs[0])
        Nx1,Ny1,Ix1,Iy1,typ1,finesse1 = probemapfile(args.inputs[1])

        print Ix0,Nx0,'...',Iy0,Ny0
        print Ix1,Nx1,'...',Iy1,Ny1

        Ix = min(Ix0,Ix1); Iy = min(Iy0,Iy1)
        Nx = max(Nx0,Nx1); Ny = max(Ny0,Ny1)

        if typ0 != typ1:
            print 'intrnl errrr',typ0,typ1
            sys.exit(1)
        typ = typ0

        u0 = loadmapfile(Ix0,Iy0,Nx,Ny,finesse0,args.inputs[0])
        u1 = loadmapfile(Ix1,Iy1,Nx,Ny,finesse1,args.inputs[1])
     
        if finesse1:
            Ixt = Ix0; Iyt = Iy0; Nxt = Nx0; Nyt = Ny0
            uc = u0
            uf = u1
        else:
            Ixt = Ix1; Iyt = Iy1; Nxt = Nx1; Nyt = Ny1
            uc = u1
            uf = u0
     
        if interp:
            uc = interpolate(Nxt,Nyt,Ixt,Iyt,uc)
     
        f = open(args.output,'w')
     
        for i in xrange(0,Nx-Ix):
            for j in xrange(0,Ny-Iy):
     
                if np.isnan(uf[i,j,0]):
     
                    if typ == 'scalar':
                        print >> f, '%d %d %f' % (j+Iy, i+Ix, uc[i,j,0])
                    elif typ == 'vector':
                        print >> f, '%d %d %f %f %f %f' % (j+Iy, i+Ix, uc[i,j,0], uc[i,j,1], uc[i,j,2], uc[i,j,3])
     
                else:
     
                    if typ == 'scalar':
                        print >> f, '%d %d %f' % (j+Iy, i+Ix, uf[i,j,0])
                    elif typ == 'vector':
                        print >> f, '%d %d %f %f %f %f' % (j+Iy, i+Ix, uf[i,j,0], uf[i,j,1], uf[i,j,2], uf[i,j,3])
     
class MergeMapTri():

    def __init__(self,args):

        import matplotlib.pyplot as plt
        import matplotlib.tri as tri
        from matplotlib.mlab import griddata

        """
        Nx0,Ny0,Ix0,Iy0,typ0,finesse0 = probemapfile(args.inputs[0])
        Nx1,Ny1,Ix1,Iy1,typ1,finesse1 = probemapfile(args.inputs[1])

        print Ix0,Nx0,'...',Iy0,Ny0
        print Ix1,Nx1,'...',Iy1,Ny1

        Ix = min(Ix0,Ix1); Iy = min(Iy0,Iy1)
        Nx = max(Nx0,Nx1); Ny = max(Ny0,Ny1)

        if typ0 != typ1:
            print 'intrnl errrr',typ0,typ1
            sys.exit(1)
        typ = typ0

        print Ix0,Nx0,'...',Iy0,Ny0
        print Ix1,Nx1,'...',Iy1,Ny1

        """

        Ix = +10000; Iy = +10000
        Nx = -10000; Ny = -10000
        for inpts in args.inputs:

            Nx0,Ny0,Ix0,Iy0,typ0,finesse0 = probemapfile(inpts)
            Ix = min(Ix0,Ix); Iy = min(Iy0,Iy)
            Nx = max(Nx0,Nx); Ny = max(Ny0,Ny)

        print 'Ix,Nx,Iy,Ny:', Ix,Nx,Iy,Ny

        ux = {}
        uy = {}
        uz = {}
        uu = {}
        # u[:,:] = vinit

        for inpts in args.inputs:

            f = open(inpts,'r')
            for line in f.readlines():
                line = line.split()
                if not line: continue # Skip empty lines
                if line[0][0] == '#': continue
                ii = int(float(line[1]))
                jj = int(float(line[0]))
                iii = ii-Ix0
                jjj = jj-Iy0
                if typ0 == 'scalar':
                    uu[iii,jjj] = float(line[2])
                elif typ0 == 'vector':
                    ux[iii,jjj] = float(line[2])
                    uy[iii,jjj] = float(line[3])
                    uz[iii,jjj] = float(line[4])
                    uu[iii,jjj] = float(line[5])
            f.close()

        xx = []; yy = []
        for rr in uu.keys():
            xx.append(rr[0])
            yy.append(rr[1])
        x = np.asarray(xx).flatten()
        y = np.asarray(yy).flatten()

        if typ0 == 'vector':
            Ux = np.asarray(ux.values()).flatten()
            Uy = np.asarray(uy.values()).flatten()
            Uz = np.asarray(uz.values()).flatten()

        U = np.asarray(uu.values()).flatten()

        # Create the Delaunay triangulation
        # triang = tri.Triangulation(x, y)
        triang = tri.Triangulation(y, x)

        # Mask off unwanted triangles.
        # xmid = x[triang.triangles].mean(axis=1)
        # ymid = y[triang.triangles].mean(axis=1)
        # mask = np.where(xmid*xmid + 2.*ymid*ymid < 100., 1, 0)
        # triang.set_mask(mask)

        # U = U / (U.max() - U.min())
        low = U.min(); upp = U.max()
        print 'low,upp:',low,upp,typ0

        if typ0 == 'vector':
            upp = 224000

        levels = np.arange(low, upp, 1.e-2 * (upp - low))

        xi,yi = np.mgrid[0:Nx:complex(0,Nx), 0:Ny:complex(0,Ny)]

        if typ0 == 'vector':
            uxi = griddata(x,y,Ux,xi,yi,interp='linear')
            uyi = griddata(x,y,Uy,xi,yi,interp='linear')
            uzi = griddata(x,y,Uz,xi,yi,interp='linear')
        uui = griddata(x,y,U,xi,yi,interp='linear')

        # fig0, ax0 = plt.subplots()
        fig1, (ax1, ax2) = plt.subplots(ncols=2)

        if typ0 == 'vector':
            strm = ax1.streamplot(yi, xi, uxi, uzi, density=5.0, color=uui, linewidth=1)

        # plt.gca().set_aspect('equal')

        cset = plt.tricontourf(triang, U, levels=levels)
        plt.colorbar(cset)

        # plt.tricontour(triang, U, colors='k')
        # plt.title('Prottt')

        plt.show()

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
     
    parser.add_argument('-i', '--inputs',  nargs = '+', help='input map files')
    parser.add_argument('-o', '--output',               help='output map file')
    args = parser.parse_args()

    # MergeMap(args)
    MergeMapTri(args)

