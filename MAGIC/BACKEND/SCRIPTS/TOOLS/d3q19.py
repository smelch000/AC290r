#!/usr/bin/env python

import sys
import numpy

class D3q19:
    
    def __init__(self,msh):

        self.icx = numpy.zeros(19, dtype='i')
        self.icy = numpy.zeros(19, dtype='i')
        self.icz = numpy.zeros(19, dtype='i')

        self.icx[0] = 0; self.icy[0] = 0; self.icz[0] = 0
        self.icx[1] = 0; self.icy[1] = 0; self.icz[1] = 1
        self.icx[2] = 1; self.icy[2] = 0; self.icz[2] = 1
        self.icx[3] =-1; self.icy[3] = 0; self.icz[3] = 1
        self.icx[4] = 0; self.icy[4] = 1; self.icz[4] = 1
        self.icx[5] = 0; self.icy[5] =-1; self.icz[5] = 1
        self.icx[6] = 0; self.icy[6] = 0; self.icz[6] =-1
        self.icx[7] =-1; self.icy[7] = 0; self.icz[7] =-1
        self.icx[8] = 1; self.icy[8] = 0; self.icz[8] =-1
        self.icx[9] = 0; self.icy[9] =-1; self.icz[9] =-1
        self.icx[10]= 0; self.icy[10]= 1; self.icz[10]=-1
        self.icx[11]= 1; self.icy[11]= 0; self.icz[11]= 0
        self.icx[12]=-1; self.icy[12]= 0; self.icz[12]= 0
        self.icx[13]= 0; self.icy[13]= 1; self.icz[13]= 0
        self.icx[14]= 0; self.icy[14]=-1; self.icz[14]= 0
        self.icx[15]= 1; self.icy[15]= 1; self.icz[15]= 0
        self.icx[16]=-1; self.icy[16]= 1; self.icz[16]= 0
        self.icx[17]=-1; self.icy[17]=-1; self.icz[17]= 0
        self.icx[18]= 1; self.icy[18]=-1; self.icz[18]= 0

        self.gridspacing = 1

        self.msh = msh

    # returns a list of all d3q19 neighbors of the i,j,k point. 
    # Neighbors can be fluid, wall, inlet or outlet nodes
    def getneighbors(self,pnt,gridspacing=False):

        if not gridspacing: self.gridspacing=self.msh.gridspacing
        i,j,k = pnt[0],pnt[1],pnt[2]
        newneigh=[]
        for n in xrange(1,19):

            ii = i + self.gridspacing * self.icx[n]
            jj = j + self.gridspacing * self.icy[n]
            kk = k + self.gridspacing * self.icz[n]

            if ii<1 or ii>self.msh.nx: continue
            if jj<1 or jj>self.msh.ny: continue
            if kk<1 or kk>self.msh.nz: continue

            i4 = self.msh.i4back(ii,jj,kk)

            ir = self.msh.binsearch(self.msh.nfluid, self.msh.itppp_f, i4)
            foundfluid = self.msh.nfluid>0 and (self.msh.itppp_f[ir]==i4)

            if not foundfluid: 
                ir = self.msh.binsearch(self.msh.nwall, self.msh.itppp_w, i4)
                foundwall = self.msh.nwall>0 and (self.msh.itppp_w[ir]==i4)

                if not foundwall: 
                    ir = self.msh.binsearch(self.msh.ninlet, self.msh.itppp_i, i4)
                    foundinlet = self.msh.ninlet>0 and (self.msh.itppp_i[ir]==i4)
     
                    if not foundinlet: 
                        ir = self.msh.binsearch(self.msh.noutlet, self.msh.itppp_o, i4)
                        foundoutlet = self.msh.noutlet>0 and (self.msh.itppp_o[ir]==i4)
     
                        if not foundoutlet:
                            newneigh.append(i4)

        return newneigh
        
    # returns a list of all d3q19 neighbors of the i,j,k point of type fluid
    def getfluidneighbors(self,pnt,gridspacing=False):

        if not gridspacing: self.gridspacing=self.msh.gridspacing
        i,j,k = pnt[0],pnt[1],pnt[2]
        newneigh=[]
        for n in xrange(1,19):

            ii = i + self.gridspacing * self.icx[n]
            jj = j + self.gridspacing * self.icy[n]
            kk = k + self.gridspacing * self.icz[n]

            if ii<1 or ii>self.msh.nx: continue
            if jj<1 or jj>self.msh.ny: continue
            if kk<1 or kk>self.msh.nz: continue

            i4 = self.msh.i4back(ii,jj,kk)

            ir = self.msh.binsearch(self.msh.nfluid, self.msh.itppp_f, i4)
            foundfluid = self.msh.nfluid>0 and (self.msh.itppp_f[ir]==i4)

        return newneigh
        
    # returns a list of all d3q19 neighbors of the i,j,k point. 
    # Neighbors can be fluid, wall, inlet or outlet nodes.
    # important points:
    # * the function is accelerated via numpy.searchsorted
    # * it is **incremental** for the newneigh array (appending data)
    def getneighbors_fast(self,pnt,newneigh=[],gridspacing=False):

        if not gridspacing: self.gridspacing=self.msh.gridspacing
        i,j,k = pnt[0],pnt[1],pnt[2]
        newneigh=[]

        if i<1 or i>self.msh.nx or j<1 or j>self.msh.ny or k<1 or k>self.msh.nz :
            return newneigh

        nfl = len(self.msh.itppp_f)
        nwl = len(self.msh.itppp_w)
        nin = len(self.msh.itppp_i)
        nou = len(self.msh.itppp_o)

        for n in range(1,19):

            ii = i + self.gridspacing * self.icx[n]
            jj = j + self.gridspacing * self.icy[n]
            kk = k + self.gridspacing * self.icz[n]

            if ii<1 or ii>self.msh.nx: continue
            if jj<1 or jj>self.msh.ny: continue
            if kk<1 or kk>self.msh.nz: continue

            i4 = self.msh.i4back(ii,jj,kk)

            ir = numpy.searchsorted(self.msh.itppp_f, i4)
            foundfluid = ir<nfl and (self.msh.itppp_f[ir]==i4)

            if not foundfluid: 
                ir = numpy.searchsorted(self.msh.itppp_w, i4)
                foundwall = ir<nwl and (self.msh.itppp_w[ir]==i4)

                if not foundwall: 
                    ir = numpy.searchsorted(self.msh.itppp_i, i4)
                    foundinlet = ir<nin and (self.msh.itppp_i[ir]==i4)
     
                    if not foundinlet: 
                        ir = numpy.searchsorted(self.msh.itppp_o, i4)
                        foundoutlet = ir<nou and (self.msh.itppp_o[ir]==i4)
     
                        if not foundoutlet:
                            newneigh.append(i4)

        return newneigh

    def getneighbors_monstre(self,pnt,newneigh=[],gridspacing=False):

        from scipy import ndimage

        DEAD_NODE=0
        FLUID_NODE=1
        INLET_NODE=3
        OUTLET_NODE=4

        nodes = np.full([nx/GRIDSPACING+1,ny/GRIDSPACING+1,nz/GRIDSPACING+1],DEAD_NODE, np.uint)
        dead_nodes=(nodes == DEAD_NODE)

        kernel_mat = ndimage.generate_binary_structure(3, 2)
        assert kernel_mat.sum() == 19
        fluid_nodes_dilated = ndimage.morphology.binary_dilation(nodes==FLUID_NODE,kernel_mat)

        wall_nodes= fluid_nodes_dilated & dead_nodes

        return wall_nodes


        
    # find neighbors only orthogonal to a given direction
    def getneighbors_aniso(self,pnt,dir):

        i,j,k = pnt[0],pnt[1],pnt[2]
        newneigh=[]
        for n in range(1,19):

            if not self.icx[n]*dir[0] + self.icy[n]*dir[1] + self.icz[n]*dir[2] == 0: continue

            ii,jj,kk = i+self.gridspacing*self.icx[n],j+self.gridspacing*self.icy[n],k+self.gridspacing*self.icz[n]
            i4 = self.msh.i4back(ii,jj,kk)

            ir = self.msh.binsearch(self.msh.nfluid, self.msh.itppp_f, i4)
            foundfluid = (self.msh.itppp_f[ir]==i4)

            if not foundfluid: 
                ir = self.msh.binsearch(self.msh.nwall, self.msh.itppp_w, i4)
                foundwall = (self.msh.itppp_w[ir]==i4)

                if not foundwall: 
                    ir = self.msh.binsearch(self.msh.ninlet, self.msh.itppp_i, i4)
                    foundinlet = (self.msh.itppp_i[ir]==i4)
     
                    if not foundinlet: 
                        ir = self.msh.binsearch(self.msh.noutlet, self.msh.itppp_o, i4)
                        foundoutlet = (self.msh.itppp_o[ir]==i4)
     
                        if not foundoutlet:
                            newneigh.append(i4)

        return newneigh
        
    def get_fluid_insature(self,fluidpts,gridspacing=False):

        if not gridspacing: self.gridspacing=self.msh.gridspacing

        insature = []
        nfl = len(fluidpts)
        for i,j,k in fluidpts:

            l_insature = False
            for n in range(1,19):

                ii = i + self.gridspacing * self.icx[n]
                jj = j + self.gridspacing * self.icy[n]
                kk = k + self.gridspacing * self.icz[n]

                i4 = self.msh.i4back(ii,jj,kk)

                ir = numpy.searchsorted(self.msh.itppp_f, i4)
                foundfluid = ir<nfl and (self.msh.itppp_f[ir]==i4)

                if not foundfluid: 
                    l_insature = True
                    break

            if l_insature:
                insature.append([i,j,k])

        return insature
        
