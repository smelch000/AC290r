#!/usr/bin/env python

import sys,os
import csv
import numpy as np
import ctypes
from myvtk import *
from math import *
from inletoutlet import *
import inletoutlet as IO

from ctypes import *
from muphywrapper import M
        
class Mesh:
    
    def __init__(self,nx=0,ny=0,nz=0,gridspacing=1):
        
        self.ID_FLUIDNODE = 1
        self.ID_WALLNODE = 2
        self.ID_INLETNODE = 3
        self.ID_OUTLETNODE = 4

        self.nfluid = 0;  self.fluid=[]
        self.nwall = 0;   self.wall=[]
        self.ninlet = 0;  self.inlet=[]
        self.noutlet = 0; self.outlet=[]

        self.nv = 0; self.nxy2 = 0; self.nx2 = 0

        if nx>0 and ny>0 and nz>0:
            self.nx = nx; self.ny = ny; self.nz = nz
            self.box_strides()

        self.itppp_f = np.zeros(0,dtype='i8')
        self.itppp_w = np.zeros(0,dtype='i8')
        self.itppp_i = np.zeros(0,dtype='i8')
        self.itppp_o = np.zeros(0,dtype='i8')

        self.i_id = np.zeros(0,dtype='i8')
        self.o_id = np.zeros(0,dtype='i8')

        self.ijk_f = np.zeros(0,dtype='i')
        self.ijk_w = np.zeros(0,dtype='i')
        self.ijk_i = np.zeros(0,dtype='i')
        self.ijk_o = np.zeros(0,dtype='i')
        self.gridspacing=gridspacing
        # M.fill_muphy_object(version='2', xpu='cpu')
        self.i_globs = None
        self.o_globs = None
        # M.MOEBIUS.muphy2wrapper_init(mycomm=0)

    def assembleUnstructuredGridPointLike(self):

        def makeugrid(itppp_):

            ugrid = vtkUnstructuredGrid()

            points = vtkPoints()
            vert = vtkVertex() # only points/vertexes - no cells
            cells = vtkCellArray()
            for n, i4 in np.ndenumerate(itppp_):
                i, j, k = self.ijk(i4)
                points.InsertNextPoint(i, j, k)
                vert.GetPointIds().SetId(0, n[0])
                cells.InsertNextCell(vert)
            ugrid.SetPoints(points)
            ugrid.SetCells(VTK_VERTEX, cells)
            return ugrid

        ug_f = makeugrid(self.itppp_f)
        ug_w = makeugrid(self.itppp_w)
        ug_i = makeugrid(self.itppp_i)
        ug_o = makeugrid(self.itppp_o)

        #if self.itppp_f.size>0: self.addarray_unstructuredgrid(ug_f, self.ID_FLUIDNODE)
        #if self.itppp_w.size>0: self.addarray_unstructuredgrid(ug_w, self.ID_WALLNODE)
        #if self.itppp_i.size>0: self.addarray_unstructuredgrid(ug_i, self.ID_INLETNODE)
        #if self.itppp_o.size>0: self.addarray_unstructuredgrid(ug_o, self.ID_OUTLETNODE)

        return ug_f, ug_w, ug_i, ug_o

    #
    # function nodes as lists and stores internally as numpy arrays
    def specifyNodes(self,itp_f,itp_w,itp_i,itp_o,i_id=None,o_id=None,i_globs=None,o_globs=None):

        self.itppp_f = np.asarray( itp_f )
        self.itppp_w = np.asarray( itp_w )
        self.itppp_i = np.asarray( itp_i )
        self.itppp_o = np.asarray( itp_o )

        self.itppp_i_id = None
        self.itppp_o_id = None

        if i_id != None: self.itppp_i_id = np.asarray( i_id )
        if o_id != None: self.itppp_o_id = np.asarray( o_id )

        self.nfluid = self.itppp_f.size
        self.nwall = self.itppp_w.size
        self.ninlet = self.itppp_i.size
        self.noutlet = self.itppp_o.size

        if i_globs: self.i_globs = i_globs
        if o_globs: self.o_globs = o_globs

        return

    def box_strides(self):
        """
        computes stride as (nx+2) and (nx+2)*(ny+2) 
        """

        self.nxy2 = long(long(self.nx+2)*long(self.ny+2))
        self.nx2 = long(self.nx+2)

        # array of points
        self.nv = long(self.nx+1)*long(self.ny+1)*long(self.nz+1)

    def loadMOEBIUSinput(self,filehdr,filedat):
        """
        function to read a .hdr/.dat file
        """

        print 'reading mesh header file:',filehdr

        h = open(filehdr,'r')
        line = h.readline().split()

        # self.nx,self.ny,self.nz = int(line[0]),int(line[1]),int(line[2])
        self.nx,self.ny,self.nz = long(line[0]),long(line[1]),long(line[2])
        line = h.readline()
        try:
          line = h.readline()
        except:
          line="1"
          
        self.gridspacing = int(line)
        h.close()

        ##################
        print 'reading mesh dat file:',filedat

        self.box_strides()

        d = open(filedat,'r')

        self.nfluid = 0
        self.nwall = 0
        self.ninlet = 0
        self.noutlet = 0
        for line in d.readlines():
            line = line.split()
            i,j,k,flg = int(line[0]),int(line[1]),int(line[2]),int(line[3])
            i4 = long(long(k)*self.nxy2 + long(j)*self.nx2 + long(i))
            if flg == 1:
                self.nfluid += 1

            elif flg == 2:
                self.nwall += 1

            elif flg == 3:
                self.ninlet += 1

            elif flg == 4:
                self.noutlet += 1

        d.close()

        print '\nnfluid: %d, wall: %d, inlet: %d, noutlet: %d' % (self.nfluid,self.nwall,self.ninlet,self.noutlet)

        self.itppp_f = np.zeros(self.nfluid,  dtype='i8')
        self.itppp_w = np.zeros(self.nwall,   dtype='i8')
        self.itppp_i = np.zeros(self.ninlet,  dtype='i8')
        self.itppp_o = np.zeros(self.noutlet, dtype='i8')

        d = open(filedat,'r')

        self.nfluid = 0
        self.nwall = 0
        self.ninlet = 0
        self.noutlet = 0
        for line in d.readlines():
            line = line.split()
            i,j,k,flg = int(line[0]),int(line[1]),int(line[2]),int(line[3])

            i4 = long(long(k)*self.nxy2 + long(j)*self.nx2 + long(i))
            if flg == 1:
                self.fluid.append([i,j,k])
                self.itppp_f[self.nfluid] = i4
                self.nfluid += 1

            elif flg == 2:
                self.wall.append([i,j,k])
                self.itppp_w[self.nwall] = i4
                self.nwall += 1

            elif flg == 3:
                self.inlet.append([i,j,k])
                self.itppp_i[self.ninlet] = i4
                self.ninlet += 1

            elif flg == 4:
                self.outlet.append([i,j,k])
                self.itppp_o[self.noutlet] = i4
                self.noutlet += 1

        d.close()

    #
    # function to load a unstructured grid file
    def loadfile_unstructuredgrid_VTK(self,fname):
        """
        load unstructured grid as .vtk file and return data
        """

        grid = vtkUnstructuredGridReader()
        grid.SetFileName(fname)
        grid.Update()

        return grid.GetOutput()

    def loadfile_unstructuredgrid_VTU(self,fname):
        """
        load unstructured grid as .vtu file and return data
        """

        grid = vtkXMLUnstructuredGridReader()
        grid.SetFileName(fname)
        grid.Update()

        return grid.GetOutput()

    def writefile_unstructuredgrid_VTK(self,fname,grid):
        """
        write unstructured grid as .vtk file
        """

        writer = vtkUnstructuredGridWriter()
        writer.SetFileName(fname)
        writer.SetInputData(grid)
        writer.Write()

    def writefile_unstructuredgrid_VTU(self,fname,grid):
        """
        write unstructured grid as .vtu file
        """

        writer = vtkXMLUnstructuredGridWriter()
        writer.SetFileName(fname)
        writer.SetInputData(grid)
        writer.Write()

    def loadfile_csv(self,fname):
        """
        load a csv file containinkg i,j,k mesh points (and optionally a pre-index)
        """

        with open(fname,'r') as csvfile:
            reader = csv.reader(csvfile,delimiter=',')
            ii = 0
            for r in reader:
                if ii==0: 
                    NFIELDS = len(r)

        nsz = reader.line_num-1
        csvfile.close()

        print 'reading csv...',fname,'...nsz:',nsz

        ijk = np.zeros([nsz,NFIELDS],dtype='int64')

        with open(fname,'r') as csvfile:
            reader = csv.reader(csvfile,delimiter=',')
            ii = 0
            for r in reader:
                if ii>0: 
                    #
                    # note how content is organized : the first 3 elements are always ijk
                    # optionally the 4th is the usued for multiblock to know 
                    # the original cell index, for subtraction purposes
                    #
                    if NFIELDS == 3:
                        ijk[ii-1,:] = int(r[0]), int(r[1]), int(r[2])

                    elif NFIELDS == 4:
                        ijk[ii-1,:] = int(r[1]), int(r[2]), int(r[3]), int(r[0]) 
                ii+=1

        return ijk

    def loadtransformfile(self,filename):
        """
        read a mesh_transform.inp file (legacy)
        """

        print 'reading mesh transform file:',filename

        d = open(filename)

        line = d.readline()
        line = line.split(':')
        line = line[1].split()
        scale = float(line[0])

        line = d.readline()
        line = line.split(':')
        line = line[1].split(',')
        translate = [float(line[0]),float(line[1]),float(line[2])]
        # print 'scale:',scale,'translate:',translate
        d.close()

        return scale,translate

    def binsearch(self,nlst, lst, value):
        """
        binary search in sorted array
        """
        lo, hi = 0, nlst-1
        while lo < hi:
            mid = (lo + hi) / 2
            if lst[mid] < value:
                lo = mid + 1
            else:
                hi = mid
        return lo

    def binsearch_fast(self, nlst, lst, value):
        """
        binary search in sorted array - fast version
        """

        # numpy version
        out = np.searchsorted(lst, value)

        """
        # muphy version DOES NOT WORK ON LARGE NUMBERS....
        c_int_p = ctypes.POINTER(ctypes.c_int)
        
        lst_p = lst.ctypes.data_as(c_int_p)

        out = M.MOEBIUS.binsearch(byref(lst_p), c_int(nlst), c_int8(value))
        """

        return out

    def ijkdir(self,idir,i4):
        """
        extract i or j or k from a i4 number
        """
        
        iz = int(i4/ self.nxy2)
        if idir==3: return iz
        
        iy = int((i4 - iz*self.nxy2)/ self.nx2)
        if idir==2: return iy
        
        return int(i4 - iz*self.nxy2 - iy*self.nx2)

    def ijk(self,i4):
        """
        extract i,j,k from a i4 number
        """
        
        d = [-1,-1,-1]

        d[2] = int(i4/ self.nxy2)
        d[1] = int((i4 - d[2]*self.nxy2)/ self.nx2)
        d[0] = int(i4 - d[2]*self.nxy2 - d[1]*self.nx2)
        return d

    def ijklist(self,i4list):
        """
        extract i,j,k from a i4 number
        """
        
        dlist = []
        # for i in xrange(i4list.size):
        #    i4 = i4list[i]
        for i4 in i4list:
            dlist.append( self.ijk(i4) )

        return dlist

    def i4back(self,i,j,k):
        """
        extract i4 from i,j,k
        """
        return long(long(k)*self.nxy2 + long(j)*self.nx2 + long(i))

    def writexyz(self,fname):
        """
        write .xyz file
        """

        f = open(fname,'w')
        f.write('%d \n\n' % (self.nfluid + self.nwall + self.ninlet + self.noutlet))

        for n in xrange(self.nfluid):
            i4 = self.itppp_f[n]
            i,j,k = self.ijk(i4)
            f.write('C %d %d %d\n' % (i,j,k))

        for n in xrange(self.nwall):
            i4 = self.itppp_w[n]
            i,j,k = self.ijk(i4)
            f.write('H %d %d %d\n' % (i,j,k))

        for n in xrange(self.ninlet):
            i4 = self.itppp_i[n]
            i,j,k = self.ijk(i4)
            f.write('O %d %d %d\n' % (i,j,k))

        for n in xrange(self.noutlet):
            i4 = self.itppp_o[n]
            i,j,k = self.ijk(i4)
            f.write('N %d %d %d\n' % (i,j,k))

        f.close()

    def VTU_to_MOEBIUSmesh(self,UG):
        """
        # write .hdr/.dat files starting from a unstructured grid
        # This should consider the cell centers !!!
        """

        print '\nin VTU_to_MOEBIUSmesh:',UG.GetNumberOfPoints()

        for pid in xrange(UG.GetNumberOfPoints()):

            coo = UG.GetPoint(pid)
            i,j,k = int(coo[0]),int(coo[1]),int(coo[2])

            self.nx = max(self.nx,i)
            self.ny = max(self.ny,j)
            self.nz = max(self.nz,k)

        self.nx += 4
        self.ny += 4
        self.nz += 4
            
        self.box_strides()

        self.nfluid = 0
        self.nwall = 0
        self.ninlet = 0
        self.noutlet = 0
        for pid in xrange(UG.GetNumberOfPoints()):

            #coo = UG.GetPoint(pid)
            #i,j,k = int(coo[0]),int(coo[1]),int(coo[2])

            flg = UG.GetPointData().GetArray('muphynode').GetValue(pid)

            # flg = self.ID_FLUIDNODE # only fluids for now !

            if flg == self.ID_FLUIDNODE: self.nfluid += 1
            elif flg == self.ID_WALLNODE: self.nwall += 1
            elif flg == self.ID_INLETNODE: self.ninlet += 1
            elif flg == self.ID_OUTLETNODE: self.noutlet += 1

        print >> sys.stderr, 'StatCount....FWIO:',self.nfluid,self.nwall,self.ninlet,self.noutlet

        self.itppp_f = np.zeros(self.nfluid,  dtype='i8')
        self.itppp_w = np.zeros(self.nwall,   dtype='i8')
        self.itppp_i = np.zeros(self.ninlet,  dtype='i8')
        self.itppp_o = np.zeros(self.noutlet, dtype='i8')

        self.nfluid = 0
        self.nwall = 0
        self.ninlet = 0
        self.noutlet = 0
        for pid in xrange(UG.GetNumberOfPoints()):

            coo = UG.GetPoint(pid)
            i,j,k = int(coo[0]),int(coo[1]),int(coo[2])

            flg = UG.GetPointData().GetArray('muphynode').GetValue(pid)

            # flg = self.ID_FLUIDNODE # only fluids for now !

            i4 = long(long(k)*self.nxy2 + long(j)*self.nx2 + long(i))
            if flg == self.ID_FLUIDNODE:
                self.fluid.append([i,j,k])
                self.itppp_f[self.nfluid] = i4
                self.nfluid += 1

            elif flg == self.ID_WALLNODE:
                self.wall.append([i,j,k])
                self.itppp_w[self.nwall] = i4
                self.nwall += 1

            elif flg == self.ID_INLETNODE:
                self.inlet.append([i,j,k])
                self.itppp_i[self.ninlet] = i4
                self.ninlet += 1

            elif flg == self.ID_OUTLETNODE:
                self.outlet.append([i,j,k])
                self.itppp_o[self.noutlet] = i4
                self.noutlet += 1

    def writeMOEBIUSinput(self,fileout,iohead=None):
        """
        write .hdr/.dat files
        """

        filehdr = fileout+'.hdr'
        filedat = fileout+'.dat'
        fileios = fileout+'.ios'

        print 'writing new mesh header file ',filehdr

        d = open(filehdr,'w')
        d.write('%d %d %d\n' % (self.nx,self.ny,self.nz))
        d.write('5 %d %d %d %d\n' % (self.nfluid,self.nwall,self.ninlet,self.noutlet))
        d.write('%d\n' % (self.gridspacing))
        d.close()

        if len(self.itppp_f) + len(self.itppp_w) + len(self.itppp_i) + len(self.itppp_o) == 0:

            if self.ijk_f.shape[0] + self.ijk_w.shape[0] + self.ijk_i.shape[0] + self.ijk_o.shape[0] == 0:
                print 'Error: both itppp and ijk arrays are empty in writeVTM !'

            if self.nx + self.ny + self.nz == 0:

                self.nx = max(np.amax(self.ijk_f[:,0]), np.amax(self.ijk_w[:,0]), np.amax(self.ijk_i[:,0]), np.amax(self.ijk_o[:,0]))
                self.ny = max(np.amax(self.ijk_f[:,1]), np.amax(self.ijk_w[:,1]), np.amax(self.ijk_i[:,1]), np.amax(self.ijk_o[:,1]))
                self.nz = max(np.amax(self.ijk_f[:,2]), np.amax(self.ijk_w[:,2]), np.amax(self.ijk_i[:,2]), np.amax(self.ijk_o[:,2]))

                # bad hack
                self.nxy2 = long(long(self.nx)*long(self.ny))
                self.nx2 = long(self.nx)

            self.itppp_f = np.zeros(self.ijk_f.shape[0], dtype='i8')
            n = 0
            for i,j,k in self.ijk_f:
                self.itppp_f[n] = self.i4back(i,j,k)
                n += 1
            
            self.itppp_w = np.zeros(self.ijk_w.shape[0], dtype='i8')
            n = 0
            for i,j,k in self.ijk_w:
                self.itppp_w[n] = self.i4back(i,j,k)
                n += 1
            
            self.itppp_i = np.zeros(self.ijk_i.shape[0], dtype='i8')
            n = 0
            for i,j,k in self.ijk_i:
                self.itppp_i[n] = self.i4back(i,j,k)
                n += 1
            
            self.itppp_o = np.zeros(self.ijk_o.shape[0], dtype='i8')
            n = 0
            for i,j,k in self.ijk_o:
                self.itppp_o[n] = self.i4back(i,j,k)
                n += 1
            
        jtppp = []
        for i in xrange(self.nfluid):  jtppp.append( [self.itppp_f[i], self.ID_FLUIDNODE, -99] )
        for i in xrange(self.nwall):   jtppp.append( [self.itppp_w[i], self.ID_WALLNODE, -99] )
        for i in xrange(self.ninlet):  jtppp.append( [self.itppp_i[i], self.ID_INLETNODE, self.itppp_i_id[i]] )
        for i in xrange(self.noutlet): jtppp.append( [self.itppp_o[i], self.ID_OUTLETNODE, self.itppp_o_id[i]] )

        jtppp.sort()

        print 'writing new mesh dat file ',filedat
        d = open(filedat,'w')

        for i4,flg,io_id in jtppp:
           i,j,k = self.ijk(i4)
           d.write('%d %d %d %d\n' % (i,j,k,flg) )

        d.close()

        if self.ninlet + self.noutlet == 0:
            return True
        else:
            io = Inletoutlet(self)
            io.writeMOEBIUSinput(fileios,iohead)
            return False

    def addarray_unstructuredgrid(self,ug,val):
        """
        add value to the unstructured vtk grid
        """

        idxarr = vtk.vtkUnsignedShortArray()
        # idxarr = vtk.vtkIntArray()
        idxarr.SetNumberOfComponents(1)
        idxarr.SetName('muphynode')

        for npt in xrange(ug.GetNumberOfPoints()):
            idxarr.InsertNextValue(val)

        ug.GetPointData().AddArray(idxarr)

    def make_unstructuredgrid(self,ijk):
        """
        create an unstructured vtk grid from a ijk 
        NB the conversion should consider the conversion from center to enclosing cell
        """

        classname = self.__class__.__name__
        funcname = self.make_unstructuredgrid.__name__
        # print '....in class:',classname, '....in function:',funcname

        grid = vtkUnstructuredGrid()

        points = vtkPoints()

        NN = len(ijk[:,0])

        if NN<=0:
            print 'Error: zero size of ijk array in ',funcname
            sys.exit(1)

        # vtkCellArray is a supporting object that explicitly represents cell connectivity.
        # The cell array structure is a raw integer list of the form:
        # (n,id1,id2,...,idn, n,id1,id2,...,idn, ...) where n is the number of points in
        # the cell, and id is a zero-offset index into an associated point list.
        voxel = vtkVoxel()

        for i,j,k in ijk:
            points.InsertNextPoint(i,j,k)

        grid.SetPoints(points)

        Lib = ctypes.CDLL( os.path.join( os.getenv("MOEBIUS_ROOT"), "BACKEND/SCRIPTS/TOOLS/points2cells.so") )

        ijk1 = (ctypes.c_int*NN)()
        ijk2 = (ctypes.c_int*NN)()
        ijk3 = (ctypes.c_int*NN)()
        ncell = (ctypes.c_int*1)()
        icell = (ctypes.c_int*(8*NN))()

        for index in xrange(NN):
            ijk1[index] = ijk[index][0]
            ijk2[index] = ijk[index][1]
            ijk3[index] = ijk[index][2]

        Lib.points2cells(ctypes.c_int(NN),ctypes.c_int(self.nx), ctypes.c_int(self.ny), ctypes.c_int(self.nz), 
                         ctypes.byref(ijk1), ctypes.byref(ijk2), ctypes.byref(ijk3), 
                         ctypes.byref(ncell), ctypes.byref(icell) )

        for i in xrange(ncell[0]):

            index  = [icell[8*i    ],
                      icell[8*i + 1],
                      icell[8*i + 2],
                      icell[8*i + 3],
                      icell[8*i + 4],
                      icell[8*i + 5],
                      icell[8*i + 6],
                      icell[8*i + 7]]

            for j in range(8):
                voxel.GetPointIds().SetId(j, index[j])

            if i%500000==0: print 'Inserting cell at ',index,' / ncell:',ncell[0]

            grid.InsertNextCell(voxel.GetCellType(), voxel.GetPointIds())

        return grid

    def clean_unstructuredgrid(self,input):
        """
        # hacked from paraview core
        """

        output= vtkUnstructuredGrid()

        Locator = vtkMergePoints()

        if input.GetNumberOfCells() == 0:
            # set up a ugrid with same data arrays as input, but no points, cells or data.
            output.Allocate(1);
            output.GetPointData().CopyAllocate(input.GetPointData(), VTK_CELL_SIZE)
            output.GetCellData().CopyAllocate(input.GetCellData(), 1)
            pts = vtkPoints()
            output.SetPoints(pts)
            pts.Delete()
            return output
     
        output.GetPointData().CopyAllocate(input.GetPointData());
        output.GetCellData().PassData(input.GetCellData());
     
        # First, create a new points array that eliminate duplicate points.
        # Also create a mapping from the old point id to the new.
        newPts = vtkPoints()
        num = input.GetNumberOfPoints()
        # id = vtkIdType()
        # newId = vtkIdType()
        newId = [0]
        # newId = vtkIdList()
        # newId.SetNumberOfIds(1)
        ptMap = vtkIdList()
        ptMap.SetNumberOfIds(num)
     
        Locator.InitPointInsertion(newPts, input.GetBounds(), num)
     
        for id in xrange(num):
            # newId.SetId(id, id)
            if Locator.InsertUniquePoint( input.GetPoint(id), newId):
                output.GetPointData().CopyData(input.GetPointData(),id,newId)
            ptMap[id] = newId
     
        output.SetPoints(newPts)
     
        # Now copy the cells.
        cellPoints = vtkIdList()
        num = input.GetNumberOfCells()
        output.Allocate(num)
        for id in xrange(num):
            # special handling for polyhedron cells
            if input.GetCellType(id) == VTK_POLYHEDRON:
                input.GetFaceStream(id, cellPoints)
                ConvertFaceStreamPointIds(cellPoints, ptMap)
            else:
                input.GetCellPoints(id, cellPoints)
                for i in xrange(cellPoints.GetNumberOfIds()):
                    cellPtId = cellPoints.GetId(i)
                    newId = ptMap[cellPtId]
                    cellPoints.SetId(i, newId)
     
            output.InsertNextCell(input.GetCellType(id), cellPoints)
     
        output.Squeeze()
     
        return output

    def bool_unstructuredgrid(self,grid1,grid2,what):
        """
        boolean of two unstructured grids made of voxels.
        (can be generalized to other type of cells)
        """

        classname = self.__class__.__name__
        funcname = self.bool_unstructuredgrid.__name__
        # print '....in class:',classname, '....in function:',funcname

        ugrid = vtkUnstructuredGrid()

        voxel = vtkVoxel()
        NPVOX = voxel.GetNumberOfPoints() # for a voxel == 8 !!

        cellPointIds = vtkIdList()

        NC1 = grid1.GetNumberOfCells()
        NC2 = grid2.GetNumberOfCells()

        setv1 = set()
        setv2 = set()

        dicta = dict()
        try:
            arr1 = grid1.GetPointData().GetArray('muphynode')
        except:
            arr1 = None
        try:
            arr2 = grid2.GetPointData().GetArray('muphynode')
        except:
            arr2 = None

        for i in xrange(NC1):

            grid1.GetCellPoints(i, cellPointIds)

            voxel = []
            for ip in xrange(NPVOX):
                pid = cellPointIds.GetId(ip)
                coo = grid1.GetPoint(pid)
                voxel.append(coo)

                if arr1:
                    dicta[coo] = arr1.GetValue(pid)

            setv1.add( tuple(voxel) )

        for i in xrange(NC2):

            grid2.GetCellPoints(i, cellPointIds)

            voxel = []
            for ip in xrange(NPVOX):
                pid = cellPointIds.GetId(ip) 
                coo = grid2.GetPoint(pid)
                voxel.append(coo)

                if arr2:
                    dicta[coo] = arr2.GetValue(pid)

            setv2.add( tuple(voxel) )

        if what=='union':
            setv = setv1 | setv2

        elif what=='intersect':
            setv = setv1 & setv2

        elif what=='difference':
            setv = setv1 - setv2

        elif what=='symmetricdifference':
            setv = setv1 ^ setv2

        else:
            print 'something wrong:',what
            sys.exit(1)

        previouspointsset = set()
        previouspoints = dict()

        voxel = vtkVoxel()
        points = vtkPoints()

        array = vtkUnsignedShortArray()
        array.SetName('muphynode')
        array.SetNumberOfComponents(1)

        ipp = 0 # total number of points
        ip = 0 # tmp counter
        for pvox in setv:

            ip = 0
            for coo in pvox:

                if not coo in previouspointsset:
                # if not coo in previouspoints.keys():

                    points.InsertNextPoint( coo )
                    array.InsertNextValue( dicta[coo] )
                    voxel.GetPointIds().SetId(ip, ipp)

                    previouspoints[ coo ] = ipp
                    previouspointsset.add(coo)

                    ipp += 1

                else:
                    ipo = previouspoints[coo]
                    voxel.GetPointIds().SetId(ip, ipo)

                ip  += 1

            ugrid.InsertNextCell(voxel.GetCellType(), voxel.GetPointIds())

        ugrid.SetPoints(points)
        ugrid.GetPointData().AddArray(array)

        # to be fixed
        # ugrid = self.clean_unstructuredgrid(ugrid)

        return ugrid


    def bool_unstructuredgrid2(self,grid1,grid2,what):
        """
        boolean of two unstructured grids made of voxels.
        (can be generalized to other type of cells)
        """

        classname = self.__class__.__name__
        funcname = self.bool_unstructuredgrid2.__name__
        # print '....in class:',classname, '....in function:',funcname

        ugrid = vtkUnstructuredGrid()

        voxel = vtkVoxel()
        NPVOX = voxel.GetNumberOfPoints() # for a voxel == 8 !!

        cellPointIds = vtkIdList()

        NC1 = grid1.GetNumberOfCells()
        NC2 = grid2.GetNumberOfCells()

        setv1 = set()
        setv2 = set()
        setp1 = set()
        setp2 = set()

        dicta = dict()
        ddd = dict()
        try:
            arr1 = grid1.GetPointData().GetArray('muphynode')
        except:
            arr1 = None
        try:
            arr2 = grid2.GetPointData().GetArray('muphynode')
        except:
            arr2 = None

        for i in xrange(NC1):

            grid1.GetCellPoints(i, cellPointIds)

            voxel = []
            for ip in xrange(NPVOX):
                pid = cellPointIds.GetId(ip)
                coo = grid1.GetPoint(pid)
                voxel.append(coo)
                setp1.add(coo)
                if arr1:
                    dicta[coo] = arr1.GetValue(ip)

            setv1.add( tuple(voxel) )

        for i in xrange(NC2):

            grid2.GetCellPoints(i, cellPointIds)

            voxel = []
            for ip in xrange(NPVOX):
                pid = cellPointIds.GetId(ip) 
                coo = grid2.GetPoint(pid)
                voxel.append(coo)
                setp2.add(coo)
                if arr2:
                    dicta[coo] = arr2.GetValue(ip)
                    ddd[coo] = arr2.GetValue(ip)
                    # print 'Subtracting....',coo,':',dicta[coo]

            setv2.add( tuple(voxel) )

        if what=='union':
            setv = setv1 | setv2
            setp = setp1 | setp2

        elif what=='intersect':
            setv = setv1 & setv2
            setp = setp1 & setp2

        elif what=='difference':
            setv = setv1 - setv2
            setp = setp1 # | setp2

        elif what=='symmetricdifference':
            setv = setv1 ^ setv2
            setp = setp1 | setp2

        else:
            print 'something wrong:',what
            sys.exit(1)

        voxel = vtkVoxel()
        points = vtkPoints()

        # the following method creates far too many point!!!
        # needs to be cleanup...still not working
        # and array values are screwed up !

        ipp = 0 # total number of points
        ip = 0 # tmp counter
        for pvox in setv:

            ip = 0
            for coo in pvox:
                if coo in setp:

                    if True:
                        points.InsertNextPoint( coo )
                        voxel.GetPointIds().SetId(ip, ipp)
                        ipp += 1
                        ip  += 1
                else:
                    print 'error here !'
                    sys.exit(1)

            ugrid.InsertNextCell(voxel.GetCellType(), voxel.GetPointIds())

        ugrid.SetPoints(points)

        if arr1 and arr2:

            array = vtkUnsignedShortArray()
            array.SetName('muphynode')
            array.SetNumberOfComponents(1)
            array.SetNumberOfValues(ipp)

            for ip in xrange(ipp):
                coo = ugrid.GetPoint(ip)
                array.SetValue(ip, dicta[coo])
                # if coo in ddd.keys():
                #     print 'Subs:',coo,':',ddd[coo]

            ugrid.GetPointData().AddArray(array)

        # to be fixed
        # ugrid = self.clean_unstructuredgrid(ugrid)

        return ugrid

    def bool_unstructuredgrid_nonfluid(self,grid1,grid2,what):
        """
        boolean of two unstructured grids made of vertices only. it neglects info about cell identity
        """

        classname = self.__class__.__name__
        funcname = self.bool_unstructuredgrid_nonfluid.__name__
        # print '....in class:',classname, '....in function:',funcname

        ugrid = vtkUnstructuredGrid()

        cellPointIds = vtkIdList()
        
        setp1 = set()
        setp2 = set()

        for pid in xrange(grid1.GetNumberOfCells()):
            setp1.add( grid1.GetPoint(pid) )

        for pid in xrange(grid1.GetNumberOfCells()):
            setp2.add( grid2.GetPoint(pid) )

        if what=='union':
            setp = setp1 | setp2

        elif what=='intersect':
            setp = setp1 & setp2

        elif what=='difference':
            setp = setp1 # | setp2

        elif what=='symmetricdifference':
            setp = setp1 | setp2

        else:
            print 'something wrong:',what
            sys.exit(1)

        points = vtkPoints()
        for coo in setp:
            points.InsertNextPoint( coo )

        ugrid.SetPoints(points)

        return ugrid



    """
    def union_unstructuredgrid_OLD(self,grid1,grid2):

        classname = self.__class__.__name__
        funcname = self.union_unstructuredgrid_OLD.__name__
        # print '....in class:',classname, '....in function:',funcname

        grid = vtkUnstructuredGrid()

        points = vtkPoints()

        # voxel = vtkCell()
        # voxel = vtkCell3D()
        voxel = vtkVoxel()
        NPVOX = voxel.GetNumberOfPoints() # for a voxel == 8 !!

        cellPointIds = vtkIdList()
        
        NP1 = grid1.GetNumberOfPoints()
        NP2 = grid2.GetNumberOfPoints()
        NC1 = grid1.GetNumberOfCells()
        NC2 = grid2.GetNumberOfCells()

        points = vtkPoints()

        for pid in xrange(NP1):
            coo = grid1.GetPoint(pid)
            points.InsertNextPoint(coo)

        for pid in xrange(NP2):
            coo = grid2.GetPoint(pid)
            points.InsertNextPoint(coo)

        grid.SetPoints(points)

        for i in xrange(NC1):

            cell = grid1.GetCell(i)
            grid1.GetCellPoints(i, cellPointIds)

            for ip in xrange(NPVOX):
                pid = cellPointIds.GetId(ip)
                coo = grid1.GetPoint(pid)
                voxel.GetPointIds().SetId(ip, pid)
                # print ip,pid,coo

            grid.InsertNextCell(voxel.GetCellType(), voxel.GetPointIds())

        for i in xrange(NC2):

            cell = grid2.GetCell(i)
            grid2.GetCellPoints(i, cellPointIds)

            for ip in xrange(NPVOX):
                pid = cellPointIds.GetId(ip)
                coo = grid2.GetPoint(pid)
                voxel.GetPointIds().SetId(ip, pid + NP1)
                # print ip,pid+NP1,coo

            grid.InsertNextCell(voxel.GetCellType(), voxel.GetPointIds())

        #print 'NP1,NP2,NPfin:',NP1,NP2,grid.GetNumberOfPoints()
        #print 'NC1,NC2,NCfin:',NC1,NC2,grid.GetNumberOfCells()

        return grid

    #
    # boolean intersect two unstructured grids made of voxels.
    # can be generalized to other type of cells
    def bool_unstructuredgrid_OLD(self,grid1,grid2,what):

        classname = self.__class__.__name__
        funcname = self.bool_unstructuredgrid_OLD.__name__
        # print '....in class:',classname, '....in function:',funcname

        grid = vtkUnstructuredGrid()

        points = vtkPoints()

        # voxel = vtkCell()
        # voxel = vtkCell3D()
        voxel = vtkVoxel()
        NPVOX = voxel.GetNumberOfPoints() # for a voxel == 8 !!

        cellPointIds = vtkIdList()
        
        NP1 = grid1.GetNumberOfPoints()
        NP2 = grid2.GetNumberOfPoints()
        NC1 = grid1.GetNumberOfCells()
        NC2 = grid2.GetNumberOfCells()

        points = vtkPoints()

        for pid in xrange(NP1):
            coo = grid1.GetPoint(pid)
            points.InsertNextPoint(coo)

        for pid in xrange(NP2):
            coo = grid2.GetPoint(pid)
            points.InsertNextPoint(coo)

        grid.SetPoints(points)

        gp = [-1,-1,-1,-1,-1,-1,-1,-1]

        gp1 = []
        for i in xrange(NC1):
            cell = grid1.GetCell(i)
            grid1.GetCellPoints(i, cellPointIds)
            gp = []
            for ip in xrange(NPVOX):
                pid = cellPointIds.GetId(ip)
                coo = grid1.GetPoint(pid)
                #voxel.GetPointIds().SetId(ip, pid)
                # gp[ip] = coo
                gp.append(coo)
            gp1.append(gp)

        gp2 = []
        for i in xrange(NC2):
            cell = grid2.GetCell(i)
            grid2.GetCellPoints(i, cellPointIds)
            gp = []
            for ip in xrange(NPVOX):
                pid = cellPointIds.GetId(ip)
                coo = grid2.GetPoint(pid)
                #voxel.GetPointIds().SetId(ip, pid)
                # gp[ip] = coo
                gp.append(coo)
            gp2.append(gp)

        for i in xrange(NC1):

            cell = grid1.GetCell(i)

            grid1.GetCellPoints(i, cellPointIds)

            gp = []
            for ip in xrange(NPVOX): 
                pid = cellPointIds.GetId(ip)
                coo = grid1.GetPoint(pid)
                # gp[ip] = coo
                gp.append(coo)

            if what=='xor':
                if gp in gp2: 
                    continue
            elif what=='intersect':
                if not gp in gp2: 
                    continue
            elif what=='exclude':
                if gp in gp2: 
                    continue

            for ip in xrange(NPVOX):
                pid = cellPointIds.GetId(ip)
                coo = grid1.GetPoint(pid)
                voxel.GetPointIds().SetId(ip, pid)
                # print ip,pid,coo

            grid.InsertNextCell(voxel.GetCellType(), voxel.GetPointIds())

        if what=='exclude':
            return grid

        for i in xrange(NC2):

            cell = grid2.GetCell(i)

            grid2.GetCellPoints(i, cellPointIds)

            gp = []
            for ip in xrange(NPVOX): 
                pid = cellPointIds.GetId(ip)
                coo = grid2.GetPoint(pid)
                # gp[ip] = coo
                gp.append(coo)

            if what=='xor':
                if gp in gp1: 
                    continue
            elif what=='intersect':
                if not gp in gp1: 
                    continue


            for ip in xrange(NPVOX):
                pid = cellPointIds.GetId(ip)
                coo = grid2.GetPoint(pid)
                voxel.GetPointIds().SetId(ip, pid + NP1)
                # print ip,pid+NP1,coo

            grid.InsertNextCell(voxel.GetCellType(), voxel.GetPointIds())

        #print 'NP1,NP2,NPfin:',NP1,NP2,grid.GetNumberOfPoints()
        #print 'NC1,NC2,NCfin:',NC1,NC2,grid.GetNumberOfCells()

        return grid
    """

    def make_simple_unstructuredgrid(self,ijk):
        """
        create and return unstructured vtk grid from ijk (NON-numpy): only points and no cells
        """

        classname = self.__class__.__name__
        funcname = self.make_simple_unstructuredgrid.__name__

        ugrid = vtkUnstructuredGrid()

        points = vtkPoints()

        NN = len(ijk[:,0])
     
        if NN<=0:
            print 'Error: zero size of ijk array in class:',classname
            sys.exit(1)
     
        # old style, only points/vertexes - no cells
        vert = vtkVertex()
        cells = vtkCellArray()

        for n in xrange(NN):

            i,j,k = ijk[n]
            points.InsertNextPoint(i,j,k)
            vert.GetPointIds().SetId(0,n)
            cells.InsertNextCell(vert)

        ugrid.SetPoints(points)

        ugrid.SetCells(VTK_VERTEX,cells)
     
        return ugrid

    def hashcols(self):
        """
        associate node index to scalar values
        """

        # unfortunately th lookuptable only works in the 0-1 range
        # self.node2hashval = {1: 0.2, 2: 0.4, 3: 0.6, 4: 0.8}
        # self.node2hashval = {1: 1.0, 2: 2.0, 3: 3.0, 4: 4.0}
        self.node2hashval = {1: 1, 2: 2, 3: 3, 4: 4}

        """
        # commented out: the default lookuptable is nice enough
        self.node2colors =  {1 : [ 0.00, 1.00, 0.00, 1.00 ], # left-most floats:RGB right-most:alpha channel
                             2 : [ 0.00, 0.00, 1.00, 0.50 ],
                             3 : [ 1.00, 1.00, 0.00, 1.00 ],
                             4 : [ 0.68, 0.42, 0.51, 1.00 ]}
        """

    def writeVTK(self, fname, ijkf, fieldval=None):
        """
        write legacy .vtk file from i,j,k,flag quadruplet (numpy array)
        """
     
        outf = open(fname, 'w')
     
        NN = ijkf.shape[0]

        print >> sys.stderr, "Writing VTK file... : ",fname
     
        print >> outf, '# vtk DataFile Version 3.0'
        print >> outf, 'output mesh_data Mauro 1.0'
        print >> outf, 'ASCII'
        print >> outf, 'DATASET UNSTRUCTURED_GRID'
        print >> outf, 'POINTS', NN, "float"
     
        for index in xrange(NN):
            print >> outf, '{0} {1} {2}'.format(ijkf[index,0], ijkf[index,1], ijkf[index,2])
     
        print >> outf, ''
        print >> outf, 'CELLS', NN, NN*2
        for index in xrange(NN):
            print >> outf, '1 {0}'.format(index)
     
        print >> outf, ''
        print >> outf, 'CELL_TYPES', NN
        for index in xrange(NN):
            print >> outf, vtk.VTK_VERTEX
     
        """
        # unneeded...everything comes from the verteces
        print >> outf, ''
        print >> outf, 'CELL_DATA', NN
        # print >> outf, 'SCALARS muphynode short 1'
        print >> outf, 'SCALARS muphynode float 1'

        print >> outf, 'LOOKUP_TABLE node_colors'
        # print >> outf, 'LOOKUP_TABLE default'
        for index in xrange(NN):
            print >> outf, '{0:1.1f}'.format(self.node2hashval[ijkf[index,3]]/4.)
        """

        print >> outf, ''
        print >> outf, 'POINT_DATA', NN
        print >> outf, 'SCALARS muphynode short'
        # print >> outf, 'SCALARS muphynode int 1'
        # print >> outf, 'SCALARS muphynode float 1'

        # The lookup table is mapped on data values between 0 and 1. So if you want 
        # to use a lookup table, you have to scale your data to the range [0,1]. The 
        # entries of the lookup table are spread equally over this range; so if you 
        # have three entries, they will be mapped to the data values 0, 0.5 and 1. 
        # Data values smaller than 0.25 are displayed in the first color, values bigger 
        # than 0.75 are displayed in the third color, and the values in between in 
        # the second color.
        #
        # Since we use an hash table for colors we scale node-type integer value
        # to float in [0,1].
        # print >> outf, 'LOOKUP_TABLE node_colors'
        print >> outf, 'LOOKUP_TABLE default'
        for index in xrange(NN):
            # print >> outf, '{0:.2}'.format(self.node2hashval[ijkf[index,3]])
            print >> outf, '{0:1}'.format(self.node2hashval[ijkf[index,3]])
            # print >> outf, '{0:1.1f}'.format(self.node2hashval[ijkf[index,3]])

        """
        print >> outf, 'LOOKUP_TABLE node_colors', len(self.node2colors)
        for index in xrange(1, len(self.node2colors)+1):
            print >> outf, '{0:.2f} {1:.2f} {2:.2f} {3:.2f}'.format(self.node2colors[index][0]/4.0, \
                                                                    self.node2colors[index][1]/4.0, \
                                                                    self.node2colors[index][2]/4.0, \
                                                                    self.node2colors[index][3]/4.0)
        """

        if fieldval != None:
            print >> outf, 'SCALARS fieldnode short'
            print >> outf, 'LOOKUP_TABLE default'
            for index in xrange(NN):
                print >> outf, fieldval
                # print >> outf, '{0:1.1f}'.format(self.node2hashval[ijkf[index,3]])

        outf.close()

    def assembleVTP(self, ijkf, spacing=1, fieldval=None):
  
        NN = ijkf.shape[0]

        points = vtkPoints()

        pnodetype = vtkUnsignedShortArray()
        pnodetype.SetName("muphynode")
        pnodetype.SetNumberOfComponents(1)

        if fieldval != None:
            fval = vtkUnsignedShortArray()
            fval.SetName("fieldval")
            fval.SetNumberOfComponents(1)

        for index in range(NN):
            i, j, k, flg = ijkf[index,0], ijkf[index,1], ijkf[index,2], ijkf[index,3]
            points.InsertNextPoint (i, j, k)
            pnodetype.InsertNextValue( flg )
            if fieldval != None:
                fval.InsertNextValue( fieldval )

        polydata = vtkPolyData()
        polydata.SetPoints(points)
        polydata.GetPointData().AddArray(pnodetype)
        if fieldval != None:
            polydata.GetPointData().AddArray(fval)

        return polydata

    def assembleVTU(self, ijkf, spacing=1):
        """
        assemble the numpy ijk array as an unstructured mesh
        """

        if self.nx + self.ny + self.nz == 0:
            print 'Error: zero nx,ny,nz in writeVTU!'
            sys.exit(1)

        ugrid = vtkUnstructuredGrid()

        points = vtkPoints()
        # vtkCellArray is a supporting object that explicitly represents cell connectivity.
        # The cell array structure is a raw integer list of the form:
        # (n,id1,id2,...,idn, n,id1,id2,...,idn, ...) where n is the number of points in
        # the cell, and id is a zero-offset index into an associated point list.
        voxel = vtkVoxel()

        pnodetype = vtkUnsignedShortArray()
        cnodetype = vtkUnsignedShortArray()
        # cnodetype = vtkIntArray()
        # cnodetype = vtkFloatArray()
     
        NN = ijkf.shape[0]

        pnodetype.SetNumberOfComponents(1)
        pnodetype.SetName("muphynode");
     
        TMP = False # True

        if TMP:
            vert = vtkVertex()
            cells = vtkCellArray()

        for index in xrange(NN):

            i,j,k,f = ijkf[index]

            points.InsertNextPoint(i,j,k)

            pnodetype.InsertNextValue(f)

            if TMP:
                vert.GetPointIds().SetId(0,index)
                cells.InsertNextCell(vert)

        ugrid.SetPoints(points)

        cnodetype.SetNumberOfComponents(1)
        cnodetype.SetName("muphynode");

        if TMP:
            ugrid.SetCells(VTK_VERTEX,cells)

        Lib = ctypes.CDLL( os.path.join( os.getenv("MOEBIUS_BIN"), "TOOLS/points2cells.so") )
     
        ijk1 = (ctypes.c_int*NN)()
        ijk2 = (ctypes.c_int*NN)()
        ijk3 = (ctypes.c_int*NN)()
        ncell = (ctypes.c_int*1)()
        icell = (ctypes.c_int*(8*NN))()
     
        nodedict = {}
        for index in xrange(NN):

            ijk1[index], ijk2[index], ijk3[index], f = ijkf[index]

            # nodedict[index] = f
            nodedict[ tuple( [ijk1[index], ijk2[index], ijk3[index]]) ] = f
     
        # construct cell according to 0,+1 directions (fowrad only)
        Lib.points2cells(ctypes.c_int(NN),ctypes.c_int(self.nx), ctypes.c_int(self.ny), ctypes.c_int(self.nz), 
                         ctypes.c_int(spacing),
                         ctypes.byref(ijk1), ctypes.byref(ijk2), ctypes.byref(ijk3), 
                         ctypes.byref(ncell), ctypes.byref(icell) )
     
        for i in xrange(ncell[0]):
     
            index  = [ icell[8*i    ], icell[8*i + 1], icell[8*i + 2], icell[8*i + 3],
                       icell[8*i + 4], icell[8*i + 5], icell[8*i + 6], icell[8*i + 7] ]

            ii,jj,kk,f = ijkf[ icell[8*i] ]

            # this order must follow the same as in the ctype lib
            ijkcell  = [ tuple([ii  ,       jj  ,       kk  ]), \
                         tuple([ii+spacing, jj  ,       kk  ]), \
                         tuple([ii  ,       jj+spacing, kk  ]), \
                         tuple([ii  ,       jj  ,       kk+spacing]), \
                         tuple([ii+spacing, jj+spacing, kk  ]), \
                         tuple([ii+spacing, jj  ,       kk+spacing]), \
                         tuple([ii  ,       jj+spacing, kk+spacing]), \
                         tuple([ii+spacing, jj+spacing, kk+spacing]) ]

            # paint cell according to zero-th vertex
            # cnodetype.InsertNextValue( nodedict[icell[8*i]] )

            col = -1
            allfluid = True
            # paint cell as fluid if all verts are fluid
            for j in xrange(8):
                if not nodedict[ tuple(ijkcell[j]) ] == self.ID_FLUIDNODE: 
                    allfluid = False
                    break
            if allfluid:
                col = self.ID_FLUIDNODE

            if col<0:
                # paint as inlet if any vert is inlet
                for j in xrange(8):
                    if not nodedict[ tuple(ijkcell[j]) ] == self.ID_INLETNODE: 
                        col = self.ID_INLETNODE
                        break
            if col<0:
                # paint as outlet if any vert is outlet
                for j in xrange(8):
                    if not nodedict[ tuple(ijkcell[j]) ] == self.ID_OUTLETNODE: 
                        col = self.ID_OUTLETNODE
                        break
            if col<0:
                # paint as wall if any vert is wall
                for j in xrange(8):
                    if not nodedict[ tuple(ijkcell[j]) ] == self.ID_WALLNODE: 
                        col = self.ID_WALLNODE
                        break

            if col<0:
                print 'Color unset !!!! cell topology error'
                for j in range(8):
                    print j,'index:',index[j], 'id:',nodedict[index[j]]
                print
                sys.exit(1)

            cnodetype.InsertNextValue( col )

     
            #idx = index[0]
            #points.InsertNextPoint( ijk1[idx], ijk2[idx], ijk3[idx] )
     
            #for idx in index:
            #    points.InsertNextPoint( ijk1[idx], ijk2[idx], ijk3[idx] )
     
            for j in range(8):
                voxel.GetPointIds().SetId(j, index[j])
     
            if i%500000==0: print 'Inserting cell at ',index,' / ncell:',ncell[0]
     
            ugrid.InsertNextCell(voxel.GetCellType(), voxel.GetPointIds())
     
        # nodetype scalar on cells
        ugrid.GetCellData().AddArray(cnodetype);

        # nodetype scalar on points
        ugrid.GetPointData().AddArray(pnodetype);

        # ugrid.GetFieldData().AddArray(cnodetype);

        #cleanFilter = vtkCleanPolyData()
        ## cleanFilter.SetInputConnection(ugrid.GetOutputPort());
        #cleanFilter.SetInput(ugrid)
        #cleanFilter.Update()

        return ugrid
     
    def meshfileToPolyData(self, MESHFNAME, spacing=1, fieldval=None):

        # open the .dat file
        try:
            meshf = open(MESHFNAME + '.dat','r')
            print >> sys.stderr, 'Reading mesh file...', MESHFNAME + '.dat'
            # start_time = time.time()
        except IOError:
            print >> sys.stderr, "Cannot find file", MESHFNAME+'.dat'
            sys.exit(1)

        # msh = Mesh()

        self.hashcols()

        meshdata = list()
        meshdataWIO = list()

        # start_time = time.time()


        for line in meshf.readlines():

            fields = tuple(line.split())

            meshdata.append( (int(fields[0]), int(fields[1]), int(fields[2]), int(fields[3])) )

            if int(fields[3]) == 1: continue # exclude fluid nodes
            # if int(fields[3]) == 3: continue # exclude inlet nodes
            # if int(fields[3]) == 4: continue # exclude outlet nodes

            meshdataWIO.append( (int(fields[0]), int(fields[1]), int(fields[2]), int(fields[3])) )

        meshf.close()
        # print >> sys.stderr, "done", time.time()-start_time, "secs"

        # write_time = time.time()
        # print >> sys.stderr, "done", time.time()-write_time, "secs"
        # print >> sys.stderr, "Total time:", time.time()-start_time, "secs"

        ijkf = np.zeros([len(meshdataWIO),4],dtype='int')
        n = 0
        for i,j,k,f in meshdataWIO:
            ijkf[n,:] = i,j,k,f
            n += 1

        pd = self.assembleVTP(ijkf, spacing=spacing, fieldval=fieldval)

        return pd

    def writeCombinedVTP(self, vtpfile, pds):

        appendFilter = vtkAppendPolyData()
        for pd in pds:
            appendFilter.AddInputData(pd)
        appendFilter.Update();

        writer = vtkXMLPolyDataWriter()
        writer.SetFileName(vtpfile)
        writer.SetDataModeToAscii()
        writer.SetInputData(appendFilter.GetOutput())
        writer.Write()

    def writeVTP(self, vtpfile, ijkf, spacing=1, fieldval=None):
        """
        write a numpy ijk array as a flat vtk unstructured mesh
        """

        polydata = self.assembleVTP(ijkf, spacing=spacing, fieldval=fieldval)

        writer = vtkXMLPolyDataWriter()
        writer.SetFileName(vtpfile)
        writer.SetDataModeToAscii()
        writer.SetInputData(polydata)
        writer.Write()

    def writeVTU(self,vtufile,ijkf,spacing=1):
        """
        write a numpy ijk array as a flat vtk unstructured mesh
        """

        ugrid = self.assembleVTU(ijkf, spacing)

        writer = vtkXMLUnstructuredGridWriter()
        writer.SetFileName(vtufile)
        writer.SetDataModeToAscii()
        writer.SetInputData(ugrid)
        writer.Write()

    def assembleVTM(self):
        """
        assemble a multiblock with fluids, wall, inlets, outlets and arrays to identify them
        """

        print 'assemblying multiblock... '

        if len(self.ijk_f) + len(self.ijk_w) + len(self.ijk_i) + len(self.ijk_o) == 0:

            if len(self.itppp_f) + len(self.itppp_w) + len(self.itppp_i) + len(self.itppp_o) == 0:
                print 'Error: both itppp and ijk arrays are empty in _writeVTM !'
                sys.exit(1)

            # convert itppp into ijk
            self.ijk_f = np.zeros( [self.nfluid, 3], dtype='int' )
            self.ijk_w = np.zeros( [self.nwall,  3], dtype='int' )
            self.ijk_i = np.zeros( [self.ninlet, 3], dtype='int' )
            self.ijk_o = np.zeros( [self.noutlet,3], dtype='int' )

            n = 0
            for i4 in self.itppp_f:
                self.ijk_f[n,:] = self.ijk(i4)
                n += 1
            if not n==self.nfluid: print 'Warning on fluid count',n,self.nfluid

            n = 0
            for i4 in self.itppp_w:
                self.ijk_w[n,:] = self.ijk(i4)
                n += 1
            if not n==self.nwall: print 'Warning on wall count',n,self.nwall

            n = 0
            for i4 in self.itppp_i:
                self.ijk_i[n,:] = self.ijk(i4)
                n += 1
            if not n==self.ninlet: print 'Warning on inlet count',n,self.ninlet

            n = 0
            for i4 in self.itppp_o:
                self.ijk_o[n,:] = self.ijk(i4)
                n += 1
            if not n==self.noutlet: print 'Warning on outlet count',n,self.noutlet
                

        if self.nx + self.ny + self.nz == 0:
            self.nx = max(np.amax(self.ijk_f[:,0]), np.amax(self.ijk_w[:,0]), np.amax(self.ijk_i[:,0]), np.amax(self.ijk_o[:,0]))
            self.ny = max(np.amax(self.ijk_f[:,1]), np.amax(self.ijk_w[:,1]), np.amax(self.ijk_i[:,1]), np.amax(self.ijk_o[:,1]))
            self.nz = max(np.amax(self.ijk_f[:,2]), np.amax(self.ijk_w[:,2]), np.amax(self.ijk_i[:,2]), np.amax(self.ijk_o[:,2]))

        print '... nx,ny,nz',self.nx,self.ny,self.nz

        # for fluid show the whole cells
        if self.nfluid>0:  
            # fluid  = self.make_simple_unstructuredgrid(self.ijk_f)
            fluid  = self.make_unstructuredgrid(self.ijk_f) # fluids : cell-like
            self.addarray_unstructuredgrid(fluid,self.ID_FLUIDNODE)

        # for walls show only points
        if self.nwall>0:   
            # wall   = self.make_unstructuredgrid(self.ijk_w)
            wall   = self.make_simple_unstructuredgrid(self.ijk_w) # nonfluids : point-like
            self.addarray_unstructuredgrid(wall,self.ID_WALLNODE)

        # for inlets show only points
        if self.ninlet>0:  
            # inlet  = self.make_unstructuredgrid(self.ijk_i)
            inlet  = self.make_simple_unstructuredgrid(self.ijk_i) # nonfluids : point-like
            self.addarray_unstructuredgrid(inlet,self.ID_INLETNODE)

        # for outlets show only points
        if self.noutlet>0: 
            # outlet = self.make_unstructuredgrid(self.ijk_o)
            outlet = self.make_simple_unstructuredgrid(self.ijk_o) # nonfluids : point-like
            self.addarray_unstructuredgrid(outlet,self.ID_OUTLETNODE)
     
        # gather all into the multiblock structure
        group = vtkMultiBlockDataGroupFilter()
     
        if self.nfluid>0:  group.AddInputData(fluid)
        if self.nwall>0:   group.AddInputData(wall)
        if self.ninlet>0:  group.AddInputData(inlet)
        if self.noutlet>0: group.AddInputData(outlet)
        # group.AddInputConnection(fluid.GetOutputPort())

        group.Update()

        return group
     
    def writeVTM(self,fileout):
        """
        write a multiblock .vtm file
        """

        print 'writing multiblock .vtm file ... '
        group = self.assembleVTM()

        writer = vtkXMLMultiBlockDataWriter()
        writer.SetFileName(fileout)
        writer.SetInputConnection(group.GetOutputPort())
        writer.Update()

    def convert_ijk_2_itppp(self,ijk):
        """
        """

        itppp = np.zeros(ijk.shape[0],dtype='i8')
        n = 0
        for n in ijk.shape[0]:
            itppp[n] = self.i4back( ijk[n] )
            n += 1
        return itppp

    def convert_itppp_2_ijk(self,itppp):
        """
        """

        n = itppp.shape[0]
        ijk = np.zeros((n,3),dtype='i')
        n = 0
        for i4 in itppp:
            ijk[n] = self.ijk(i4)
            n += 1
        return ijk


if __name__ == '__main__':

    msh = Mesh()

    grid1 = msh.loadfile_unstructuredgrid_VTK('A.vtk')
    grid2 = msh.loadfile_unstructuredgrid_VTK('B.vtk')

    grid = msh.bool_unstructuredgrid(grid1,grid2,'union')
    msh.writefile_unstructuredgrid_VTK('AB.vtk',grid)

    grid = msh.bool_unstructuredgrid(grid1,grid2,'intersect')
    msh.writefile_unstructuredgrid_VTK('ABi.vtk',grid)

    grid = msh.bool_unstructuredgrid(grid1,grid2,'difference')
    msh.writefile_unstructuredgrid_VTK('ABe.vtk',grid)

    grid = msh.bool_unstructuredgrid(grid1,grid2,'symmetricdifference')
    msh.writefile_unstructuredgrid_VTK('ABx.vtk',grid)


    """
    # for the fun of it....remove faces from cubes

    gf1 = vtkGeometryFilter()
    gf2 = vtkGeometryFilter()
    if VTK_MAJOR_VERSION <= 5:
        gf1.SetInput(grid1)
        gf2.SetInput(grid2)
    else:
        gf1.SetInputData(grid1)
        gf2.SetInputData(grid2)
    gf1.Update(); 
    gf2.Update(); 

    poly1 = gf1.GetOutput()
    poly2 = gf2.GetOutput()

    poly1.BuildLinks()
    poly1.DeleteCell(1)
    poly1.DeleteCell(2)
    poly1.DeleteCell(3)
    poly1.DeleteCell(4)
    poly1.RemoveDeletedCells()

    grid2 = vtkUnstructuredGrid()
    grid2.ShallowCopy(poly1)
    grid2.Update()

    # m.writefile_unstructuredgrid_VTK('AB2.vtk',grid2)

    writer =  vtkXMLPolyDataWriter()
    writer.SetFileName("A.vtp");
    if VTK_MAJOR_VERSION <= 5:
        writer.SetInput(poly1);
    else:
        writer.SetInputData(pol1);
    #writer.SetDataModeToBinary();
    #writer.SetDataModeToAscii();
    writer.Write();

    """
