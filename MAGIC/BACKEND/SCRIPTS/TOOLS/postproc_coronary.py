#!/usr/bin/env python

import sys,os,argparse
import numpy as np
from myvtk import *
from math import *
import ConfigParser
# import quadrature
import mygrace

from muphy2wrapper import *

VOXCM = 200         # voxels per cm
VCH = 0.2           # characteristic velocity in m/s # optional change
MACH = 0.025        # characteristic Mach # optional change

RHO = 0.22 # lattice density - set as in powerflow
MASSDENSITY = 1000. # blood mass density in MKS
VISCOSITY = 4.e-6   # blood viscosity in MKS
PCH = 100.           # characteristic pressure in mmHg

SLICE_ANALYSIS_FREQ = 10

# CLIP_IT = True # slow ... but will be slightly better at bifurcations
CLIP_IT = False # fast

WRITE_SINGLE_SLICE_FILES_WITH_OBS = False
WRITE_CLINE_SLICE_FILES_WITH_OBS = True
 
#######################################################

def convert(ug,ARROLD,ARRNEW,factor,compo):

    arr = vtkFloatArray()
    arr.SetNumberOfComponents(compo)
    arr.SetName(ARRNEW)

    # print ug.GetCellData().GetArray(ARROLD)
    for cid in xrange(ug.GetNumberOfCells()):

        # qty = ug.GetCellData().GetArray(ARROLD).GetValue(cid)
        qty = ug.GetCellData().GetArray(ARROLD).GetTuple(cid)

        for q in qty:
            arr.InsertNextValue( q * factor )

        # if cid<30000 and cid%5000==0:
        #     print cid,ARROLD,':',qty

    return arr

#######################################################

def convert_cpmapping(ug, um, name):

    arr = vtkFloatArray()
    arr.SetNumberOfComponents(1)
    arr.SetName(name)

    PaTommHg = 0.007500616827042

    PCH = um.getCharacteristicPressure_mmHg()

    fact = PaTommHg * (um.PhysicalMassDensity * um.Vch**2) / (RHO * MACH**2)

    # Pressure coefficient cp mapping 
    for cid in xrange(ug.GetNumberOfCells()):

        rho = ug.GetCellData().GetArray('density').GetValue(cid)

        # p = PCH + fact * (rho - RHOIN)
        p = PCH + fact * (rho - RHO)

        arr.InsertNextValue( p )

    return arr


#######################################################

def SliceIt(xf,ugn,transform,kk,polyslice):

    print '...slicing along:',xf

    # the XSection is only used to take the centroid and compute the normal
    xreader = vtkPolyDataReader()
    xreader.SetFileName(xf)
    xreader.Update()

    xsect = vtkTransformPolyDataFilter()
    xsect.SetInputConnection(xreader.GetOutputPort())
    xsect.SetTransform(transform)
    xsect.Update()

    # first point of XSections is reserved for centerline point (centroid)
    npts_contour = xsect.GetOutput().GetNumberOfPoints()

    coo0 = xsect.GetOutput().GetPoint( 0 )
    coo1 = xsect.GetOutput().GetPoint( 1 )
    # coo2 = xsect.GetOutput().GetPoint( 2 )
    coo2 = xsect.GetOutput().GetPoint( 3 )

    AA = coo1[0]-coo0[0], coo1[1]-coo0[1], coo1[2]-coo0[2]
    BB = coo2[0]-coo0[0], coo2[1]-coo0[1], coo2[2]-coo0[2]

    dx = AA[1]*BB[2] - AA[2]*BB[1]
    dy = AA[2]*BB[0] - AA[0]*BB[2]
    dz = AA[0]*BB[1] - AA[1]*BB[0]

    dd = sqrt(dx**2 + dy**2 + dz**2)

    if dd<1.e-6:
        print 'kk:',kk,'dx,dy,dz',dx,dy,dz; print 'coo0:',coo0; print 'coo1:',coo1; print 'coo2:',coo2

    norm = [dx/dd, dy/dd, dz/dd]

    # pts = vtkPoints()
    # for pid in xrange(npts_contour):
    #     coo = xsect.GetOutput().GetPoint( pid )
    #     pts.InsertNextPoint(coo)

    plane = vtkPlane()
    coo = xsect.GetOutput().GetPoint( 0 )
    plane.SetOrigin(coo)
    plane.SetNormal(norm)

    # cube = vtkCylinder()
    cube = vtkSphere()
    cube.SetCenter(coo)
    # cube.SetRadius(10) # clip the data inside a cube. Should be large enough in lattice units
    # cube.SetRadius(20) # clip the data inside a cube. Should be large enough in lattice units
    cube.SetRadius(30) # clip the data inside a cube. Should be large enough in lattice units
    # cube.SetRadius(100) # clip the data inside a cube. Should be large enough in lattice units

    if CLIP_IT:
        clipper = vtkClipDataSet()
        clipper.SetClipFunction(cube)
        if VTK_MAJOR_VERSION <= 5:
            clipper.SetInput(ugn)
        else:
            clipper.SetInputData(ugn)
        clipper.InsideOutOn()
        clipper.Update()

    cutter = vtkCutter()
    cutter.SetCutFunction(plane)

    if CLIP_IT:
        cutter.SetInputConnection(clipper.GetOutputPort())
    else:
        cutter.SetInputData(ugn)

    cutter.Update()

    cellPointIds = vtkIdList()

    arr_dens = cutter.GetOutput().GetCellData().GetArray('density_MKS')
    arr_press = cutter.GetOutput().GetCellData().GetArray('pressure_mmHg')
    arr_vel = cutter.GetOutput().GetCellData().GetArray('velocity_MKS')

    area,dens,press,flowrate = 0,0,0,0
    for i in xrange(cutter.GetOutput().GetNumberOfCells()):

        cell = cutter.GetOutput().GetCell(i)
        cutter.GetOutput().GetCellPoints(i, cellPointIds)
        if cell.GetNumberOfPoints() != 3: 
            print 'Not A Triangle'
            sys.exit(1)

        V0 = cutter.GetOutput().GetPoint( cellPointIds.GetId(0) )
        V1 = cutter.GetOutput().GetPoint( cellPointIds.GetId(1) )
        V2 = cutter.GetOutput().GetPoint( cellPointIds.GetId(2) )

        AA = V1[0]-V0[0], V1[1]-V0[1], V1[2]-V0[2]
        BB = V2[0]-V0[0], V2[1]-V0[1], V2[2]-V0[2]

        dx = AA[1]*BB[2] - AA[2]*BB[1]
        dy = AA[2]*BB[0] - AA[0]*BB[2]
        dz = AA[0]*BB[1] - AA[1]*BB[0]

        dA = 0.5 * sqrt( dx**2 + dy**2 + dz**2 )
        dA *= (1.e-3/scale)**2 # transform to m

        area += dA

        dens += arr_dens.GetValue(i) * dA

        press += arr_press.GetValue(i) * dA

        flowrate += (norm[0]*arr_vel.GetComponent(i,0) + \
                     norm[1]*arr_vel.GetComponent(i,1) + \
                     norm[2]*arr_vel.GetComponent(i,2)) * dA

    if area<=0: return None, None, None, None, None

    dens /= area
    press /= area
    # flowrate /= area
    avevel = flowrate / area

    # area in mm, dens in Kg/m3, press in mmHg, flowrate in ml/min, avevel in m/s

    if vtk.VTK_MAJOR_VERSION <= 5:
        # polyslice.AddInputConnection(cutter.GetProducerPort())
        polyslice.AddInputConnection(cutter.GetOutputPort())
        # polyslice.AddInput(cutter.GetOuput())
    else:
        polyslice.AddInputConnection(cutter.GetOutputPort())
        # polyslice.AddInput(cutter.GetOutput())
        # polyslice.AddInputData(cutter)
        # polyslice.AddInputData(cutter.GetOuput())

    polyslice.Update()

    if WRITE_SINGLE_SLICE_FILES_WITH_OBS:

        # write the slice_ files
        writer = vtkPolyDataWriter()
        writer.SetFileName( os.path.join(cdr, 'slice_'+ os.path.split(xf)[1]) )

        if VTK_MAJOR_VERSION <= 5:
            writer.SetInput(cutter.GetOutput())
        else:
            writer.SetInputData(cutter.GetOutput())
        writer.Write()

    """

    # computation of cross-section integrals via quadratures

    ntriquad = 0
    area_xsection_integral = 0; rho_xsection_integral = 0; flux_xsection_integral = 0

    # triangles = vtkCellArray()
    # polytriangle = vtkPolyData()

    # Loop over triangles in a cross section
    for pid in xrange(1,npts_contour+1):

        tri = vtkTriangle()
        tri.GetPointIds().SetId(0,0)
        tri.GetPointIds().SetId(1,(pid-1) % (npts_contour-1) + 1)
        tri.GetPointIds().SetId(2,(pid) % (npts_contour-1) + 1)
        # triangles.InsertNextCell(tri)

        T = quadrature.Triangle()

        coo = xsect.GetOutput().GetPoint( 0 )
        T.x1 = coo[0]; T.y1 = coo[1]; T.z1 = coo[2]

        coo = xsect.GetOutput().GetPoint( (pid-1) % (npts_contour-1) + 1 )
        T.x2 = coo[0]; T.y2 = coo[1]; T.z2 = coo[2]

        coo = xsect.GetOutput().GetPoint( (pid) % (npts_contour-1) + 1 )
        T.x3 = coo[0]; T.y3 = coo[1]; T.z3 = coo[2]

        probePoints = vtkPoints()

        for q in T.GetQuadraturePoints():
            probePoints.InsertNextPoint(q)

        probePolyData = vtkPolyData()
        probePolyData.SetPoints(probePoints)

        xprobe = vtkProbeFilter()
        if VTK_MAJOR_VERSION <= 5:
            xprobe.SetInput(probePolyData)
            xprobe.SetSource(ugn)
        else:
            xprobe.SetInputData(probePolyData)
            xprobe.SetSourceData(ugn)
        xprobe.PassPointArraysOn()
        xprobe.Update()

        valids = xprobe.GetValidPoints()

        ntup = valids.GetNumberOfTuples()

        if ntup != 4: continue

        arr_dens = xprobe.GetOutput().GetPointData().GetArray('density_MKS')
        arr_vel = xprobe.GetOutput().GetPointData().GetArray('velocity_MKS')
        dens = []
        vel = []
        for pid in xrange(ntup):
            dens.append( arr_dens.GetValue(pid) )
            vel.append( [arr_vel.GetComponent(pid,0),
                         arr_vel.GetComponent(pid,1),
                         arr_vel.GetComponent(pid,2)] )

        ntriquad += 1

        if ntup != 4: continue

        area_tri, rho_tri, flux_tri = T.GetQuadrature(dens,vel)

        area_xsection_integral += area_tri
        rho_xsection_integral += rho_tri
        flux_xsection_integral += flux_tri

    # print ' no. of sampled triangles:',ntriquad, 'no. of triangles:',npts_contour-2

    if ntriquad != npts_contour-2: continue

    # print '   INTs:',rho_xsection_integral,flux_xsection_integral

    print kk, area_xsection_integral, rho_xsection_integral, flux_xsection_integral
    print >> fa, kk, area_xsection_integral, rho_xsection_integral, flux_xsection_integral
    fa.flush()
    """

    return area*1.e+6, dens, press, 60.e+6*flowrate, avevel

#######################################################
def read_ugrid(filename):

    filenametrunc, extension = os.path.splitext(filename)

    print '\nReading...',filename

    if extension == '.vtk':

       reader = vtkUnstructuredGridReader()
       reader.SetFileName(filename)
       reader.Update()
       ug = reader.GetOutput()

    elif extension == '.pvtu':

       reader = vtkXMLGenericDataObjectReader()
       reader.SetFileName(filename)
       reader.Update()

       p2c = vtkPointDataToCellData()
       p2c.SetInputConnection(reader.GetOutputPort())
       # p2c.SetInput(reader.GetOutput())
       p2c.PassPointDataOn()
       p2c.Update()

       ug = p2c.GetOutput()

       print 'UG n.points:',ug.GetNumberOfPoints(), ' n.cells:',ug.GetNumberOfCells()

    else:
      print 'unknown input file extension:',extension
      sys.exit(1)

    print '...done reading'

    return ug

####################################################
def read_mesh_transform(filename):

    d = open(filename)
    line = d.readline()
    line = line.split(':')
    line = line[1].split()
    scale = float(line[0])

    line = d.readline()
    line = line.split(':')
    line = line[1].split(',')
    translate = [float(line[0]),float(line[1]),float(line[2])]
    d.close()

    return scale,translate

####################################################
def write_volume_to_surface_ugrid(ugn, filename):

    # general-purpose filter to extract geometry (and associated data) from any type of dataset
    # Geometry is obtained as follows: all 0D, 1D, and 2D cells are extracted. All 2D faces that are used by only one 3D cell (i.e., boundary faces) are extracted
    geomFilt = vtkGeometryFilter()
    if VTK_MAJOR_VERSION <= 5:
        geomFilt.SetInput(ugn)
    else:
        geomFilt.SetInputData(ugn)
    geomFilt.Update()

    # Generates triangles from input polygons and triangle strips. 
    # It also generates line segments from polylines unless PassLines is off,
    # and generates individual vertex cells from vtkVertex point lists unless PassVerts is off
    triFilt = vtkTriangleFilter()
    if VTK_MAJOR_VERSION <= 5:
        triFilt.SetInput(geomFilt.GetOutput())
    else:
        triFilt.SetInputData(geomFilt.GetOutput())
    triFilt.Update()

    # reduce the number of triangles in the triangle mesh
    quadClust = vtkQuadricClustering()
    if VTK_MAJOR_VERSION <= 5:
        quadClust.SetInput(triFilt.GetOutput())
    else:
        quadClust.SetInputData(triFilt.GetOutput())
    quadClust.SetNumberOfDivisions(350,350,350)
    quadClust.CopyCellDataOn()
    quadClust.Update()

    # transform point coordinates
    """
    pdtrans = vtkTransformPolyDataFilter()
    if VTK_MAJOR_VERSION <= 5:
        pdtrans.SetInput(quadClust.GetOutput())
    else:
        pdtrans.SetInputData(quadClust.GetOutput())
    pdtrans.SetTransform(transform)
    pdtrans.Update()
    """

    # adjust coordinates using a windowed sinc function interpolation kernel to relax the mesh
    # making the cells better shaped and the vertices more evenly distributed 
    smooth = vtkWindowedSincPolyDataFilter()
    if VTK_MAJOR_VERSION <= 5:
        # smooth.SetInput(pdtrans.GetOutput())
        smooth.SetInput(quadClust.GetOutput())
    else:
        # smooth.SetInputData(pdtrans.GetOutput())
        smooth.SetInputData(quadClust.GetOutput())
    smooth.SetNumberOfIterations(30)
    smooth.SetPassBand(.15)
    smooth.Update()

    # merge duplicate points, and/or remove unused points and/or remove degenerate cells
    cleaned = vtkCleanPolyData()
    if VTK_MAJOR_VERSION <= 5:
        cleaned.SetInput( smooth.GetOutput() )
    else:
        cleaned.SetInputData( smooth.GetOutput() )
    cleaned.Update()

    # generates triangles from input polygons and triangle strips. It also generates line segments 
    # from polylines unless PassLines is off, and generates individual vertex cells from vtkVertex point lists unless PassVerts is off
    tri = vtkTriangleFilter()
    if VTK_MAJOR_VERSION <= 5:
        tri.SetInput(cleaned.GetOutput())
    else:
        tri.SetInputData(cleaned.GetOutput())
    tri.PassVertsOff()
    tri.Update()
    a = tri.GetOutput()

    ctop = vtkCellDataToPointData()
    if VTK_MAJOR_VERSION <= 5:
        ctop.SetInput(a)
    else:
        ctop.SetInputData(a)
    ctop.Update()
    geom   = ctop.GetOutput()

    """
    nCells = geom.GetNumberOfCells()
    pnts   = geom.GetPoints()
    scals  = geom.GetPointData().GetArray('pressure_mmHg')

    print 'scals:',scals
    pntCrds = []
    pntScals   = []
    for i in xrange(nCells):
        a = vtkIdList()
        geom.GetCellPoints(i,a)
        nIds = a.GetNumberOfIds()
        for idd in xrange(nIds):
                id = a.GetId(idd)
                pntCrds.append(pnts.GetPoint(id)[0])
                pntCrds.append(pnts.GetPoint(id)[1])
                pntCrds.append(pnts.GetPoint(id)[2])
                pntScals.append(scals.GetTuple(id)[0])
    """

    if filename:

        writer = vtkPolyDataWriter()
        writer.SetFileName(filename)

        if VTK_MAJOR_VERSION <= 5:
            writer.SetInput(geom)
        else:
            writer.SetInputData(geom)
        writer.Write()

    print '...done!'
#######################################################
#######################################################

if __name__ == '__main__':

    print 'Postprocessing Coronary ... [ VTK v.',vtk.VTK_MAJOR_VERSION,']'

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--coronary_config', required=True,  help='coronary configuration (skeleton.conf)')
    parser.add_argument('-i', '--input_unstruct_mesh', required=True,  help='mesh file (.vtu)')
    parser.add_argument('-o', '--output_surface',      required=True,  help='surface file (.vtp)')
    parser.add_argument('-u', '--output_ugrid',        required=False, help='mesh file (.vtu)')
    args = parser.parse_args()

    CFILE = args.coronary_config
    IFILE = args.input_unstruct_mesh
    OUGRIDFILE = args.output_ugrid
    OSURFFILE = args.output_surface

    class skel:
        pass

    skel.CFILES = {} # centerlines
    skel.XFILES = {} # crosssections

    if not os.path.isfile(CFILE):
        print >> sys.stderr, "Cannot find mesh file", IFILE
        sys.exit(1)

    if not os.path.isfile(IFILE):
        print >> sys.stderr, "Cannot find mesh file", IFILE
        sys.exit(1)

    # configuration to get the artery clines/xsections
    cfg = ConfigParser.SafeConfigParser()

    try:
        cfg.read(CFILE)
    except:
        print 'Config file:', CFILE,' missing !!'
        sys.exit(1)

    try:
        skel.branches = []

        units = cfg.get('general', 'skeletonunits')

        FTRANSF = cfg.get('general', 'mesh_transform')

        if not os.path.isfile(FTRANSF):
            print >> sys.stderr, "Cannot find file",FTRANSF
            sys.exit(1)

        side = cfg.get('general', 'side')

        if side != 'LCA' and side != 'RCA':
            print >> sys.stderr, 'coronary can only be LCA or RCA...got: ',side
            sys.exit(1)

        inletarea = float(cfg.get('general', 'inletarea'))

        if inletarea < 0:
            print 'Warning: Negative inlet area, resetting to 1.0'
            inletarea = 1.0

        for pname in cfg.sections():
            if pname == 'general': continue
            bb = {}
            bb['folder'] = cfg.get(pname, 'folder')
            outletarea = float(cfg.get(pname, 'outletarea'))
            if outletarea < 0:
                print 'Warning: Negative outlet area, resetting to 1.0'
                outletarea = 1.0
            bb['outletarea'] = outletarea
            skel.branches.append(bb)

    except:
        print 'Config file:', CFILE,' internal error !!'
        sys.exit(1)

    # for RCA, the order of vessels is : RCA, others
    if side=='RCA':
        br2 = []
        while len(br2)<len(skel.branches):
            for bb in skel.branches:
                if bb['folder'].find('RCA')>=0: 
                    br2.append(bb)
                    break
            for bb in skel.branches:
                if bb['folder'].find('LCX')<0 and bb['folder'].find('LCX')<0: 
                    br2.append(bb)
    # for LCA, the order of vessels is : LAD, LCX, others
    elif side=='LCA':
        br2 = []
        while len(br2)<len(skel.branches):
            for bb in skel.branches:
                if bb['folder'].find('LAD')>=0: 
                    br2.append(bb)
                    break
            for bb in skel.branches:
                if bb['folder'].find('LCX')>=0: 
                    br2.append(bb)
                    break
            for bb in skel.branches:
                if bb['folder'].find('LCX')<0 and bb['folder'].find('LCX')<0: 
                    br2.append(bb)
    skel.branches = br2

    CLINEDIRS = []
    for b in skel.branches:
        CLINEDIRS.append(b['folder'])

    for cf in CLINEDIRS:

        # if not os.path.isfile(cf):
        if not os.path.isdir(cf):
            # print >> sys.stderr, "Cannot find cline file", cf
            print >> sys.stderr, "Cannot find cline dir", cf
        else:
            fl = cf.split('CLINE_')
            cfl = os.path.join(cf, fl[1] + '.vtk')

            if not os.path.isfile(cfl):
                print >> sys.stderr, "Cannot find cline file", cfl
                sys.exit(1)

            skel.CFILES[cf] = cfl

            xf = []
            for f in os.listdir(cf):
                dummy, extension = os.path.splitext(f)
                # print dummy
                cfx = os.path.join(cf,f)
                if f[0:9]=='crossSect' and os.path.isfile(cfx) and extension=='.vtk':
                    xf.append(cfx)
            skel.XFILES[cf] = xf
            # print xf

    um = UnitsMapping()

    um.setResolution(VOXCM) #vox/cm
    um.setPhysicalMassDensity(MASSDENSITY) #Kg/m3
    um.setPhysicalViscosity(VISCOSITY) #m2/s
    um.setCharacteristicVelocity(VCH) #m/s
    um.setCharacteristicPressure_mmHg(PCH) # mmHg

    um.setCharacteristicMach(MACH) # best practice for resolution 300

    um.addInlet(id=0, velocity=0., area=inletarea)#

    ic = 1
    for bb in skel.branches:
        um.addOutlet(id=ic, name=bb['folder'], area=bb['outletarea'])#
        ic += 1

    um.setOutletMethodOutflow()
    
    um.workout()

    # read old unstruct grid with data in lattice units
    ug = read_ugrid(IFILE)

    # make new unstruct grid with data converted to physical units
    ugn = vtkUnstructuredGrid()

    points = vtkPoints()

    # cells = vtkCells()
    # cells = vtkCellArray()
    # vert = vtkVertex()
    # voxel = vtkVoxel()

    for pid in xrange(ug.GetNumberOfPoints()):

        coo = ug.GetPoint(pid)
        i,j,k = int(coo[0]),int(coo[1]),int(coo[2])
        points.InsertNextPoint(i,j,k)

        """
        ug.InsertNextCell(voxel.GetCellType(), voxel.GetPointIds())
        vert.GetPointIds().SetId(0,n)
        cells.InsertNextCell(vert)
        """

    ugn.SetPoints(points)

    cellPointIds = vtkIdList()
    cell = vtkVoxel()
    # cell = vtkVertex()

    for cid in xrange(ug.GetNumberOfCells()):

        ug.GetCellPoints(cid, cellPointIds)

        for ip in xrange( cell.GetNumberOfPoints() ):

            pid = cellPointIds.GetId(ip)
            cell.GetPointIds().SetId(ip, pid)

        ugn.InsertNextCell(cell.GetCellType(), cell.GetPointIds())

    for i in range( ug.GetCellData().GetNumberOfArrays() ):
        print i,'______Arrays:',ug.GetCellData().GetArrayName(i) # ,ug.GetCellData().GetArray(i).GetName()

    arr0 = convert_cpmapping(ug, um, 'pressure_mmHg')
    arr1 = convert(ug, 'density',  'density_MKS',  1000.0/RHO, 1)
    arr2 = convert(ug, 'velocity', 'velocity_MKS', um.Dx/um.Dt, 3)

    ugn.GetCellData().AddArray(arr0)
    ugn.GetCellData().AddArray(arr1)
    ugn.GetCellData().AddArray(arr2)

    if OUGRIDFILE:

        print 'writing new unstruct grid file for debugging ...'
        writer = vtkUnstructuredGridWriter()
        writer.SetFileName(OUGRIDFILE)
        if VTK_MAJOR_VERSION <= 5:
            writer.SetInput(ugn)
        else:
            writer.SetInputData(ugn)
        writer.Write()

    # read the mesh_transform.inp
    scale,translate = read_mesh_transform(FTRANSF)

    transform = vtkTransform()
    transform.Translate(translate)
    transform.Scale([scale,scale,scale])

    class o:
        pass

    slicetime = 0.
    # loop over centerlines and 
    for cdr in CLINEDIRS:

        cf = skel.CFILES[cdr]
        slicefile = 'slices_'+os.path.split( cf.split('CLINE_')[1] )[1]

        print '\nProbing wdir:', os.path.split(cf)[0], \
              ' CLineFile:',os.path.split(cf)[1], \
              ' SliceFile:', slicefile

        clinereader = vtkPolyDataReader()
        clinereader.SetFileName(cf)
        clinereader.Update()

        # Cline = clinereader.GetOutput()

        Cline = vtkTransformPolyDataFilter()
        Cline.SetInputConnection(clinereader.GetOutputPort())
        Cline.SetTransform(transform)
        Cline.Update()

        probe = vtkProbeFilter()
        if VTK_MAJOR_VERSION <= 5:
            probe.SetInput(Cline.GetOutput())
            probe.SetSource(ugn)
        else:
            probe.SetInputData(Cline.GetOutput())
            probe.SetSourceData(ugn)
        probe.Update()

        froot = cf.split('.vtk')[0]

        fp = open(froot+'_pressure_mmHg.dat','w')
        fd = open(froot+'_density_MKS.dat','w')
        fv = open(froot+'_velocity_MKS.dat','w')

        # print 'Valid Points:',probe.GetOutput().GetNumberOfPoints(),probe.GetValidPoints().GetNumberOfTuples()

        valids = probe.GetValidPoints()
        valid_ids = []
        for pid in xrange(valids.GetNumberOfTuples()):
            valid_ids.append( int(valids.GetTuple1(pid)) )

        # for pid in xrange(probe.GetOutput().GetNumberOfPoints()):
        for pid in valid_ids:

            print >> fp, pid,probe.GetOutput().GetPointData().GetArray('pressure_mmHg').GetValue(pid)
            print >> fd, pid,probe.GetOutput().GetPointData().GetArray('density_MKS').GetValue(pid)

            vx = probe.GetOutput().GetPointData().GetArray('velocity_MKS').GetComponent(pid,0)
            vy = probe.GetOutput().GetPointData().GetArray('velocity_MKS').GetComponent(pid,1)
            vz = probe.GetOutput().GetPointData().GetArray('velocity_MKS').GetComponent(pid,2)

            print >> fv, pid,math.sqrt(vx**2 + vy**2 + vz**2),vx,vy,vz

        fp.close()
        fd.close()
        fv.close()

        polyslice = vtkAppendPolyData()

        agrobs = [o(), o(), o(), o(), o()]

        agrobs[0].name = 'Area [mm\S2\N]'
        agrobs[1].name = 'Density [Kg/m\S3\N]'
        agrobs[2].name = 'Pressure [mmHg]'
        agrobs[3].name = 'Flowrate [ml/min]'
        agrobs[4].name = '<Vel> [m/s]'
        agrobs[0].dct = {}
        agrobs[1].dct = {}
        agrobs[2].dct = {}
        agrobs[3].dct = {}
        agrobs[4].dct = {}

        start = time.time()
        # loop over slices and compute flowrates
        kk = -1
        for xf in skel.XFILES[cdr]:

            kk += 1

            if kk % SLICE_ANALYSIS_FREQ != 0: continue

            areamm2,dens,press,flowratemlmin,avevel = SliceIt(xf,ugn,transform,kk,polyslice)

            if areamm2 != None:
                agrobs[0].dct[kk] = areamm2
                agrobs[1].dct[kk] = dens
                agrobs[2].dct[kk] = press
                agrobs[3].dct[kk] = flowratemlmin
                agrobs[4].dct[kk] = avevel

        stop = time.time()
        slicetime += stop - start
        mygrace.writevalues(fname = froot+'_xsectionals.dat', 
                            title = cf, 
                            xname = 'Cross-Section Id', 
                            yname = 'Cross-Section Values', 
                            obs = agrobs)

        if WRITE_CLINE_SLICE_FILES_WITH_OBS:

            # Remove any duplicate points from polyslice
            cleanFilter = vtkCleanPolyData()
            cleanFilter.SetInputConnection(polyslice.GetOutputPort())
            cleanFilter.Update()

            # write the slice_ files
            writer = vtkPolyDataWriter()
            writer.SetFileName( os.path.join(cdr, slicefile ) )

            if vtk.VTK_MAJOR_VERSION <= 5:
                writer.SetInput(cleanFilter.GetOutput())
            else:
                writer.SetInputData(cleanFilter.GetOutput())
            writer.Write()

    # from the volumetric ugrid create the surface version and write on file
    write_volume_to_surface_ugrid(ugn,OSURFFILE)

    print 'Slicing Time:',slicetime

