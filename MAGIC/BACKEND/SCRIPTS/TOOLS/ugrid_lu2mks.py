#!/usr/bin/env python

import sys,os,argparse,re
import numpy as np
from myvtk import *
from math import *
import prepare
# import quadrature

from muphy2wrapper import *

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

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

        if cid<3:
            print cid,ARROLD,':',qty

    return arr

def convert_cpmapping(ug, um, name):

    arr = vtkFloatArray()
    arr.SetNumberOfComponents(1)
    arr.SetName(name)

    PaTommHg = 0.007500616827042

    PCH = um.getCharacteristicPressure_mmHg()

    RHOIN = prepare.RHO
    UIN = prepare.um.getInletLatticeVelocity(name='Inlet')
    MACH = prepare.um.getCharacteristicMach()

    fact = PaTommHg * (um.PhysicalMassDensity * um.Vch**2) / (RHOIN * MACH**2)

    # Pressure coefficient cp mapping 
    for cid in xrange(ug.GetNumberOfCells()):

        rho = ug.GetCellData().GetArray('density').GetValue(cid)

        p = PCH + fact * (rho - RHOIN)

        arr.InsertNextValue( p )

    return arr


def dist2(v1,v2):
    return (v1[0] - v2[0])**2 + (v1[1] - v2[1])**2 + (v1[2] - v2[2])**2

if __name__ == '__main__':

    print 'Vtk Version:',vtk.VTK_MAJOR_VERSION

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_unstruct_mesh', required=True,  help='mesh file (.vtu)')
    parser.add_argument('-m', '--mesh_transform',      required=True,  help='mesh_transform (.inp)')
    # parser.add_argument('-c', '--input_clines',        required=False, help='cline files (.vtk)')
    parser.add_argument('-d', '--input_clinedirs',        required=False, help='cline dirs (dir1,dir2,...)')
    parser.add_argument('-o', '--output_surface',      required=True,  help='surface file (.vtp)')
    parser.add_argument('-u', '--output_ugrid',        required=False, help='mesh file (.vtu)')
    args = parser.parse_args()

    IFILE = args.input_unstruct_mesh
    FTRANSF = args.mesh_transform
    CFILES = {} # centerlines
    XFILES = {} # crosssections
    if args.input_clinedirs:

    # if args.input_clines:
        CLINEDIRS = args.input_clinedirs.split(',')
    else:
        CLINEDIRS = []
    OUGRIDFILE = args.output_ugrid
    OSURFFILE = args.output_surface
 
    if not os.path.isfile(IFILE):
        print >> sys.stderr, "Cannot find mesh file", IFILE
        sys.exit(1)

    if not os.path.isfile(FTRANSF):
        print >> sys.stderr, "Cannot find file",FTRANSF
        sys.exit(1)

    for cf in CLINEDIRS:

        # if not os.path.isfile(cf):
        if not os.path.isdir(cf):
            # print >> sys.stderr, "Cannot find cline file", cf
            print >> sys.stderr, "Cannot find cline dir", cf
            sys.exit(1)
        else:
            fl = cf.split('CLINE_')
            cfl = os.path.join(cf, fl[1] + '.vtk')

            if not os.path.isfile(cfl):
                print >> sys.stderr, "Cannot find cline file", cfl
                sys.exit(1)

            CFILES[cf] = cfl

            ordered_files = natural_sort( os.listdir(cf) )

            xf = []
            for f in ordered_files:
                dummy, extension = os.path.splitext(f)
                cfx = os.path.join(cf,f)
                if f[0:9]=='crossSect' and os.path.isfile(cfx) and extension=='.vtk':
                    print 'XSECT:',cfx
                    xf.append(cfx)
            XFILES[cf] = xf

            if XFILES[cf]==[]:
                print 'WARNING .... no crossSect* files....crosssection analysis (flowrate) will not be done'


    prepare.units()
    print 'Char MACH is :', prepare.um.getCharacteristicMach()

    print 'reading...'

    filename, extension = os.path.splitext(IFILE)

    if extension == '.vtk':
       reader = vtkUnstructuredGridReader()
       reader.SetFileName(IFILE)
       reader.Update()
       ug = reader.GetOutput()

    elif extension == '.pvtu':
       reader = vtkXMLGenericDataObjectReader()
       reader.SetFileName(IFILE)
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

        # dens = ug.GetPointData().GetArray('density').GetValue(pid)
        # press = 1./3. * dens
        # arr.InsertNextValue(press)

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

    arr0 = convert_cpmapping(ug, prepare.um, 'pressure_mmHg')
    arr1 = convert(ug, 'density',  'density_MKS',  1000.0/prepare.RHO, 1)
    arr2 = convert(ug, 'velocity', 'velocity_MKS', prepare.um.Dx/prepare.um.Dt, 3)

    ugn.GetCellData().AddArray(arr0)
    ugn.GetCellData().AddArray(arr1)
    ugn.GetCellData().AddArray(arr2)

    # print ugn

    print 'writing...'

    if OUGRIDFILE:
        writer = vtkUnstructuredGridWriter()
        writer.SetFileName(OUGRIDFILE)
        if VTK_MAJOR_VERSION <= 5:
            writer.SetInput(ugn)
        else:
            writer.SetInputData(ugn)
        writer.Write()

    d = open(FTRANSF)
    line = d.readline()
    line = line.split(':')
    line = line[1].split()
    scale = float(line[0])

    # radiusSqr = (scale * 20.)**2
    radiusSqr = (scale * 10.)**2

    line = d.readline()
    line = line.split(':')
    line = line[1].split(',')
    translate = [float(line[0]),float(line[1]),float(line[2])]
    d.close()

    transform = vtkTransform()
    transform.Translate(translate)
    transform.Scale([scale,scale,scale])

    for cdr in CLINEDIRS:

        cf = CFILES[cdr]

        print 'probing cline...',cf

        clinereader = vtkPolyDataReader()
        clinereader.SetFileName(cf)
        clinereader.Update()

        # cline = clinereader.GetOutput()

        cline = vtkTransformPolyDataFilter()
        cline.SetInputConnection(clinereader.GetOutputPort())
        cline.SetTransform(transform)
        cline.Update()

        probe = vtkProbeFilter()
        if VTK_MAJOR_VERSION <= 5:
            probe.SetInput(cline.GetOutput())
            probe.SetSource(ugn)
        else:
            probe.SetInputData(cline.GetOutput())
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

            print >> fv, pid,sqrt(vx**2 + vy**2 + vz**2),vx,vy,vz

        fp.close()
        fd.close()
        fv.close()

        # fa = open( os.path.join(cdr, 'xsectionals.dat'),'wb')
        fa = open( froot+'_xsectionals.dat','wb')

        print >> fa,'@ title "',cf,'"'
        print >> fa,'@ yaxis  label "Cross-Sectional Values"'
        print >> fa,'@ xaxis  label "Cross-Section Id"'
        print >> fa,'@    xaxis  label char size 1.200000'
        print >> fa,'@    yaxis  label char size 1.200000'

        print >> fa,'@    s0 legend  "Area [mm\S2\N]"'
        print >> fa,'@    s1 legend  "Density [Kg/m\S3\N]"'
        print >> fa,'@    s2 legend  "Pressure [mmHg]"'
        # print >> fa,'@    s3 legend  "Flowrate [m\S3\N/s]"'
        print >> fa,'@    s3 legend  "Flowrate [ml/min]"'
        print >> fa,'@    s4 legend  "<Vel> [m/s]"'

        kk = -1
        for xf in XFILES[cdr]:

            kk += 1

            print 'probing xf:',xf
            xreader = vtkPolyDataReader()
            xreader.SetFileName(xf)
            xreader.Update()

            xsect = vtkTransformPolyDataFilter()
            xsect.SetInputConnection(xreader.GetOutputPort())
            xsect.SetTransform(transform)
            xsect.Update()

            # first point of XSections is reserved for centerline point
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
 
            # implicit function (not a source) for clipping
            # clipvol = vtkBox()
            # clipvol.SetBounds(coo[0]-15,coo[0]+15,coo[1]-15,coo[1]+15,coo[2]-15,coo[2]+15)

            radius = 20
            # clipvol = vtkCylinder()
            clipvol = vtkSphere()
            # clipvol.SetCenter([0.,0.,0.])
            clipvol.SetCenter(coo)
            clipvol.SetRadius(radius) # clip the data inside a clipvol. Should be large enough in lattice units
            # clipvol.SetRadius(20) # clip the data inside a clipvol. Should be large enough in lattice units
            # clipvol.SetRadius(30) # clip the data inside a clipvol. Should be large enough in lattice units
            # clipvol.SetRadius(100) # clip the data inside a clipvol. Should be large enough in lattice units
            transform2 = vtkTransform()
            transform2.Identity()

            def AlmostEquals(a,b):
                return abs(a-b)<1.e-6

            def L2Norm(a):
                return sqrt(a[0]**2 + a[1]**2 + a[2]**2)

            if not AlmostEquals(norm[0], 0.0) or not AlmostEquals(norm[1], 1.0) or not AlmostEquals(norm[2], 0.0):

                rotationAxis = [-norm[2], 0.0, -norm[0]]
                rotationAxisMagnitude = L2Norm(rotationAxis)
                rotationAxis = [val / rotationAxisMagnitude for val in rotationAxis]
                angle = (180./3.1415) * (acos(norm[1] / L2Norm(norm))) # Find the rotation angle
                transform2.RotateWXYZ(angle, rotationAxis[0], rotationAxis[1], rotationAxis[2])

            # transform2.Translate(coo[0], coo[1], coo[2])

            # clipvol.SetTransform(transform2)

            """
            clipper = vtkClipDataSet()
            if VTK_MAJOR_VERSION <= 5:
                clipper.SetInput(ugn)
            else:
                clipper.SetInputData(ugn)
            clipper.InsideOutOn()
            clipper.SetClipFunction(clipvol)
            # clipper.SetClipFunction(plane)
            clipper.Update()
            """

            cutter = vtkCutter()
            cutter.SetCutFunction(plane)
            """
            cutter.SetInputConnection(clipper.GetOutputPort())
            """
            cutter.SetInputData(ugn)
            # cutter.SetInputData(ugn)
            cutter.Update()

            cellPointIds = vtkIdList()

            # breakdown to avoid unconnected pieces of the cutter
            conn = vtkPolyDataConnectivityFilter()

            if VTK_MAJOR_VERSION<=5:
                conn.SetInput(cutter.GetOutput())
            else:
                conn.SetInputData(cutter.GetOutput())

            conn.SetExtractionModeToClosestPointRegion()
            conn.SetClosestPoint(coo)
            conn.Update()

            Slice = conn.GetOutput()

            cf = vtkCenterOfMass()
            if VTK_MAJOR_VERSION <= 5:
                cf.SetInput(Slice)
            else:
                cf.SetInputData(Slice)
            cf.SetUseScalarsAsWeights(False)
            cf.Update()
 
            # print 'Coms:', cf.GetCenter(),'Ctr:',coo

            #arr_dens = cutter.GetOutput().GetCellData().GetArray('density_MKS')
            #arr_press = cutter.GetOutput().GetCellData().GetArray('pressure_mmHg')
            #arr_vel = cutter.GetOutput().GetCellData().GetArray('velocity_MKS')
            arr_dens = Slice.GetCellData().GetArray('density_MKS')
            arr_press = Slice.GetCellData().GetArray('pressure_mmHg')
            arr_vel = Slice.GetCellData().GetArray('velocity_MKS')

            area,dens,press,flowrate = 0,0,0,0
            # for i in xrange(cutter.GetOutput().GetNumberOfCells()):
            for i in xrange(Slice.GetNumberOfCells()):

                # cell = cutter.GetOutput().GetCell(i)
                # cutter.GetOutput().GetCellPoints(i, cellPointIds)
                cell = Slice.GetCell(i)
                Slice.GetCellPoints(i, cellPointIds)
                if cell.GetNumberOfPoints() != 3: 
                    print 'Not A Triangle'
                    sys.exit(1)

                # V0 = cutter.GetOutput().GetPoint( cellPointIds.GetId(0) )
                # V1 = cutter.GetOutput().GetPoint( cellPointIds.GetId(1) )
                # V2 = cutter.GetOutput().GetPoint( cellPointIds.GetId(2) )
                V0 = Slice.GetPoint( cellPointIds.GetId(0) )
                V1 = Slice.GetPoint( cellPointIds.GetId(1) )
                V2 = Slice.GetPoint( cellPointIds.GetId(2) )

                """
                AA = V1[0]-V0[0], V1[1]-V0[1], V1[2]-V0[2]
                BB = V2[0]-V0[0], V2[1]-V0[1], V2[2]-V0[2]

                dx = AA[1]*BB[2] - AA[2]*BB[1]
                dy = AA[2]*BB[0] - AA[0]*BB[2]
                dz = AA[0]*BB[1] - AA[1]*BB[0]

                dA = 0.5 * sqrt( dx**2 + dy**2 + dz**2 )
                dA *= (1.e-3/scale)**2 # transform to m2

                area += dA

                dens += arr_dens.GetValue(i) * dA

                press += arr_press.GetValue(i) * dA

                flowrate += (norm[0]*arr_vel.GetComponent(i,0) + \
                             norm[1]*arr_vel.GetComponent(i,1) + \
                             norm[2]*arr_vel.GetComponent(i,2)) * dA
                """

                # if dist2(coo, V0)<radiusSqr and dist2(coo, V1)<radiusSqr and dist2(coo, V2)<radiusSqr:
                if dist2(coo, V0)<radiusSqr or dist2(coo, V1)<radiusSqr or dist2(coo, V2)<radiusSqr:

                    AA = V1[0]-V0[0], V1[1]-V0[1], V1[2]-V0[2]
                    BB = V2[0]-V0[0], V2[1]-V0[1], V2[2]-V0[2]

                    dx = AA[1]*BB[2] - AA[2]*BB[1]
                    dy = AA[2]*BB[0] - AA[0]*BB[2]
                    dz = AA[0]*BB[1] - AA[1]*BB[0]

                    dA = 0.5 * math.sqrt( dx**2 + dy**2 + dz**2 )
                    dA *= (1.e-3/scale)**2 # transform to m2

                    area += dA

                    dens += arr_dens.GetValue(i) * dA

                    press += arr_press.GetValue(i) * dA

                    flowrate += (norm[0]*arr_vel.GetComponent(i,0) + \
                                 norm[1]*arr_vel.GetComponent(i,1) + \
                                 norm[2]*arr_vel.GetComponent(i,2)) * dA

                    #ringLine = str(0.001*(nrings-idx)*cl.stepLen)+' '+str(flowrate*-1.0E6*60.0)+'\n' #1.0E6*60.0 to go from m^3 / sec to mlLiters/min
                    #outFile.write(ringLine)

            avevel = 0
            if area>0:
                dens /= area
                press /= area
                # flowrate /= area
                avevel = flowrate / area

                # area in mm2, dens in Kg/m3, press in mmHg, flowrate in ml/min, avevel in m/s
                print >> fa, kk, area * 1.e+6, dens, press, 60.e+6 * flowrate, avevel

                fa.flush()

            """
            # if kk%12==0:
            if kk >= 85:

                sample = vtkSampleFunction()
                sample.SetImplicitFunction(clipvol);
                # sample.SetModelBounds(xmin, xmax, ymin, ymax, zmin, zmax)
                sample.SetModelBounds(coo[0]-radius, coo[0]+radius, coo[1]-radius, coo[1]+radius, coo[2]-radius, coo[2]+radius)
                sample.SetSampleDimensions(10,10,10)
                sample.ComputeNormalsOff()
                sample.Update()

                writer = vtkXMLImageDataWriter()
                # writer.SetFileName( 'clip'+str(kk)+'.vtk' )
                writer.SetFileName( 'clip'+str(kk)+'.vti' )
                if VTK_MAJOR_VERSION <= 5:
                    writer.SetInput(sample.GetOutput())
                else:
                    writer.SetInputData(sample.GetOutput())
                writer.Write()

                writer = vtkPolyDataWriter()
                writer.SetFileName('cut'+str(kk)+'.vtk')
                if VTK_MAJOR_VERSION <= 5:
                    # writer.SetInput(cutter.GetOutput())
                    writer.SetInput(Slice)
                else:
                    # writer.SetInputData(cutter.GetOutput())
                    writer.SetInputData(Slice)
                writer.Write()

            if kk==95: sys.exit(1)
            """

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

        fa.close()

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

    if OSURFFILE:
        writer = vtkPolyDataWriter()
        writer.SetFileName(OSURFFILE)
        if VTK_MAJOR_VERSION <= 5:
            writer.SetInput(geom)
        else:
            writer.SetInputData(geom)
        writer.Write()

    print '...done!'
