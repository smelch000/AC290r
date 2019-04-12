#!/usr/bin/env python

#
# script to fill a closed, water-tight surface with mesh points:
#   check if the points fall into a single inlet or *single* outlet box
#   use a scale factor to alter mesh resolution
#
# after this phase, wrap_by_wall.py should be run to wrap with wall nodes
#

from myvtk import *
import sys
import argparse
import numpy as np
import wrap_by_wall
import mesh
import d3q19

parser = argparse.ArgumentParser()

parser.add_argument('-s', '--surface',    required=True, help='input stl file')
parser.add_argument('-C', '--clip',   required=False, default=False, help='clip region stl file')
parser.add_argument('-A', '--aclip',   required=False, default=False, help='anticlip region stl file')
parser.add_argument('-i', '--inlet',      required=False, default=False, help='inlet stl file')
parser.add_argument('-o', '--outlet',     required=False, default=False, help='outlet stl file')
parser.add_argument('-O', '--output',     required=True, default='bgkflag', help='output')
parser.add_argument('-S', '--SCALE',      required=False, default=1.0, help='scale factor')
parser.add_argument('-G', '--GRIDSPACING', required=False, default=1, help='grid spacing')
parser.add_argument('-w', '--wrapbywall', required=False, action='store_true', help='wrap by wall at the end')
parser.add_argument('-V', '--visual', required=False, default=False, help='3d viz')
parser.add_argument('-Ib', '--Ib', required=False, default=False, help='begin of i')
parser.add_argument('-Jb', '--Jb', required=False, default=False, help='begin of j')
parser.add_argument('-Kb', '--Kb', required=False, default=False, help='begin of k')
parser.add_argument('-Ie', '--Ie', required=False, default=False, help='end of i')
parser.add_argument('-Je', '--Je', required=False, default=False, help='end of j')
parser.add_argument('-Ke', '--Ke', required=False, default=False, help='end of k')

args = parser.parse_args()

readsurf = vtkSTLReader(); readsurf.SetFileName( args.surface ); readsurf.Update()
filtsurf = readsurf

SCALE = float(args.SCALE)
print 
print 
print ' *** ',sys.argv[0],' ***'
print 
print sys.argv[1:]
print 'Scaling factor:',args.SCALE
print 'Mesh Grid :',args.GRIDSPACING
print 

if args.inlet:
    print 'Inlet from file:',args.inlet
    readin = vtkSTLReader(); readin.SetFileName( args.inlet ); readin.Update()
    filtin = readin

if args.outlet:
    print 'Outlet from file:',args.outlet
    readout = vtkSTLReader(); readout.SetFileName( args.outlet ); readout.Update()
    filtout = readout

if args.clip:
    print 'Clip region from file:',args.clip
    readclip = vtkSTLReader(); readclip.SetFileName( args.clip ); readclip.Update()
    filtclip = readclip

if args.aclip:
    print 'Aclip region from file:',args.aclip
    readaclip = vtkSTLReader(); readaclip.SetFileName( args.aclip ); readaclip.Update()
    filtaclip = readaclip

bnd = filtsurf.GetOutput().GetBounds()
o = (bnd[0], bnd[2], bnd[4])
print 'Surface origin: %8.3f %8.3f %8.3f '%o, 'Bounds: %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f '%bnd
 

def scaleGeom(filtgeom,SCALE):
    t = vtkTransform()
    t.Scale(SCALE, SCALE, SCALE)
    tF = vtkTransformPolyDataFilter()
    tF.SetInputConnection(filtgeom.GetOutputPort())
    tF.SetTransform(t)
    tF.Update()
    return tF

def translateGeom(filtgeom,Tx,Ty,Tz):
    t = vtkTransform()
    t.Translate(Tx, Ty, Tz)
    tF = vtkTransformPolyDataFilter()
    tF.SetInputConnection(filtgeom.GetOutputPort())
    tF.SetTransform(t)
    tF.Update()
    return tF

# now TRANSLATING
if bnd[0]<1 or bnd[2]<1 or bnd[4]<1: # translation to positive octant

    filtsurf = translateGeom(filtsurf, -bnd[0]+10, -bnd[2]+10, -bnd[4]+10)
    if args.inlet:  filtin = translateGeom(filtin, -bnd[0]+10, -bnd[2]+10, -bnd[4]+10)
    if args.outlet: filtout = translateGeom(filtout, -bnd[0]+10, -bnd[2]+10, -bnd[4]+10)
    if args.clip: filtclip = translateGeom(filtclip, -bnd[0]+10, -bnd[2]+10, -bnd[4]+10)
    if args.aclip: filtaclip = translateGeom(filtaclip, -bnd[0]+10, -bnd[2]+10, -bnd[4]+10)

bnd = filtsurf.GetOutput().GetBounds()
o = (bnd[0], bnd[2], bnd[4])
print 'Surface origin: %8.3f %8.3f %8.3f'%o, 'surface bounds: %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f '%bnd
 
# now SCALING
filtsurf = scaleGeom(filtsurf, SCALE)
if args.inlet:  filtin = scaleGeom(filtin, SCALE)
if args.outlet: filtout = scaleGeom(filtout, SCALE)
if args.clip: filtclip = scaleGeom(filtclip, SCALE)
if args.aclip: filtaclip = scaleGeom(filtaclip, SCALE)
bnd = filtsurf.GetOutput().GetBounds()
o = (bnd[0], bnd[2], bnd[4])
print 'Surface origin: %8.3f %8.3f %8.3f'%o, 'surface bounds: %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f '%bnd
 


surface = filtsurf.GetOutput()
if args.inlet:  inlet = filtin.GetOutput()
if args.outlet: outlet = filtout.GetOutput()
if args.clip: clip = filtclip.GetOutput()
if args.aclip: aclip = filtaclip.GetOutput()

#XR = int(bnd[0]*SCALE)-2, int(bnd[1]*SCALE)+2
#YR = int(bnd[2]*SCALE)-2, int(bnd[3]*SCALE)+2
#ZR = int(bnd[4]*SCALE)-2, int(bnd[5]*SCALE)+2
XR = int(bnd[0])-2, int(bnd[1])+2
YR = int(bnd[2])-2, int(bnd[3])+2
ZR = int(bnd[4])-2, int(bnd[5])+2
print '\nInitial discrete bounds:',XR,YR,ZR

points = vtkPoints()

if args.clip:
    clipbnd = clip.GetBounds()
    XR = int(max(XR[0], clipbnd[0]-2)), int(min(XR[1],clipbnd[1]+2))
    YR = int(max(YR[0], clipbnd[2]-2)), int(min(YR[1],clipbnd[3]+2))
    ZR = int(max(ZR[0], clipbnd[4]-2)), int(min(ZR[1],clipbnd[5]+2))
    print 'resetting to clip region bounds:', XR,YR,ZR

GRIDSPACING = int(args.GRIDSPACING)

ibegin = (XR[0] / GRIDSPACING) * GRIDSPACING
jbegin = (YR[0] / GRIDSPACING) * GRIDSPACING
kbegin = (ZR[0] / GRIDSPACING) * GRIDSPACING

iend = XR[1]
jend = YR[1]
kend = ZR[1]

if args.Kb: kbegin = int(args.Kb)
if args.Jb: jbegin = int(args.Jb)
if args.Ib: ibegin = int(args.Ib)
if args.Ke: kend = int(args.Ke)
if args.Je: jend = int(args.Je)
if args.Ie: iend = int(args.Ie)

print '\nDiscrete bounds:',ibegin,':',iend,'  ',jbegin,':',jend,'  ',kbegin,':',kend
print

print 'Initial points filling...'
ntry = 0
for k in xrange(kbegin, kend, GRIDSPACING):

    for j in xrange(jbegin, jend, GRIDSPACING):

        for i in xrange(ibegin, iend, GRIDSPACING):

            points.InsertNextPoint(i, j, k)
            ntry += 1
print '... done with size',ntry
 
pointsPolydata = vtkPolyData()
pointsPolydata.SetPoints(points);
 
# Points inside test
select_surf = vtkSelectEnclosedPoints()

print 'default tolerance:', select_surf.GetTolerance(),
select_surf.SetTolerance(0.)
print '    new tolerance:', select_surf.GetTolerance()

select_surf.SetInputData(pointsPolydata);
select_surf.SetSurfaceData(surface);
select_surf.Update();

if args.inlet:  
    select_inlet = vtkSelectEnclosedPoints()
    select_inlet.SetTolerance(0.)
    select_inlet.SetInputData(pointsPolydata);
    select_inlet.SetSurfaceData(inlet);
    select_inlet.Update();
if args.outlet: 
    select_outlet = vtkSelectEnclosedPoints()
    select_outlet.SetTolerance(0.)
    select_outlet.SetSurfaceData(outlet);
    select_outlet.SetInputData(pointsPolydata);
    select_outlet.Update();
if args.clip: 
    select_clip = vtkSelectEnclosedPoints()
    select_clip.SetTolerance(0.)
    select_clip.SetInputData(pointsPolydata);
    select_clip.SetSurfaceData(clip);
    select_clip.Update();
if args.aclip: 
    select_aclip = vtkSelectEnclosedPoints()
    select_aclip.SetTolerance(0.)
    select_aclip.SetInputData(pointsPolydata);
    select_aclip.SetSurfaceData(aclip);
    select_aclip.Update();
 
FLUID = 1
WALL = 2
INLET = 3
OUTLET = 4

xmn = 1e5; xmx =-1e5
ymn = 1e5; ymx =-1e5
zmn = 1e5; zmx =-1e5
#itp_f,itp_w,itp_i,itp_o = [],[],[],[]

nfluid,nwall,ninlet,noutlet = 0,0,0,0

for n in xrange( pointsPolydata.GetNumberOfPoints() ):
    if select_surf.IsInside(n):

        p = pointsPolydata.GetPoint(n)
        P = int(p[0]), int(p[1]), int(p[2])

        xmn = min(xmn,P[0]); xmx = max(xmx,P[0])
        ymn = min(ymn,P[1]); ymx = max(ymx,P[1])
        zmn = min(zmn,P[2]); zmx = max(zmx,P[2])

        if args.clip and not select_clip.IsInside(n): continue
        if args.aclip and    select_aclip.IsInside(n): continue

        if args.inlet and select_inlet.IsInside(n):
            ninlet += 1
        elif args.outlet and select_outlet.IsInside(n):
            noutlet += 1
        else:
            nfluid += 1

nx = xmx + 2*GRIDSPACING
ny = ymx + 2*GRIDSPACING
nz = zmx + 2*GRIDSPACING

nxy2 = long(nx+2)*long(ny+2)
nx2 = long(nx+2)
msh = mesh.Mesh(nx,ny,nz)

msh.itppp_f = np.zeros(nfluid,  dtype='i8')

msh.nfluid,msh.nwall,msh.ninlet,msh.noutlet = 0,0,0,0

for n in xrange( pointsPolydata.GetNumberOfPoints() ):
    if select_surf.IsInside(n):

        p = pointsPolydata.GetPoint(n)

        i,j,k = int(p[0]), int(p[1]), int(p[2])
        i4 =  k*nxy2 + j*nx2 + i

        if args.clip and not select_clip.IsInside(n): continue
        if args.aclip and    select_aclip.IsInside(n): continue

        if args.inlet and select_inlet.IsInside(n):
            #itp_i.append(i4)
            msh.inlet.append([i,j,k])
            msh.itppp_i[msh.ninlet] = i4
            msh.ninlet += 1
        elif args.outlet and select_outlet.IsInside(n):
            #itp_o.append(i4)
            msh.outlet.append([i,j,k])
            msh.itppp_o[msh.noutlet] = i4
            msh.noutlet += 1
        else:
            #itp_f.append(i4)
            msh.fluid.append([i,j,k])
            msh.itppp_f[msh.nfluid] = i4
            msh.nfluid += 1

print 'Number of Fluid,Wall,Inlet,Outlet: ',msh.nfluid,msh.nwall,msh.ninlet,msh.noutlet

#msh.itppp_f = np.asarray( itp_f )
#msh.itppp_w = np.asarray( itp_w )
#msh.itppp_i = np.asarray( itp_i )
#msh.itppp_o = np.asarray( itp_o )

print

if args.wrapbywall:
    # wrap_by_wall.wrapit(infile='bgkflag',outfile='FWIO')

    print '\n  wall wrapping'

    lat = d3q19.D3q19(msh)

     # create the neighbors layer around all nodes
    newcomers = []

    itppp = np.copy(msh.fluid)
    for pnt in itppp:

        # find all neighbors of point (fluid, inlet or outlet)
        newneigh = lat.getneighbors(pnt, gridspacing=GRIDSPACING)

        for i4 in newneigh:
            newcomers.append(i4)

    print '  wall with duplicates:   ',len(newcomers)

    # remove duplicates
    newcomers = list(set(newcomers))

    print '  wall without duplicates:',len(newcomers)

    print '  wall initial points filling...'
    wpoints = vtkPoints()
    for i4 in newcomers:

        i,j,k = msh.ijkdir(1,i4), msh.ijkdir(2,i4), msh.ijkdir(3,i4)
        wpoints.InsertNextPoint(i, j, k)
 
    wpointsPolydata = vtkPolyData()
    wpointsPolydata.SetPoints(wpoints);
 
    # Points inside test
    select_surf = vtkSelectEnclosedPoints()
    select_surf.SetTolerance(0.)
    select_surf.SetInputData(wpointsPolydata);
    select_surf.SetSurfaceData(surface);
    select_surf.Update();

    # take wall nodes outside the surface
    newcomersinside = []
    for n in xrange( wpointsPolydata.GetNumberOfPoints() ):

        if select_surf.IsInside(n): continue # exclude internal nodes

        p = wpointsPolydata.GetPoint(n)
        i,j,k = int(p[0]), int(p[1]), int(p[2])
        newcomersinside.append( k*nxy2 + j*nx2 + i )

    msh.itppp_w = np.zeros(len(newcomersinside), dtype='i8')

    msh.nwall = 0
    for i4 in newcomersinside:

        msh.itppp_w[msh.nwall] = i4
        msh.nwall += 1
    print '  wall final number of outer nodes',msh.nwall
    print
 
    msh.itppp_w.sort()

print 'writing mesh on files',args.output + '*'

msh.writeMOEBIUSinput(args.output)




if not args.visual:   sys.exit(1)   ############################ STOP HERE JUST IN CASE




#insideArray = select_surf.GetOutput().GetPointData().GetArray("SelectedPoints")
#for i in xrange(insideArray.GetNumberOfTuples()) :
#    print i," selected point as array value : ", insideArray.GetComponent(i,0)
 
cubeMapper = vtkPolyDataMapper()
cubeMapper.SetInputConnection(readsurf.GetOutputPort());
 
cubeActor = vtkActor()
cubeActor.SetMapper(cubeMapper);
cubeActor.GetProperty().SetOpacity(0.5);
 
#Points mapper, actor
#First, apply vtkVertexGlyphFilter to make cells around points, vtk only render cells.

vertexGlyphFilter = vtkVertexGlyphFilter()
if VTK_MAJOR_VERSION <= 5:
  vertexGlyphFilter.AddInput(pointsPolydata);
else:
  vertexGlyphFilter.AddInputData(pointsPolydata);

vertexGlyphFilter.Update();
 
pointsMapper = vtkPolyDataMapper()
pointsMapper.SetInputConnection(vertexGlyphFilter.GetOutputPort());
 
pointsActor = vtkActor()
pointsActor.SetMapper(pointsMapper);
pointsActor.GetProperty().SetPointSize(2);
pointsActor.GetProperty().SetColor(0.0,0.0,1);
 
#Create a renderer, render window, and interactor
renderer = vtkRenderer()
renderWindow = vtkRenderWindow()
renderWindow.AddRenderer(renderer);
renderWindowInteractor = vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)
 
# Add the actor to the scene
renderer.AddActor(cubeActor);
renderer.AddActor(pointsActor)
renderer.SetBackground(.0, 1,.0);
 
# Render and interact
renderWindow.Render();
renderWindowInteractor.Start();
