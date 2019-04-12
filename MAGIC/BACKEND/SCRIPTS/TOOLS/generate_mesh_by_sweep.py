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
import wrap_by_wall

parser = argparse.ArgumentParser()

parser.add_argument('-s', '--surface',    required=True, help='input stl file')
parser.add_argument('-b', '--boundregion',   required=False, default=False, help='bounding region stl file')
parser.add_argument('-i', '--inlet',      required=False, default=False, help='inlet stl file')
parser.add_argument('-o', '--outlet',     required=False, default=False, help='outlet stl file')
parser.add_argument('-S', '--SCALE',      required=False, default=1.0, help='scale factor')
parser.add_argument('-w', '--wrapbywall', required=False, default=False, help='wrap by wall at the end')
parser.add_argument('-V', '--visual', required=False, default=False, help='3d viz')

args = parser.parse_args()

readsurf = vtkSTLReader(); readsurf.SetFileName( args.surface ); readsurf.Update()

if args.inlet:
    print 'Reading inlet from file:',args.inlet
    readin = vtkSTLReader(); readin.SetFileName( args.inlet ); readin.Update()

if args.outlet:
    print 'Reading outlet from file:',args.outlet
    readout = vtkSTLReader(); readout.SetFileName( args.outlet ); readout.Update()

if args.boundregion:
    print 'Reading bounding region from file:',args.boundregion
    readbound = vtkSTLReader(); readbound.SetFileName( args.boundregion ); readbound.Update()

SCALE = float(args.SCALE)
print 
print 
print ' ***** Running : ',sys.argv[0],' ****** '
print 
print 'Scaling/Mesh resolution factor:',SCALE

# cubeSource = vtkCubeSource()
# # cubeSource.SetXLength(2.0);
# cubeSource.Update();
 
# surface = cubeSource.GetOutput();
surface = readsurf.GetOutput();
if args.inlet:  inlet = readin.GetOutput();
if args.outlet: outlet = readout.GetOutput();
if args.boundregion: bound = readbound.GetOutput();
 
points = vtkPoints()

bnd = surface.GetBounds()
# o = (0,0,0)
o = (bnd[0]-1, bnd[2]-1, bnd[4]-1)
XR = int(bnd[0]*SCALE)-2, int(bnd[1]*SCALE)+2
YR = int(bnd[2]*SCALE)-2, int(bnd[3]*SCALE)+2
ZR = int(bnd[4]*SCALE)-2, int(bnd[5]*SCALE)+2
print 'surface bounds:',XR,YR,ZR

if args.boundregion:
    boundbnd = bound.GetBounds()
    XR = int(max(XR[0], boundbnd[0]*SCALE-2)), int(min(XR[1],boundbnd[1]*SCALE+2))
    YR = int(max(YR[0], boundbnd[2]*SCALE-2)), int(min(YR[1],boundbnd[3]*SCALE+2))
    ZR = int(max(ZR[0], boundbnd[4]*SCALE-2)), int(min(ZR[1],boundbnd[5]*SCALE+2))
    #XR = int(max(XR[0], boundbnd[0]-2)), int(min(XR[1],boundbnd[1]+2))
    #YR = int(max(YR[0], boundbnd[2]-2)), int(min(YR[1],boundbnd[3]+2))
    #ZR = int(max(ZR[0], boundbnd[4]-2)), int(min(ZR[1],boundbnd[5]+2))
    print 'resetting to boundregion region bounds:', XR,YR,ZR

print 'Initial points filling...'
for k in xrange(ZR[0], ZR[1], 1):
    for j in xrange(YR[0], YR[1], 1):
        for i in xrange(XR[0], XR[1], 1):

            points.InsertNextPoint(i/SCALE, j/SCALE, k/SCALE)

print '... done'
 
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
if args.boundregion: 
    select_bound = vtkSelectEnclosedPoints()
    select_bound.SetTolerance(0.)
    select_bound.SetInputData(pointsPolydata);
    select_bound.SetSurfaceData(bound);
    select_bound.Update();
 
print 'writing points on files bgkflag* ...'

FLUID = 1
WALL = 2
INLET = 3
OUTLET = 4

b = open('bgkflag.dat','w')
h = open('bgkflag.hdr','w')
if args.inlet or args.outlet: 
    ios = open('bgkflag.ios','w')
# f = open('innerpoints.xyz','w')
nfl = 0
nwl = 0
nin = 0
nout = 0
xmn = 1e5; xmx =-1e5
ymn = 1e5; ymx =-1e5
zmn = 1e5; zmx =-1e5
ins = []
outs = []
for i in xrange( pointsPolydata.GetNumberOfPoints() ):
    if select_surf.IsInside(i):

        p = pointsPolydata.GetPoint(i)

        # P = int(SCALE*(p[0]-o[0])), int(SCALE*(p[1]-o[1])), int(SCALE*(p[2]-o[2]))
        P = int(p[0]-o[0]), int(p[1]-o[1]), int(p[2]-o[2])

        if args.boundregion and not select_bound.IsInside(i):
            continue

        if args.inlet and select_inlet.IsInside(i):
            print >> b, '%d %d %d %d' % ( P[0], P[1], P[2], INLET )
            #print >> f, 'H %.1f %.1f %.1f' % P
            ins.append( P )
            nin += 1
        elif args.outlet and select_outlet.IsInside(i):
            print >> b, '%d %d %d %d' % ( P[0], P[1], P[2], OUTLET )
            #print >> f, 'O %.1f %.1f %.1f' % P
            outs.append( P )
            nout += 1
        else:
            print >> b, '%d %d %d %d' % ( P[0], P[1], P[2], FLUID )
            #print >> f, 'C %.1f %.1f %.1f' % P
            nfl += 1
        xmn = min(xmn,P[0]); xmx = max(xmx,P[0])
        ymn = min(ymn,P[1]); ymx = max(ymx,P[1])
        zmn = min(zmn,P[2]); zmx = max(zmx,P[2])

print >> h, xmx+2, ymx+2, zmx+2
print >> h, '5 ',nfl,nwl,nin,nout

print 'Number of Fluid,Wall,Inlet,Outlet: ',nfl,nwl,nin,nout

b.close()
h.close()

if args.inlet or args.outlet: 
    print >> ios, 2 
    print >> ios, '1     inlet     pressure   0  -1   0     -0.000000     1.000000    -0.000000      0.000000 inlet_0_IN.stl'
    print >> ios, '2    outlet         flow   0   0  -1      0.000000     0.000000    -1.000000      0.000000 outlet_0_OUT.stl'

    print >> ios, len(ins)
    for i,j,k in ins:
        print >> ios, i,j,k,1

    print >> ios, len(outs)
    for i,j,k in outs:
        print >> ios, i,j,k,2

    ios.close()

if args.wrapbywall:
    wrap_by_wall.wrapit(infile='bgkflag',outfile='FWIO')


print 

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
