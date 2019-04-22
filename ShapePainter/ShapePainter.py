#!/usr/bin/env python

# from myvtk import *
from vtk import *
import random
import math
import subprocess, shlex, os ,sys
import EXTRAS.genmeshx as GX
import EXTRAS.bgkflag2xyz as bgkflag2xyz

def SP(NX, NY, NZ, 
       RADIUS, 
       INBC, OUTBC, 
       INVAL, OUTVAL, 
       OUTSTL,
       VISUAL, 
       PERIODIC = None):
    """
    NX, NY, NZ: size of box in lattice units (l.u.)
    RADIUS: radius of cylinder in unstenotic region (l.u.)
    INBC, OUTBC: keywords for inlet/outlet BC ('flow' or 'pressure')
    INVAL, OUTVAL: values for boundary conditions (l.u.)
    OUTSTL: name of output file for triangular mesh
    VISUAL: flag to interactive visualization
    PERIODIC: periodicity along x,y,z
    """

    HIZ = NZ/10.

    OFF = 0

    # paint pixels on the canvas
    volume = vtkImageData()
    volume.SetExtent(OFF, OFF+NX-1, OFF, OFF+NY-1, OFF, OFF+NZ-1)
    volume.AllocateScalars(VTK_DOUBLE,1)

    dxh = OFF + NX/2.; dyh = OFF + NY/2.; dzh = OFF + NZ/2.

    Scalars = vtkFloatArray()

    for k in xrange(OFF, OFF+NZ):
        for j in xrange(OFF, OFF+NY):
            for i in xrange(OFF, OFF+NX):

                d = 0

                r = math.sqrt((i - dxh)**2 + (j - dyh)**2)
                zz = abs(k - NZ/2 - OFF)

                if zz > HIZ/2.:
                    if r < RADIUS : d = 1
                else:
                    if r < RADIUS/2.: d = 1

                if k == OFF or k == OFF + NZ-1: d = 0

                Scalars.InsertNextValue( d )

    volume.GetPointData().SetScalars(Scalars)

    # run the marching cubes
    Marching = vtkMarchingCubes()
    Marching.SetInputData(volume)
    Marching.SetValue(0, 0.5)
    Marching.ComputeNormalsOn()
    Marching.Update()

    # transform into a triangular mesh
    triangles = vtkTriangleFilter()
    triangles.SetInputConnection(Marching.GetOutputPort())

    # write the mesh onfile
    writer = vtkSTLWriter()
    writer.SetInputConnection(triangles.GetOutputPort())
    writer.SetFileName(OUTSTL)
    writer.SetFileTypeToASCII()
    writer.Write()

    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(triangles.GetOutputPort())
 
    actor = vtkActor()
    actor.SetMapper(mapper)

    i_globs, o_globs = [], []
    if PERIODIC == None:
        # define the inlet box to select the inlet region
        IN = vtkCubeSource()
        IN.SetCenter(dxh, dyh, OFF+1)
        IN.SetXLength(1.8*dxh); IN.SetYLength(1.8*dyh); IN.SetZLength(10)
        tIN = vtkTriangleFilter()
        tIN.SetInputConnection(IN.GetOutputPort())
        mIN = vtkPolyDataMapper()
        mIN.SetInputConnection(IN.GetOutputPort())
        aIN = vtkActor()
        aIN.SetMapper(mIN)
        aIN.GetProperty().SetColor(1.,0.,0.)
        aIN.GetProperty().SetOpacity(0.5)
        writer = vtkSTLWriter()
        writer.SetInputConnection(tIN.GetOutputPort())
        writer.SetFileName('IN.stl')
        writer.SetFileTypeToASCII()
        writer.Write()

        # define the outlet box to select the inlet region
        OUT = vtkCubeSource()
        OUT.SetCenter(dxh, dyh, OFF+NZ)
        OUT.SetXLength(1.8*dxh); OUT.SetYLength(1.8*dyh); OUT.SetZLength(10)
        tOUT = vtkTriangleFilter()
        tOUT.SetInputConnection(OUT.GetOutputPort())
        mOUT = vtkPolyDataMapper()
        mOUT.SetInputConnection(OUT.GetOutputPort())
        aOUT = vtkActor()
        aOUT.SetMapper(mOUT)
        aOUT.GetProperty().SetColor(0.,0.,1.)
        aOUT.GetProperty().SetOpacity(0.5)
        writer = vtkSTLWriter()
        writer.SetInputConnection(tOUT.GetOutputPort())
        writer.SetFileName('OUT.stl')
        writer.SetFileTypeToASCII()
        writer.Write()

        # build mesh nodes inside the generated surfaces
        mesh = GX.buildinside(triangles, None, None, [tIN], [tOUT], True, False, 1)

        i_globs.append( {'id':1, 'bctype':INBC, 'iodir':(0,0,-1), 'ioval':INVAL} )
        o_globs.append( {'id':2, 'bctype':OUTBC, 'iodir':(0,0,+1), 'ioval':OUTVAL} )

    else:

        # build mesh nodes inside the generated surfaces
        mesh = GX.buildinside(triangles, None, None, [], [], True, False, 1, PERIODIC)

    mesh.specifyNodes(mesh.itppp_f, mesh.itppp_w, mesh.itppp_i, mesh.itppp_o, 
                      i_id = mesh.itppp_i_id, o_id = mesh.itppp_o_id,
                      i_globs = i_globs, o_globs = o_globs)

    # write mesh on file. We provide the bounding box here.
    mesh.writeMOEBIUSinput('bgkflag')
 

    # write file xyz
    class E(): pass
    e = E()
    e.infiles = ['bgkflag']
    e.outfile = 'bgkflag'
    e.outwall = 'wall'

    bgkflag2xyz.main(e)

    if VISUAL:

        renderer = vtkRenderer()
        renderer.AddActor(actor)
        renderer.AddActor(aIN)
        renderer.AddActor(aOUT)
        renderer.SetBackground(0,0,0)
        renderer.ResetCamera()
 
        renderWindow = vtkRenderWindow()
        renderWindow.AddRenderer(renderer)
 
        iren = vtkRenderWindowInteractor()
        iren.SetRenderWindow(renderWindow)
        iren.Initialize()
        iren.Start()

if __name__ == '__main__' :

    RADIUS = 10
    LENGTH = 20 * RADIUS

    SP( NX = 2*RADIUS + 4, 
        NY = 2*RADIUS + 4, 
        NZ = LENGTH, 
        RADIUS = RADIUS,
        INBC = 'flow', OUTBC = 'pressure',
        INVAL = '0.01', OUTVAL = '0.0',
        OUTSTL = 'SP.stl',
        VISUAL = False)

