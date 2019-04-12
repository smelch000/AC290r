#!/usr/bin/env python

import sys,os,argparse
import numpy as np
from myvtk import *
from math import *
import ConfigParser
# import quadrature

from muphy2wrapper import *

RES = 50

#######################################################

def CenterOfPoly(poly):
    """ 
    return center of mass of polydata
    """
    centerOfMassFilter = vtkCenterOfMass()
    if VTK_MAJOR_VERSION <= 5:
        centerOfMassFilter.SetInput(poly)
    else:
        centerOfMassFilter.SetInputData(poly)
    centerOfMassFilter.SetUseScalarsAsWeights(False)
    centerOfMassFilter.Update()
    return centerOfMassFilter.GetCenter()

def CleanPolyData(poly):
    Clean = vtkCleanPolyData()
    Clean.SetInputData(poly)
    # Clean.SetTolerance(1.0) # this will slow down by far
    Clean.Update()
    return Clean
    

def TaubinFilter(poly):

    Smooth = vtkWindowedSincPolyDataFilter()
    Smooth.SetInputData(poly)
    Smooth.SetNumberOfIterations(200)
    Smooth.SetPassBand(0.1)
    Smooth.SetBoundarySmoothing(0)
    # Smooth.SetFeatureEdgeSmoothing(0)
    # Smooth.SetNonManifoldSmoothing(1)
    # Smooth.SetNormalizeCoordinates(1)
    Smooth.Update()
    return Smooth

if __name__ == '__main__':

    print 'Sock and Smooth ... VTK v.',vtk.VTK_MAJOR_VERSION

    parser = argparse.ArgumentParser()

    parser.add_argument('-s', '--stl', required=True,  help='lumen (.stl)')
    parser.add_argument('-c', '--cline', required=True,  help='cline (.vtk)')
    args = parser.parse_args()

    STLFILE = args.stl
    CLINEFILE = args.cline
 
    if not os.path.isfile(STLFILE):
        print >> sys.stderr, "Cannot find stl file", STLFILE
        sys.exit(1)

    if not os.path.isfile(CLINEFILE):
        print >> sys.stderr, "Cannot find cline file", CLINEFILE
        sys.exit(1)

    readerSTL = vtkSTLReader()
    readerSTL.SetFileName(STLFILE)
    readerSTL.Update()

    readerCLINE = vtkGenericDataObjectReader()
    readerCLINE.SetFileName(CLINEFILE)
    readerCLINE.Update()

    mapSTL = vtkDataSetMapper()
    if vtk.VTK_MAJOR_VERSION <= 5:
        mapSTL.SetInput(readerSTL.GetOutput())
    else:
        mapSTL.SetInputData(readerSTL.GetOutput())
    actSTL = vtkActor()
    actSTL.SetMapper(mapSTL)
    actSTL.GetProperty().SetOpacity(1.0)
    actSTL.GetProperty().SetColor( [1.,0.,0.] )
 
    mapCLINE = vtkDataSetMapper()
    if vtk.VTK_MAJOR_VERSION <= 5:
        mapCLINE.SetInput(readerCLINE.GetOutput())
    else:
        mapCLINE.SetInputData(readerCLINE.GetOutput())
    actCLINE = vtkActor()
    actCLINE.SetMapper(mapCLINE)

    """
    Wrap = vtkSphereSource()
    Wrap.SetThetaResolution( RES )
    Wrap.SetPhiResolution( RES )
    Wrap.SetCenter( CenterOfPoly(readerCLINE.GetOutput()) )
    # bnds = readerCLINE.GetOutput().GetBounds()
    bnds = readerSTL.GetOutput().GetBounds()
    Wrap.SetRadius( max(bnds[1]-bnds[0], bnds[3]-bnds[2], bnds[5]-bnds[4]) )
    """

    Wrap1 = vtkTubeFilter()
    Wrap1.SetInputConnection( readerCLINE.GetOutputPort())
    Wrap1.SetRadius(1.)
    Wrap1.SetNumberOfSides(50)

    Wrap1.Update()
    # Wrap1 = CleanPolyData(Wrap1.GetOutput())

    mapWrap1 =vtkPolyDataMapper();
    mapWrap1.SetInputConnection(Wrap1.GetOutputPort());
    actWrap1 =vtkActor();
    actWrap1.SetMapper(mapWrap1);
    actWrap1.GetProperty().SetOpacity(0.3)
 
    Wrap2 = vtkTubeFilter()
    Wrap2.SetInputConnection( readerCLINE.GetOutputPort())
    Wrap2.SetRadius(1.5)
    Wrap2.SetNumberOfSides(200)

    Wrap2.Update()
    # Wrap2 = CleanPolyData(Wrap2.GetOutput())

    mapWrap2 =vtkPolyDataMapper();
    mapWrap2.SetInputConnection(Wrap2.GetOutputPort());
    actWrap2 =vtkActor();
    actWrap2.SetMapper(mapWrap2);
    actWrap2.GetProperty().SetOpacity(0.3)
 
    Smooth1 = vtkSmoothPolyDataFilter()
    Smooth1.SetInputConnection(0, Wrap1.GetOutputPort())
    Smooth1.SetInputData(1, readerSTL.GetOutput())
    Smooth1.SetNumberOfIterations(100)
    # Smooth1.SetNumberOfIterations(1000)
    # Smooth1.SetRelaxationFactor(0.0001)
    Smooth1.SetRelaxationFactor(0.1)
    Smooth1.FeatureEdgeSmoothingOff()
    # Smooth1.FeatureEdgeSmoothingOn()
    Smooth1.BoundarySmoothingOn()
    Smooth1.Update()
    mapSmooth1 = vtkPolyDataMapper();
    mapSmooth1.SetInputConnection(Smooth1.GetOutputPort());
    actSmooth1 = vtkActor();
    actSmooth1.SetMapper(mapSmooth1);
    actSmooth1.GetProperty().SetOpacity(0.3)
    actSmooth1.GetProperty().SetColor( [0.,1.,0.] )

    if True:
        Smooth2 = vtkSmoothPolyDataFilter()
        Smooth2.SetInputConnection(0, Wrap2.GetOutputPort())
        Smooth2.SetInputData(1, readerSTL.GetOutput())
        Smooth2.SetNumberOfIterations(10)
        # Smooth2.SetNumberOfIterations(100)
        # Smooth2.SetNumberOfIterations(1000)
        Smooth2.SetRelaxationFactor(0.1)
        # Smooth2.SetRelaxationFactor(0.5)
        # Smooth2.FeatureEdgeSmoothingOff()
        # Smooth2.SetEdgeAngle(1.)
        # Smooth2.SetFeatureAngle(90.)
        Smooth2.FeatureEdgeSmoothingOn()
        Smooth2.BoundarySmoothingOn()
        # Smooth2.BoundarySmoothingOff()
        Smooth2.Update()

        # Smooth2 = TaubinFilter(Smooth2.GetOutput())
    else:
        Smooth2 = TaubinFilter(Wrap2.GetOutput())

    mapSmooth2 = vtkPolyDataMapper();
    mapSmooth2.SetInputConnection(Smooth2.GetOutputPort());
    actSmooth2 = vtkActor();
    actSmooth2.SetMapper(mapSmooth2);
    actSmooth2.GetProperty().SetOpacity(0.3)
    actSmooth2.GetProperty().SetColor( [0.,1.,1.] )

    renderer = vtkRenderer()

    renderer.SetBackground(0.1, 0.2, 0.4)
 
    renderer.AddActor(actCLINE)
    renderer.ResetCamera()
    renderer.AddActor(actSTL)
    # renderer.AddActor(actWrap1)
    # renderer.AddActor(actSmooth1)
    # renderer.AddActor(actWrap2)
    renderer.AddActor(actSmooth2)

    render_window = vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(600, 600)
 
    interactor = vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)
 
    print 'STL     n. of cells:',readerSTL.GetOutput().GetNumberOfCells()
    print 'smooth2 n. of cells:',Smooth2.GetOutput().GetNumberOfCells()

    interactor.Initialize()
    interactor.Start()
