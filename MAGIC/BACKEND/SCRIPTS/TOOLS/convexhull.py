#!/usr/bin/env python

from myvtk import *
import sys

def fakeConvexHull(poly):

    sphereSource = vtkSphereSource()
    sphereSource.SetRadius(500)
    sphereSource.SetPhiResolution(300)
    sphereSource.SetThetaResolution(300)
    sphereSource.Update()
 
    points = poly.GetPoints()
 
    smoothFilter = vtkSmoothPolyDataFilter()

    smoothFilter.SetInputConnection(0, sphereSource.GetOutputPort())
    smoothFilter.SetInputData(1, poly)

    # smoothFilter.SetNumberOfIterations(15)
    # smoothFilter.SetNumberOfIterations(100)
    # smoothFilter.SetNumberOfIterations(200)
    smoothFilter.SetNumberOfIterations(1000)

    # smoothFilter.SetRelaxationFactor(0.1)
    # smoothFilter.SetRelaxationFactor(0.5)
    smoothFilter.SetRelaxationFactor(1.0)
    # smoothFilter.FeatureEdgeSmoothingOff()

    smoothFilter.FeatureEdgeSmoothingOn()
    smoothFilter.BoundarySmoothingOn()

    smoothFilter.Update()
 
    return smoothFilter # .GetOutputPort()

filein = sys.argv[1]
fileout = sys.argv[2]

reader = vtkSTLReader()
reader.SetFileName(filein)
reader.Update()

cvx = fakeConvexHull(reader.GetOutput())

writer = vtkSTLWriter()
writer.SetInputConnection(cvx.GetOutputPort())
writer.SetFileTypeToASCII()
writer.SetFileName(fileout)
writer.Update()

