#!/usr/bin/env python

import sys
from myvtk import *


filein = sys.argv[1]
fileout = sys.argv[2]

reader = vtkSTLReader()
reader.SetFileName(filein)
reader.Update()

normals = vtkPolyDataNormals()
normals.SetInputConnection(reader.GetOutputPort())
normals.ConsistencyOn()
normals.SplittingOff()
normals.Update()

# input = reader.GetOutput()
# normals.GetOutput().GetPointData().SetNormals(input.GetPointData().GetNormals())

reverseSense = vtkReverseSense()
reverseSense.SetInputConnection(normals.GetOutputPort());
# reverseSense.SetInputConnection(reader.GetOutputPort());
# reverseSense.SetInput(reader.GetOutput());
reverseSense.ReverseNormalsOn();
reverseSense.Update();

# map = vtk.vtkPolyData()
# map.SetInputConnection(reverseSense.GetOutputPort())
 
normals = reverseSense.GetOutput().GetPointData().GetNormals()

writer = vtkSTLWriter()
# writer.SetInputConnection(normals.GetOutputPort())
writer.SetInputConnection(reverseSense.GetOutputPort())
# writer.SetInputData(map.GetOutput())
writer.SetFileTypeToASCII()
writer.SetFileName(fileout)
writer.Update()
