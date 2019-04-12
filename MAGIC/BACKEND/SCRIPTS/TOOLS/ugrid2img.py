#!/usr/bin/env python

import sys,os,argparse,re
from vtk import *

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',  required=True,  help='input (vtk or ensight) file')
parser.add_argument('-s', '--spacing',  required=True,  help='new spacing')
args = parser.parse_args()

ext = os.path.splitext(args.input)[-1].lower()

if ext == '.vtk':

    reader = vtkUnstructuredGridReader()
    reader.SetFileName( args.input  )
    reader.Update()

    ug = reader.GetOutput()

elif ext == '.vtu':

    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName( args.input  )
    reader.Update()

    ug = reader.GetOutput()

elif ext == '.case':

    reader = vtkEnSightGoldReader()
    reader.SetCaseFileName( args.input  )
    # reader.SetTimeValue( 1 )
    reader.Update()
    geom = vtkGeometryFilter()
    geom.SetInputConnection( reader.GetOutputPort() )
    geom.GetOutput()
    ug = geom.GetOutput()

# reader.SetFileName( 'bifurca.vtk'  )


bounds = 6*[0]
ug.GetBounds(bounds)
print 'ug bounds:',bounds

# spacing = [0.5, 0.5, 0.5] 
# spacing = [1.0, 1.0, 1.0]
# spacing = [2.0, 2.0, 2.0]
spacing = [float(args.spacing), float(args.spacing), float(args.spacing)]

Image = vtkImageData()
Image.SetSpacing(spacing)
 
dim = [int((bounds[1] - bounds[0]) / spacing[0]),
       int((bounds[3] - bounds[2]) / spacing[1]),
       int((bounds[5] - bounds[4]) / spacing[2]) ]

origin = [bounds[0] + spacing[0]/2, 
          bounds[2] + spacing[1]/2, 
          bounds[4] + spacing[2]/2]

Image.SetDimensions(dim)
Image.SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1)
Image.SetOrigin(origin)
 
if VTK_MAJOR_VERSION <= 5:
  Image.SetScalarTypeToUnsignedChar()
  Image.AllocateScalars()
else:
  # Image.AllocateScalars(VTK_UNSIGNED_CHAR,1)
  Image.AllocateScalars(VTK_FLOAT, 2)

xprobe = vtkProbeFilter()
if VTK_MAJOR_VERSION <= 5:
    xprobe.SetInput(Image)
    xprobe.SetSource(ug)
else:
    xprobe.SetInputData(Image)
    xprobe.SetSourceConnection(reader.GetOutputPort())
xprobe.PassPointArraysOn()
# xprobe.ComputeToleranceOn()
# xprobe.SetTolerance(1.e-3)
# xprobe.SpatialMatchOn()
# xprobe.SpatialMatchOff()
xprobe.Update()
print 'probing done, tolerance',xprobe.GetTolerance()

probeout_typecast = vtkImageCast()
probeout_typecast.SetInputData(xprobe.GetOutput())
probeout_typecast.SetOutputScalarTypeToUnsignedShort()
probeout_typecast.Update()

writer = vtkXMLImageDataWriter()
writer.SetFileName('pippo2.vti')
if VTK_MAJOR_VERSION <= 5:
    # writer.SetInputConnection(Image.GetProducerPort())
    writer.SetInputConnection(probeout_typecase.GetProducerPort())
else:
    # writer.SetInputData(Image)
    writer.SetInputData(probeout_typecast.GetOutput())
writer.Write()

