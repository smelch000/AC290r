#!/usr/bin/env python

import sys

# generates a sphere (closed surface, vtkPolyData) and converts it into vtkImageData
# where the foreground voxels are 1 and the background voxels 0
# Internally vtkPolyDataToImageStencil is utilized. 

RESOLUTION = 1
# RESOLUTION = 2
# RESOLUTION = 4

from vtk import *
import math,sys

# read the STL file
reader = vtkSTLReader()
reader.SetFileName( sys.argv[1] )
reader.Update()

whiteImage = vtkImageData()
bounds = 6*[0];
reader.GetOutput().GetBounds(bounds);

spacing = [1./RESOLUTION, 1./RESOLUTION, 1./RESOLUTION]

whiteImage.SetSpacing(spacing);
 
# compute dimensions
dim = 3*[0]
for i in xrange(3):
    # dim[i] = math.ceil(((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]))
    dim[i] = int(((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]))

print dim

whiteImage.SetDimensions(dim);
whiteImage.SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);
 
origin = [bounds[0] + spacing[0]/2, bounds[2] + spacing[1]/2, bounds[4] + spacing[2]/2]

whiteImage.SetOrigin(origin);
 
whiteImage.AllocateScalars(VTK_UNSIGNED_CHAR,1);

INVALUE = 255
MIDVALUE = 100
OUTVALUE = 0

# fill the image with foreground voxels:
for i in xrange( whiteImage.GetNumberOfPoints() ):
    whiteImage.GetPointData().GetScalars().SetTuple1(i, INVALUE);
 
# polygonal data - image stencil:
pol2stenc = vtkPolyDataToImageStencil()
pol2stenc.SetInputData(reader.GetOutput())

pol2stenc.SetOutputOrigin(origin)
pol2stenc.SetOutputSpacing(spacing)
pol2stenc.SetOutputWholeExtent(whiteImage.GetExtent())
pol2stenc.Update()

# cut the corresponding white image and set the background:
imgstenc = vtkImageStencil()
imgstenc.SetInputData(whiteImage)
imgstenc.SetStencilConnection(pol2stenc.GetOutputPort())
imgstenc.ReverseStencilOff()
imgstenc.SetBackgroundValue(OUTVALUE)
imgstenc.Update()
 
# thresh = vtkThreshold()
# thresh.SetInputConnection(imgstenc.GetOutputPort())
# thresh.ThresholdBetween(0, 256)

thresh = vtkImageThreshold()
thresh.SetInputConnection(imgstenc.GetOutputPort())
thresh.ThresholdByUpper(MIDVALUE)
thresh.SetInValue(INVALUE)
thresh.SetOutValue(OUTVALUE)
thresh.ReleaseDataFlagOff()
thresh.Update()

img = thresh.GetOutput()
# img = imgstenc.GetOutput()

# print img
# print img.GetPointData()
x_dim, y_dim, z_dim = img.GetDimensions()

hf = open('AAA.hdr','w')
df = open('AAA.dat','w')
dm = open('AAA.xyz','w')

print >> hf, x_dim, y_dim, z_dim

# print >> df, x_dim*y_dim*z_dim
# print >> df
nlo = 0; nhi = 0
for x in xrange(x_dim):
    for y in xrange(y_dim):
        for z in xrange(z_dim):
            val = img.GetScalarComponentAsFloat(x,y,z,0)
            if val < MIDVALUE:
                nlo += 1
            else:
                print >> dm,'C',x,y,z
                nhi += 1

print 'nodes:',nhi

# write = vtkXMLUnstructuredGridWriter()
# write.SetInputData( thresh.GetOutput() )
# write.SetFileName('thre.vtu')
# write.Update()

# writer = vtkMetaImageWriter()
# writer.SetFileName("SphereVolume.mhd");
# if VTK_MAJOR_VERSION <= 5:
#   writer.SetInput(imgstenc.GetOutput());
# else:
#   writer.SetInputData(imgstenc.GetOutput());
# writer.Write();  
 
