#!/usr/bin/env python

from myvtk import *
import argparse
import os

def downsampleImageFilter(imagefilter, dlevel):
    """ 
    Downsampling filter
    """

    filt = vtkImageShrink3D()

    if vtk.VTK_MAJOR_VERSION <= 5:
        filt.SetInput(imagefilter.GetOutput())
    else:
        filt.SetInputConnection(imagefilter.GetOutputPort())

    filt.SetShrinkFactors(dlevel,dlevel,dlevel)

    filt.Update()

    return filt

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Downsample image.')

    parser.add_argument('-i', '--input',   required=True,  help='input VTI file')
    parser.add_argument('-d', '--downsample',   required=True,  help='downsampling frequency 2,3,4...')
    parser.add_argument('-o', '--output',  required=True,  help='output VTI file')

    args = parser.parse_args()

    FILEIN = args.input
    DOWNSAMPLE = int(args.downsample)
    FILEOUT = args.output

    # Read the file

    froot,fext = os.path.splitext(FILEIN)

    if fext == '.vti':
        origreader = vtkXMLImageDataReader()
    elif fext =='.mhd' or fext =='.mha':
        origreader = vtkMetaImageReader()
    else:
        print 'input img format unrecognized'
        sys.exit(1)

    origreader.SetFileName(FILEIN)
    origreader.Update()

    shrinkFilter = downsampleImageFilter(origreader, DOWNSAMPLE)

    """
    # Downsampling filter
    shrinkFilter = vtkImageShrink3D()
    if vtk.VTK_MAJOR_VERSION <= 5:
        shrinkFilter.SetInput(origreader.GetOutput());
    else:
        shrinkFilter.SetInputConnection(origreader.GetOutputPort());
    shrinkFilter.SetShrinkFactors(DOWNSAMPLE,DOWNSAMPLE,DOWNSAMPLE)
    """

    fname, ext = os.path.splitext(FILEOUT)
    # FILEOUT = fname +'_downsampled_'+str(DOWNSAMPLE)+'.vti'

    if fext == '.vti':
        writer = vtkXMLImageDataWriter()
        writer.SetFileName( FILEOUT )
    elif fext =='.mhd' or fext =='.mha':
        writer = vtkMetaImageWriter()
        if fext =='.mha': print 'Warning ... renaming file as .mhd/.raw pair (vtk bug)'
        writer.SetFileName( fname+'.mhd' )
        writer.SetRAWFileName( fname+'.raw' )
    else:
        print 'output img format unrecognized'
        sys.exit(1)

    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput( shrinkFilter.GetOutput() )
    else:
        writer.SetInputConnection( shrinkFilter.GetOutputPort() )

    froot,fext = os.path.splitext(FILEOUT)

    writer.Update()
    writer.Write()

    print '...done writing files'
