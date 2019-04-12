#!/usr/bin/env python

import getopt,sys,os,argparse
from myvtk import *

def convert_dicom_meta_image(datadir, outfile):

	# read dicom dir and convert to vti single file 
	dicom_reader = vtk.vtkDICOMImageReader()
	dicom_reader.SetDirectoryName(datadir)
	dicom_reader.Update()
	print 'reading dicom...:',datadir

	writer = vtk.vtkMetaImageWriter()
	writer.SetInput(dicom_reader.GetOutput())
	if vtk.VTK_MAJOR_VERSION <= 5:
		writer.SetInput(dicom_reader.GetOutput())
	else:
		writer.SetInputData(dicom_reader.GetOutput())
	writer.SetFileName(outfile)
	writer.Write()
	print 'written MetaImage file into ',outfile

def main(args):

    datadir = args.dicomdir
    outfile = args.outfile

    # read dicom dir and convert to vti single file 
    dicom_reader = vtk.vtkDICOMImageReader()
    dicom_reader.SetDirectoryName(datadir)
    dicom_reader.Update()
    print 'reading dicom...:',datadir

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetInputData(dicom_reader.GetOutput())
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(dicom_reader.GetOutput())
    else:
        writer.SetInputData(dicom_reader.GetOutput())
    writer.SetFileName(outfile)
    writer.Write()
    print 'written vti file into ',outfile

###############################
def iniparser(parser):

    parser.add_argument('-d', '--dicomdir', default='bgkflag',  help='dicom folder')
    parser.add_argument('-o', '--outfile', default='bgkflag',  help='output image file (.vti)')

    return parser

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = iniparser(parser)
    args = parser.parse_args()

    main(args)
