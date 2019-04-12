#!/usr/bin/env python

import sys,os,argparse
from myvtk import *
# import dicom
# from dicom import *

def main(args):

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(args.infile)
    reader.Update()
    print 'reading vti file ',args.infile

    #dicom_writer = vtk.vtkDICOMImageWriter()
    #dicom_writer.SetDirectoryName(args.datadir)
    #dicom_writer.Update()
    #print 'writing dicom...:',args.datadir

    # meta = vtk.vtkDICOMMetaData()
    # meta = vtkDICOMMetaData()
    # meta.SetAttributeValue( PatientName , " Doe John ")
    # meta.SetAttributeValue( ScanningSequence , " GR ")
    # meta.SetAttributeValue( SequenceVatiant , " SP ")
    # meta.SetAttributeValue( ScanOptions , "")
    # meta.SetAttributeValue( MRAcquisitionType , "2D")

    # Plug the generator and meta data into the writer.
    # writer = vtk.vtkDICOMWriter()
    writer = vtk.vtkDICOMWriter()
    writer.SetInputConnection( reader.GetOutputPort() )
    # writer.SetMetaData( meta )
    writer.SetGenerator( generator )
    # Set the output filename format as a printf-style string
    writer.SetFilePattern ("%s/IM -0001 -%04.4 d. dcm ")
    writer.SetFilePrefix(args.datadir) # Set the directory to write the files into.
    writer.Write()

##############################
def iniparser(parser):

    parser.add_argument('-i', '--infile', default=None,  help='input image file (.vti)')
    parser.add_argument('-d', '--datadir', default=None,  help='output dicom directory')

    return parser

###############################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = iniparser(parser)
    args = parser.parse_args()

    main(args)
