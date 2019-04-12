#!/usr/bin/env python

import argparse,sys,os
from myvtk import *

# read vtk and convert to xml
def main(vtkfile, xmlfile):

    vtk_reader = vtk.vtkGenericDataObjectReader()
    vtk_reader.SetFileName(args.infile)
    vtk_reader.Update()
    print 'reading vtk...:',args.infile

    # passA = vtk.vtkPassArrays()
    # passA.SetInputData( vtk_reader.GetOutput() )
    # passA.SetInput( vtk_reader.GetOutput() )
    # passA.SetRemoveArrays(True)
    # # passA.RemoveArraysOn()
    # passA.ClearArrays() # clear all arrays
    # passA.ClearPointDataArrays()
    # passA.ClearCellDataArrays()
    # passA.ClearFieldDataArrays()
    # passA.ClearFieldTypes() # clear all
    # passA.Update()

    # print vtk_reader.GetOutput()

    writer = vtk.vtkXMLUnstructuredGridWriter()
    # writer.SetDataModeToAscii()
    writer.SetDataModeToBinary()
    writer.SetInputConnection(vtk_reader.GetOutputPort())
    # writer.SetInputConnection(passA.GetOutputPort())
    writer.SetFileName(args.outfile)
    writer.Write()
    print 'written xml file into ',args.outfile

##############################
def iniparser(parser):

    parser.add_argument('-i', '--infile', default=None,  help='input image file (.vti)')
    parser.add_argument('-o', '--outfile', default=None,  help='output image file (.xml)')

    return parser

###############################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = iniparser(parser)
    args = parser.parse_args()

    main(args)
