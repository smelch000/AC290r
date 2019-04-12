#!/usr/bin/env python

import argparse,sys,os
from myvtk import *

# read vtk and convert to ply
def main(args):

	vtk_reader = vtk.vtkGenericDataObjectReader()
	vtk_reader.SetFileName(args.infile)
	vtk_reader.Update()
	# print 'reading vtk...:',args.infile

	writer = vtk.vtkPLYWriter()
	writer.SetFileTypeToASCII()
	writer.SetInput(vtk_reader.GetOutput())
	writer.SetFileName(args.outfile)
	writer.Write()
	# print 'written ply file into ',args.outfile


###############################
def iniparser(parser):

    parser.add_argument('-i', '--infile', default=None,  help='input mesh file (.vtk)')
    parser.add_argument('-o', '--outfile', default=None,  help='output mesh file (.ply)')

    return parser

###############################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = iniparser(parser)
    args = parser.parse_args()

    main(args)
