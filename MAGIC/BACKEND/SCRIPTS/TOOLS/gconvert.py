#!/usr/bin/env python
#
import sys,os,argparse
import bgkflag2vtk
import bgkflag2vtu
import bgkflag2xyz
import bgkflag_inp2hdr
import csv2bgkflag
import dicom2vti, vti2dicom
import vtk2xml, vtk2ply
import nas2stl


###################
def main(args):

    if args.intype == 'mesh' or args.intype == 'bgkflag':

        if args.outtype == 'vtk':

            bgkflag2vtk.main(args)

        elif args.outtype == 'vtu':

            bgkflag2vtu.main(args)

        elif args.outtype == 'xyz':

            bgkflag2xyz.main(args)

        else:
            print 'unknown type of output file'

    elif args.intype == 'vtk':
        pass

    elif args.intype == 'vti':
        pass

    elif args.intype == 'csv':

        if args.outtype == 'mesh' or args.outtype == 'bgkflag':

            csv2bgkflag.main(args)

        else:
            print 'unknown type of output file'

    else:
        print 'intype not recognized', args.intype


###############################
def iniparser(parser):

    parser.add_argument('-t', '--intype',   required=True, default=None, help='input file type')
    parser.add_argument('-T', '--outtype',  required=True, default=None, help='output file type')
    parser.add_argument('-i', '--infiles', required=True, nargs='+',     default=None, help='input files prefix')
    parser.add_argument('-o', '--outfile', required=True,                default=None, help='output file prefix')
    parser.add_argument('-w', '--outwall', required=False, default=None,               help='output wall file prefix')

    return parser

###################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = iniparser(parser)
    args = parser.parse_args()

    main(args)
