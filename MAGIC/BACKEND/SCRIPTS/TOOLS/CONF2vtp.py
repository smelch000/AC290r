#!/usr/bin/env python
  
import sys
import argparse
from vtk import *

###############################
def writevtp(x,y,z,filenm):

    points = vtkPoints()
 
    for i in xrange(len(x)):
        points.InsertNextPoint( x[i],y[i],z[i] )
 
    # Create a polydata object and add the points to it.
    polydata = vtkPolyData()
    polydata.SetPoints(points)
 
    writer =  vtkXMLPolyDataWriter()
    writer.SetFileName(filenm)
    writer.SetInputData(polydata)
 
    # Optional - set the mode. The default is binary.
    writer.SetDataModeToBinary();
    #writer.SetDataModeToAscii();
 
    writer.Write();

###############################
def main(args):
    print args

    f = open(args.infile,'r')

    for itime in xrange(100000):

        try:
            natms = int(f.readline())
            print 'natms:',natms
            f.readline()

            x,y,z = [],[],[]
            for i in xrange(natms):
                l = f.readline().split()
                x.append(float(l[1]))
                y.append(float(l[2]))
                z.append(float(l[3]))

            rt = '%d'%itime
            rt = rt.zfill(3)
            writevtp(x,y,z,'CONF_'+rt+'.vtp')
        except:
            break
        

###############################
def iniparser(parser):

    parser.add_argument('-i', '--infile', default='CONF.xyz',  help='input file (.xyz)')
    parser.add_argument('-o', '--outfile', default=None,  help='output file (.vtp)')

    return parser

###############################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = iniparser(parser)
    args = parser.parse_args()

    main(args)


