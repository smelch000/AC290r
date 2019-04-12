#!/usr/bin/env python

import getopt, sys, os
from myvtk import *

###################
def usage():
    print >> sys.stderr, "Usage:",sys.argv[0],"-i vtk_file -o vtk_envelope_file -s stl_envelope_file"
    print >> sys.stderr, "Options:"
    print >> sys.stderr, ""
    print >> sys.stderr, "\t-i vtk_file"
    print >> sys.stderr, "\t--input vtk_file"
    print >> sys.stderr, "\t\tSpecifies volume-mesh file (bgkflag originated)."
    print >> sys.stderr, ""
    print >> sys.stderr, "\t-o vtk_envelope_file"
    print >> sys.stderr, "\t--output vtk_envelope_file."
    print >> sys.stderr, "\t\tIf unset stdout is used."
    print >> sys.stderr, ""
    print >> sys.stderr, "\t-s stl_envelope_file"
    print >> sys.stderr, "\t--output_stl stl_envelope_file."
    


MESHFNAME=""
OUTFNAME=""
STLFNAME=""

opts, args = getopt.getopt(sys.argv[1:], "i:o:s:", ["input=","output=","output_stl="])
for o, a in opts:
    if o in ("-i", "--input"):
        if (a == "-"):
            MESHFNAME = "/dev/stdin"
        else:
            MESHFNAME = a
    elif o in ("-o", "--output"):
        OUTFNAME = a
    elif o in ("-s", "--output_stl"):
        STLFNAME = a
    else:
        usage()
        sys.exit(1)

if (MESHFNAME == ""):
    usage()
    sys.exit(1)

if (OUTFNAME == ""):
    outf = sys.stdout
else:
    try:
        outf = open(OUTFNAME, 'w')
    except IOError:
        print >> sys.stderr, "Cannot write file", OUTFNAME
        sys.exit(1)

if (STLFNAME == ""):
    usage()
    sys.exit(1)


print 'reading vtk...'
reader = vtk.vtkDataSetReader()
reader.SetFileName(MESHFNAME)
reader.Update()

# reader = vtk.vtkPointSource()
# reader.SetNumberOfPoints(250)

print 'converting unstruct mesh to polydata ...'
gF =  vtk.vtkGeometryFilter()
gF.SetInputConnection( reader.GetOutputPort() )
gF.Update()
        
print 'cleaning data...'
cleaner = vtk.vtkCleanPolyData()
cleaner.SetInputConnection (gF.GetOutputPort())
cleaner.Update()

print 'delaunay envelope...'
delaunay3D = vtk.vtkDelaunay3D()
delaunay3D.SetInputConnection (cleaner.GetOutputPort())
delaunay3D.SetAlpha(1.0)

# print 'converting unstruct mesh to polydata ...'
gF =  vtk.vtkGeometryFilter()
gF.SetInputConnection( delaunay3D.GetOutputPort() )
gF.Update()
        
print 'write vtk file --->',OUTFNAME
writer = vtk.vtkDataSetWriter()
# writer.SetInputConnection ( gF.GetOutputPort() )
writer.SetInputConnection ( delaunay3D.GetOutputPort() )
writer.SetFileName(OUTFNAME)
writer.Write()

print 'write stl file --->',STLFNAME
writer = vtk.vtkSTLWriter()
writer.SetInputConnection ( gF.GetOutputPort() )
writer.SetFileName(STLFNAME)
writer.Write()

print 'Done!'

