import sys
from myvtk import *

if len(sys.argv) != 3:
    print "Usage:", sys.argv[0]," -i infile -o outfile"
    sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]

pdr = vtk.vtkPolyDataReader()
pdr.SetFileName(infile)
pdr.Update()

f = open(outfile,'w')
for i in xrange(pdr.GetOutput().GetNumberOfPoints()):
    x,y,z = pdr.GetOutput().GetPoint(i)
    print >> f, x,y,z

