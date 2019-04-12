#!/usr/bin/python
# read an STL file and convert it to a INP file in which
# each "atom" is a triangle vertex to be used for bouzidi BC
import sys
from atomhandler import *

if (__name__ == "__main__"):
  if len(sys.argv) < 2:
    print "\tUsage: stl2inp filebasename\n\tConvert <filebasename>.stl into <filebasename>.inp format.\n"
    sys.exit(1)
  a=AtomHandler(sys.argv[1]+'.stl')
  a.loadAtoms()
  a.setFileName(sys.argv[1]+'.inp')
  a.saveAtoms()
  a.setFileName(sys.argv[1]+'.xyz')
  a.saveAtoms("XYZ")
  
  
  
  
