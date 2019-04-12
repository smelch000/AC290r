#!/usr/bin/python
# convert an "atom.inp" file into a XYZ file
from atomhandler import *
import sys

if __name__ == '__main__':

    if len( sys.argv ) != 3:
        print 'Usage :',sys.argv[0],' inpfile xyzfile'
        sys.exit(1)

    a=AtomHandler(sys.argv[1])

    a.loadAtoms()

    a.setFileName(sys.argv[2])
    a.saveAtoms(formato='XYZ')


