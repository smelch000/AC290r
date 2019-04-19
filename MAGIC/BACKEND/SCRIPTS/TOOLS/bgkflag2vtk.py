#!/usr/bin/env python

import sys,os
import argparse
import ctypes

from TOOLS.mesh import *

def iniparser(parser):

    parser.add_argument('-i', '--infile', required=False,  help='mesh file (.dat)')
    parser.add_argument('-o', '--outfile', required=True,  help='output complete vtp file')
    parser.add_argument('-w', '--outwall', required=False, default=None, help='output wall vtp file')
    parser.add_argument('-id', '--identity', required=False, default=None, help='node identity')

    return parser

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = iniparser(parser)
    args = parser.parse_args()

    msh = Mesh()

    pd = msh.meshfileToPolyData(args.infile)
    msh.writePolyData(pd, msh, args.outfile)
