#!/usr/bin/env python

import sys
import numpy as np
import argparse
from TOOLS.bgkflag_IO import *
from TOOLS.mesh import *

def main(args):

    nx,ny,nz,iskip,nfluid,nwall,ninlet,noutlet = bgkflag_hdr_read(args.infiles[0]+'.hdr')

    WASTE = False

    if WASTE:
        iflg,nfluid,nwall,ninlet,noutlet = bgkflag_dat_memwaste_read(args.infiles[0]+'.dat',nx,ny,nz)
    else:
        ijkf,nfluid,nwall,ninlet,noutlet = bgkflag_dat_read(args.infiles[0]+'.dat')

    print 'nfluid,nwall,ninlet,noutlet',nfluid,nwall,ninlet,noutlet

    print 'writing file ...',args.outfile+'.xyz'

    ga = open(args.outfile+'.xyz','w')

    if args.outwall != None:
        print 'and writing file ',args.outwall+'.xyz'
        gw = open(args.outwall+'.xyz','w')
        gw.write('%d\n\n' % (nwall) )
    ga.write('%d\n\n' % (nfluid + nwall + ninlet + noutlet) )
    # ga.write('%d\n\n' % (ninlet + noutlet) )

    if WASTE:
        for k in xrange(nz):
            for j in xrange(ny):
                for i in xrange(nx):

                    fl = iflg[i,j,k]
                    if fl == Mesh().ID_FLUIDNODE: ga.write('F %d %d %d\n' % (i+1,j+1,k+1))

                    elif fl == Mesh().ID_WALLNODE:
                        ga.write('H %d %d %d\n' % (i+1,j+1,k+1))
                        if args.outwall != None: gw.write('H %d %d %d\n' % (i+1,j+1,k+1))

                    elif fl == Mesh().ID_INLETNODE: ga.write('C %d %d %d\n' % (i+1,j+1,k+1))

                    elif fl == Mesh().ID_OUTLETNODE: ga.write('O %d %d %d\n' % (i+1,j+1,k+1))
    else:
        for n in xrange(len(ijkf)):
            i,j,k,fl = ijkf[n,:]
            if fl == Mesh().ID_FLUIDNODE: ga.write('F %d %d %d\n' % (i+1,j+1,k+1))

            elif fl == Mesh().ID_WALLNODE: 
                ga.write('H %d %d %d\n' % (i+1,j+1,k+1))
                if args.outwall != None:
                    gw.write('H %d %d %d\n' % (i+1,j+1,k+1))

            elif fl == Mesh().ID_INLETNODE: ga.write('C %d %d %d\n' % (i+1,j+1,k+1))

            elif fl == Mesh().ID_OUTLETNODE: ga.write('O %d %d %d\n' % (i+1,j+1,k+1))

    ga.close()

    if args.outwall != None:
        gw.close()

###############################
def iniparser(parser):

    parser.add_argument('-i', '--infiles', nargs = '+', default='bgkflag', required=True, help='input mesh file (.dat)')
    parser.add_argument('-o', '--outfile', default='bgkflag',  required=False, help='output mesh file (.xyz)')
    parser.add_argument('-w', '--outwall', default=None,  help='output wall file (.xyz)')

    return parser

###############################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = iniparser(parser)
    args = parser.parse_args()

    print '\n *** ',sys.argv[0],' ***\n'
    
    main(args)
