#!/usr/bin/env python

import os, sys

def get_settings(outfile):

    print (outfile)
    f = open(outfile,'w')
    print >> f, "Host name:"
    f.close()
    os.system( 'hostname >> ' + outfile)


    f = open(outfile, 'a')
    print >> f, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print >> f, "  PROCESSOR"
    print >> f

    print >> f, "Byte order: %s" % sys.byteorder
    try:
        print >> f, "Float info: %s" % sys.float_info
    except AttributeError:
        print >> f, "Float info not available"
    print >> f, "Platform: %s" % sys.platform

    print >> f, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print >> f, " PATHS"
    print >> f

    print >> f, "System Path ($PATH):"
    f.close()
    os.system( 'echo $PATH >> ' + outfile)
    f = open(outfile, 'a')
    print >> f, "Python Path (sys.path):"
    print >> f, sys.path

    print >> f, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print >> f, " PACKAGE VERSIONS"
    print >> f

    print >> f, "Python:",
    print >> f, sys.version
    print >> f, sys.version_info, '\n'

    print >> f, "Numpy:",
    try:
        import numpy
        print >> f, numpy.version.version
    except ImportError:
        print >> f, "Not available"

    print >> f, "Scipy:",
    try:
        import scipy
        print >> f, scipy.version.version
    except ImportError:
        print >> f, "Not available"

    print >> f, "Matplotlib:",
    try:
        import matplotlib
        print >> f, matplotlib.__version__
    except ImportError:
        print >> f, "Not available"

    print >> f, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print >> f, " ENVIRONMENT"
    print >> f
    for k,v in os.environ.items():
        print >> f, ( "%s: %s" % (k,v) )

    f.close()


if __name__ == '__main__':

    outfile = 'system_settings.txt'
    print ('writing sys settings on file: ', outfile)

    get_settings(outfile)
