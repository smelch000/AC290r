#!/usr/bin/env python

import sys,os
import numpy

class Inletoutlet():

    def __init__(self,msh):

        # self.itppp_i_id = numpy.zeros(msh.nv, dtype='i8')
        # self.itppp_o_id = numpy.zeros(msh.nv, dtype='i8')
        self.itppp_i_id = msh.itppp_i_id
        self.itppp_o_id = msh.itppp_o_id
        self.i_globs = msh.i_globs
        self.o_globs = msh.o_globs
        self.nin_id = 0
        self.nout_id = 0
        self.type = {}
        self.dir = {}
        self.iohead = []
        self.msh = msh

        # if self.msh.nx==0 and self.msh.ny==0 and self.msh.nz==0:
        #     print 'Called inletoutlet initializer with empty mesh data'
        #     sys.exit(1)

    def loadMOEBIUSinput(self,filename):

        if not os.path.isfile(filename):
            return False

        io = open(filename,'r')
        
        niohead = int(io.readline())
        for i in xrange(niohead):
        
            line = io.readline()
            self.iohead.append(line)
        
            line = line.split()
            id = int(line[0])
            self.dir[id-1] = [int(line[3]),int(line[4]),int(line[5])]
            self.type[id-1] = line[1]
        
            if self.type[id-1]=='inlet': 
                self.nin_id += 1
            elif self.type[id-1]=='outlet': 
                self.nout_id += 1
        
        self.iodict = []
        for i in xrange(self.nin_id + self.nout_id):
            self.iodict.append({})
        
        nin = int( io.readline().split()[0] )
        for i in xrange(nin):
            line = io.readline().split()
            i,j,k,id = int(line[0]),int(line[1]),int(line[2]),int(line[3])
        
            i4 = self.msh.i4back(i,j,k)
        
            # check
            ir = self.msh.binsearch(self.msh.ninlet, self.msh.itppp_i, i4)
            if not self.msh.itppp_i[ir]==i4 : print 'inlet point not matching mesh',ir,i4
        
            self.itppp_i_id[ir] = id
            self.iodict[id-1][i4] = [ir,self.type[id-1],id,self.dir[id-1]]
        
        nout = int( io.readline().split()[0] )
        for i in xrange(nout):
            line = io.readline().split()
            i,j,k,id = int(line[0]),int(line[1]),int(line[2]),int(line[3])
        
            i4 = self.msh.i4back(i,j,k)
        
            # check
            ir = self.msh.binsearch(self.msh.noutlet, self.msh.itppp_o, i4)
            if not self.msh.itppp_o[ir]==i4 : print 'outlet point not matching mesh',ir,i4
        
            self.itppp_o_id[ir] = id
            self.iodict[id-1][i4] = ir
        
        io.close()

        return True

    def writeMOEBIUSinput(self,filename,iohead=None):

        # standard extension for inlet/outlet
        if not filename.endswith('.ios'):
            filename = filename + '.ios'

        print >> sys.stderr, 'writing new inlet/outlet file ',filename

        if self.msh.itppp_i.shape + self.msh.itppp_o.shape == 0:

            self.msh.itppp_i = self.msh.convert_ijk_2_itppp(self.msh.ijk_i)
            self.msh.itppp_o = self.msh.convert_ijk_2_itppp(self.msh.ijk_o)

        elif self.msh.ijk_i.shape[0] + self.msh.ijk_o.shape[0] == 0:

            self.msh.ijk_i = self.msh.convert_itppp_2_ijk(self.msh.itppp_i)
            self.msh.ijk_o = self.msh.convert_itppp_2_ijk(self.msh.itppp_o)

        else:
            print >> sys.stderr, 'no itppp nor ijk for inlet/outlet....'
            sys.exit(1)

        d = open(filename,'w')

        if iohead:
            self.iohead = iohead

            # write header in inlet/outlet file
            d.write('%d\n' % (len(self.iohead)))

            for l in self.iohead:
                d.write('%s\n' % (l))
        else:

            iMAX = 0
            if len(self.msh.itppp_i_id) > 0:
                iMAX = self.msh.itppp_i_id.max()

            oMAX = 0
            if len(self.msh.itppp_o_id) > 0:
                oMAX = self.msh.itppp_o_id.max()

            IOM = max(iMAX, oMAX)

            print >> d, IOM

            # for i in xrange(1, self.msh.itppp_i_id.max()+1 ):
            if self.i_globs == None:
                for i in xrange(1, iMAX+1 ):
                    print >> d, i, ' inlet pressure   0 0 +1   .0 .0 .0     0.000'
            else:
                for glb in self.i_globs:
                    print >> d, glb['id'], 'inlet', glb['bctype'], glb['iodir'], '     ', glb['iodir'], '     ', glb['ioval']

            # for i in xrange(self.msh.itppp_o_id.max(), self.msh.itppp_o_id.max()+1 ):
            # for i in xrange( iMAX, oMAX+1 ):
            if self.o_globs == None:
                for i in xrange( 1, oMAX+1 ):
                    print >> d, i, 'outlet pressure   0 0 +1   .0 .0 .0     0.000'
            else:
                for glb in self.o_globs:
                    print >> d, glb['id'], 'outlet', glb['bctype'], glb['iodir'], '     ', glb['iodir'], '     ', glb['ioval']

        # write inlet data
        if self.msh.inlet>0: d.write('%d\n' % (self.msh.ninlet))

        for n in xrange(self.msh.ninlet):
            i,j,k = self.msh.ijk( self.msh.itppp_i[n] )
            id = self.itppp_i_id[n]
            d.write('%d %d %d %d\n' % (i,j,k,id))

        # write outlet data
        if self.msh.inlet>0: d.write('%d\n' % (self.msh.noutlet))

        for n in xrange(self.msh.noutlet):
            i,j,k = self.msh.ijk( self.msh.itppp_o[n] )
            id = self.itppp_o_id[n]
            d.write('%d %d %d %d\n' % (i,j,k,id))

        d.close()
