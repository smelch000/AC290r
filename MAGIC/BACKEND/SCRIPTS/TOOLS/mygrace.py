#/usr/bin/env python
import sys
import operator

# Example of usage:
#
# class o:
#    pass
#obs = [o(),o()]
#obs[0].name = 'CAZ'
#obs[0].dct = {1:10, 2:20, 3:18}
#obs[1].name = 'TAZ'
#obs[1].dct = {1:9, 2:21, 3:15}
#
#mygrace.writevalues('pippo.agr','TITOLO','XVALS','YVALS',obs)

#
def writevalues(fname,title,xname,yname,obs):

    fa = open(fname, 'wb')

    print >> fa,'@ title "' + title + '"'
    print >> fa,'@ yaxis  label "' + yname + '"'
    print >> fa,'@ xaxis  label "' + xname + '"'
    print >> fa,'@ xaxis  label char size 1.200000'
    print >> fa,'@ yaxis  label char size 1.200000'

    inc = 0
    for o in obs:
        print >> fa, '@ s' + str(inc) + ' legend  "' + o.name + '"'
        inc += 1

    for o in obs:
        for k in sorted(o.dct):
            print >> fa, k, o.dct[k]
        print >> fa
        
    
