#!/usr/bin/env python

def main(args):

    f = open(args.infile+'.inp','r')

    h = open(args.outfile+'.hdr','w')
    d = open(args.outfile+'.dat','w')

    lh = f.readline()
    h.write(lh)
    lh = f.readline()

    nf = 0
    nw = 0
    ni = 0
    no = 0
    for l in f.readlines():

        d.write(l)

        l=l.split()

        i,j,k,t=int(l[0]),int(l[1]),int(l[2]),int(l[3])

        if t==1: nf+=1
        if t==2: nw+=1
        if t==3: ni+=1
        if t==4: no+=1

    h.write(lh[0]+' %d %d %d %d\n' % (nf,nw,ni,no))

###############################
def iniparser(parser):

    parser.add_argument('-i', '--infile', default='bgkflag',  help='input file prefix (.inp)')
    parser.add_argument('-o', '--outfile', default='bgkflag',  help='output file prefix (.dat/hdr)')

    return parser

###############################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = iniparser(parser)
    args = parser.parse_args()

    main(args)
