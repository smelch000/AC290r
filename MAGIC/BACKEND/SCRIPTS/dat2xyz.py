#!/usr/bin/python

import sys

if (__name__ == "__main__"):
  if len(sys.argv) < 2:
    print "\tUsage: dat2xyz filebasename\n\tConvert <filebasename>.dat into <filebasename>.xyz format.\n"
    sys.exit(1)
  linee=[]
  try:
    nodi=open(sys.argv[1]+'.dat','r').readlines()
  except:
    sys.stderr.write("Error. Problem in reading "+sys.argv[1]+".dat\n")
    sys.exit(1)
  for nodo in nodi:
    n=nodo.split()
    linee.append(n[3]+" "+n[0]+" "+n[1]+" "+n[2])

  f=open(sys.argv[1]+".xyz",'w')
  f.write(str(len(linee))+"\n1\n")
  for linea in linee:
    f.write(linea+"\n")

  f.close()

