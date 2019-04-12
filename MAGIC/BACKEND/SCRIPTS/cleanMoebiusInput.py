#!/usr/bin/env python
from moebiusBox import *
import sys

if (__name__ == "__main__"):
  if len(sys.argv) < 3:
      print "Usage: cleanMoebiusInput.py basename maxStepping [inletheader] [outletheader]"
      print "\tRead <basename>n.[hdr|dat] files, where n=1,2,4,8,...,<maxStepping> "
      print "\tand clean bad nodes in the inter-resolution interfaces, plus"
      print "\toperating a walls and holes repairing."
      print "\tWrite resulting meshes into corresponding 'c<basename>n' files."
      print "\tDefault inlet string: \"inlet velocity      0  0   -1   0  0  1  0.002 \""
      print "\tDefault outlet string: \"outlet pressure  0  0   1  0  0  -1   0.0\""
      print "\tWARNING: Only the given inlet and outlet string are put in IOS files."
      sys.exit(0)

  MAXSTEPPING=int(sys.argv[2])
  BASEFILE=sys.argv[1]
  #SCALE=8

  steps=[]
  for i in range(0,21):
    steps=[pow(2,i)]+steps
    if pow(2,i)==MAXSTEPPING: break

  if i==21: raise Exception("Error: <maxStepping> must be a power of 2 ")
  #steps=[64,32,16,8,4,2,1]
  #steps=[8,4,2,1]

    
  """
  if SCALE==8:
    nx,ny,nz=832,832,4256 #scale 8
    #nx,ny,nz=896,896,4352 #scale 8
  elif SCALE==4:
    nx,ny,nz=512,512,2240 #scale 4
  elif SCALE==2:
    nx,ny,nz=256,256,1120

  """
  instring=" inlet velocity      0  0   -1   0  0  1  0.002 " if len(sys.argv)<4 else sys.argv[3]
  outstring=" outlet pressure  0  0   1  0  0  -1   0.0" if len(sys.argv)<5 else sys.argv[4]
  f=[]
  w=[]
  i=[]
  o=[]
  for s in steps:
    ff=open(BASEFILE+str(s)+".hdr",'r')
    [nx,ny,nz]=ff.readline().split()
    nx,ny,nz=int(nx),int(ny),int(nz)
    [defnode,nfl,nwl,nin,nout]=ff.readline().split()
    defnode,nfl,nwl,nin,nout=int(defnode),int(nfl),int(nwl),int(nin),int(nout)
    ff.close()
    f.append(moebiusBox('FLUID',nx,ny,nz,s))
    w.append(moebiusBox('WALL',nx,ny,nz,s) if nwl>0 else False) 
    i.append(moebiusBox('INLET',nx,ny,nz,s,instring) if nin>0 else False) 
    o.append(moebiusBox('OUTLET',nx,ny,nz,s,outstring) if nout>0 else False) 
    

  for j,s in enumerate(steps):
    f[j].load(BASEFILE+str(s),'FLUID')
    if w[j]:
       w[j].load(BASEFILE+str(s),'WALL')
       f[j].add([w[j]])
    if i[j]:
       i[j].load(BASEFILE+str(s),'INLET')
       f[j].add([i[j]])
    if o[j]:
       o[j].load(BASEFILE+str(s),'OUTLET')
       f[j].add([o[j]])
        
  print "Nodes repaired: ",f[0].repairFluidHoles(f[1]), "on mesh ",f[0].spacing
  for j in range(1,len(f)-1):
    print "Nodes repaired: ",f[j].repairFluidHoles(f[j-1],f[j+1]), "on mesh ",f[j].spacing
  print "Nodes repaired: ",f[-1].repairFluidHoles(f[-2]), "on mesh ",f[-1].spacing

  print "Walls added: ",f[0].repairWalls(f[1]), "on mesh ",f[0].spacing
  for j in range(1,len(f)-1):
    print "Walls added: ",f[j].repairWalls(f[j-1],f[j+1]), "on mesh ",f[j].spacing
  print "Walls added: ",f[-1].repairWalls(f[-2]), "on mesh ",f[-1].spacing


  for j in range(1,len(f)):
    print "Nodes cleaned:", f[j].removeQ27(f[j-1]), "on mesh ",f[j].spacing

  for j,s in enumerate(steps):
    f[j].write('c'+BASEFILE+str(s),'DAT')
    f[j].write('c'+BASEFILE+str(s),'XYZ')
    if i[j] or o[j]: f[j].write('c'+BASEFILE+str(s),'IOS')
