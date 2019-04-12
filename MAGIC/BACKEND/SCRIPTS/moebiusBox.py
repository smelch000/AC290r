#!/usr/bin/env python
# 6-Jun-2018
# P.M.

""" usage example:

from moebiusBox import *
f2=moebiusBox('FLUID',200,48,48,2)    
f1=moebiusBox('FLUID',200,48,48,1)    

w=moebiusBox('WALL',200,48,48,2)    
i=moebiusBox('INLET',200,48,48,2,"inlet flow      -1  0   0   1  0  0  0.002 ")    
o=moebiusBox('OUTLET',200,48,48,2,"outlet pressure  +1  0   0  -1  0  0   0.0")    

f2.doBox([4,4,4],[198,46,46])

i.doLayerYZ([2,2],[48,48],2)
o.doLayerYZ([2,2],[48,48],200)

w.doBox([2,2,2],[200,48,48])
w.doBox([4,4,4],[198,46,46],remove=True)

f1.doBox([78,10,10],[122,38,38])
f2.doBox([82,24,24],[118,26,26],remove=True) #create a "hole" in low-res mesh

f2.add([w,i,o])
f2.write('bgkflag2','XYZ')
f2.write('bgkflag2','DAT')
f2.write('bgkflag2','IOS')

f1.write('bgkflag1','DAT')
f1.write('bgkflag1','XYZ')
"""
import math
import sys

class moebiusBox:
  #     0 1 2 3  4  5  6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
  dirx=[0,1,0,0,-1, 0, 0,1,1,0,-1,-1, 0, 1, 1, 0,-1,-1, 0,-1, 1,-1, 1,-1, 1,-1, 1]
  diry=[0,0,1,0, 0,-1, 0,1,0,1, 1, 0,-1,-1, 0, 1,-1, 0,-1,-1,-1, 1, 1,-1,-1, 1, 1]
  dirz=[0,0,0,1, 0, 0,-1,0,1,1, 0, 1, 1, 0,-1,-1, 0,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1]
  NDIR=19
  NDIREXT=27

  def __init__(self,kind,NX,NY,NZ,spacing,header=""):
    #header is to be put on the top of ios file for I/O nodes
    # e.g. header="inlet flow      -1  0   0   1  0  0  0.002"
    self.kind=self.__getType(kind)
    if kind == 'INLET' or kind == 'OUTLET':
      self.header=header
      
    self.spacing=spacing
    if (NX % spacing) + (NY % spacing) + (NZ % spacing) != 0:
      sys.stderr.write("WARNING: bounding box are not multiple of spacing!\n")
      
    self.NX,self.NY,self.NZ=int(NX/spacing),int(NY/spacing),int(NZ/spacing)
    self.nodes={}
    self.mynodes={}
    self.objlist=[self] #add itself in the added object list

  def doBox(self,fromVertex,toVertex,remove=False):
    fromVertex=self.__normalise(fromVertex)
    toVertex=self.__normalise(toVertex)
    if remove:
      if self.nodes:
        for k in range(fromVertex[2],toVertex[2]+1): 
          for j in range(fromVertex[1],toVertex[1]+1): 
            for i in range(fromVertex[0],toVertex[0]+1):
                try:
                  del self.nodes[self.__roll(i,j,k)]
                except:
                  pass
              
    else:
      #if len(self.nodes) >0: self.node.clear()
      for k in range(fromVertex[2],toVertex[2]+1): 
        for j in range(fromVertex[1],toVertex[1]+1): 
          for i in range(fromVertex[0],toVertex[0]+1):
              self.__addnode(i,j,k)
            
   
  def doLayerXY(self,fromXY,toXY,Z,remove=False):
    self.doBox([fromXY[0],fromXY[1],Z],[toXY[0],toXY[1],Z],remove)
  def doLayerXZ(self,fromXZ,toXZ,Y,remove=False):
    self.doBox([fromXZ[0],Y,fromXZ[1]],[toXZ[0],Y,toXZ[1]],remove)
  def doLayerYZ(self,fromYZ,toYZ,X,remove=False):
    self.doBox([X,fromYZ[0],fromYZ[1]],[X,toYZ[0],toYZ[1]],remove)
   
  def doCylinder(self,baseCenter,otherBaseCenter,baseRadius,remove=False):
    b=self.__normalise(baseCenter,True)
    d=self.__normalise(otherBaseCenter,True)
    d=[d[0]-b[0],d[1]-b[1],d[2]-b[2]]
    height=math.sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2])
    r=float(baseRadius)/float(self.spacing)    
    if remove:
      if self.nodes:
        for k in range(1,self.NZ+1): 
          for j in range(1,self.NY+1): 
            for i in range(1,self.NX+1):
              p=[float(i),float(j),float(k)]
              if self.__isInCylinder(p[0],p[1],p[2],b[0],b[1],b[2],d[0],d[1],d[2],r,height):
                try:
                  del self.nodes[self.__roll(i,j,k)]
                except:
                  pass
    else:         
      #if len(self.nodes) >0: self.node.clear()
      for k in range(1,self.NZ+1): 
        for j in range(1,self.NY+1): 
          for i in range(1,self.NX+1):
            p=[float(i),float(j),float(k)]
            if self.__isInCylinder(p[0],p[1],p[2],b[0],b[1],b[2],d[0],d[1],d[2],r,height):
              self.__addnode(i,j,k)

  #return True if (px,py,pz) is within the cylinder whose center base is (qx,qy,qz), radius is r,
  # axis is oriented along (dx,dy,dz) direction and height is h. 
  def __isInCylinder(self,px,py,pz,qx,qy,qz,dx,dy,dz,r,h):
    ppx=px-qx; ppy=py-qy; ppz=pz-qz;
    # do pp curl d
    curlx=ppy*dz-ppz*dy
    curly=ppz*dx-ppx*dz
    curlz=ppx*dy-ppy*dx
    n2=dx*dx+dy*dy+dz*dz
    
    pdot=ppx*dx+ppy*dy+ppz*dz
    
    return (curlx*curlx+curly*curly+curlz*curlz )<=r*r*n2 and pdot >=0. and pdot*pdot <= h*h*n2

  # make the union of all the objects of same spacing AND bounding box
  # by adding their nodes into self.nodes.
  # NOTE: self.mynodes is unchanged
  def add(self,objlist):
    
    for o in objlist:
      if o.spacing != self.spacing or o.NX != self.NX or o.NY != self.NY or o.NZ != self.NZ:
        sys.stderr.write("WARNING: an object has different spacing or bounding box. It'll be ignored!")
        continue
      self.objlist.append(o)      
      self.nodes.update(o.nodes)
      
  def load(self,filebasename,kind):
    tipodacaricare=self.__getType(kind)

    hdr=open(filebasename+'.hdr','r').readlines()
    if int(hdr[2].split()[0]) != self.spacing:
      raise Exception("IncompatibleSpacing")

    campi=hdr[0].split()
    n = self.__normalise([int(campi[0]),int(campi[1]),int(campi[2])])
    if n != [self.NX,self.NY,self.NZ]:
      raise Exception("IncompatibleBoxSize")

    campi=hdr[1].split()
    nnodes=int(campi[tipodacaricare])
    nodes = open(filebasename+'.dat','r').readlines()
    for linea in nodes:
      campi=linea.split()
      if int(campi[3]) !=tipodacaricare: continue
      [i,j,k]=self.__normalise([int(campi[0]),int(campi[1]),int(campi[2])])
      self.__addnode(i,j,k)
      
  def loadFromSTL(self,filebasename,z1,z2):
    import vtk

    reader = vtk.vtkSTLReader()
    reader.SetFileName(filebasename+'.stl')
    
    surfacePoly = vtk.vtkPolyData()
    #if vtk.VTK_MAJOR_VERSION <= 5:
    surfacePoly.SetInput(reader.GetOutput())
    #else:
    #  surfacePoly.SetInputConnection(reader.GetOutputPort())      
    points=vtk.vtkPoints()
    for x in range(1,self.NX+1):
      for y in range(1,self.NY+1):
        for z in range(z1,z2+1):
          points.InsertNextPoint([x,y,z])
    pointsPoly=vtk.vtkPolyData()
    pointsPoly.SetPoints(points);

    selectEnclosed = vtk.vtkSelectEnclosedPoints()
    selectEnclosed.SetInputData(pointsPly)
    selectEnclosed.SetSurfaceData(surfacePoly)
    selectEnclosed.Update()  
    
    for i in points.GetNumberOfPoints():
      if selectEnclosed.IsInside(n): self.__addnode(points.GetPoint(n))
      
      
  def write(self,filebasename,formato="DAT"):
      
    f=open(filebasename+'.'+formato.lower(),'w')
    if formato.upper() == "DAT":

      #create and write header
      h=open(filebasename+'.hdr','w')
      h.write(str(self.NX*self.spacing)+" "+str(self.NY*self.spacing)+" "+str(self.NZ*self.spacing)+"\n")
      nfl,nwall,ninl,nout=0,0,0,0
      for k in self.nodes:
        tipo=self.nodes[k][3].kind
        if tipo==1:
          nfl = nfl + 1
        elif tipo==2:
          nwall = nwall + 1
        elif tipo==3:
          ninl = ninl + 1
        elif tipo==4:
          nout = nout + 1
      h.write("5 "+str(nfl)+" "+str(nwall)+" "+str(ninl)+" "+str(nout)+"\n"+str(self.spacing)+'\n')
      h.close()
      
      #write DAT file    
      n=self.nodes #alias
      keys=n.keys()
      keys.sort() # the following loop must be ordered in x,y,z
      for k in keys:
        m=[]
        for i in [0,1,2]: m.append(str(n[k][i]*self.spacing))
        f.write(m[0]+" "+m[1]+" "+m[2]+" "+str(n[k][3].kind)+"\n")
    
    elif formato.upper() == "IOS":
   
      i=0
      for o in self.objlist:
        if o.kind != 3 and o.kind != 4: continue
        i = i +1

      if i==0: raise Exception("IOSRequiresIONodes")
    
      f.write(str(i)+'\n')
      i=1
      indice={}
      for o in self.objlist:
        if o.kind != 3 and o.kind != 4: continue
        if o.header == "": raise Exception("EmptyHeaderInIONodes")
        indice[o]=str(i)        
        f.write(str(i)+" "+o.header+'\n')
        i = i +1

      for tipo in [3,4]: #first inlet then outlet
        noone=True
        for o in indice:
          if o.kind != tipo: continue
          noone=False
          f.write(str(len(o.mynodes))+'\n')
          keys=list(o.mynodes.keys())
          keys.sort()
          for k in keys:
            n=o.mynodes[k]
            m=[]
            for j in [0,1,2]: m.append(str(n[j]*self.spacing))
            f.write(m[0]+" "+m[1]+" "+m[2]+" "+indice[o]+"\n")
        if noone:
          f.write('0\n')
           
    elif formato.upper() == "XYZ":
      nome=['','F','W','I','O','D']
      n=self.nodes #alias
      f.write(str(len(n))+"\n1\n")
      for k in n:
        m=[]
        for i in [0,1,2]: m.append(str(n[k][i]*self.spacing))
        f.write(nome[n[k][3].kind]+" "+m[0]+" "+m[1]+" "+m[2]+"\n")
    else:
      raise Exception("FileFormatUnknown")
    
    f.close()
    
  def __addnode(self,i,j,k):
    key=self.__roll(i,j,k)
    self.nodes[key]=[i,j,k,self]          
    if self.kind == 3 or self.kind == 4 :
      self.mynodes[key]=[i,j,k,self] #needed to save self nodes into IOS file

  def __roll(self,i,j,k):
    return i+(j-1)*self.NX+(k-1)*self.NX*self.NY
  
  def __normalise(self,v,convert=False):
    v1=[]
    for i in range(len(v)):
      if convert:
        v1.append(float(v[i])/float(self.spacing))
      else:
        v1.append(v[i]/self.spacing)
    return v1
  
  def __getType(self,kind):
    if kind == 'FLUID':
      return 1
    elif kind == 'WALL':
      return 2
    elif kind == 'INLET':
      return 3
    elif kind == 'OUTLET':
      return 4
    elif kind == 'DEAD':
      return 5
    else:
      raise Exception('NodeTypeUnknown')

  def compenetrate2(self,othermesh,inverse=False,dirs=[]):
    if self.spacing != 2*othermesh.spacing and  2*self.spacing != othermesh.spacing:
      raise Exception('Compenetration must involve meshes of double resolution')
    otherIsFiner= (self.spacing> othermesh.spacing)
    newnodes={}
    for q in self.nodes:
      if self.nodes[q][3].kind != 1: continue
      if not self.__isInsature(self.nodes[q][:3],othermesh): continue
      [i0,j0,k0]=self.nodes[q][:3]
      connections=0
      if otherIsFiner:
        for p in range(1,self.NDIR):
          i=i0*2+self.dirx[p]
          j=j0*2+self.diry[p]
          k=k0*2+self.dirz[p]
          if othermesh.__isOutOfBounds([i,j,k]): continue
          if othermesh.__roll(i,j,k) in othermesh.nodes:
            if othermesh.nodes[othermesh.__roll(i,j,k)][3].kind == 1:
              connections = connections + 1           
          
        if othermesh.__roll(i0*2,j0*2,k0*2) in othermesh.nodes:
          if othermesh.nodes[othermesh.__roll(i0*2,j0*2,k0*2)][3].kind != 1:
            connections = 0 #there MUST be a coinciding node           
        else:
          connections = 0 #there MUST be a coinciding node 
      else:
        if (i0%2,j0%2,k0%2)==(0,0,0):
          i=i0/2    
          j=j0/2    
          k=k0/2    
          if othermesh.__roll(i,j,k) in othermesh.nodes:
            if othermesh.nodes[othermesh.__roll(i,j,k)][3].kind == 1:
              connections = connections + 2
        for p in range(1,self.NDIR):
          i=i0+self.dirx[p]
          j=j0+self.diry[p]    
          k=k0+self.dirz[p]    
          if (i%2,j%2,k%2) == (0,0,0):
            i=i/2
            j=j/2
            k=k/2
            if othermesh.__isOutOfBounds([i,j,k]): continue
            if othermesh.__roll(i,j,k) in othermesh.nodes:
              if othermesh.nodes[othermesh.__roll(i,j,k)][3].kind == 1:
                connections = connections + 1           
      if connections<2: #insaturable node, populate its neighboring
        if inverse:
          newnodes[q]=1
        else:
          pdirs= range(1,self.NDIR) if len(dirs)==0 else dirs
          for p in pdirs:
            i=i0+self.dirx[p]      
            j=j0+self.diry[p]      
            k=k0+self.dirz[p]      
            kk=self.__roll(i,j,k)
            if self.__isOutOfBounds([i,j,k]) or kk in self.nodes: continue
            if kk not in newnodes: newnodes[kk]=[i,j,k,self]
            
            
            
    if inverse:
      for q in newnodes: del self.nodes[q]
    else:
      self.nodes.update(newnodes)
    return len(newnodes)   

  def removeQ27(self,coarsemesh):
    if 2*self.spacing != coarsemesh.spacing:
      raise Exception('Lonely node Q27 removing must involve a mesh of half resolution')
    newnodes={}
    for q in self.nodes:
      if self.nodes[q][3].kind != 1: continue
      [i0,j0,k0]=self.nodes[q][:3]
      if i0%2==0 or j0%2==0 or k0%2==0: continue
      if not self.__isInsature([i0,j0,k0],coarsemesh): continue
      connected=True
      for p in range(self.NDIR,self.NDIREXT):
        i=i0+self.dirx[p]      
        j=j0+self.diry[p]      
        k=k0+self.dirz[p]      
        if self.__isOutOfBounds([i,j,k]): continue
        if coarsemesh.__roll(i/2,j/2,k/2) not in coarsemesh.nodes:
          connected=False
          break        
      if connected: newnodes[self.__roll(i0,j0,k0)]=1
    for q in newnodes: del self.nodes[q]
    return len(newnodes)   


  def compenetrate(self,finermesh):
    if self.spacing != 2*finermesh.spacing:
      raise Exception('1st step compenetration must involve a finer mesh')
    newnodes={}
    for q in self.nodes:
      if self.nodes[q][3].kind != 1: continue
      if not self.__isInsature(self.nodes[q][:3],finermesh): continue
      [i0,j0,k0]=self.nodes[q][:3]
      for p in range(1,self.NDIREXT):
        i=i0+self.dirx[p]      
        j=j0+self.diry[p]      
        k=k0+self.dirz[p]      
        if self.__isOutOfBounds([i,j,k]): continue
        kk=self.__roll(i,j,k)
        if kk in self.nodes: continue
        if kk not in newnodes: newnodes[kk]=[i,j,k,self]
    self.nodes.update(newnodes)
    return len(newnodes)   
  
  def coating(self,finer,coarser):
    newnodes={}
    walls=moebiusBox('WALL',self.NX*self.spacing,self.NY*self.spacing,self.NZ*self.spacing,self.spacing)
    for q in self.nodes:
      if self.nodes[q][3].kind != 1: continue
      ii=self.nodes[q][:3]
      if not self.__isInsature(ii,finer): continue
      if coarser:
        if not self.__isInsature(ii,coarser): continue
      for p in range(1,self.NDIR):
        i=ii[0]+self.dirx[p]      
        j=ii[1]+self.diry[p]      
        k=ii[2]+self.dirz[p]      
        if self.__isOutOfBounds([i,j,k]): continue
        kk=self.__roll(i,j,k)
        if kk in self.nodes: continue
        if kk not in newnodes: newnodes[kk]=[i,j,k,walls]
    self.nodes.update(newnodes)
    return len(newnodes)   
        
       
  def __isInsature(self,ii, othermesh):
    [i0,j0,k0]=ii
    if othermesh.spacing != self.spacing:
      connected=False
      if othermesh.spacing > self.spacing:
        for p in range(1,self.NDIREXT):
          i=i0+self.dirx[p]
          i = i + i%2
          j=j0+self.diry[p]      
          j = j + j%2
          k=k0+self.dirz[p]      
          k = k + k%2
          if self.__isOutOfBounds([i,j,k]): continue
          if othermesh.__roll(i/2,j/2,k/2) in othermesh.nodes:
            connected=True
            break
        """
        if not connected:
          for p in range(1,self.NDIREXT):
            i=i0+3*self.dirx[p]
            i = i + i%2
            j=j0+3*self.diry[p]      
            j = j + j%2
            k=k0+3*self.dirz[p]      
            k = k + k%2
            if self.__isOutOfBounds([i,j,k]): continue
            if othermesh.__roll(i/2,j/2,k/2) in othermesh.nodes:
              connected=True
              break
        """
      else:
        for p in range(1,self.NDIREXT):
          i=i0+self.dirx[p]
          j=j0+self.diry[p]      
          k=k0+self.dirz[p]      
          if self.__isOutOfBounds([i,j,k]): continue
          if othermesh.__roll(i*2,j*2,k*2) in othermesh.nodes:
            connected=True
            break
        """
        if not connected:
          for p in range(1,self.NDIREXT):
            i=i0+2*self.dirx[p]
            j=j0+2*self.diry[p]      
            k=k0+2*self.dirz[p]      
            if self.__isOutOfBounds([i,j,k]): continue
            if othermesh.__roll(i*2,j*2,k*2) in othermesh.nodes:
              connected=True
              break
        """
          
          
      if not connected: return False
      
    for p in range(1,self.NDIR):
      i=i0+self.dirx[p]      
      j=j0+self.diry[p]      
      k=k0+self.dirz[p]      
      if self.__isOutOfBounds([i,j,k]): continue
      if self.__roll(i,j,k) not in self.nodes: return True
    return False
  
  def repairFluidHoles(self,othermesh1,othermesh2=False):
    
    newnodes={}
    for q in self.nodes:
      if self.nodes[q][3].kind != 1: continue
      ii=self.nodes[q][:3]
      [i0,j0,k0]=ii
   
      if othermesh2:
        connected= self.__isInsature(ii,othermesh1) or self.__isInsature(ii,othermesh2)
      else:
        connected= self.__isInsature(ii,othermesh1)
      if connected: continue
    
      for p in range(1,self.NDIR):
        i=i0+self.dirx[p]      
        j=j0+self.diry[p]      
        k=k0+self.dirz[p]      
        if self.__isOutOfBounds([i,j,k]): continue
        kk=self.__roll(i,j,k)
        if kk not in self.nodes:
          if kk not in newnodes: newnodes[kk]=[i,j,k,self]
                
    self.nodes.update(newnodes)
    return len(newnodes)   
    
  def repairWalls(self,othermesh1,othermesh2=False):
    
    newnodes={}
    walls=moebiusBox('WALL',self.NX*self.spacing,self.NY*self.spacing,self.NZ*self.spacing,self.spacing)
   
    for q in self.nodes:
      if self.nodes[q][3].kind != 1: continue
      ii=self.nodes[q][:3]
      [i0,j0,k0]=ii
   
      if othermesh2:
        connected= self.__isInsature(ii,othermesh1) or self.__isInsature(ii,othermesh2)
      else:
        connected= self.__isInsature(ii,othermesh1)
      if connected: continue
    
      for p in range(1,self.NDIR):
        i=i0+self.dirx[p]      
        j=j0+self.diry[p]      
        k=k0+self.dirz[p]      
        if self.__isOutOfBounds([i,j,k]): continue
        kk=self.__roll(i,j,k)
        if kk not in self.nodes:
          if kk not in newnodes: newnodes[kk]=[i,j,k,walls]
                
    self.nodes.update(newnodes)
    return len(newnodes)   
    

  def __isOutOfBounds(self,ii):
    [i,j,k]=ii
    return i <1 or i > self.NX or j <1 or j > self.NY or k <1 or k > self.NZ 
