#!/usr/bin/env python
# 13-Dec-2018
# P.M.

""" usage example:

from atomhandler import *

a1=AtomHandler("atom1.inp")
a1.loadAtoms()
a2=AtomHandler("atom2.inp")
a2.loadAtoms()
#database to resolve lacking vdw interactions
database=[
"X Y 12-6   19551.4      15.2404      0.00000      0.00000      10.0000",    
"X Z 12-6   98854.7      34.2694      0.00000      0.00000      10.0000",    
"X XX 12-6   260913.      58.2708      0.00000      0.00000      10.0000"    
]
a1.add(a2,"new",True,database)
#or
#a1.add(a2,"new",True, open('mydatabase','r').readlines() )
a1.replaceSection("MODEL","molecular",["MODEL\n","Azzzz"])
a1.translate(10.,10.,10)
a1.setFileName("atom_new.inp")
a1.saveAtoms()

"""

import sys
import math
class AtomHandler:
  PDBATOMX=30
  PDBATOMY=38
  PDBATOMZ=46
  ##construct the AtomHandler object associated to a given inp file
  # @param[in] fname The associated inp filename (extension included).
  #              It can be changed afterward with setFileName
  # @param[in] commentchars The string containig one or more comment chars
  def __init__(self,fname,commentchars='!'):
    self.commentchars=commentchars
    self.lines=[]
    self.setFileName(fname)
    
  ##change the filename for i/o operations
  # @param[in] fname The associated inp filename (extension included).
  def setFileName(self,fname):
    self.filename=fname
  ## Change the string of comment chars
  #
  # This takes effect only at a subsequent reading with loadAtoms
  # @param[in] commentchars The string containig one or more comment chars
  def setCommentchars(self,commentchars):
    self.commentchars=commentchars

  def __dist(self,p1,p2):
      return (pow(p1[0]-p2[0],2)+pow(p1[1]-p2[1],2)+pow(p1[2]-p2[2],2))

  ## Load all infos from the associated inp file
  #
  # WARNING: lines containing '\n' chars (they will be treated as SINGLE lines).
  # should be avoided!
  # WARNING: scalex,scaley,scalez are only used when loading STL files
  
  def loadAtoms(self,scalex=1.,scaley=1.,scalez=1.):
    
    if self.filename[-4:].upper()==".STL":
      linee=open(self.filename,'r').readlines()
      it=iter(linee)
      atoms={}
      while True:
        linea=self.__jumpTo("vertex",it)
        if not linea: break
        coord=linea.split()
        atoms[(float(coord[1])*scalex,float(coord[2])*scaley,float(coord[3])*scalez)]=1
        for i in [1,2]: #just to repeat twice
          linea=it.next()
          coord=linea.split()
          if coord[0] != "vertex":
            sys.stderr.write("Error: STL file "+ self.filename+" corrupted\n")
            sys.exit(1)
          atoms[(float(coord[1])*scalex,float(coord[2])*scaley,float(coord[3])*scalez)]=1
        
      atomlist=atoms.keys()
      
      triangles=[]
      it=iter(linee)
      maxdist=0.
      while True:
        linea=self.__jumpTo("vertex",it)
        if not linea: break
        v=[]
        coord=linea.split()
        v.append(atomlist.index((float(coord[1])*scalex,float(coord[2])*scaley,float(coord[3])*scalez))+1)

        for i in [1,2]: #just to repeat twice
          linea=it.next()
          coord=linea.split()
          v.append(atomlist.index((float(coord[1])*scalex,float(coord[2])*scaley,float(coord[3])*scalez))+1)

        triangles.append(v)
        
      self.lines=self.lines+["MODEL","# title","molecular types    1","_Nodes","nummols    1"]
      self.lines.append("atoms "+str(len(atomlist)))
      for i in range(0,len(atomlist)):
        self.lines.append("P  1.  0.  1  1  1. 1. 1. 1   ! name, mass, chge, rept, frzt, inertia x/y/z, frzr")
      self.lines.append("angles "+str(len(triangles)))
      for i in range(0,len(triangles)):
        t=triangles[i]
        self.lines.append("bcon "+str(t[0])+" "+str(t[1])+" "+" "+str(t[2])+" 1. 45.")
        p1=[atomlist[t[0]-1][0],atomlist[t[0]-1][1],atomlist[t[0]-1][2]]
        p2=[atomlist[t[1]-1][0],atomlist[t[1]-1][1],atomlist[t[1]-1][2]]
        p3=[atomlist[t[2]-1][0],atomlist[t[2]-1][1],atomlist[t[2]-1][2]]
        maxdist=max(maxdist,self.__dist(p1,p2),self.__dist(p1,p3),self.__dist(p2,p3))
      self.lines=self.lines+["hydro",\
      ".true. 0.0   1.0  1.0   1.0  0.0      ! passive scalar, tumbling coeff, visc_enhancer, csi, oblate, smooth",\
      "0.010000   0.0  ! gammas","finish","vdw    1",\
      "P P lj  0   2.0  1.0  1.0  "+str(math.sqrt(maxdist)*1.01),"CLOSEMODEL"]
      self.lines=self.lines+["CONFIG",str(len(atomlist))]
      for i in range(0,len(atomlist)):
        a=atomlist[i]
        self.lines.append(str(a[0])+" "+str(a[1])+" "+str(a[2])+" 0 0 0 "+str(i+1))
      self.lines.append("CLOSECONFIG")
        
      #self.filename=self.filename[:-4]+".inp"
    else:
    
      self.lines=self.cleanLines(open(self.filename,'r').readlines())
       
      
  ## Save all data to the associated inp, xyz or pdb file
  #
  # WARNING: there must not be empty lines after keyword 'finish'. The first one
  # is assumed to contain the molecule name.
  # @param[in] formato Can be INP|XYZ in the XYZ case it saves atom names and coordinates only
  def saveAtoms(self,formato="INP"):
    f=open(self.filename,'w')
    if formato.upper()=="INP":
      for linea in self.lines:
          f.write(linea+"\n")
          """
          if moreReadable:
            l=linea.split()[0].upper()
            if 'CLOSE' in l or l in ['MOLECULAR']: f.write('\n')
          """
    elif formato.upper()=="XYZ":
      confs,order= self.getMolecularTypeNamesAndConf()
      confs_keys=sorted(confs,key=order.__getitem__)
      n=0
      for mol in confs_keys:
        n += len(confs[mol])
                 
      f.write(str(n)+"\n"+" 1\n")
      
      i=1
      for mol in confs_keys:
        for c in confs[mol]:
          cc=c.split()
          #f.write(" ".join([mol,cc[0],cc[1],cc[2],str(i+1)])+"\n")
          f.write(" ".join([mol,cc[0],cc[1],cc[2]])+"\n")
          i += 1 

    elif formato.upper()=="PDB":
      confs,order= self.getMolecularTypeNamesAndConf()
      confs_keys=sorted(confs,key=order.__getitem__)
      i=1
      for mol in confs_keys:
        atoms=self.getMoleculeAtomList(mol)
        n=len(atoms)
        if n<1: continue
        for j,c in enumerate(confs[mol]):
          cc=c.split()
          #f.write(" ".join([mol,cc[0],cc[1],cc[2],str(i+1)])+"\n")
          f.write('%6s%5d %5s             %8.3f%8.3f%8.3f\n' % ("ATOM  ", i,atoms[j%n],float(cc[0]),float(cc[1]),float(cc[2])))
          i += 1 


    f.close()
    
  ## Advance the iterator/enumerator <it> to the 1st line with <stringa> in the first field
  #
  # @returns the line it jumped to (and the line number if <it> is enumerate(iter()) object )
  # @param[in] stringa The string to be matched at the first field (CASE INSENSITIVE)
  # @it[in,out] The iterator (or enumerate(iterator) )  that has to be advanced 
  def __jumpTo(self,stringa,it):
    s=stringa.upper()
    if isinstance(it,enumerate):
        while True:
            try:
                i,ll=it.next()
            except StopIteration:
                #sys.stderr.write("Warning: "+stringa+" not found in "+self.filename+"\n")
                return -1,""
            if ll.split()[0].upper() == s: break
          
        return i,ll    
    else:  
        while True:
            try:
                ll=it.next()
            except StopIteration:
                #sys.stderr.write("Warning: "+stringa+" not found in "+self.filename+"\n")
                return ""
            if ll.split()[0].upper() == s: break
          
        return ll    

  ## Extract and returns contiguous lines from the INP infos.
  #
  # Extract and returns contiguous lines from the INP file previously loaded
  # with LoadAtoms.
  # NOTICE: Blank line or line starting with a char in self.commentchars are ignored
  # in this loading.
  # The lines extracted are in form of a list of strings that are COPIED from self.lines.
  # The lines go from that starting with <fromstring> (included) up to the one starting with
  # <tostring> (included only if <tostringIncluded>=True) or to the end if the latter is not given.
  # If <readNumLines> is given as True then at most N lines are returned with N
  # taken from the second field of the line starting with <fromstring>
  #
  # @param[in] fromstring The string to be matched in the 1st field of the 1st line extracted.
  # @param[in] tostring The string to be matched in the 1st field of the line just following the last
  # one extracted (or matching the 1st field of the last line extracted if <tostringIncluded>=True.
  # @param[in] tostringIncluded If True the last line extracted has <tostring> in the 1st field.
  # @param[in] readNumLines If True N lines are read with N being read at the 2nd field of the 1st
  # line extracted.
  
  def extractSection(self,fromstring,tostring="",readNumLines=False,tostringIncluded=False,occurrence=1):
    linee=self.lines
    if not linee: return []
    il=iter(linee)

    ##jump to the <occurrence> occurrence of <fromstring>
    nocc=0
    while nocc<occurrence: 
      ll=self.__jumpTo(fromstring,il)
      ll=ll.split()
      if not ll: return []
      nocc=nocc+1
      
    if readNumLines and len(ll)>1:
        n=int(ll[1])
    else:
        n=len(linee)
    out=[" ".join(ll)] ## include fromstring line
    n += 1 ##because fromstring is included
    while len(out) < n:
        ll=il.next()
        if ll.split()[0].upper() == tostring.upper():
          if tostringIncluded: out.append(ll)
          break
        out.append(ll)

    return out
  
  ## Replace lines from the <occurence>-th <fromstring> to <tostring> (both excluded)
  #
  # Remove all lines enclosed between the <occurence>-th occurrence of <fromstring>,
  # and the following <tostring> (both excluded) or to the end if the latter is empty.
  # Replace them with <sectionbody>.
  # @returns success Returns False (and no repace is done) if either <fromstring>
  # or <tostring> was not found   
  def replaceSection(self,fromstring,tostring,sectionbody,occurrence=1):
    linee=self.lines
    if not linee: return False
    il=enumerate(iter(linee))

    ##jump to the <occurrence> occurrence of <fromstring>
    for k in range(0,occurrence):
      i,ll=self.__jumpTo(fromstring,il)
    if not ll:
      sys.stderr.write("Warning: no replace done. "+fromstring+" not found in "+self.filename+"\n")
      return False
    newlines=linee[:(i+1)]+sectionbody
    if tostring:
      ##jump to the 1st occurrence of <tostring>
      j,ll=self.__jumpTo(tostring,il)
      if not ll:
        sys.stderr.write("Warning: no replace done. "+tostring+" not found in "+self.filename+"\n")
        return False
      newlines=newlines+linee[j:]

    self.lines=newlines
    return True

  ## Return the {names of molecular types: CONFIG} dict, {order in which they appear in self.filename} dict.
  #
  # Return the tuple with the names of all molecular types as the keys of a dictionary
  # in which the value is a list containing the configuration lines
  # of the corresponding atoms in CONFIG section (it is EMPTY
<<<<<<< HEAD
  # for _Wall 'molecule') and the order in which he molecular
=======
  # for _Wall 'molecule') and the order in which the molecular
>>>>>>> master
  # types are read from self.filename.
  # @returns {molecule name: [atom configuration]},[order of names]
  def getMolecularTypeNamesAndConf(self):
    linee = self.lines
    if not linee: return {},{}
    il=iter(linee)
    ##jump to the 1st occurrence of "molecular"
    ll=self.__jumpTo("molecular",il)
    ll=ll.split()
    if len(ll)<3: return {},{}
    n=int(ll[2]) ##get the no. of molecular types

    order={}
    iorder=1
    conffrom={}
    confto={}
    counter=0
    nummols = 1
    nextIsName=True
    thereIsWall=False
    while len(confto)<n:
        try:
            ll = il.next().split()
        except StopIteration:
            break ##end of file
        if nextIsName:
            ##this now must be the name of the molecule
            key=ll[0]
            isWall= self.__isWall(key) 
            if isWall: thereIsWall = True
            ##store the first line index in atoms configuration for this molecule
            conffrom[key]=counter
            order[key]=iorder
            iorder = iorder + 1
            nextIsName=False ##next line cannot be mol name
            continue
       
        if ll[0].lower() == "nummols":
            ## "nummols" must be written before "atoms"
            nummols=int(ll[1])
        elif ll[0].lower() == "atoms": 
            ##notice that counter is not increased if the mol is a _Wall
            ##because that doesn't have atoms in CONFIG
            if not isWall: counter += int(ll[1])*nummols
            ##store the last line index (plus one) in atoms configuration for the latest found molecule
            confto[key]=counter ##at the end, this is the tot. no. of atoms (apart from W)
        elif ll[0].lower() == "finish":
            nextIsName= True ##next line must contain mol name

    if n > len(conffrom):
      sys.stderr.write("Warning: found less molecules than the given types number in "+self.filename+"\n")
    
    config=self.getConfiguration()
    natoms = len(config)
    #if thereIsWall: counter -= 1 ##because wall atom does not have a configuration
    if natoms != counter:
        sys.stderr.write("Error: the tot. no. of atoms in CONFIG doesn't match that in the MODEL section in "+ self.filename+"\n")
        sys.exit(1)
      
    configurations={}
    for mol in conffrom:
        configurations[mol]=config[conffrom[mol]:confto[mol]]
        
    return configurations,order
  
  #return header section ("MODEL" included).
  def getHead(self):
      return self.extractSection("MODEL","molecular") 

  ## Return atom species' name of a given molecule type
  #
  # @param molName The name of the molecule type to read
  # @returns[out] atomnames i.e. <atomnames>:
  # the atom names in the order they appear in the molecular model.
  def getMoleculeAtomList(self,molName):
      
      model= self.extractSection(molName,"finish")
      il=iter(model)
      #jump to the 1st occurrence of "atoms"
      ll=self.__jumpTo("atoms",il)
      ll=ll.split()
      if not ll or len(ll)<2:
          #sys.stderr.write("Warning: atoms section absent/corrupted in "+molName+" in "+self.filename+"\n")
          return []
    
      n=int(ll[1])
  
      atoms=[]
      i=0
      while i<n:
          try:
              ll=il.next()
          except StopIteration:
              sys.stderr.write("Error: found less atoms than the given atoms number for "+molName+" in "+self.filename+"\n")
              sys.exit(1)
              
          atoms.append(ll.split()[0])
          i += 1
      return atoms
    


  ## Return model section and atom species' name of a given molecule type
  #
  # @param molName The name of the molecule type to read
  # @returns[out] (model,atomnames) i.e. <model>: the model section associated
  # to the <molName> molecule specie name, ('finish' keyword included), <atomnames>:
  # the unordered *set* of atom names.
  def getMoleculeModel(self,molName):
      
      model= self.extractSection(molName,"finish",tostringIncluded=True)
      il=iter(model)
      #jump to the 1st occurrence of "atoms"
      ll=self.__jumpTo("atoms",il)
      ll=ll.split()
      if not ll or len(ll)<2:
          #sys.stderr.write("Warning: atoms section absent/corrupted in "+molName+" in "+self.filename+"\n")
          return [],set()
    
      n=int(ll[1])
  
      atoms=set()
      i=0
      while i<n:
          try:
              ll=il.next()
          except StopIteration:
              sys.stderr.write("Error: found less atoms than the given atoms number for "+molName+" in "+self.filename+"\n")
              sys.exit(1)
              
          atoms.add(ll.split()[0].upper())
          i += 1
      return model,atoms
    
  ##return vdw section (non-bonding interactions).
  def getVdw(self):
      return self.extractSection("vdw","",readNumLines= True)
  ##return anything beyond "CLOSECONFIG" (included)
  def getTail(self):
      return self.extractSection("CLOSECONFIG")
  ##return atoms configuration section
  def getConfiguration(self):
      out=self.extractSection("CONFIG","CLOSECONFIG")
      return out[2:] #remove the first two lines, i.e. "CONFIG" and the tot. atom number

  ## Return "clean" lines
  #
  # Return the list of lines after removing commented (i.e.
  # with a comment char as a first non-blank character) and empty lines.
  # It can be used to "import" external lines (e.g. a vdw database)
  # that need to be cleaned.
  # @param[in] lines If given, the "external" lines to be cleaned.
  # Otherwise self.lines will be cleaned.
  # @returns the list of clean lines
  #
  def cleanLines(self,lines):
      """ too slow!!
      for i,l in enumerate(lines):
          if '\n' in l:
            l=l.split('\n')
            newlines=newlines+l
          else:
            newlines.append(l)
      """      
      ##exclude commented and empty lines
      newlines=filter(lambda l: l and l != "\n" and l.split()[0][0] not in self.commentchars,lines)
      ##remove leading/trailing spaces
      return map(lambda l: l.strip(),newlines)
   
  ## Add 3 fields to the atoms configuration line
  #
  # Add "counter molCounter 'molName'" fields to each atom configuration line (in config[molName])
  # the first counter is just global to the entire section, while molCounter runs within
  # the molecule type whose name is molName.
  # @param[in,out] config The dict such that config[molName] is the list of molName's atoms configuration
  # @param[in] molName The molecule name
  # @param[in] counter The global index of config lines (initially should be 1)
  # @param[in] molCounter The starting molecular index in config lines (initially
  # should be 1 for each molecule type)
  # @returns index,molIndex updated values 
  def __addInfoInConfig(self,config,molName,counter,molCounter):
      jj=counter
      for j,c in enumerate(config[molName]):
        cc=c.split()[0:6] ## take only the 1st "basical" 6 fields
        cc=cc+[str(counter),str(j+molCounter),molName] ##add 3 further fields
        config[molName][j]=" ".join(cc)
        counter += 1
      return counter,counter-jj+molCounter
  
  def __isWall(self,molName): return (molName.upper()=="_WALL" or molName.upper()=="WALL") 

  ## Add and combine another inp infos to the self one
  #
  # @param[in] otherAtoms The reference to another AtomHandler object
  # @param[in] newheader The string giving the header of the resulting inp
  # @param[in] addMolInfoInConfig If True then add "index molIndex 'molName'"
  # fields to each atom configuration line (default: False).
  # @param[in] listDataBaseVdw The list of vdw lines that can be used to specify the 
  # non-bonding interactions between atom species of self and otherAtoms that
  # are not listed into the vdw sections of both of these.
  # If these interactions lacks even in this list, then a ">>>MISSING<<<" keyword will
  # be put into the resulting vdw section.
  
  def add(self,otherAtoms,newheader,addMolInfoInConfig=False,listDataBaseVdw=[]):
    
      configurations,order=self.getMolecularTypeNamesAndConf()
      configurations_keys= sorted(configurations,key=order.__getitem__)

      otherConfigurations,otherOrder=otherAtoms.getMolecularTypeNamesAndConf()
      otherConfigurations_keys= sorted(otherConfigurations,key=otherOrder.__getitem__)
      
      conftot=[]
      i=1
      for mol in configurations_keys:
          if addMolInfoInConfig:
            i,j=self.__addInfoInConfig(configurations,mol,i,1)
          conftot=conftot+configurations[mol]
            
          if mol in otherConfigurations:
              sys.stdout.write("Molecule type "+mol+" found in both "+self.filename+" and "+otherAtoms.filename+"\n")
              if addMolInfoInConfig:
                i,j=self.__addInfoInConfig(otherConfigurations,mol,i,j)
              #append in conftot list the conf list of the other inp
              conftot=conftot+otherConfigurations[mol]
      for mol in otherConfigurations_keys:
          if mol not in configurations:
              if addMolInfoInConfig:
                i,j=self.__addInfoInConfig(otherConfigurations,mol,i,1)
              conftot=conftot+otherConfigurations[mol]

      molecules={}
      atoms=set()
      nummols={}

      for mol in configurations_keys:
          molecules[mol],atoms0=self.getMoleculeModel(mol)
          atoms.update(atoms0) #accumulate atom names of all molecules
          nummols[mol]=int(molecules[mol][1].split()[1]) #get nummols from the "nummols n" line 
          
      for mol in otherConfigurations_keys:
        if mol not in configurations:
          molecules[mol],atoms0=otherAtoms.getMoleculeModel(mol)
          if not molecules[mol]:
            sys.stderr.write("Warning: "+mol+" atoms section absent/corrupted in "+otherAtoms.filename+"\n")
          atoms.update(atoms0) #accumulate atom names of all molecules
          nummols[mol]=int(molecules[mol][1].split()[1]) #get nummols from the "nummols n" line 
        elif not self.__isWall(mol):
          #it's a molecule (not Wall) already in self inp, so update nummols
          nummols[mol] += 1
          molecules[mol][1] = "nummols "+str(nummols[mol])
      ##WARNING: any '\n' causes the line to be splitted (empty lines will be removed)
      newlines=[]

      for mol in configurations_keys:
          #append this mol model if not wall
          if not self.__isWall(mol):
              newlines=newlines+molecules[mol]
      for mol in otherConfigurations_keys:
        if mol not in configurations:
          #append this mol model if not wall
          if not self.__isWall(mol):
              newlines=newlines+molecules[mol]
      ## append _Wall molecule at the last position <<<<<<
      for mol in molecules:
        if self.__isWall(mol):
           newlines=newlines+molecules[mol]
           break
         
      newlines= ["MODEL",newheader,"molecular types "+str(len(molecules))] + newlines

      vdws=self.getVdw()[1:] #exclude "vdw n" line
      othervdws=otherAtoms.getVdw()[1:] #exclude "vdw n" line
      pairs=[]
      pvdw=[]#will contain the vdw line associated to the unordered set {atom1,atom2}
      for vdw in vdws:
        v=vdw.split()
        s={v[0],v[1]}
        if s not in pairs:
          pairs.append(s)
          pvdw.append(vdw)
      #append vdws of the other inp that are not present in self inp
      for vdw in othervdws:
        v=vdw.split()
        s={v[0],v[1]}
        if s not in pairs:
          pairs.append(s)
          pvdw.append(vdw)
          
      #convert the given database (list) of vdws into
      #[(atom1,atom2),...] list and [vdw1,vdw2,...] list
      dataBaseAtoms=[]
      dataBaseVdw=[]
      if listDataBaseVdw:
        listDataBaseVdw=self.cleanLines(listDataBaseVdw)#clean lines
        for vdw in listDataBaseVdw:
          v=vdw.split()
          s={v[0],v[1]}
          if s not in dataBaseAtoms:
            dataBaseAtoms.append(s)
            dataBaseVdw.append(vdw)
          
       
      lackingVdw=" >>>MISSING<<<"
      atoms=list(atoms)
      missing=False
      for atom1 in atoms:
        for atom2 in atoms:
          s={atom1,atom2}
          if s in pairs: continue
          if s not in dataBaseAtoms:
            pairs.append(s)
            pvdw.append(atom1+" "+atom2+" "+lackingVdw)
            missing=True
          else:
            pairs.append(s)
            pvdw.append(dataBaseVdw[dataBaseAtoms.index(s)])
      
      if missing: sys.stdout.write("Warning: there are missing non-bonding pair interactions in "+self.filename+"\n")
      
      newlines.append("vdw "+str(len(pairs)))
      for v in pvdw:
        newlines.append(v)
      newlines = newlines + ["CLOSEMODEL","CONFIG",str(len(conftot))]
      newlines = newlines + conftot + ["CLOSECONFIG"]
      ##add TABULA section, if any.
      tab=self.extractSection("TABULA","CLOSETABULA",tostringIncluded=True)
      if not tab:
        tab=otherAtoms.extractSection("TABULA","CLOSETABULA",tostringIncluded=True)
        
      self.lines=newlines+tab
      
  ## Get atoms bounding box limits
  #
  #@returns atoms bounding box limits
  def getBoundingBox(self):
      conf=self.getConfiguration()
      a=conf[0].split()
      xM=xm=float(a[0])
      yM=ym=float(a[1])
      zM=zm=float(a[2])
      for atom in conf:
          a=atom.split()
          xM,yM,zM=max(xM,float(a[0])),max(yM,float(a[1])),max(zM,float(a[2]))
          xm,ym,zm=min(xm,float(a[0])),min(ym,float(a[1])),min(zm,float(a[2]))
      return (xm,ym,zm,xM,yM,zM)
    
  def translate(self,x,y,z):
      conf=self.getConfiguration()
      #newconf=list(conf)
      for i,atom in enumerate(conf):
          a=atom.split()
          xa,ya,za=float(a[0])+x,float(a[1])+y,float(a[2])+z
          conf[i]=" ".join([str(xa),str(ya),str(za)]+a[3:])
          
      conf = [str(len(conf))]+conf #prepend the no. of atoms
      self.replaceSection("CONFIG","CLOSECONFIG",conf)

  # WARNING: this ONLY scales coords in CONFIG section! Length parameters (e.g. cutoff) are untouched
  def scale(self,sx,sy,sz):
      if sx==1. and sy==1. and sz==1.: return
      conf=self.getConfiguration()
      #newconf=list(conf)
      for i,atom in enumerate(conf):
          a=atom.split()
          xa,ya,za=float(a[0])*sx,float(a[1])*sy,float(a[2])*sz
          conf[i]=" ".join([str(xa),str(ya),str(za)]+a[3:])
          
      conf = [str(len(conf))]+conf #prepend the no. of atoms
      self.replaceSection("CONFIG","CLOSECONFIG",conf)
  ##Translate atom coordinates to match a given pdb file
  #
  # Translate all atoms coordinates such that the 1st atom will match
  # the 1st atom position in PDBfile (if FirstCalpha==False, otherwise find the first CA atom)
  # @param[in] PDBfile The name of the PDB file
  # @param[in] FirstCalpha True if one wants to refer to the 1st CA atom
  # @returns x,y,z as the applied traslation
  def translateToPDBFile(self,PDBfile,FirstCalpha=False):
    il=iter(open(PDBfile,'r').readlines())
    if not FirstCalpha:
      l=self.__jumpTo('ATOM',il) #.split()
<<<<<<< HEAD
      if l and len(l)>47:
        x_o,y_o,z_o= float(l[31:(31+8)]),float(l[39:(39+8)]),float(l[47:(47+8)])
=======
      if l and len(l)>PDBATOMZ:
        x_o,y_o,z_o= float(l[PDBATOMX:(PDBATOMX+8)]),float(l[PDBATOMY:(PDBATOMY+8)]),float(l[PDBATOMZ:(PDBATOMZ+8)])
>>>>>>> master
      else:
        sys.stderr.write("Error: unable to read ATOM coordinates from "+PDBfile+"\n")
        sys.exit(1)
    else:
      l=self.__jumpTo('ATOM',il) #.split()
      while True:
        if l.split()[2]=="CA":
<<<<<<< HEAD
          if l and len(l)>47:
            x_o,y_o,z_o= float(l[31:(31+8)]),float(l[39:(39+8)]),float(l[47:(47+8)])
=======
          if l and len(l)>PDBATOMZ:
            x_o,y_o,z_o= float(l[PDBATOMX:(PDBATOMX+8)]),float(l[PDBATOMY:(PDBATOMY+8)]),float(l[PDBATOMZ:(PDBATOMZ+8)])
>>>>>>> master
          else:
            sys.stderr.write("Error: unable to read ATOM coordinates from "+PDBfile+"\n")
            sys.exit(1)
          break          
        try:
            l=il.next() #.split()
        except StopIteration:
          sys.stderr.write("Error: unable to read ATOM coordinates from "+PDBfile+"\n")
          sys.exit(1)
          
      
    ## get the first atom (must be CAlpha atom if FirstCalpha == True)
    a=self.getConfiguration()[0].split()
    x,y,z=x_o-float(a[0]),y_o-float(a[1]),z_o-float(a[2])
    self.translate(x,y,z)
    sys.stdout.write("Atoms in "+self.filename+" translated by "+str(x)+","+str(y)+","+str(z)+"\n")
    return x,y,z
  
  ##Eliminate distant bonds
  #
  # Eliminate from bonds section all interactions involving atoms whose distance is more than
  # <maxdist>
  # @param[in] maxdist The maximum admitted distance for bonds
  # WARNING: To evaluate distances it considers only the atom coordinates of the
  # FIRST molecule of each kind in CONFIG section
  def pruneBonds(self,maxdist):
    confs,order=self.getMolecularTypeNamesAndConf()
    confs_keys=sorted(confs,key=order.__getitem__)
    nm=1
    il=enumerate(iter(self.lines))      
    for mol in confs_keys:
      bonds=self.extractSection("bonds","angles",nm)
      #natmsPerMol=int((sqrt((len(bonds)-1)*8+1)+1)/2)
      
      newbonds=[]
      for ind in range(1,len(bonds)):
        b=bonds[ind].split()
        if b[0]=="opep" or b[0]=="opem" or b[0]=="12-6":
          i=int(b[1])-1
          c=confs[mol][i].split()
          xi,yi,zi=float(c[0]),float(c[1]),float(c[2])
          j=int(b[2])-1
          c=confs[mol][j].split()
          xj,yj,zj=float(c[0]),float(c[1]),float(c[2])
          if self.__dist([xi,yi,zi],[xj,yj,zj])<maxdist*maxdist:
            newbonds.append(bonds[ind])
        else:
          newbonds.append(bonds[ind])
      self.replaceSection("bonds","angles",newbonds,nm)

      ##jump to the nm-th occurrence of "bonds"
      k,ll=self.__jumpTo("bonds",il)
      self.lines[k]= "bonds "+str(len(newbonds))+" "+" ".join(bonds[0].split()[2:])
      nm = nm + 1
      
    return len(newbonds)
  

  ## Translate bounding box 
  #
  # Translate all atoms coord.s such that the "origin" of the bounding box
  # is at the given point (default: the origin)
  # @param[in] ox,oy,oz The new bounding box "origin" coords.
  # @returns xt,yt,zt the performed translation
  def translateBoundingBoxTo(self,ox=0.,oy=0,oz=0.):
    x0,y0,z0,x,y,z=self.getBoundingBox()
    self.translate(ox-x0,oy-y0,oz-z0)
    return ox-x0,oy-y0,oz-z0
  

  ## Change the name of a molecule
  #
  # Return True if the replacement was done
  # @param[in] oldname The molecule name to be changed
  # @param[in] newname The new molecule name to replace <oldname>
  # @returns True if one name change is done
  # WARNING: only one replacement is done
  def changeMolName(self,oldname,newname):
    linee = self.lines
    if not linee: return False
    il=enumerate(iter(linee))
    
    ##jump to the 1st occurrence of "molecular"
    i,ll=self.__jumpTo("molecular",il)
    ll=ll.split()
    if len(ll)<3: return False
    n=int(ll[2]) ##get the no. of molecular types
    counter=0
    nextIsName=True
    while counter<n:
        try:
            i,ll = il.next()
            ll = ll.split()
        except StopIteration:
            break ##end of file
        if nextIsName:
            ##this now must be the name of the molecule
            key=ll[0]
            if key == oldname:
              self.lines[i]=newname
              return True
            counter = counter +1
            nextIsName = False
            continue
       
        if ll[0].lower() == "finish":
            nextIsName= True ##next line must contain mol name

    return False
  

<<<<<<< HEAD
=======

>>>>>>> master
