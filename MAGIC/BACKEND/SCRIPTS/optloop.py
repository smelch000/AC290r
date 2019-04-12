#!/usr/bin/python
#6/04/18

from modeller import *
from modeller.automodel import *    # Load the automodel class
import sys
import os

if (__name__ == "__main__"):
  if len(sys.argv) < 2:
      print "Usage: optloop.py pdbfile "
      print "\tGet the ATOMs from the <pdbfile>.pdb, the alignment from <pdbfile>.ali file,"
      print "\tthe missing loops infos from <pdbfile>.txt (both created by generate_alignment.py)"
      print "\tand write to <pdbfile>_fill_*.pdb file the complete modeled and refined sequence."
      sys.exit(0)

  log.verbose()
  env = environ()
  # directories for input atom files
  env.io.atom_files_directory = ['.','../NMDA']


  class MyLoopModel(loopmodel):
      
    def select_loop_atoms(self):
      #return selection(self)
      return selection(self.residue_range(RESID_FROM,RESID_TO))

    def special_restraints(self, aln):
      rsr = self.restraints
      at = self.atoms

    #CHAIN = "A"
    #for CHAIN in ["A","B","C","D"]:
  RESID_FROM=""
  RESID_TO=""



  start=1

  a = automodel(env, alnfile = sys.argv[1]+'.ali',
                knowns = sys.argv[1],sequence = sys.argv[1]+"_fill")
  #              knowns = sys.argv[1]+'.pdb', sequence = sys.argv[2])
  a.starting_model= start
  a.ending_model  = start+9

  a.very_fast()
  a.make()

  minmodpdf=1e10
  for model in a.outputs:
    if model['failure'] is not None: continue
    if model['molpdf'] <minmodpdf:
      minmodpdf=model['molpdf']
      best_model=model

  #OFFSET= the very first resid number in the template pdb file
  #(even if it is a miising one) MINUS ONE.
  #This because in the resulting models residues are re-numbered from 1

  inloops=open(sys.argv[1]+'.txt',"r").readlines()
  for i in range(0,len(inloops)):
    if inloops[i].find("#OFFSET")>-1:
      OFFSET=int(inloops[i+1])
      break
    
  res_from=[]
  res_to=[]
  found=False
  for linea in inloops:
    if not found:
      if linea.find("#from to")==-1:
        continue
      else:  
        found =True
        continue
    res_from.append(int(linea.split()[0]))
    res_to.append(int(linea.split()[1]))
  """  
  if CHAIN is "A":
    OFFSET=22
    res_from=[23,831,53,802,95 ,657, 442,617,583]
    res_to  =[24,847,57,808,102,663,449 ,622,604]
  elif CHAIN is "B":
    OFFSET=26
    res_from=[27,842,440,803,540,570,616]
    #res_to  =[29,852,450,820,550,601,629]
    res_to  =[29,846,450,820,550,601,629]
  elif CHAIN is "C":
    OFFSET=22
    res_from=[23,834,95 ,617, 442,617,583]
    #res_to  =[31,852,451,806,549,601,629]
    res_to  =[31,846,451,806,549,601,629]
  elif CHAIN is "D":
    OFFSET=26
    res_from=[27,839,440,803,541,568,615]
    res_to  =[31,852,451,806,549,601,629]
  """
  start=11

  print "\n\n============== OPTIMIZATION STARTED ==================== \n\n"

  #numbering offset is applied because here the starting pdb file is a MODEL result (coming
  #from the previous very_fast automodel)
  for i in range(0,len(res_to)):
    RESID_FROM=str(res_from[i]-OFFSET)
    RESID_TO=str(res_to[i]-OFFSET)

    m = MyLoopModel(env, inimodel = best_model['name'],sequence = sys.argv[1]+"_fill")

    m.loop.starting_model = start
    m.loop.ending_model   = start+9
    m.loop.md_level       = refine.very_slow

    m.make()
    
    minmodpdf=1e10

    for model in m.loop.outputs:
      if model['failure'] is not None: continue
      if model['molpdf'] <minmodpdf:
        minmodpdf=model['molpdf']
        best_model=model
        
    start=start+10
    
    print "\n\n>>>>>>  MODEL "+best_model['name']+" DONE <<<<<<<<<\n\n"

  print "\n\n>>>>>>  END: BEST FINAL MODEL= "+best_model['name']+" <<<<<<<<<\n\n"
  
  os.rename(best_model['name'], best_model['name'][:-4]+"_BEST.pdb")
