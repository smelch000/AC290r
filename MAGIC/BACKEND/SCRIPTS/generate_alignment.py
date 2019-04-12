#!/home/miocchi/anaconda2/bin/python
#6/04/18
# Get the ATOMs from the <pdbfile>.pdb pertaining to the given"
# chain <chain> and write to <outfile>.pdb file the chain atoms, to <outfile>.ali"
# the complete and uncomplete sequences and to <outfile>.txt infos about missing loops."
# by P. Miocchi


from modeller import *
import textwrap
import sys
if (__name__ == "__main__"):
  if len(sys.argv) < 4:
      print "Usage: generate_alignment.py pdbfile outfile chain"
      print "\tGet the ATOMs from the <pdbfile>.pdb pertaining to the given"
      print "\tchain <chain> and write to <outfile>.pdb file the chain atoms, to <outfile>.ali"
      print "\tthe complete and uncomplete sequences and to <outfile>.txt infos about missing loops."
      print "\tWARNING: the sequence to be constructed must be called <outfile>_fill."
      sys.exit(0)


  chain=sys.argv[3]
  amin={}
  amin["ALA"]="A"
  amin["ARG"]="R"
  amin["ASN"]="N"
  amin["ASP"]="D"
  amin["CYS"]="C"
  amin["GLU"]="E"
  amin["GLN"]="Q"
  amin["GLY"]="G"
  amin["HIS"]="H"
  amin["ILE"]="I"
  amin["LEU"]="L"
  amin["LYS"]="K"
  amin["MET"]="M"
  amin["PHE"]="F"
  amin["PRO"]="P"
  amin["SER"]="S"
  amin["THR"]="T"
  amin["TRP"]="W"
  amin["TYR"]="Y"
  amin["VAL"]="V"


  env = environ()
  env.io.atom_files_directory = "./"

  # Create a new empty alignment and model:
  aln = alignment(env)
  mdl = model(env)

  # Read the whole atom file
  mdl.read(file=sys.argv[1], model_segment=('FIRST:'+sys.argv[3], 'END:'+sys.argv[3]))

  # Add the model sequence to the alignment
  aln.append_model(mdl, align_codes=sys.argv[1], atom_files=sys.argv[1])

  aln.write(file=sys.argv[2]+'.ali')
  #read the first three lines from the ali (PIR format) file produced
  #by modeller
  ali=open(sys.argv[2]+'.ali','r')
  alifirstlines=ali.readline()+ali.readline()+ali.readline()
  ali.close()
  alifirstlines=alifirstlines.replace(sys.argv[1],sys.argv[2])

  pdb = open(sys.argv[1]+".pdb","r").readlines()

  #find all the missing sequences reading from 6 lines
  #beyond "MISSING RESIDUES"
  i=0
  for linea in pdb:
    if linea.find("MISSING RESIDUES")>-1: break
    i=i+1

  i=i+6
  #build the list of resid index that are KNOWN as missed
  missed={}
  while True:
    riga=pdb[i].split()
    i=i+1
    if len(riga)<5: break
    #missed=missed+amin(riga[2])
    if riga[3] == chain:
      missed[int(riga[4])]=riga[2]
    
  #build the sequence of detected aminoacids for the specified chain
  #and write them into outpdbchain.pdb file
  detected={}
  previous=""
  outpdb=open(sys.argv[2]+".pdb","w")
  for linea in pdb:
    riga=linea.split()
    if riga[0] != "ATOM": continue
    if riga[4] != chain: continue
    if riga[5] != previous: detected[int(riga[5])]=riga[3]
    previous=riga[5]
    outpdb.write(linea)

  outpdb.close()
  
  if not detected:
    print "Error: no chain "+chain+"'s atom found!"
    sys.exit(1)
    
  print ">>> Chain "+chain+" ATOMS written on file "+sys.argv[2]+".pdb\n"

  #find widest seq. index range
  posmin=min(min(detected.keys()),min(missed.keys()))
  posmax=max(max(detected.keys()),max(missed.keys()))

    
  #print alignments to std output and to outalifile.ali file
  #and build the sequence of resid that are missing both in the ATOM list
  #and in the MISSING RESIDUES list. This is a rare but not impossible situation.
  diff=""
  sequenza=""
  reallymissed={}

  j=1
  miss=True
  missedresid_from=[]
  missedresid_to=[]
  for i in range(posmin,posmax+1):
    if i in detected:
      miss=True
      diff=diff+amin[detected[i]]
      sequenza=sequenza+amin[detected[i]]
    elif i in missed:
      if miss:
        missedresid_from.append(j)
        missedresid_to.append(j)
        miss=False
        
      missedresid_to[-1]=j
      diff=diff+'-'
      sequenza=sequenza+amin[missed[i]]
    else:
      miss=True
      reallymissed[i]=1

    j += 1
      
  print "CHAIN: ",chain,"\n"
  print "TOTALLY MISSED= ",len(reallymissed),"\n"
  out=textwrap.wrap(diff+"*",75,break_on_hyphens=False,break_long_words=True)
  ali=open(sys.argv[2]+'.ali','w')
  ali.write(alifirstlines)
  for linea in out:
    print linea
    ali.write(linea+"\n")

  ali.write(">P1;"+sys.argv[2]+"_fill"+"\nsequence:::::::::\n") 

  print "\n________COMPLETE SEQUENCE________________\n"

  out=textwrap.wrap(sequenza+"*",75,break_on_hyphens=False,break_long_words=True)
  for linea in out:
    print linea
    ali.write(linea+"\n")
      
  ali.close()
  print ">>> Alignment written on file "+sys.argv[2]+".ali\n"

  #write infos to <outloops> file useful to optloop.py to model and refine the missing loops  
  loops=open(sys.argv[2]+'.txt','w')
  print "\n________MISSING LOOPS________________\n"
  loops.write("#MISSED LOOPS= \n"+str(len(missedresid_from))+"\n")
  loops.write("#MISSED RESID= \n"+str(len(missed))+"\n")
  loops.write("#TOTALLY MISSED RESID= \n"+str(len(reallymissed))+"\n")
  print "OFFSET= ",posmin-1
  loops.write("#OFFSET= \n"+str(posmin-1)+"\n")
  loops.write("#from to (in original numbering)\n")
  for i in range(0,len(missedresid_from)):
    loops.write(str(missedresid_from[i]+posmin-1)+" "+str(missedresid_to[i]+posmin-1)+"\n")
    print missedresid_from[i]+posmin-1,missedresid_to[i]+posmin-1
    
  loops.close()
  print ">>> Missed loops infos written on file "+sys.argv[2]+".txt\n"
  
