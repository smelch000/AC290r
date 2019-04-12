#!/bin/bash

# for each receptor chain convert : OPEP ===> atom_A.inp.FLEXI, atom_B.inp.FLEXI, ...
#                            plus : OPEP ===> atom_A.inp.EN,    atom_B.inp.EN,    ...
# Needs xxxx.tgz dowloaded from opep.galaxy.ibpc.fr

if [ $# -le 0 ]; then 
  echo "Usage: opep2muphy.sh chain1 chain2 ..."
  exit
fi
#path to converter.f90 executable and scsc.dat file
CONVERTERPATH=./
BOX=300
RESCALEMASS=.true.

#files='t7vrh wypzs qe2x8 izyn1'

# do the Fully Flexible (FLEXI) OPEP -> MUPHY conversion
cnt=0
for file in $*; do
   tar -zxf $file.tgz ichain.dat 
   tar -zxf $file.tgz scale.dat
   F=`tar -ztf $file.tgz|grep list` 
   tar -zxf $file.tgz $F
   mv $F $file.list
   F=`tar -ztf $file.tgz|grep top` 
   tar -zxf $file.tgz $F
   mv $F $file.top
   F=`tar -ztf $file.tgz|grep pdb` 
   tar -zxf $file.tgz $F
   mv $F $file.pdb
   echo '------- PROCESSING FLEXI -------------' $file
   echo '------- PROCESSING FLEXI -------------' $file
   echo '------- PROCESSING FLEXI -------------' $file
   $CONVERTERPATH/converter << EOF
$file
$BOX $BOX $BOX
1.8
.false.
$CONVERTERPATH
.false.
.false.
$RESCALEMASS
EOF
   mv atom.inp atom_$file.inp.FLEXI
done

# do the Elastic Network (EN) OPEP -> MUPHY conversion
for file in $*; do
   tar -zxf $file.tgz ichain.dat 
   tar -zxf $file.tgz scale.dat
   F=`tar -ztf $file.tgz|grep list` 
   tar -zxf $file.tgz $F
   mv $F $file.list
   F=`tar -ztf $file.tgz|grep top` 
   tar -zxf $file.tgz $F
   mv $F $file.top
   F=`tar -ztf $file.tgz|grep pdb` 
   tar -zxf $file.tgz $F
   mv $F $file.pdb
   echo '------- PROCESSING EN -------------' $file
   echo '------- PROCESSING EN -------------' $file
   echo '------- PROCESSING EN -------------' $file
   $CONVERTERPATH/converter << EOF
$file
$BOX $BOX $BOX
1.8
.false.
$CONVERTERPATH
.false.
.true.
.false.
6
EOF
   mv atom.inp atom_$file.inp.EN
done

