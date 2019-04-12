#!/usr/bin/env python

import subprocess,shlex,os,sys
import ConfigParser
from utils.pyutils import my_create_dir, my_copy_file
from ctypes import *


RES = os.getenv('RES')
print 'Mesh Resolution:',RES

exe = os.getenv('MOEBIUS_BIN')+'BACKEND/SCRIPT/TOOLS/generate_mesh'

cfg = ConfigParser.SafeConfigParser()
cfg.read('mesh.conf')

model = cfg.get('general', 'file_model')
print 'model:',model

plaques = cfg.get('general', 'file_plaques')
print 'plaques:',plaques

if plaques == 'None':
    plaques=None

num_nodes = cfg.get('general', 'num_nodes')

ios = ''
for i in range(int(num_nodes)):

    io = cfg.get('node_'+str(i), 'file_name')

    if io:
        if len(ios)==0:
            ios += '-i '+io+' -b pressure '
        else:
            ios += '-o '+io+' -B flow '

tgtdir = 'mesh_'+RES+'0voxcm'

if not os.path.isdir(tgtdir):
    my_create_dir(tgtdir)

if plaques:
    # if not my_create_dir(tgtdir): print 'Cannot create directory:', tgtdir
    lumdir = tgtdir+'/preplaques'
    plqdir = tgtdir+'/plaques'
    if not my_create_dir( lumdir ): print 'Cannot create lum directory :',lumdir
    if not my_create_dir( plqdir ): print 'Cannot create plq directory :',plqdir

    cmd = exe + ' -f '+model+' -s '+RES+' -m preplaques_bgkflag '+ios
    subprocess.call( shlex.split(cmd) )
    subprocess.call( shlex.split('mv preplaques_bgkflag.hdr '+lumdir) )
    subprocess.call( shlex.split('mv preplaques_bgkflag.dat '+lumdir) )
    subprocess.call( shlex.split('mv preplaques_bgkflag.ios '+lumdir) )
    subprocess.call( shlex.split('cp mesh_transform.inp '+lumdir) )

    subprocess.call( shlex.split('mv mesh_transform.inp '+tgtdir) )

    # cmd = exe + ' -f '+plaques+' -s '+RES+' -y '+tgtdir+'/mesh_transform.inp'+' -m plaques_bgkflag '
    cmd = exe + ' -f '+plaques+' -s '+RES+' -y '+tgtdir+'/mesh_transform.inp'
    subprocess.call( shlex.split(cmd) )
    subprocess.call( shlex.split('mv bgkflag.hdr '+plqdir+'/plaques.hdr') )
    subprocess.call( shlex.split('mv bgkflag.dat '+plqdir+'/plaques.dat') )
    subprocess.call( shlex.split('mv bgkflag.ios '+plqdir+'/plaques.ios') )


    Lib = CDLL(os.getenv("MOEBIUS_BIN")+"/BACKEND/SCRIPTS/TOOLS/combine_lumen_plaques.so")
            
    name1 = os.path.join(lumdir,'preplaques_bgkflag')
    name2 = os.path.join(plqdir,'plaques')
    name3 = os.path.join(tgtdir,'bgkflag')

    print 'Combining lumen and plaque meshes...'

    Lib.combine(c_int(len(name1)), name1, c_int(len(name2)), name2, c_int(len(name3)), name3)

else:
    # if not my_create_dir(tgtdir): print 'Cannot create directory:', tgtdir
    cmd = exe + ' -f '+model+' -s '+RES+' -m bgkflag '+ios
    subprocess.call( shlex.split(cmd) )

    subprocess.call( shlex.split('cp mesh_transform.inp '+tgtdir) )
    subprocess.call( shlex.split('mv bgkflag.hdr '+tgtdir) )
    subprocess.call( shlex.split('mv bgkflag.dat '+tgtdir) )
    subprocess.call( shlex.split('mv bgkflag.ios '+tgtdir) )


print '...done!'

    

