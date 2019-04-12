#!/usr/bin/env python
#
# This script couples the geometry in an STL file with the shear 
# stress data produced by the MOEBIUS simulator (for that geometry).
# The output is gives as a VTK file.
#
# Written by Mauro Bisson (mauro.bis@gmail.com)
#
import sys
import time
import getopt

def dist2(n, m):
	return ((n[0]-m[0])*(n[0]-m[0]) + (n[1]-m[1])*(n[1]-m[1]) + (n[2]-m[2])*(n[2]-m[2]))

def usage():
	print >> sys.stderr, "Usage:",sys.argv[0],"-g stl_file -m mesh_file -o out_file",
	print >> sys.stderr, "[-r transform_file] [-s scale_factor] [-t translation_vector]"
	print >> sys.stderr, "Options:"
	print >> sys.stderr, "\t-g stl_file"
	print >> sys.stderr, "\t--geometry stl_file"
	print >> sys.stderr, "\t\tSpecifies the STL file."
	print >> sys.stderr, ""
	print >> sys.stderr, "\t-m mesh_file"
	print >> sys.stderr, "\t--mesh mesh_file"
	print >> sys.stderr, "\t\tSpecifies the surface mesh file containing the shear stress values."
	print >> sys.stderr, "\t\tIf mesh_file is '-' then shear stress values are read from stdin."
	print >> sys.stderr, ""
	print >> sys.stderr, "\t-o out_file"
	print >> sys.stderr, "\t--output out_file"
	print >> sys.stderr, "\t\tSpecifies the output file."
	print >> sys.stderr, "\t\tIf this option is not set the standard output is used."
	print >> sys.stderr, ""
	print >> sys.stderr, "\t-r transform_file"
	print >> sys.stderr, "\t--readtransf transform_file"
	print >> sys.stderr, "\t\tSpecifies the file containing the transformation applied to the mesh."
	print >> sys.stderr, ""
	print >> sys.stderr, "\t-s scale_factor"
	print >> sys.stderr, "\t--scale scale_factor"
	print >> sys.stderr, "\t\tSpecifies the scale factor applied to the STL geometry."
	print >> sys.stderr, "\t\tIf not specified the scale factor defaults to 1.0."
	print >> sys.stderr, "\t\tIf -r option is also used, the scale factor specified in the"
	print >> sys.stderr, "\t\ttransform_file takes precedence."
	print >> sys.stderr, ""
	print >> sys.stderr, "\t-t translation_vector"
	print >> sys.stderr, "\t--transl translation_vector"
	print >> sys.stderr, "\t\tSpecifies the translation applied to the STL geometry. The"
	print >> sys.stderr, "\t\ttranslation vector is specified as 'Va,Vb,Vc' (without quotes)."
	print >> sys.stderr, "\t\tIf not specified the translation vector defaults to (0.0,0.0,0.0)."
	print >> sys.stderr, "\t\tIf -r option is also used, the translation vector specified in the"
	print >> sys.stderr, "\t\ttransform_file takes precedence."

def readneline(f):
	l = f.readline()
	while not l.strip():
		l = f.readline()
	return l

MSTEP = 2
TRANSL = [0.0, 0.0, 0.0]
SCALE = 1.0
MESHFNAME = ""
STLFNAME = ""
TRANSFNAME = ""
OUTFNAME = ""

opts, args = getopt.getopt(sys.argv[1:], "g:m:s:t:o:r:", ["geometry=","mesh=","scale=","transl=","output=","readtransf="])
for o, a in opts:
        if o in ("-g", "--geometry"):
		STLFNAME = a
	elif o in ("-m", "--mesh"):
		if (a == "-"):
			MESHFNAME = "/dev/stdin"
		else:
			MESHFNAME = a
	elif o in ("-s", "--scale"):
		SCALE = float(a)
	elif o in ("-t", "--transl"):
		TRANSL = [float(x) for x in a.split(',')]
	elif o in ("-r", "--readtransf"):
		TRANSFNAME = a
	elif o in ("-o", "--output"):
		OUTFNAME = a
	else:
		usage()
		sys.exit(1)

if (MESHFNAME == "" or STLFNAME == ""):
	usage()
	sys.exit(1)

try:
	meshf = open(MESHFNAME,'r')
except IOError:
	print >> sys.stderr, "Cannot find file", MESHFNAME
	sys.exit(1)

try:
	stlf = open(STLFNAME, 'r')
except IOError:
	print >> sys.stderr, "Cannot find file", STLFNAME
	sys.exit(1)
	
if (OUTFNAME == ""):
	outf = sys.stdout
else:
	try:
		outf = open(OUTFNAME, 'w')
	except IOError:
		print >> sys.stderr, "Cannot write file", OUTFNAME
		sys.exit(1)

if (TRANSFNAME != ""):
	try:
		transff = open(TRANSFNAME, 'r')
	except IOError:
		print >> sys.stderr, "Cannot find file", TRANSFNAME
		sys.exit(1)
	for line in transff:
		fields = tuple(line.split())
		if (fields[0] == "scale_factor:"):
			SCALE = float(fields[1])
		elif (fields[0] == "translation_vector:"):
			TRANSL = [float(x) for x in fields[1].split(',')]
	transff.close()

mesh_stress = dict()
mesh_pressure = dict()

mesh_stress_min = mesh_pressure_min = sys.float_info.max
mesh_stress_max = mesh_pressure_max = sys.float_info.min

print >> sys.stderr, "Reading SHEAR STRESS and PRESSURE file...",
start_time = time.time()
for line in meshf:
        fields = tuple(line.split())
	p = (int(fields[0]), int(fields[1]), int(fields[2]))
	mesh_stress[p] = float(fields[3])
	mesh_pressure[p] = float(fields[4])

	mesh_stress_min = min(mesh_stress_min, mesh_stress[p])
	mesh_stress_max = max(mesh_stress_max, mesh_stress[p])
	
	mesh_pressure_min = min(mesh_pressure_min, mesh_pressure[p])
	mesh_pressure_max = max(mesh_pressure_max, mesh_pressure[p])

meshf.close()
print >> sys.stderr, "done", time.time()-start_time, "secs"

#print >> sys.stderr, 'pressure raange:', mesh_pressure_min, mesh_pressure_max
#print >> sys.stderr, 'stress raange:', mesh_stress_min, mesh_stress_max

stl_points = dict()
stl_points_rev = dict()
stl_tris = set()
count = 0

# Here we need readline so we cannot use file iterations...
print >> sys.stderr, "Reading STL file...",
stlread_time = time.time()
line = readneline(stlf) #stlf.readline()
while True:
	line = readneline(stlf) #stlf.readline() # facet normal
	fields = tuple(line.split())
	if (fields[0] == "endsolid"):
		break
	#fields = tuple(line.split())
	# ignore for now
	#p = (float(fields[2]), float(fields[3]), float(fields[4]))
	#print 'NORMAL {0} {1} {2}'.format(p[0], p[1], p[2])
	
	line = readneline(stlf) #stlf.readline() # reads "outer loop ..."
	
	line = readneline(stlf) #stlf.readline() # vertex 1
	fields = tuple(line.split())
	x = float(fields[1])*SCALE+TRANSL[0]
	y = float(fields[2])*SCALE+TRANSL[1]
	z = float(fields[3])*SCALE+TRANSL[2]
	p1 = (x,y,z)
	if (not p1 in stl_points) :
		stl_points[p1] = count
		stl_points_rev[count] = p1
		count = count+1

	line = readneline(stlf) #stlf.readline() # vertex 2
	fields = tuple(line.split())
	x = float(fields[1])*SCALE+TRANSL[0]
	y = float(fields[2])*SCALE+TRANSL[1]
	z = float(fields[3])*SCALE+TRANSL[2]
	p2 = (x,y,z)
	if (not p2 in stl_points) :
		stl_points[p2] = count
		stl_points_rev[count] = p2
		count = count+1
	
	line = readneline(stlf) #stlf.readline() # vertex 3
	fields = tuple(line.split())
	x = float(fields[1])*SCALE+TRANSL[0]
	y = float(fields[2])*SCALE+TRANSL[1]
	z = float(fields[3])*SCALE+TRANSL[2]
	p3 = (x,y,z)
	if (not p3 in stl_points) :
		stl_points[p3] = count
		stl_points_rev[count] = p3
		count = count+1

	stl_tris.add((stl_points[p1], stl_points[p2], stl_points[p3]))
	line = readneline(stlf) #stlf.readline() # endloop
	line = readneline(stlf) #stlf.readline() # endfacet
stlf.close()
print >> sys.stderr, "done", time.time()-stlread_time, "secs"

print >> sys.stderr, "Merging shear stress and model data...",
merge_time = time.time()
stl_stress = dict()
stl_pressure = dict()
rng = range(-MSTEP,MSTEP+1)
rng = [(x,y,z) for x in rng for y in rng for z in rng]
for n in range(count):
	stl_p = stl_points_rev[n]
	istl_p = (int(round(stl_p[0])), int(round(stl_p[1])), int(round(stl_p[2])))
	stl_stress[n] = 0.0
	stl_pressure[n] = 0.0
#	d = dist2((0,0,0),(MSTEP,MSTEP,MSTEP))
	
	st_ncount = 0
	pr_ncount = 0
	
	for (i,j,k) in rng:
		neigh = (istl_p[0]+i, istl_p[1]+j, istl_p[2]+k)
		
		if (not neigh in mesh_stress): continue
		stl_stress[n] = stl_stress[n] + mesh_stress[neigh]
		st_ncount = st_ncount + 1

		if (not neigh in mesh_pressure): continue
		stl_pressure[n] = stl_pressure[n] + mesh_pressure[neigh]
		pr_ncount = pr_ncount + 1

	# we set min value for triangle points with no mesh-point neighbors
	if (st_ncount != 0):
		stl_stress[n] = stl_stress[n] / st_ncount
	else:
		stl_stress[n] = mesh_stress_min
	
	if (pr_ncount != 0):
		stl_pressure[n] = stl_pressure[n] / pr_ncount
	else:
		stl_pressure[n] = mesh_pressure_min

print >> sys.stderr, "done", time.time()-merge_time, "secs"

print >> sys.stderr, "Writing VTK output file...",
vtkwrite_time = time.time()
print >> outf, '# vtk DataFile Version 3.0'
print >> outf, 'output sstress Mauro 1.0'
print >> outf, 'ASCII'
print >> outf, 'DATASET UNSTRUCTURED_GRID'
print >> outf, 'POINTS', count, "float"

for index in range(count):
	print >> outf, '{0} {1} {2}'.format(stl_points_rev[index][0], stl_points_rev[index][1], stl_points_rev[index][2])

print >> outf, ''
print >> outf, 'CELLS', len(stl_tris), len(stl_tris)*4
for t in stl_tris:
	print >> outf, '3 {0} {1} {2}'.format(t[0], t[1], t[2])

print >> outf, ''
print >> outf, 'CELL_TYPES', len(stl_tris)
for t in stl_tris:
	print >> outf, '5' # VTK_TRIANGLE

print >> outf, ''
print >> outf, 'POINT_DATA', count
print >> outf, 'SCALARS ESS[Pa] float'
print >> outf, 'LOOKUP_TABLE default'
for index in range(count):
	print >> outf, '{0:.4}'.format(stl_stress[index])

print >> outf, ''
print >> outf, 'SCALARS PRESSURE[Pa] float'
print >> outf, 'LOOKUP_TABLE default'
for index in range(count):
	print >> outf, '{0:.4}'.format(stl_pressure[index])

outf.close()
print >> sys.stderr, "done", time.time()-vtkwrite_time, "secs"
print >> sys.stderr, "Total time:", time.time()-start_time, "secs"
sys.exit(0)
