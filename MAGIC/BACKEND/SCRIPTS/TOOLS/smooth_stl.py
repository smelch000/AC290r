#!/usr/bin/env python
#
# This script smooths the geometry in an STL file by applying the
# lamda-mu algorithm. Optionally, the geometry may be interpolated
# before applying the smoothing.
#
# Written by Mauro Bisson (mauro.bis@gmail.com)
#
import sys
import getopt
import math
import time

def dist2(n, m):
	return ((n[0]-m[0])*(n[0]-m[0]) + (n[1]-m[1])*(n[1]-m[1]) + (n[2]-m[2])*(n[2]-m[2]))

def usage():
	print >> sys.stderr, "Usage:",sys.argv[0],"-g stl_file",
	print >> sys.stderr, "Options:"
	print >> sys.stderr, "\t-g stl_file"
	print >> sys.stderr, "\t--geometry stl_file"
	print >> sys.stderr, "\t\tSpecifies the STL file to smooth."
	print >> sys.stderr, ""
	print >> sys.stderr, "\t-i"
	print >> sys.stderr, "\t--interpolate"
	print >> sys.stderr, "\t\tPerforms triangle interpolation on input mesh (default: no)."
	print >> sys.stderr, ""
	print >> sys.stderr, "\t-s scale_factor"
	print >> sys.stderr, "\t--scale scale_factor"
	print >> sys.stderr, "\t\tSpecifies the Laplacian smoothing scale factor (default: 0.3)."
	print >> sys.stderr, ""
	print >> sys.stderr, "\t-n itarations"
	print >> sys.stderr, "\t--numiter itarations"
	print >> sys.stderr, "\t\tSpecifies the number of times Laplacian smoothing has to be applied (default: 5)."
	print >> sys.stderr, ""
	print >> sys.stderr, "\t-o out_file"
        print >> sys.stderr, "\t--output out_file"
        print >> sys.stderr, "\t\tSpecifies the output file."
        print >> sys.stderr, "\t\tIf this option is not set the standard output is used."

def readneline(f):
	l = f.readline()
	while not l.strip():
		l = f.readline()
	return l

def dot_prod(a,b):
	return sum(n*m for n,m in zip(a,b))

def vec_prod(a,b):
	return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def len_vect(a):
	return math.sqrt(sum(n*n for n in a))

def dis_vect(a,b):
	return math.fabs(len_vect(a[0]-b[0], a[1]-b[1], a[2]-b[2]))

STLFNAME = ""
OUTFNAME = ""
INTERPOLATE = 0

# we will use the Taubin lambda/mu smoothing algorithm with an optimal k_PB=0.1
K_PB = 0.1
SCALE_POS = 0.3
SCALE_NEG = SCALE_POS/(K_PB*SCALE_POS-1.0)
ITERATIONS = 10

opts, args = getopt.getopt(sys.argv[1:], "g:s:n:o:i", ["geometry=","scale=","numiter=","output=","interpolate"])
for o, a in opts:
        if o in ("-g", "--geometry"):
		STLFNAME = a
        elif o in ("-s", "--scale"):
		SCALE_POS = float(a)
		SCALE_NEG = SCALE_POS/(K_PB*SCALE_POS-1.0)
        elif o in ("-n", "--numiter"):
		ITERATIONS = int(a)
        elif o in ("-i", "--interpolate"):
		INTERPOLATE = 1
	elif o in ("-o", "--output"):
		OUTFNAME = a
	else:
		usage()
		sys.exit(1)

#print >> sys.stderr, SCALE_POS ,SCALE_NEG

if (STLFNAME == ""):
	usage()
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

vertex2id = dict()
id2vertex = dict()
triangles = set()
count = 0

print >> sys.stderr, "Reading STL file...",
start_time = time.time()
# Here we need readline so we cannot use file iterations...
line = readneline(stlf) #stlf.readline()
while True:
	line = readneline(stlf) #stlf.readline() # facet normal
	fields = tuple(line.split())

	if (fields[0] == "endsolid"):
		break
	
	n = (float(fields[2]), float(fields[3]), float(fields[4]))
	
	line = readneline(stlf) #stlf.readline() # reads "outer loop ..."
	
	line = readneline(stlf) #stlf.readline() # vertex 1
	fields = tuple(line.split())
	x = float(fields[1])
	y = float(fields[2])
	z = float(fields[3])
	p1 = (x,y,z)
	if (not p1 in vertex2id) :
		vertex2id[p1] = count
		id2vertex[count] = p1
		count = count+1

	line = readneline(stlf) #stlf.readline() # vertex 2
	fields = tuple(line.split())
	x = float(fields[1])
	y = float(fields[2])
	z = float(fields[3])
	p2 = (x,y,z)
	if (not p2 in vertex2id) :
		vertex2id[p2] = count
		id2vertex[count] = p2
		count = count+1
	
	line = readneline(stlf) #stlf.readline() # vertex 3
	fields = tuple(line.split())
	x = float(fields[1])
	y = float(fields[2])
	z = float(fields[3])
	p3 = (x,y,z)
	if (not p3 in vertex2id) :
		vertex2id[p3] = count
		id2vertex[count] = p3
		count = count+1

	if INTERPOLATE == 1:
		# splice each triangle in 4 by using middle points in each side
		p12 = ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2, (p1[2]+p2[2])/2)
		if (not p12 in vertex2id) :
			vertex2id[p12] = count
			id2vertex[count] = p12
			count = count+1
	
		p23 = ((p2[0]+p3[0])/2, (p2[1]+p3[1])/2, (p2[2]+p3[2])/2)
		if (not p23 in vertex2id) :
			vertex2id[p23] = count
			id2vertex[count] = p23
			count = count+1
		
		p31 = ((p1[0]+p3[0])/2, (p1[1]+p3[1])/2, (p1[2]+p3[2])/2)
		if (not p31 in vertex2id) :
			vertex2id[p31] = count
			id2vertex[count] = p31
			count = count+1

	# insert trinagles (p1,p2,p3) such that (p2-p1)X(p3-p1) is directed
	# as the normal n
	side1 = (p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])
	side2 = (p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2])
	if (dot_prod(vec_prod(side1, side2), n) > 0.0):
		if INTERPOLATE == 0:
			triangles.add((vertex2id[p1], vertex2id[p2], vertex2id[p3]))
		else:
			triangles.add(( vertex2id[p1], vertex2id[p12], vertex2id[p31]))
			triangles.add((vertex2id[p12],  vertex2id[p2], vertex2id[p23]))
			triangles.add((vertex2id[p12], vertex2id[p23], vertex2id[p31]))
			triangles.add((vertex2id[p23],  vertex2id[p3], vertex2id[p31]))
	else:
		if INTERPOLATE == 0:
			triangles.add((vertex2id[p2], vertex2id[p1], vertex2id[p3]))
		else:
			triangles.add(( vertex2id[p12],  vertex2id[p1], vertex2id[p31]))
			triangles.add((  vertex2id[p2], vertex2id[p12], vertex2id[p23]))
			triangles.add(( vertex2id[p23], vertex2id[p12], vertex2id[p31]))
			triangles.add((  vertex2id[p3], vertex2id[p23], vertex2id[p31]))


	line = readneline(stlf) #stlf.readline() # endloop
	line = readneline(stlf) #stlf.readline() # endfacet
stlf.close()
print >> sys.stderr, "done", time.time()-start_time, "secs"
print >> sys.stderr, "Using", len(triangles), "triangles", count, "vertices."

# we no longer need vertex2id...
del(vertex2id)

degree = dict([(n,0) for n in range(count)])
stream = dict()
print >> sys.stderr, "Generating geometry stream...",
stream_time = time.time()
# for every vertex v do 
#	stream[v] = set{w | {v,w} in EDGES and v < w}
#	degree[v] = |{w | {v,w} in EDGES}|
# NB: verices connected only to neighbors with higer ids won't appear in stream!
for t in triangles:

	# t = (a, b, c)
	tmp = list(t)
	tmp.sort()

	if not tmp[0] in stream:
		stream[tmp[0]] = set()
	if not tmp[1] in stream[tmp[0]]:
		stream[tmp[0]].add(tmp[1])
		degree[tmp[0]] = degree[tmp[0]] + 1
		degree[tmp[1]] = degree[tmp[1]] + 1
	if not tmp[2] in stream[tmp[0]]:
		stream[tmp[0]].add(tmp[2])
		degree[tmp[0]] = degree[tmp[0]] + 1
                degree[tmp[2]] = degree[tmp[2]] + 1

	if not tmp[1] in stream:
		stream[tmp[1]] = set()
	if not tmp[2] in stream[tmp[1]]:
	        stream[tmp[1]].add(tmp[2])
		degree[tmp[1]] = degree[tmp[1]] + 1
                degree[tmp[2]] = degree[tmp[2]] + 1

print >> sys.stderr, "done", time.time()-stream_time,"secs"

print >> sys.stderr, "Smoothing..."
delta = dict([(id,[0.0,0.0,0.0]) for id in range(count)])
smooth_time = time.time()
for it in range(ITERATIONS*2):

	if it%2 == 0: 
		it_time = time.time()
	for id in range(count):

		# vertices with id higer than those of all neighs
		# don't appear in stream
		if not id in stream:
			continue
		
		v = id2vertex[id]
		for nid in stream[id]:
	
			nv = id2vertex[nid]
			diff = (nv[0]-v[0], nv[1]-v[1], nv[2]-v[2])

			delta[id][0] = delta[id][0] + diff[0]
			delta[id][1] = delta[id][1] + diff[1]
			delta[id][2] = delta[id][2] + diff[2]
		
			delta[nid][0] = delta[nid][0] - diff[0]
			delta[nid][1] = delta[nid][1] - diff[1]
			delta[nid][2] = delta[nid][2] - diff[2]

	for id in range(count):

		v = id2vertex[id]

		delta[id][0] = delta[id][0] / degree[id]
		delta[id][1] = delta[id][1] / degree[id]
		delta[id][2] = delta[id][2] / degree[id]

		if it%2 == 0:
			scale = SCALE_POS
		else:
			scale = SCALE_NEG

		id2vertex[id] = (v[0]+delta[id][0]*scale, v[1]+delta[id][1]*scale, v[2]+delta[id][2]*scale)
		delta[id] = [0.0, 0.0, 0.0]
	
	if it%2 == 1:
		print >> sys.stderr, "\titeration", it/2, "running time:",time.time()-it_time,"secs"

print >> sys.stderr, "...done", time.time()-smooth_time, "secs"

print >> sys.stderr, "Writing STL output file...",
write_time = time.time()
print >> outf, 'solid ascii'
for t in triangles:
	v0 = id2vertex[t[0]]
	v1 = id2vertex[t[1]]
	v2 = id2vertex[t[2]]

	# compute normal using vertices order, directed as (v1-v0)X(v2-v0)
	side1 = (v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2])
        side2 = (v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2])
	n = vec_prod(side1, side2)

	n_len = len_vect(n)

	if n_len == 0.0: continue

	print >> outf, ' facet normal {0} {1} {2}'.format(n[0]/n_len, n[1]/n_len, n[2]/n_len)
	print >> outf, '  outer loop'
	print >> outf, '   vertex {0} {1} {2}'.format(v0[0], v0[1], v0[2])
	print >> outf, '   vertex {0} {1} {2}'.format(v1[0], v1[1], v1[2])
	print >> outf, '   vertex {0} {1} {2}'.format(v2[0], v2[1], v2[2])
	print >> outf, '  endloop'
	print >> outf, ' endfacet'

print >> outf, 'endsolid'
outf.close()
print >> sys.stderr, "done", time.time()-write_time, "secs"
print >> sys.stderr, "Total time:", time.time()-start_time, "secs"
sys.exit(0)
