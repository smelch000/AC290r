#!/usr/bin/env python
#
# This script filters the data in an STL file by removing the
# connected components whose volume is less than a given threshold
# with respect of the component with the maximum volume.
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
	print >> sys.stderr, 'Usage:',sys.argv[0],'-g stl_file',
	print >> sys.stderr, 'Options:'
	print >> sys.stderr, '\t-g stl_file'
	print >> sys.stderr, '\t--geometry stl_file'
	print >> sys.stderr, '\t\tSpecifies the STL file to filter.'
	print >> sys.stderr, ''
	print >> sys.stderr, '\t-o out_file'
        print >> sys.stderr, '\t--output out_file'
        print >> sys.stderr, '\t\tSpecifies the output file.'
        print >> sys.stderr, '\t\tIf this option is not set the standard output is used.'

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

def norm_tri(a,b,c):

	norm = vec_prod( (b[0]-a[0], b[1]-a[1], b[2]-a[2]),
			 (c[0]-a[0], c[1]-a[1], c[2]-a[2]) )
	l = len_vect(norm)
        print 'V1',a,'V2',b,'V3',c,'norm:',norm,'len_vect(norm):',l
	return tuple([c/l for c in norm])

def dSur(v1, v2, v3):

	v12 = (v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2])
	v13 = (v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2])

	vprod = vec_prod(v12, v13)
	return len_vect(vprod)/2.0

def dVol(v1, v2, v3, n):

        bz = (v1[2] + v2[2] + v3[2])/3.0;

        v1proj = (v2[0] - v1[0], v2[1] - v1[1], 0.0)
        v2proj = (v3[0] - v1[0], v3[1] - v1[1], 0.0)

        vprod = vec_prod(v1proj, v2proj)
        tri_proj_area = math.fabs(0.5*vprod[2])*bz;

        return tri_proj_area*math.copysign(1.0, n[2]);

STLFNAME = ''
OUTFNAME = ''

opts, args = getopt.getopt(sys.argv[1:], 'g:o:', ['geometry=','output='])
for o, a in opts:
        if o in ('-g', '--geometry'):
		STLFNAME = a
	elif o in ('-o', '--output'):
		OUTFNAME = a
	else:
		usage()
		sys.exit(1)

if (STLFNAME == ''):
	usage()
	sys.exit(1)

try:
	stlf = open(STLFNAME, 'r')
except IOError:
	print >> sys.stderr, 'Cannot find file', STLFNAME
	sys.exit(1)
	
if (OUTFNAME == ''):
	outf = sys.stdout
else:
	try:
		outf = open(OUTFNAME, 'w')
	except IOError:
		print >> sys.stderr, 'Cannot write file', OUTFNAME
		sys.exit(1)

# temp dicts
points2id = dict()
edges2id = dict()
triangles2id = dict()
# permanent dicts
edges2tris = dict()
id2normals = dict()

dupTris = 0

print >> sys.stderr, 'Reading STL file...',
start_time = time.time()
# Here we need readline so we cannot use file iterations...
line = readneline(stlf) #stlf.readline()
while True:
	line = readneline(stlf) #stlf.readline() # facet normal
	fields = tuple(line.split())

	if (fields[0] == 'endsolid'):
		break
	
	norm = (float(fields[2]), float(fields[3]), float(fields[4]))
	if len_vect(norm) > 1.0E-6:
		norm = tuple([c/len_vect(norm) for c in norm])

	line = readneline(stlf) #stlf.readline() # reads 'outer loop ...'
	line = readneline(stlf) #stlf.readline() # vertex 1
	fields = tuple(line.split())
	x = float(fields[1])
	y = float(fields[2])
	z = float(fields[3])
	p1 = (x,y,z)

	line = readneline(stlf) #stlf.readline() # vertex 2
	fields = tuple(line.split())
	x = float(fields[1])
	y = float(fields[2])
	z = float(fields[3])
	p2 = (x,y,z)
	
	line = readneline(stlf) #stlf.readline() # vertex 3
	fields = tuple(line.split())
	x = float(fields[1])
	y = float(fields[2])
	z = float(fields[3])
	p3 = (x,y,z)

	# zero-length normals are recomputed based on verices order
	mynorm = norm_tri(p1, p2, p3)
	if (math.fabs(1.0-dot_prod(norm, mynorm)) > 1.0E-3):
		norm = mynorm

	if p1 not in points2id:
		p1Id = len(points2id)
		points2id[p1] = p1Id
	else:
		p1Id = points2id[p1]

	if p2 not in points2id:
		p2Id = len(points2id)
		points2id[p2] = p2Id
	else:
		p2Id = points2id[p2]

	if p3 not in points2id:
		p3Id = len(points2id)
		points2id[p3] = p3Id
	else:
		p3Id = points2id[p3]

	# ignore zero-area triangles
	if p1Id == p2Id or p1Id == p3Id or p2Id == p3Id:
		line = readneline(stlf) #stlf.readline() # endloop
		line = readneline(stlf) #stlf.readline() # endfacet
		continue

	# build the 3 edges (pA, pB) s.t. pA < pB
	e1 = (min(p1Id,p2Id), max(p1Id,p2Id))
	e2 = (min(p1Id,p3Id), max(p1Id,p3Id))
	e3 = (min(p2Id,p3Id), max(p2Id,p3Id))

	if not e1 in edges2id:
		e1Id = len(edges2id)
		edges2id[e1] = e1Id
	else:
		e1Id = edges2id[e1]

	if not e2 in edges2id:
		e2Id = len(edges2id)
		edges2id[e2] = e2Id
	else:
		e2Id = edges2id[e2]

	if not e3 in edges2id:
		e3Id = len(edges2id)
		edges2id[e3] = e3Id
	else:
		e3Id = edges2id[e3]

	tri = tuple(sorted([e1Id, e2Id, e3Id]))

	# discards duplicated trinagles
	if not tri in triangles2id:
		triId = len(triangles2id)
		triangles2id[tri] = triId
		
		id2normals[triId] = norm
	
		if not e1Id in edges2tris:
			edges2tris[e1Id] = []
		edges2tris[e1Id].append(triId)
		
		if not e2Id in edges2tris:
			edges2tris[e2Id] = []
		edges2tris[e2Id].append(triId)
		
		if not e3Id in edges2tris:
			edges2tris[e3Id] = []
		edges2tris[e3Id].append(triId)
	else:
		dupTris += 1

	line = readneline(stlf) #stlf.readline() # endloop
	line = readneline(stlf) #stlf.readline() # endfacet
stlf.close()

# remove dictionaries no longer needed and invert the others
id2points = dict((val, key) for key, val in points2id.iteritems())
id2edges = dict((val, key) for key, val in edges2id.iteritems())
id2triangles = dict((val, key) for key, val in triangles2id.iteritems())

del(points2id)
del(edges2id)
del(triangles2id)

print >> sys.stderr, 'done', time.time()-start_time, 'secs (found',\
		     len(id2triangles), 'triangles,', len(id2edges), 'edges and', len(id2points), 'points)' 
print >> sys.stderr, 'Removed', dupTris, 'duplicated triangles'

volume = 0.0
surface = 0.0
minX = minY = minZ = sys.float_info.max
maxX = maxY = maxZ = sys.float_info.min
print >> sys.stderr, 'Computing stats...',
stats_time = time.time()
for tid, t in id2triangles.items():

	triPointIds = list(set(id2edges[t[0]] + id2edges[t[1]] + id2edges[t[2]]))

	v1 = id2points[triPointIds[0]]
	v2 = id2points[triPointIds[1]]
	v3 = id2points[triPointIds[2]]

	surface += dSur(v1, v2, v3)
	volume += dVol(v1, v2, v3, id2normals[tid])

	minX = min(v1[0], minX); minX = min(v2[0], minX); minX = min(v3[0], minX)
	maxX = max(v1[0], maxX); maxX = max(v2[0], maxX); maxX = max(v3[0], maxX)

	minY = min(v1[1], minY); minY = min(v2[1], minY); minY = min(v3[1], minY)
	maxY = max(v1[1], maxY); maxY = max(v2[1], maxY); maxY = max(v3[1], maxY)

	minZ = min(v1[2], minZ); minZ = min(v2[2], minZ); minZ = min(v3[2], minZ)
	maxZ = max(v1[2], maxZ); maxZ = max(v2[2], maxZ); maxZ = max(v3[2], maxZ)

print >> sys.stderr, 'done', time.time()-stats_time, 'secs'

if OUTFNAME:
	print >> sys.stderr, 'Writing output file...',
	write_time = time.time()

print >> outf, 'geometry_surface:', surface
print >> outf, 'geometry_volume:', volume
print >> outf, 'bounding_box_volume:', (maxX-minX)*(maxY-minY)*(maxZ-minZ)
outf.close()
if OUTFNAME:
	print >> sys.stderr, 'done', time.time()-write_time, 'secs'

print >> sys.stderr, 'Total time:', time.time()-start_time, 'secs'
sys.exit(0)
