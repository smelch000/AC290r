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
	print >> sys.stderr, 'Usage:',sys.argv[0],'-g stl_file [-t volume_threshold] [-o out_file]'
	print >> sys.stderr, 'Options:'
	print >> sys.stderr, '\t-g stl_file'
	print >> sys.stderr, '\t--geometry stl_file'
	print >> sys.stderr, '\t\tSpecifies the STL file to filter.'
	print >> sys.stderr, ''
	print >> sys.stderr, '\t-t volume_threshold'
	print >> sys.stderr, '\t--threshold volume_threshold'
	print >> sys.stderr, '\t\tSpecifies the volume threshold (between 0.0 and 1.0) below which connected components will be removed (default: 0.0).'
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

def dVol(v1, v2, v3, n):

        bz = (v1[2] + v2[2] + v3[2])/3.0;

        v1proj = (v2[0] - v1[0], v2[1] - v1[1], 0.0)
        v2proj = (v3[0] - v1[0], v3[1] - v1[1], 0.0)

        vprod = vec_prod(v1proj, v2proj)
        tri_proj_area = math.fabs(0.5*(vprod[2])) * bz;

        return math.copysign(tri_proj_area, n[2]);

def find_first_key(f, dic, notFound):

	ret = notFound
	for key, val in dic.iteritems():
		if not f(key, val): continue
		ret = key
		break
	return ret


def print_triangle(tid):

	print '\ntriangle id:', seed
	eIds = id2triangles[seed]
	print 'edge ids:', eIds
	point_set = tuple(set([pid for eid in eIds for pid in id2edges[eid]]))
	print 'point ids:', point_set
	print 'point coordinates:'
	for pid in point_set:
		print '\t', id2points[pid]


STLFNAME = ''
THRESHOLD = 0.0
OUTFNAME = ''

opts, args = getopt.getopt(sys.argv[1:], 'g:t:o:', ['geometry=','threshold=','output='])
for o, a in opts:
        if o in ('-g', '--geometry'):
		STLFNAME = a
        elif o in ('-t', '--threshold'):
		THRESHOLD = float(a)
		if not (THRESHOLD >= 0.0 and THRESHOLD <= 1.0):
			usage()
			sys.exit(1)
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

print >> sys.stderr, 'Constructing triangle adjacency list...',
graph_time = time.time()

# construct the adjacency list of triangle graph
triGraph = dict((tid, []) for tid in id2triangles)
for eid, triIds in edges2tris.iteritems():
	if len(triIds) == 1:
		# edge', eid, 'belongs to a single triangle'
		pass
	else:
		couples = [(t1, t2) for t1 in triIds for t2 in triIds if t1 != t2]
		for tid1, tid2 in couples:
			triGraph[tid1].append(tid2)

print >> sys.stderr, 'done', time.time()-graph_time, 'secs'
del(edges2tris)

# perform a Breadth-first search to iteratively remove triangles 
# with less than 3 neighboring triangles from the adjacency list
print >> sys.stderr, 'Removing loosely connected triangles...',
bfs_time = time.time()
queue = [];
while True:
	if len(queue) == 0:
		seed = find_first_key(lambda k,v: len(v) < 3, triGraph, -1)
		if seed == -1: break
		queue = [seed]

	# take curr triangle to remove
	cur = queue.pop(0)
	for neigh in triGraph[cur]:
		# remove curr from the adjacency list of its neighbors
		triGraph[neigh].remove(cur)
		if len(triGraph[neigh]) < 3:
			# if the degree of neighbor dropped below 2 AND
			# it's not already in the queue...
			if neigh not in queue: queue.append(neigh)
	# remove curr from the queue
	triGraph.pop(cur)

print >> sys.stderr, 'done', time.time()-bfs_time, 'secs'

# perform a Breadth-first search to:
# - find all the connected components and compute their volume
# - associate every tringle to the connected component it belongs to
component = dict()
tri2comp = dict((tid, -1) for tid in triGraph)

print >> sys.stderr, 'Analyzing connected components...',
bfs_time = time.time()

queue = []; compId = -1
while True:

	if len(queue) == 0:
		seed = find_first_key(lambda k,v: v == -1, tri2comp, -1)
		if seed == -1: break
		
		compId += 1
		queue = [seed]
		tri2comp[seed] = compId

		triPointIds = list(set(id2edges[id2triangles[seed][0]] +
				       id2edges[id2triangles[seed][1]] +
				       id2edges[id2triangles[seed][2]]))
		component[compId] = dVol(id2points[triPointIds[0]],
					 id2points[triPointIds[1]],
					 id2points[triPointIds[2]],
					 id2normals[seed])
	
	cur = queue.pop(0)
	for neigh in triGraph[cur]:
		if tri2comp[neigh] != -1: continue # already visited
		tri2comp[neigh] = compId

		triPointIds = list(set(id2edges[id2triangles[neigh][0]] +
				       id2edges[id2triangles[neigh][1]] +
				       id2edges[id2triangles[neigh][2]]))
		component[compId] += dVol(id2points[triPointIds[0]],
					  id2points[triPointIds[1]],
					  id2points[triPointIds[2]],
					  id2normals[neigh])
		queue.append(neigh)

#print >> sys.stderr, 'Found', len(component), 'connected components:'
#for compId, vol in component.iteritems():
#	print >> sys.stderr, '\t', compId, ' -> volume:', vol

# mark filtered-out components
if len(component.values()) == 0:
    sys.exit(0)

volCut = math.fabs(max(component.values())) * THRESHOLD
for compId, vol in component.iteritems():
	if math.fabs(vol) < volCut:
		component[compId] = 'cutted'

print >> sys.stderr, 'done', time.time()-bfs_time, 'secs (found',\
		     len(component), 'connected components,',\
		     component.values().count('cutted'), 'filtered out)'

print >> sys.stderr, 'Writing STL output file...',
write_time = time.time()
print >> outf, 'solid ascii'
for tid in triGraph:

	if component[tri2comp[tid]] == 'cutted': continue

	edges = id2triangles[tid]
	p1Id, p2Id, p3Id = list(set(id2edges[edges[0]] + id2edges[edges[1]] + id2edges[edges[2]]))

	p1 = id2points[p1Id]
	p2 = id2points[p2Id]
	p3 = id2points[p3Id]

	# compute normal using vertices order, directed as (p1-p0)X(p2-p0)
	side1 = (p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2])
        side2 = (p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2])

	# ensure triangle vertices are written in counter-clockwise order with
	# respect to triangle normal
	if (dot_prod(vec_prod(side1, side2), id2normals[tid]) < 0.0):
		tmp = p1
		p1 = p2
		p2 = tmp

	print >> outf, ' facet normal {0} {1} {2}'.format(id2normals[tid][0], id2normals[tid][1], id2normals[tid][2])
	print >> outf, '  outer loop'
	print >> outf, '   vertex {0} {1} {2}'.format(p1[0], p1[1], p1[2])
	print >> outf, '   vertex {0} {1} {2}'.format(p2[0], p2[1], p2[2])
	print >> outf, '   vertex {0} {1} {2}'.format(p3[0], p3[1], p3[2])
	print >> outf, '  endloop'
	print >> outf, ' endfacet'

print >> outf, 'endsolid'
outf.close()
print >> sys.stderr, 'done', time.time()-write_time, 'secs'
print >> sys.stderr, 'Total time:', time.time()-start_time, 'secs'
sys.exit(0)
