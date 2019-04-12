#!/usr/bin/env python

import os
import re
import sys
import time
import math
import getopt

def usage():
    print >> sys.stderr, "Usage:",sys.argv[0],"-g stl_file -c vtk_file -s px,py,pz -e px,py,pz -n numcuts [-m]"
    print >> sys.stderr, "Options:"
    print >> sys.stderr, "\t-g stl_file"
    print >> sys.stderr, "\t--geometry stl_file"
    print >> sys.stderr, "\t\tSpecifies the STL file."
    print >> sys.stderr, ""
    print >> sys.stderr, "\t-c vtk_file"
    print >> sys.stderr, "\t--cline vtk_file"
    print >> sys.stderr, "\t\tSpecifies the VTK file containing a centerline for stl_file"
    print >> sys.stderr, "\t\tspecifed with option -g."
    print >> sys.stderr, ""
    print >> sys.stderr, "\t-s px,py,pz"
    print >> sys.stderr, "\t--startp px,py,pz"
    print >> sys.stderr, "\t\tSpecifies the centerline point from which cross sections are cut"
    print >> sys.stderr, ""
    print >> sys.stderr, "\t-e px,py,pz"
    print >> sys.stderr, "\t--endp px,py,pz"
    print >> sys.stderr, "\t\tSpecifies the centerline point to which cross sections are cut"
    print >> sys.stderr, ""
    print >> sys.stderr, "\t-n numcuts"
    print >> sys.stderr, "\t--number px,py,pz"
    print >> sys.stderr, "\t\tSpecifies the number of cross sections cut between the start and"
    print >> sys.stderr, "\t\tend point"
    print >> sys.stderr, ""
    print >> sys.stderr, "\t-m"
    print >> sys.stderr, "\t--smooth"
    print >> sys.stderr, "\t\tSpecifies whether the centerline in vtk_file must be smoothed"
    print >> sys.stderr, "\t\tbefore cuttin the cross sections."

def dist2(n, m):
    return (n[0]-m[0])**2 + (n[1]-m[1])**2 + (n[2]-m[2])**2

def get_normal(n, m):
    l = math.sqrt(dist2(n,m))
    if math.fabs(l-0.0) < 1.0E-10: return None
    return tuple([(m[0]-n[0])/l,
              (m[1]-n[1])/l,
              (m[2]-n[2])/l])

def vec_prod(a,b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def len_vect(a):
    return math.sqrt(sum(n*n for n in a))

def dSur(v1, v2, v3):

    v12 = (v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2])
    v13 = (v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2])

    vprod = vec_prod(v12, v13)
    return len_vect(vprod)/2.0

def polygon_metrics(points, n):
    
    center = [0.0, 0.0, 0.0]
    for i in range(n):
        center = map(lambda x,y:x+y, center, points.GetPoint(i))
    center = map(lambda x: x/n, center)

    peri = 0.0
    area = 0.0
    for i in range(1, n):
        peri += math.sqrt(math.fabs(dist2(points.GetPoint(i-1), points.GetPoint(i))))
        area += dSur(points.GetPoint(i-1), points.GetPoint(i), center)
        
    peri += math.fabs(dist2(points.GetPoint(n-1), points.GetPoint(0)))
    area += dSur(points.GetPoint(n-1), points.GetPoint(0), center)

    return peri, area

def vect_mean(vlist):
    if not vlist: return None
    l = len(vlist)
    tmp = (0.0, 0.0, 0.0)
    for i in range(l):
        tmp = map(sum, zip(tmp, vlist[i]))

    return map(lambda x: x/l, tmp)

if __name__ == '__main__':

    STLFNAME = ''
    CLINEFNAME = ''
    STAPOINT = None
    ENDPOINT = None
    CUTNUM = 1
    SMOOTH = False
    SMOOTH_WIDTH = 5

    opts, args = getopt.getopt(sys.argv[1:], "g:c:s:e:n:m", ["geometry=","cline=","startp=","endp=","number=","smooth"])

    if not opts:
        usage()
        sys.exit(1)

    for o, a in opts:
        if o in ("-g", "--geometry"):
            STLFNAME = a
        elif o in ("-c", "--cline"):
            CLINEFNAME = a
        elif o in ("-s", "--startp", "-e", "--endp"):
            point = re.split('[ ]*,[ ]*', a.strip())
            point = filter(lambda x:x, point)
            try:
                point = tuple(float(coo) for coo in point)
            except ValueError:
                print >> sys.stderr, 'Bad point specification:', a, '\n'
                usage()
                sys.exit(1)
            if len(point) != 3:
                print >> sys.stderr, 'Bad number of coordinates for point:', a, '\n'
                usage()
                sys.exit(1)
            if o in ("-s", "--startp"):
                if STAPOINT:
                    usage()
                    sys.exit(1)
                STAPOINT = point
            else:
                if ENDPOINT:
                    usage()
                    sys.exit(1)
                ENDPOINT = point
        elif o in ("-n", "--number"):
            CUTNUM = int(a)
        elif o in ("-m", "--smooth"):
            SMOOTH = True
        else:
            usage()
            sys.exit(1)

    if not STLFNAME or not CLINEFNAME:
        print >> sys.stderr, 'Both geometry and centerline file must be specified.'
        usage()
        sys.exit(1)

    if STLFNAME:
        if not os.path.isfile(STLFNAME):
            print >> sys.stderr, 'Cannot find file', STLFNAME
            sys.exit(1)

    if CLINEFNAME:
        if not os.path.isfile(CLINEFNAME):
            print >> sys.stderr, 'Cannot find file', CLINEFNAME
            sys.exit(1)

    if CUTNUM < 1:
        print 'Number of cuts must be > 0!'
        usage()
        sys.exit(1)

    # moved here to pay import latency after parameters checking
    from myvtk import *

    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(CLINEFNAME)
    reader.Update()

    clinePolyData = reader.GetOutput()

    # we only consider the first VTK_POLY_LINE cell in the file
    line = clinePolyData.GetCell(0)
    if line.GetCellType() != vtk.VTK_POLY_LINE:
        print 'VTK_POLY_LINE expected in file', CLINEFNAME
        sys.exit(1)

    points = line.GetPoints()
    nump = line.GetNumberOfPoints()

    if nump < 2:
        print 'Too few points for center line!'
        sys.exit(1)

    # perform laplacian smooth, if requested
    if SMOOTH:
        for l in range(1,SMOOTH_WIDTH):
    
            pbuf = []
            nbuf = [points.GetPoint(i) for i in range(min(nump, l+1))]

            for pid in range(nump):
            
                tmp = vect_mean(pbuf + nbuf)

                if len(pbuf) >= l: pbuf.pop(0)
                pbuf += [points.GetPoint(pid)]
                
                nbuf.pop(0)
                if pid+l+1 < nump: nbuf += [points.GetPoint(pid+l+1)]
    
                points.SetPoint(pid, tmp)

    # find the points on the CLine that are the nearest to the specified start and end
    startId = 0
    endId   = nump-1
    if STAPOINT or ENDPOINT:
        sDist = sys.maxint
        eDist = sys.maxint
        for pid in range(nump):
            p = points.GetPoint(pid)
            if STAPOINT:
                d2 = dist2(p, STAPOINT)
                if d2 < sDist:
                    startId = pid
                    sDist = d2
            if ENDPOINT:
                d2 = dist2(p, ENDPOINT)
                if d2 < eDist:
                    endId = pid
                    eDist = d2

    # the point range in the CLine is [startId,...,endId]

    print 'IDs of starting and ending points:', startId, endId
    incr = 1 if startId < endId else -1
    pIdList = range(startId, endId+incr, incr)

    length = 0.0
    #for pid in pIdList[1:]:
    for i in range(1, len(pIdList)):
        length += math.sqrt(dist2(points.GetPoint(pIdList[i-1]),
                    points.GetPoint(pIdList[i  ])))

    print 'Length of center line section: {0: >8.3f}'.format(length)
    if CUTNUM > 1:
        stepLen = length / (CUTNUM-1)
    else:
        stepLen = 0
    print 'Cuts distance: {0: >8.3f}'.format(stepLen)

    # find cut planes positions
    cutPlanes = []
    currIdx = 0
    currp = points.GetPoint(pIdList[currIdx])

    nextIdx = 1
    while True: # it happens that clines have first point duplicated...
        n = get_normal(currp, points.GetPoint(pIdList[nextIdx]))
        if n: break
        nextIdx += 1

    cutPlanes.append([currp, n])

    for i in range(CUTNUM-1): # we always start from 0 even if there are initial duplicate points
        clen = 0.0

        while True:
            nextIdx = currIdx+1
            nextp = points.GetPoint(pIdList[nextIdx])
            d = math.sqrt(dist2(currp, nextp))

            if (clen + d) > stepLen: break
            if nextIdx == len(pIdList)-1: break
        
            clen += d
            currIdx = nextIdx
            currp = nextp
        
        dl = stepLen-clen
        ratio = dl/d
        #print '\tCurrent polyline length:', clen+dl
        #print '\tCurrent polyline segment:', pIdList[currIdx], pIdList[nextIdx]

        p = tuple([currp[0] + ratio*(nextp[0]-currp[0]),
                   currp[1] + ratio*(nextp[1]-currp[1]),
                   currp[2] + ratio*(nextp[2]-currp[2])])

        cutPlanes.append([p, get_normal(currp, nextp)])
        currp = p

    stl = vtk.vtkSTLReader()
    stl.SetFileName(STLFNAME)
    pdw = vtk.vtkPolyDataWriter()

    fsection = open('sections.dat','w')

    for i in range(len(cutPlanes)): 
        
        p, n = cutPlanes[i]
     
        #print 'Cross section {0}: position {1: >8.3f} {2: >8.3f} {3: >8.3f}'.format(i, p[0], p[1], p[2])
        #print 'Cross section {0}:   normal {1: >8.3f} {2: >8.3f} {3: >8.3f}'.format(i, n[0], n[1], n[2])
        
        plane = vtk.vtkPlane()
        plane.SetOrigin(p)
        plane.SetNormal(n)
     
        cutEdges = vtk.vtkCutter()
        cutEdges.SetInputConnection(stl.GetOutputPort())
        cutEdges.SetCutFunction(plane)
        #cutEdges.GenerateCutScalarsOn()
        #cutEdges.SetValue(0, 0.0)
     
        cutStrips = vtk.vtkStripper()
        cutStrips.SetInputConnection(cutEdges.GetOutputPort())
        cutStrips.Update()
        
        #cutPoly = vtk.vtkPolyData()
        #cutPoly.SetPoints(cutStrips.GetOutput().GetPoints())
        #cutPoly.SetLines(cutStrips.GetOutput().GetLines())
     
        cutPoly = vtk.vtkPolyData()
        cutPoly.SetPoints(cutStrips.GetOutput().GetPoints())
        cutPoly.SetLines(cutStrips.GetOutput().GetLines())
     
        conn = vtk.vtkPolyDataConnectivityFilter()
        if vtk.VTK_MAJOR_VERSION <= 5:
            conn.SetInput(cutPoly)
        else:
            conn.SetInputData(cutPoly)
        conn.SetExtractionModeToClosestPointRegion() # extract region nearest to point
        conn.SetClosestPoint(p)
        conn.Update()
     
        section = conn.GetOutput()
     
        # compute polygons metrics
        peri, area = polygon_metrics(section.GetCell(0).GetPoints(), 
                                    section.GetCell(0).GetNumberOfPoints()-1) # last is the repetition of the first!
     
        # print 'Cross section {0}: perimeter {1: >8.3f} - area {2: >8.3f} - m.diam {3: >8.3f}'.format(i, 
        #                                                 peri, 
        #                                                 area, 
        #                                                 math.sqrt(area/(math.pi))*2.0)
     
        print >> fsection, '{0} {1: >8.3f} {2: >8.3f}'.format(i, area, math.sqrt(area/(math.pi))*2.0)
     
        # write data with only azimutal lines
        # pdw.SetInput(section) 
        # pdw.SetFileName('crossSect'+str(i)+'.vtk')
        # pdw.Write()
     
        # create cells with triangular cells
        merge = vtk.vtkAppendPolyData()
        center = vtk.vtkPolyData()
        cc = vtk.vtkPoints()
        cc.InsertNextPoint(p) # insert center
        center.SetPoints(cc)
     
        if vtk.VTK_MAJOR_VERSION <= 5:
            merge.AddInput(center) 
            merge.AddInput(section) 
        else:
            merge.AddInputData(center) 
            merge.AddInputData(section) 
        merge.Update() 
     
        merge.GetOutput().DeleteCells()
        segmCells = vtk.vtkCellArray()
     
        line = vtk.vtkLine()
     
        nump = section.GetNumberOfPoints()
     
        SHOWTRIANGLES=False
        for k in range(1,nump+1):
     
            if SHOWTRIANGLES:
                t = vtk.vtkTriangle()
                t.GetPointIds().SetId(0,0)
                t.GetPointIds().SetId(1,k)
                t.GetPointIds().SetId(2,k%nump+1)
                segmCells.InsertNextCell(t)
            else:
                line.GetPointIds().SetId(0,k)
                line.GetPointIds().SetId(1,k%nump+1)
                segmCells.InsertNextCell(line)
     
        merge.GetOutput().SetLines(segmCells)
        # print '# of Cells:',merge.GetOutput().GetNumberOfCells()
     
        field = vtk.vtkFieldData()
        field.SetNumberOfTuples(3)
     
        val = vtk.vtkFloatArray()
        val.SetName("area")
        val.InsertNextValue(area)
        field.AddArray(val)
     
        val = vtk.vtkFloatArray()
        val.SetName("mean_diameter")
        val.InsertNextValue( 2.0 * math.sqrt(area/math.pi) )
        field.AddArray(val)
     
        val = vtk.vtkFloatArray()
        val.SetName("perimeter")
        val.InsertNextValue(peri)
        field.AddArray(val)
     
        merge.GetOutput().SetFieldData(field)
     
        merge.Update()
     
        if vtk.VTK_MAJOR_VERSION <= 5:
            pdw.SetInput(merge.GetOutput()) 
        else:
            pdw.SetInputData(merge.GetOutput()) 
     
        if i<10:
            pad = '00'+str(i)
        elif i<100:
            pad = '0'+str(i)
        elif i<1000:
            pad = str(i)
        else:
            pad = 'XXX'+str(i)
        pdw.SetFileName('crossSect'+pad+'.vtk')
     
        pdw.Write()


    ## uncomment the block below to write cut normals in a vtk file
    ## for review
    #wpoints = vtk.vtkPoints()
    #for p in cutPlanes:
    #    wpoints.InsertNextPoint(p[0])
    #    # to visualize normals add a point along them
    #    q = [p[0][0]+p[1][0], p[0][1]+p[1][1], p[0][2]+p[1][2]]
    #    wpoints.InsertNextPoint(q)
    #
    #polydata = vtk.vtkPolyData()
    #polydata.SetPoints(wpoints)
    #
    #segmCells = vtk.vtkCellArray()
    #for i in range(len(cutPlanes)):
    #    line = vtk.vtkLine()
    #    line.GetPointIds().SetId(0,2*i)
    #    line.GetPointIds().SetId(1,2*i+1)
    #    segmCells.InsertNextCell(line)
    #
    #polydata.SetLines(segmCells)
    #
    #pdw.SetInput(polydata) 
    #pdw.SetFileName('cut_normals.vtk')
    #pdw.Write()
