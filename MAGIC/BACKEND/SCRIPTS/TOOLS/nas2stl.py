#!/usr/bin/env python

import sys,os,time,math,logging
import argparse

from myvtk import *
try:
        import numpy
except:
        logging.info('numpy unfound.')

def getMePoints(fp):
        points = []
        for line in fp:
                if 'GRID*' in line:
                        row1 = line.split()
                        nextLine = fp.next()
                        row2 =nextLine.split()
                        x = float(row1[3])
                        y = float(row1[4].split('*')[0] )
                        z = float(row2[1])
                        points.append( (x,y,z) )

        return points

def getTri(fp):
        tris = []
        for line in fp:
                if 'CTRIA3' in line:
                        row = line.split()
                        triId = int(row[2] )
                        v1    = int(row[3] )
                        v2    = int(row[4] )
                        v3    = int(row[5] )
                        tris.append( (triId, v1, v2, v3 ) )

        return tris 

def getNameMap(fp):
        names = {}
	count = 0 
        for line in fp:
		count = count + 1
                if '$' in line and ('GRID POINTS' not in line) and ('ELEMENTS' not in line ) and ('ANSA_NAME' not in line) and count > 10:
			print 'in getNameMap', line
			name = line.replace('$','').replace('\n','')
                        nextLine = fp.next()
			if 'PSHELL' in nextLine.split()[0] :
				PID = nextLine.split()[1]
				print name, PID
                        	#PID = int(row[1] )
				names[PID] = name 
	return names

def mkpolyData(points,tris):
        vpoints = vtk.vtkPoints()
        for p in points:
                vpoints.InsertNextPoint(p[0], p[1], p[2] )

        vtris = vtk.vtkCellArray()
        ids = vtk.vtkFloatArray()
        for t in tris:
                vtri = vtk.vtkTriangle()
                vtri.GetPointIds().SetId(0, t[1]-1)
                vtri.GetPointIds().SetId(1, t[2]-1)
                vtri.GetPointIds().SetId(2, t[3]-1)
                ids.InsertNextValue(t[0])       
                vtris.InsertNextCell(vtri)
                                                          
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(vpoints)
        polydata.SetPolys(vtris)
        polydata.Modified()
        polydata.GetCellData().SetScalars(ids)
        #polydata.Update()

        return polydata

def VtkFromNasFile(infile):     
        a = open(infile,'r')
        points = getMePoints(a)
        b = open(infile,'r')
        tris = getTri(b)
        b = open(infile,'r')
	nameMap= getNameMap(b)
	print 'in vtkFromNasFIle, MAP=',nameMap
        polyData =  mkpolyData(points,tris)
        return ( nameMap, polyData ) 

def NasFromVtkPolyData(outfile,mypoly,mymap): 
        FileStr = '''$ =========================================
$ Version:
$ Created:      
$ =========================================

CEND
$ GRID POINTS
BEGIN BULK
'''
        def mkLine(cds,i):
                blank = "        "
                grid  = "GRID*   "
                myID  = "{:>8}".format(i)
                zero  = "{:>8}".format(0)
                crd1  = "{:>16}".format( "%E"%(crds[0]) )
                crd2  = "{:>16}".format( "%E"%(crds[1]) )
                end   = "*G{:<6}".format(i)
                line1 =  grid+blank+myID+blank+zero+crd1+crd2+end+'\n'  
                crd2  = "{:>16}".format( "%E"%(crds[2]) )
                start   = "*G{:<6}".format(i)
                line2=start+crd2+'\n' 
                return line1+line2


        pnts = mypoly.GetPoints()       
        npts = pnts.GetNumberOfPoints()
        logging.info("npts ="+str(npts) )
        for i in range(1,npts+1):
                crds = pnts.GetPoint(i-1)
                nextLines = mkLine(crds,i)
                FileStr = FileStr + nextLines

        FileStr = FileStr + '$ELEMENTS\n'
        poly = mypoly.GetPolys()        
        nPoly = poly.GetNumberOfCells()
        idArr = mypoly.GetCellData().GetScalars()
	haveIds = True
	if idArr is None:
		logging.info('We have no scalar Map yet!')
		haveIds = False
                if mymap == None:
                    mymap = {}
		mymap['1'] = 'wall'
        poly.InitTraversal()
        logging.info('Npoly ='+str(nPoly) )
        countIds = {}   
        for i in range(0,nPoly):
                cell = vtk.vtkIdList()
                poly.GetNextCell(cell) 
		if haveIds is True:
                	mytag = int(idArr.GetTuple(i)[0])       
		else:
			mytag = 1

                if mytag in countIds:
                        countIds[mytag] += 1
                else:
                        countIds[mytag] = 1
                        
                tri= 'CTRIA3  ' 
                myID  = "{:>8}".format(i+1)
                tag   = "{:>8}".format(mytag)
                vtx1  = "{:>8}".format( cell.GetId(0)+1 )
                vtx2  = "{:>8}".format( cell.GetId(1)+1 )
                vtx3  = "{:>8}".format( cell.GetId(2)+1 )
                nextLine = tri+myID+tag+vtx1+vtx2+vtx3+'\n'
                FileStr = FileStr + nextLine
        nLines = npts*2 + nPoly+ 7
        logging.info( 'UniquePIDs =',countIds)
        logging.info(nLines)
	print 'mymap=',mymap
	print 'countIds=',countIds
        for curId in countIds.keys():
                tagLine1 = '$ANSA_NAME;%i;PSHELL;~\n'%(curId)
                if str(curId) in mymap:
                	tagLine2 = '$%s\n'%( mymap[str(curId)] )
		else:
                	tagLine2 = '$wall\n'
                tagLine3 = 'PSHELL  %i       1                                                       +%i\n'%(curId,nLines)
                tagLine4 = '+%i\n'%nLines
                FileStr = FileStr + tagLine1 + tagLine2 + tagLine3 + tagLine4 
                nLines = nLines + 4     
        FileStr = FileStr + 'ENDDATA\n'
        open(outfile,'w').write(FileStr)



def getPoints(infile):
        f = open(infile,'r')
        crd = []
        allPts = []
        for line in f:
                if 'GRID' in line:
                        vals = line.split(' ')
                        myvals = [ x for x in vals if x != '' ]
                        if len(myvals) > 3: 
                                crd.append(float(myvals[3]) )
                                tmp = myvals[4].split('*')[0]
                                crd.append(float(tmp))
                else:
                        if '*G' in line:
                                tmp = line.split()[1]
                                crd.append(float(tmp) )
                                allPts.append(crd)
                                crd = []
        return allPts

def getTris(infile,num):
        f = open(infile,'r')
        allTris = []
        for line in f:
                if 'CTRIA' in line:
                        tmp = line.split()
                        if tmp[2] == num:
                                allTris.append( [ int(tmp[3]), int(tmp[4]), int(tmp[5]) ])
        return allTris

def getBox(infile, offset=0.0):
        allPts = getPoints(infile)

        xLow  = allPts[0][0]    
        xHigh = allPts[0][0]    
        yLow  = allPts[0][1]    
        yHigh = allPts[0][1]    
        zLow  = allPts[0][2]    
        zHigh = allPts[0][2]    
        for point in allPts:

                if xLow > point[0]:
                        xLow = point[0]
                if xHigh < point[0]:
                        xHigh = point[0]

                if yLow > point[1]:
                        yLow = point[1]
                if yHigh < point[1]:
                        yHigh = point[1]

                if zLow > point[2]:
                        zLow = point[2]
                if zHigh < point[2]:
                        zHigh = point[2]
                
        return [ (xLow - offset, xHigh+ offset ), (yLow - offset , yHigh + offset ), (zLow - offset , zHigh+ offset) ]

def getNames(infile,name):
	print infile
        f = open(str(infile),'r')
        names = []
        for line in f:
                if name in line.lower():
                        names.append(line.replace('$','').replace('\n',''))     
        return names

def getOutletNames(infile):
        return getNames(infile,'outlet')

def getInletNames(infile):
        return getNames(infile,'inlet')


def getInletArea(infile):
        f = open(infile,'r')
        crd = []
        allPts = []
        
        for line in f:
                if 'inlet' in line.lower():
                        nextLine = f.next()
                        logging.log('Inlet ID ***', nextLine.split(' ')[2] ,'***')
                        id = nextLine.split(' ')[2]

        allPts = getPoints(infile)
        tris   = getTris(infile, id)
        triTotal = 0.0
        
        for t in tris:
                pt1 = allPts[t[0]-1]
                pt2 = allPts[t[1]-1]
                pt3 = allPts[t[2]-1]
                u = [ pt1[0] - pt2[0] , 
                      pt1[1] - pt2[1] , 
                      pt1[2] - pt2[2] ]
                v = [ pt1[0] - pt3[0] , 
                      pt1[1] - pt3[1] , 
                      pt1[2] - pt3[2] ]

                uxv_i = u[1]*v[2] - u[2]*v[1] 
                uxv_j = u[2]*v[0] - u[0]*v[2] 
                uxv_k = u[0]*v[1] - u[1]*v[0]
                area = math.sqrt( uxv_i*uxv_i + uxv_j*uxv_j + uxv_k*uxv_k)
                triTotal = triTotal + area 

        return triTotal/2.0 

def getAllCentroids( infile):
        faces = getInletNames(infile)
        faces.extend( getOutletNames(infile) ) 

        centers = [ infile ]    
        for faceNm in faces:
                centers.append( [faceNm, getFaceCentroid(infile, faceNm)] )

        return centers


def getFaceCentroid(infile, faceNm):
        f = open(infile,'r')
        crd = []
        allPts = []

        
        for line in f:
                if faceNm.lower() in line.lower():
                        nextLine = f.next()
                        id = nextLine.split(' ')[2]
        
                        
        allPts = getPoints(infile)
        tris   = getTris(infile, id)
        t1, t2 , t3 = zip(*tris)
        t1 = list(t1)
        t2 = list(t2)
        t3 = list(t3)
        t1.extend(t2)
        t1.extend(t3)
        t1 = numpy.int_(t1)
        t4 = numpy.unique(t1)

        myPoints = []
        for i in t4:
                myPoints.append(allPts[i-1])

        xSum = 0.0
        ySum = 0.0
        zSum = 0.0
        for pt in myPoints:
                xSum = xSum + pt[0] 
                ySum = ySum + pt[1] 
                zSum = zSum + pt[2] 
        xAvg = xSum / len(myPoints)     
        yAvg = ySum / len(myPoints)     
        zAvg = zSum / len(myPoints)     

        return (xAvg, yAvg, zAvg) 


def main(args):

        logging.info( getAllCentroids(args.infile) )

        # vpdo = nas.VtkFromNasFile(sys.argv[1])
        vpdo = VtkFromNasFile( args.infile )

        stlw = vtk.vtkSTLWriter()
        stlw.SetFileName( args.outfile )

        if vtk.VTK_MAJOR_VERSION <= 5:
            stlw.SetInput(vpdo[1])
        else:
            stlw.SetInputData(vpdo[1])

        stlw.Write()

###############################
def iniparser(parser):

    parser.add_argument('-i', '--infile', default='bgkflag',  help='input mesh file (.nas)')
    parser.add_argument('-o', '--outfile', default='bgkflag',  help='input mesh file (.stl)')

    return parser

###############################
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = iniparser(parser)
    args = parser.parse_args()

    main(args)
