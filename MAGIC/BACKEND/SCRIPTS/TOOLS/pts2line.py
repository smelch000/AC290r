#!/usr/bin/env python

from vtk import *

# readerSTL = vtkGenericDataObjectReader()
# readerSTL.SetFileName('lumen_ascii_fillholes.stl')
# readerSTL.Update()

reader = vtkGenericDataObjectReader()
reader.SetFileName('cline_partial_points.vtk')
reader.Update()

geomFilt = vtkGeometryFilter()
if VTK_MAJOR_VERSION <= 5:
    geomFilt.SetInput(reader.GetOutput())
else:
    geomFilt.SetInputData(reader.GetOutput())
geomFilt.Update()


# for i in xrange(npts):
#     print geomFilt.GetOutput().GetPoint(i)

intSegmCells = vtkCellArray()

NSKIP = 5
npts = geomFilt.GetOutput().GetNumberOfPoints() - NSKIP
for pid in xrange(0, npts-1, NSKIP):
        
    pif = min(pid + NSKIP, npts-1)

    line = vtkLine()
    line.GetPointIds().SetId(0,pid)
    line.GetPointIds().SetId(1,pif)

    intSegmCells.InsertNextCell(line)
                    
finalSeg = vtkPolyData()
finalSeg.SetPoints(geomFilt.GetOutput().GetPoints())
finalSeg.SetLines(intSegmCells)

pdw = vtkPolyDataWriter()

pdw.SetInputData(finalSeg)

"""
# dijkstra requires points and connections (from finalSeg these are lines)
dijkstra = vtkDijkstraGraphGeodesicPath()
if vtk.VTK_MAJOR_VERSION <= 5:
    dijkstra.SetInput(finalSeg)
else:
    dijkstra.SetInputData(finalSeg)
dijkstra.SetStartVertex(0)
dijkstra.SetEndVertex(npts-1)
dijkstra.Update()
# lenArray = vtk.vtkDoubleArray()
# dijkstra.GetCumulativeWeights(lenArray)

pdw.SetInputConnection(dijkstra.GetOutputPort())
"""

pdw.SetFileName('cline_fin.vtk')
pdw.Write()

