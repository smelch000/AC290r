#!/usr/bin/env python

import os,sys

from myvtk import *

def filterShells(filterIn=None, filenameIn=None):

    if filterIn != None and filenameIn != None:
        print 'Usage: filterShells PolyDataFilter or file.stl'
        return

    if filenameIn != None:
        filterIn = vtkSTLReader()
        filterIn.SetFileName(filenameIn)
        filterIn.Update()

    if hasattr(filterIn, 'GetOutput'):
       dataIn = filterIn.GetOutput()
    else:
       dataIn = filterIn

    # print 'DATAIN:',dataIn

    #props = vtkMassProperties()
    #props.SetInputData(filterIn)
    #props.Update()
    #print 'vol:',props.GetVolume(), 'surf:',props.GetSurfaceArea()

    conn = vtkPolyDataConnectivityFilter()
    conn.SetInputData(dataIn)

    # conn.SetExtractionModeToLargestRegion()
    conn.SetFullScalarConnectivity(True)
    conn.ScalarConnectivityOn()
    # conn.SetScalarRange(0, 100)
    conn.ColorRegionsOn()
    conn.Update()

    conn.SetExtractionModeToAllRegions()

    return conn

def filterShellsRetrieve(conn, percent=100., sortmethod='volume'):

    nregions = conn.GetNumberOfExtractedRegions()
    # print 'nregions:',nregions
    # print '::::', conn.GetExtractionModeAsString()

    thresh = vtkThreshold()
    dgeoms = {}
    for i in xrange(nregions):

        thresh.SetInputConnection(conn.GetOutputPort())
        thresh.ThresholdBetween(i, i+1.001)
        thresh.SetInputArrayToProcess(1, 0, 0, 0, "regionId")
        thresh.Update()

        geom = vtkGeometryFilter()
        geom.SetInputConnection(thresh.GetOutputPort())
        geom.Update()

        props = vtkMassProperties()
        props.SetInputConnection(geom.GetOutputPort())
        props.Update()

        #writer = vtkSTLWriter()
        #writer.SetFileName('T' + str(i) + '.stl')
        #writer.SetInputConnection(geom.GetOutputPort())
        #writer.Update()

        if sortmethod == 'volume':
            dgeoms[props.GetVolume()] = geom

        elif sortmethod in ['surface', 'surfacearea', 'surfaceArea']:
            dgeoms[props.GetSurfaceArea()] = geom

        # print i,'vol:',props.GetVolume(),'surf:',props.GetSurfaceArea(),'Region Size:',conn.GetRegionSizes().GetTuple(i)[0],thresh.GetOutput().GetNumberOfPoints()

    n = len(dgeoms.keys())
    v0 = sorted(dgeoms.keys())[n-1]

    appendF = vtk.vtkAppendPolyData()
    for v in sorted(dgeoms.keys(),reverse=True):

        if v/v0 > percent/100.:

            #print v/v0, dgeoms[v].GetOutput().GetNumberOfPoints()

            appendF.AddInputData(dgeoms[v].GetOutput())
            appendF.Modified()
            appendF.Update()

    appendF.Update()

    return appendF

if __name__ == '__main__':

    # first find the statistics of solid

    reader = vtkSTLReader()
    reader.SetFileName('multiShell.stl')
    reader.Update()

    #conn = filterShells(filterIn=reader.GetOutput())

    # conn = filterShells(filterIn=reader)

    conn = filterShells(filenameIn='multiShell.stl')

    # then obtain the filter with data
    geomFilt = filterShellsRetrieve(conn, percent=3., sortmethod='volume')

    #writer = vtkSTLWriter()
    #writer.SetFileName('T.stl')
    #writer.SetInputConnection(geomFilt.GetOutputPort())
    #writer.Update()

