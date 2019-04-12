#!/usr/bin/env python


# Example of vtkBooleanOperationPolyDataFilter :
#       Computes the boundary of the union, intersection, or difference volume computed 
#       from the volumes defined by two input surfaces.
#
# The two surfaces do not need to be manifold, but if they are not, unexpected results 
# may be obtained. The resulting surface is available in the first output of the filter. 
# The second output contains a set of polylines that represent the intersection 
# between the two input surfaces.

import vtk
import sys
import argparse

class ErrorObserver:

   def __init__(self):
       self.__ErrorOccurred = False
       self.__ErrorMessage = None
       self.CallDataType = 'string0'

   def __call__(self, obj, event, message):
       self.__ErrorOccurred = True
       self.__ErrorMessage = message

   def ErrorOccurred(self):
       occ = self.__ErrorOccurred
       self.__ErrorOccurred = False
       return occ

   def ErrorMessage(self):
       return self.__ErrorMessage


def Dice(volI,volD1,volD2):
    # volI : intersection: true positive
    # volD1 : difference : false positive
    # volD2 : difference : false negative
    return ( (2.*volI) / (volD1 + volD2 + 2.*volI) )

def Union(gONE,gTWO):
    err0 = ErrorObserver()

    bOp = vtk.vtkBooleanOperationPolyDataFilter()
    bOp.AddObserver('ErrorEvent', err0)

    bOp.SetOperationToUnion()
    # bOp.SetOperationToIntersection()
    # bOp.SetOperationToDifference()
    if vtk.VTK_MAJOR_VERSION <= 5:
        bOp.SetInputConnection( 0, gONE.GetProducerPort() )
        bOp.SetInputConnection( 1, gTWO.GetProducerPort() )
    else:
        bOp.SetInputData( 0, gONE.GetOutput() )
        bOp.SetInputData( 1, gTWO.GetOutput() )

    bOp.Update()
    if err0.ErrorOccurred():
        return None

    return bOp

def Intersection(gONE,gTWO):
    err1 = ErrorObserver()
    bOp = vtk.vtkBooleanOperationPolyDataFilter()
    bOp.AddObserver('ErrorEvent', err1)

    # bOp.SetOperationToUnion()
    bOp.SetOperationToIntersection()
    # bOp.SetOperationToDifference()

    if vtk.VTK_MAJOR_VERSION <= 5:
        bOp.SetInputConnection( 0, gONE.GetProducerPort() )
        bOp.SetInputConnection( 1, gTWO.GetProducerPort() )
    else:
        bOp.SetInputData( 0, gONE.GetOutput() )
        bOp.SetInputData( 1, gTWO.GetOutput() )

    bOp.Update()
    if err1.ErrorOccurred():
        print '++++++++++++++++++++++ ERRRRRRRRR'
        return None

    return bOp

def Difference(gONE,gTWO):
    err2 = ErrorObserver()
    bOp = vtk.vtkBooleanOperationPolyDataFilter()

    # bOp.SetOperationToUnion()
    # bOp.SetOperationToIntersection()
    bOp.SetOperationToDifference()

    bOp.AddObserver('ErrorEvent', err2)

    if vtk.VTK_MAJOR_VERSION <= 5:
        bOp.SetInputConnection( 0, gONE.GetProducerPort() )
        bOp.SetInputConnection( 1, gTWO.GetProducerPort() )
    else:
        bOp.SetInputData( 0, gONE.GetOutput() )
        bOp.SetInputData( 1, gTWO.GetOutput() )

    bOp.Update()
    if err2.ErrorOccurred():
        return None

    return bOp

def getProps(bX):
    err = ErrorObserver()
    props = vtk.vtkMassProperties()
    props.AddObserver('ErrorEvent', err)
    props.SetInputConnection(bX.GetOutputPort())
    props.Update()
    if err.ErrorOccurred():
        return -100000.,-100000.
    else:
        return props.GetVolume(),props.GetSurfaceArea()

def IsEncapsulated(poly1,poly2):

    pointsPolydata = vtk.vtkPolyData()
    pointsPolydata.SetPoints( poly2.GetPoints() )

    selectEnclosedPoints = vtk.vtkSelectEnclosedPoints()
    if vtk.VTK_MAJOR_VERSION <= 5:
        selectEnclosedPoints.SetInput(pointsPolydata);
    else:
        selectEnclosedPoints.SetInputData(pointsPolydata);

    if vtk.VTK_MAJOR_VERSION <= 5:
        selectEnclosedPoints.SetSurface(poly1)
    else:
        selectEnclosedPoints.SetSurfaceData(poly1)

    selectEnclosedPoints.Update();
    return selectEnclosedPoints.IsInside(0)

fileOutputWindow = vtk.vtkFileOutputWindow()
fileOutputWindow.SetFileName( "output.txt" );
# Note that the SetInstance function is a static member of vtkOutputWindow.
vtk.vtkOutputWindow.SetInstance(fileOutputWindow)

parser = argparse.ArgumentParser()

parser.add_argument('-1', '--input1',  required=True,  help='input of 1st stl ')
parser.add_argument('-2', '--input2',  required=True,  help='input of 2nd stl')
parser.add_argument('-x', '--operation',  required=True,  help='operation (union, intersection)')
args = parser.parse_args()

gONE = vtk.vtkSTLReader()
gONE.SetFileName( args.input1 )
gONE.Update()

mONE = vtk.vtkPolyDataMapper()
if vtk.VTK_MAJOR_VERSION <= 5:
  mONE.SetInputConnection( gONE.GetProducerPort() )
else:
  # mONE.SetInputData( gONE )
  mONE.SetInputData( gONE.GetOutput() )
# mONE.ScalarVisibilityOff()
mONE.ScalarVisibilityOn()

aONE = vtk.vtkActor()
aONE.SetMapper( mONE )

gTWO = vtk.vtkSTLReader()
gTWO.SetFileName( args.input2 )
gTWO.Update()

mTWO = vtk.vtkPolyDataMapper()
if vtk.VTK_MAJOR_VERSION <= 5:
  mTWO.SetInputConnection( gTWO.GetProducerPort() )
else:
  # mTWO.SetInputData( gTWO )
  mTWO.SetInputData( gTWO.GetOutput() )
# mTWO.ScalarVisibilityOff()

aTWO = vtk.vtkActor()
aTWO.SetMapper( mTWO )

if args.operation == 'union':
    # print 'Union....'
    bOp = Union(gONE,gTWO)

elif args.operation == 'intersection':
    # print 'Intersection....'
    bOp = Intersection(gONE,gTWO)

elif args.operation == 'difference':
    # print 'Difference....'
    bOp = Difference(gONE,gTWO)

else:
    # print 'Operation unrecognized!!!'
    sys.exit(1)

bU = Union(gONE,gTWO)
bI = Intersection(gONE,gTWO)
bD1 = Difference(gONE,gTWO)
bD2 = Difference(gTWO,gONE)

volU,surfU = getProps(bU)
# print 'Union....            Volume:',volU,'  Surface',surfU
volI,surfI = getProps(bI)
# print 'Intersection....     Volume:',volI,'  Surface',surfI
volD1,surfD1 = getProps(bD1)
# print 'Difference 1-2....   Volume:',volD1,'  Surface',surfD1
volD2,surfD2 = getProps(bD2)
# print 'Difference 2-1....   Volume:',volD2,'  Surface',surfD2

TWOinONE = False
ONEinTWO = False

if volU<0:
    TWOinONE = IsEncapsulated(gONE.GetOutput(),gTWO.GetOutput())

    ONEinTWO = IsEncapsulated(gTWO.GetOutput(),gONE.GetOutput())

    if TWOinONE:
        vIN,S = getProps(gTWO)
        vOUT,S = getProps(gONE)
    elif ONEinTWO:
        vIN,S = getProps(gONE)
        vOUT,S = getProps(gTWO)

    score = Dice(vIN, abs(vIN-vOUT), 0.0)
else:
    score = Dice(volI, volD1, volD2)

print '              Global Dice Score: %f '%score

bMapper = vtk.vtkPolyDataMapper()
bMapper.SetInputConnection( bOp.GetOutputPort() );
bMapper.ScalarVisibilityOff();

bActor = vtk.vtkActor()
bActor.SetMapper( bMapper )

# The usual rendering stuff is created.
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

bActor.GetProperty().SetColor(0,0,1)
aONE.GetProperty().SetColor(1,0,0)
aTWO.GetProperty().SetColor(0,1,1)
aONE.GetProperty().SetOpacity(0.1)
aTWO.GetProperty().SetOpacity(0.1)

# Add the actors to the renderer, set the background and size
ren.AddActor(bActor);
ren.AddActor(aONE);
ren.AddActor(aTWO);

ren.SetBackground(1.,1.,1.)
renWin.SetSize(500, 500)
ren.ResetCamera()

iren.Initialize()
renWin.Render()



cONE = vtk.vtkPolyDataConnectivityFilter()
if vtk.VTK_MAJOR_VERSION<=5:
    cONE.SetInput(gONE.GetOutput())
else:
    cONE.SetInputData(gONE.GetOutput())
# cONE.SetExtractionModeToLargestRegion()
# cONE.SetExtractionModeToAllRegions()
cONE.Update()
sONE = cONE.GetRegionSizes()

cTWO = vtk.vtkPolyDataConnectivityFilter()
if vtk.VTK_MAJOR_VERSION<=5:
    cTWO.SetInput(gTWO.GetOutput())
else:
    cTWO.SetInputData(gTWO.GetOutput())
# cTWO.SetExtractionModeToLargestRegion()
# cTWO.SetExtractionModeToAllRegions()
cTWO.Update()
sTWO = cTWO.GetRegionSizes()

cONE.SetExtractionModeToAllRegions()
cONE.Update()
cTWO.SetExtractionModeToAllRegions()
cTWO.Update()

nONE = cONE.GetNumberOfExtractedRegions()
nTWO = cTWO.GetNumberOfExtractedRegions()

# print '\nNumber of individual plaques:',nONE,nTWO,'\n\n'

cONE.SetExtractionModeToSpecifiedRegions()
cTWO.SetExtractionModeToSpecifiedRegions()

score = {}
size = {}
for i in xrange(nONE):

    cONE.InitializeSpecifiedRegionList()
    cONE.AddSpecifiedRegion(i)
    cONE.Update()
    bitONE  = vtk.vtkCleanPolyData()
    bitONE.SetInputData(cONE.GetOutput())
    bitONE.Update()
    pONE = vtk.vtkPolyData()
    pONE.DeepCopy(bitONE.GetOutput() )
    fONE = vtk.vtkGeometryFilter()
    fONE.SetInputData(pONE)
    fONE.Update()

    # mmONE = vtk.vtkPolyDataMapper()
    # mmONE.SetInputData( fONE.GetOutput() )
    # aaONE = vtk.vtkActor()
    # aaONE.SetMapper( mmONE )
    # aaONE.GetProperty().SetColor(1,0,0)
    # aaONE.GetProperty().SetOpacity(0.3)
    # ren.AddActor(aaONE);
    # renWin.Render()

    for j in xrange(nTWO):

        cTWO.InitializeSpecifiedRegionList()
        cTWO.AddSpecifiedRegion(j)
        cTWO.Update()
        bitTWO  = vtk.vtkCleanPolyData()
        bitTWO.SetInputData(cTWO.GetOutput())
        bitTWO.Update()
        pTWO = vtk.vtkPolyData()
        pTWO.DeepCopy(bitTWO.GetOutput() )
        pTWO.Modified()
        fTWO = vtk.vtkGeometryFilter()
        fTWO.SetInputData(pTWO)
        fTWO.Update()

        # mmTWO = vtk.vtkPolyDataMapper()
        # mmTWO.SetInputData( fTWO.GetOutput() )
        # aaTWO = vtk.vtkActor()
        # aaTWO.SetMapper( mmTWO )
        # aaTWO.GetProperty().SetColor(0,1,0)
        # aaTWO.GetProperty().SetOpacity(0.3)
        # # ren.AddActor(aaTWO);
        # # renWin.Render()

        # vol0,surf0 = getProps2(pONE)
        vol0,surf0 = getProps(fONE)
        # vol0,surf0 = getProps(bitONE)

        bU = Union(fONE,fTWO)
        bI = Intersection(fONE,fTWO)
        bD1 = Difference(fONE,fTWO)
        bD2 = Difference(fTWO,fONE)

        volU,surfU = getProps(bU)
        volI,surfI = getProps(bI)
        volD1,surfD1 = getProps(bD1)
        volD2,surfD2 = getProps(bD2)

        vIN,S = getProps(fONE)
        vOUT,S = getProps(fTWO)

        if volU > 0 and volI > 0 and volD1 > 0 and volD2 > 0:
            # ren.AddActor(aaTWO);
            # renWin.Render()
            size[i] = volU
            score[i] = Dice(volI, volD1, volD2)
            # print i,j,'volI:',volI,'volD1:',volD1,'voldD2:',volD2,
            # print i,j,' Dice Score:', score[i],'....size:',size[i],getProps(fONE)[0],getProps(fTWO)[0],'SIZEs:',sONE.GetTuple(i)[0], sTWO.GetTuple(j)[0]

            print 'id1-id2: %4d %4d:' % (i,j), \
                  ' Dice Score: %f' % score[i], \
                  '        volumes: %f %f' % (vIN,vOUT)
        else:
            TWOinONE = IsEncapsulated(fONE.GetOutput(),fTWO.GetOutput())
            ONEinTWO = IsEncapsulated(fTWO.GetOutput(),fONE.GetOutput())

            if TWOinONE:
                vIN,S = getProps(fTWO)
                vOUT,S = getProps(fONE)
            elif ONEinTWO:
                vIN,S = getProps(fONE)
                vOUT,S = getProps(fTWO)
            else:
                continue
                # print 'BOH!'

            score[i] = Dice(vIN, abs(vIN-vOUT), 0.0)

            print 'id1,id2: %4d %4d:' % (i,j), \
                  ' Dice Score: %f' % score[i], \
                  '        volumes: %f %f' % (vIN,vOUT)
            # print i,j,'No Intersect', getProps(fONE)[0],getProps(fTWO)[0],'SIZEs:',sONE.GetTuple(i)[0], sTWO.GetTuple(j)[0]

        # ren.RemoveActor(aaTWO)
        # renWin.Render()

    # name = raw_input()

    # ren.RemoveActor(aaONE)
    # renWin.Render()

iren.Start()
