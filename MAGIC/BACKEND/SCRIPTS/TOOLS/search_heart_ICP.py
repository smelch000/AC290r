#!/usr/bin/env pythonw

import math,sys,os,random,time
import argparse
from myvtk import *
import numpy as np
from downsample_image import *

MARCHINGLEVEL =  50 # 100 # 200

# Threshold interval for detecting the heart
WLEVELS =  150, 900
# WLEVELS =  100, 900
# WLEVELS = 200, 900

# Threshold interval for probes
# PLEVELS = 150, 900
PLEVELS = 250, 900

# flag to represent the heart as :
#   True: a skin bject      (---> better if used via Delaunay3D)
#   False: a set of points falling within the Threshold window
SURFACE_PROBING = True 
# SURFACE_PROBING = False

TRAOLD = None
SOURCE = None
FEATURE = None

"""
This code finds the best location for the heart inside a CT scan.
It uses a pre-defined heart shape from manual carving (e.g. case21_opt5.0_downsampled_10.vti)
and applies translations - rotations - scaling (R-T-S) in fixed order until a best match is found

The scan and probe images are :
1. thresholded in a given HU level window (image -> polydata)
2. polydata are transformed as enclosing surfaces (featuring outer shape and inner chambers)
3. the scan surface is gaussian splattered (as a new image) in order to do the fitting
4. the probing points for the fit are the polydata points of the heart template

* the best match is given by summing all splatter values obtained via vtk probes

Known issues:
    The method typically finds the heart but at times is upside-down (apex reversed)
    Success rate could be raised by using weights on apex and/or aorta to stress those features
"""

def fakeConvexHull(poly):

    sphereSource = vtkSphereSource()
    sphereSource.SetRadius(500)
    sphereSource.SetPhiResolution(150)
    sphereSource.SetThetaResolution(150)
    sphereSource.Update()
 
    points = poly.GetPoints()
 
    smoothFilter = vtkSmoothPolyDataFilter()

    smoothFilter.SetInputConnection(0, sphereSource.GetOutputPort())
    smoothFilter.SetInputData(1, poly)

    # smoothFilter.SetNumberOfIterations(15)
    # smoothFilter.SetNumberOfIterations(100)
    # smoothFilter.SetNumberOfIterations(200)
    smoothFilter.SetNumberOfIterations(1000)

    # smoothFilter.SetRelaxationFactor(0.1)
    smoothFilter.SetRelaxationFactor(0.5)
    # smoothFilter.FeatureEdgeSmoothingOff()
    # smoothFilter.FeatureEdgeSmoothingOn()
    # smoothFilter.BoundarySmoothingOn()

    smoothFilter.Update()
 
    return smoothFilter # .GetOutputPort()

def TransformPolyActorMatrix(actor,matrix):
    """ 
    Transform actor and returns both the transformed actor and a new transformed polydata
    It works in absolute terms
    """

    TN = vtkTransform()
    TN.SetMatrix( matrix )
    TN.Update()

    actor.SetUserTransform( TN )

    ptrans = vtkTransformPolyDataFilter()

    ptrans.SetTransform(TN)

    if VTK_MAJOR_VERSION<=5:
        ptrans.SetInput(actor.GetMapper().GetInput())
    else:
        ptrans.SetInputData(actor.GetMapper().GetInput())
    ptrans.Update()

    # print CenterOfPoly(ptrans.GetOutput())
    return ptrans.GetOutput()

class Imgview():
    """ 
    ugly class to create the color map with given opacity... to be largely ameliorated
    and to render the image as a vtkVolume
    """
    def __init__(self,imagedata,renderer,uint=False):

        self.volumeColor = vtkColorTransferFunction()

        if uint:
            # self.volumeColor.AddRGBPoint(-1021, 1.0, 0.0, 1.0)
            self.volumeColor.AddRGBPoint(  0,  .0,  .0,  .0)
            self.volumeColor.AddRGBPoint( 10,  .0,  .0, 1.0)
            self.volumeColor.AddRGBPoint(100,  .0, 1.0,  .0)
            self.volumeColor.AddRGBPoint(200, 1.0,  .0,  .0)
            self.volumeColor.AddRGBPoint(256,  .0,  .0,  .0)
        else:
            ## self.volumeColor.AddRGBPoint(-3024, 1.0, 0.0, 1.0)
            self.volumeColor.AddRGBPoint(-1021, 1., 0., 1.)
            self.volumeColor.AddRGBPoint(   0,  0., 1., 0.)
            self.volumeColor.AddRGBPoint( 100,  0., 0., 1.)
            self.volumeColor.AddRGBPoint( 300,  0., 1., 0.)
            self.volumeColor.AddRGBPoint( 600,  0., 1., 1.)
            self.volumeColor.AddRGBPoint( 1000, 0., 0., 0.)
            self.volumeColor.AddRGBPoint( 1500, 0., 0., 0.)
            #self.volumeColor.AddRGBPoint(80, 0.0, 0.0, 1.0)
            #self.volumeColor.AddRGBPoint(100, 1.0, 0.0, 0.0)
            #self.volumeColor.AddRGBPoint(600, 1.0, 1.0, 1.0)


        # self.volumeColor.ClampingOn()

        #  cast noise image to unsigned short
        self.imageCast = vtkImageCast()
        # self.imageCast.SetInputConnection(imagedata)
        if VTK_MAJOR_VERSION <= 5:
            self.imageCast.SetInput(imagedata)
        else:
            self.imageCast.SetInputData(imagedata)
        self.imageCast.SetOutputScalarTypeToUnsignedShort()
        self.imageCast.Update()


        '''
        self.v16 = vtk.vtkImageChangeInformation()
        self.v16.SetInputConnection(imagedata)
        #SM self.v16.SetOutputSpacing(spacing[0], spacing[1], s.GetZSpacing())
        # self.v16.SetOutputSpacing(spacing[0], spacing[1], s.GetSpacing()[2])
        self.v16.Update()
        '''

        #self.opac=[(0.000000,0.000000),
        #            (156.101963,0.052203),
        #            (234.303331,0.435345),
        #            (325.738777,0.765805),
        #            (417.174224,0.919061),
        #            (498.984886,0.995690),
        #            (1500.000000,1.000000)]

        #self.opac=[ (0.000, 0.0),
        #            (200.0, 0.0),
        #            (300.0, 0.9),
        #            (400.0, 0.0),
        #            (600.0, 0.0),
        #            (700.0, 0.0),
        #            (1500., 0.0)]

        if uint:
            self.opac=[ (  0, 0.0),
                        ( 11, 0.01),
                        (100, 0.01),
                        (254, 0.8),
                        (255, 0.8)]
            #self.opac=[ ( 00, 0.0),
            #            ( 10, 0.8),
            #            (154, 0.8),
            #            (155, 0.0)]
            #self.opac=[ (200, 0.0),
            #            (210, 0.8),
            #            (254, 0.8),
            #            (255, 0.0)]

        else:
            #self.opac=[ (-3201, 0.0),
            #            ( 99.0, 0.0),
            #            ( 100.0, 0.1),
            #            (600.0, 0.8),
            #            (601.0, 0.0),
            #            (1500., 0.0)]
            # self.opac=[ (-200, 0.4),
            #             ( 99.0, 0.4),
            #             ( 100.0, 0.5),
            #             (600.0, 0.8),
            #             (601.0, 0.5),
            #             (1500., 0.5)]
            self.opac=[(0.000000,0.000000),
               (156.101963,0.052203),
               (234.303331,0.135345),
               (325.438777,0.165805),
               (417.474224,0.119061),
               (498.384886,0.195690),
               (1500.000000,0.000000)]

        # print 'UINT:',uint

        self.volumeScalarOpacity = vtk.vtkPiecewiseFunction()
        self.volumeScalarOpacity.ClampingOff()
        # self.volumeScalarOpacity.ClampingOn()
        for pnt in self.opac:
            # print "(%f,%f),"%(pnt[0],pnt[1]) 
            self.volumeScalarOpacity.AddPoint(pnt[0],pnt[1])

        self.volumeProperty = vtkVolumeProperty()
        self.volumeProperty.SetColor(self.volumeColor)
        self.volumeProperty.SetScalarOpacity(self.volumeScalarOpacity)
        ##self.volumeProperty.SetGradientOpacity(volumeGradientOpacity)
        self.volumeProperty.SetInterpolationTypeToLinear()
        # self.volumeProperty.ShadeOn()
        self.volumeProperty.SetAmbient(0.4)
        self.volumeProperty.SetDiffuse(0.7)
        self.volumeProperty.SetSpecular(0.2)

        if False:
            rayCastFunction = vtkVolumeRayCastCompositeFunction()
            rayCastFunction.SetCompositeMethodToInterpolateFirst()
            self.volumeMapper = vtkVolumeRayCastMapper()
            self.volumeMapper.SetInputConnection(self.imageCast.GetOutputPort())
            # self.volumeMapper.SetInputConnection(self.v16.GetOutputPort())
            # self.volumeMapper.SetInputConnection(imagedata.GetOutputPort())
            self.volumeMapper.SetVolumeRayCastFunction(rayCastFunction)
            self.volumeMapper.SetSampleDistance(0.5)
        else:
            self.volumeMapper = vtk.vtkSmartVolumeMapper()
            self.volumeMapper.SetBlendModeToComposite()
            # self.volumeMapper[mapObj].SetInputConnection(self.scaleShifter[mapObj].GetOutputPort())
            self.volumeMapper.SetInputData(self.imageCast.GetOutput())
            self.volumeMapper.SetSampleDistance(0.5)


        self.volume = vtk.vtkVolume()
        self.volume.SetMapper(self.volumeMapper)
        self.volume.SetProperty(self.volumeProperty)
        self.volume.PickableOff()	       
	        
        renderer.AddViewProp(self.volume)

def KeyPressEvent(caller,eventId):
    global TRAOLD,SOURCE,probes,FEA1

    key = iren.GetKeySym()

    if key=='Left':
        if i1 != None:
            renderer[0].AddViewProp(i1.volume)
        for probe in probes:
            # probe.AOptPoly.GetProperty().SetRepresentationToPoints()
            probe.AOptPoly.VisibilityOff()
        # FEA1.VisibilityOff()
        FEA1.GetProperty().SetOpacity(0.1)

    elif key=='Right':
        if i1 != None:
            renderer[0].RemoveViewProp(i1.volume)
        for probe in probes:
            probe.AOptPoly.VisibilityOn()
            # probe.AOptPoly.GetProperty().SetRepresentationToWireframe()
            probe.AOptPoly.GetProperty().SetRepresentationToSurface()
        # FEA1.VisibilityOn()
        FEA1.GetProperty().SetOpacity(1.0)

    elif key=='Up':
        if i2 != None:
             renderer[1].AddViewProp(i2.volume)
        if AScanSurface != None:
            AScanSurface.VisibilityOn()
        # for probe in probes:
        #     # probe.AOptPoly.GetProperty().SetRepresentationToPoints()

    elif key=='Down':
        if i2 != None:
            renderer[1].RemoveViewProp(i2.volume)
        if AScanSurface != None:
            AScanSurface.VisibilityOff()
        # for probe in probes:
        #     # probe.AOptPoly.GetProperty().SetRepresentationToSurface()

    elif key=='i':

        icp = vtk.vtkIterativeClosestPointTransform()

        icp.SetSource(SOURCE)
        # icp.SetSource(probes[0].APoly.GetMapper().GetInput())
        icp.SetTarget(target)

        """
        icp.GetLandmarkTransform().SetModeToRigidBody()
        # icp.GetLandmarkTransform().SetModeToSimilarity()
        # icp.GetLandmarkTransform().SetModeToAffine()
        # icp.DebugOn()
        # icp.SetCheckMeanDistance(1)
        icp.SetMaximumNumberOfLandmarks(100) # (10000)
        # icp.SetMaximumNumberOfIterations(10)# (100)
        icp.SetMaximumMeanDistance(0.0001)
        # icp.SetMaximumMeanDistance(0.01)
        # icp.StartByMatchingCentroidsOn()
        icp.Modified()
        icp.Update()
        """

        icp.SetCheckMeanDistance(1)
        icp.SetMaximumMeanDistance(0.001)
        icp.SetMaximumNumberOfIterations(30)
        icp.SetMaximumNumberOfLandmarks(50)
        #  actor.SetUserTransform(icp[ridx][sidx])


        # probes[0].AOptPoly.SetUserTransform(icp)

        # TRACONC = vtkTransform()

        # TRA = icp.GetLandmarkTransform()
        # TRACONC.SetMatrix(TRA.GetMatrix())

        # TRACONC.SetMatrix(icp.GetMatrix())
        # if TRAOLD != None:
        #     pass
        #     # TRACONC.Concatenate(TRAOLD)

        renderer[1].RemoveActor(probes[0].AOptPoly)

        SOURCE = TransformPolyActorMatrix(probes[0].AOptPoly, icp.GetMatrix())

        probes[0].AOptPoly, probes[0].MOptPoly = ActorGlyphVisual(SOURCE)
        probes[0].AOptPoly.GetProperty().SetColor(0,1,0)
        renderer[1].AddActor(probes[0].AOptPoly)

        # SOURCE = TransformPolyActorMatrix(probes[0].AOptPoly, icp.GetLandmarkTransform().GetMatrix())

        # probes[0].AOptPoly.SetUserTransform(TRACONC)
        # probes[0].AOptPoly.SetUserTransform(TRA)
        # probes[0].AOptPoly.SetUserTransform(icp)

        # print '0 icp'
        # print '1 icp',icp.GetMatrix()
        # print '2 icp',icp.GetLandmarkTransform().GetMatrix()
        # print '3 icp',probes[0].AOptPoly.GetUserTransform().GetMatrix()

        # TRAOLD = probes[0].AOptPoly.GetUserTransform()

    for r in renderer: r.Render()
    renderWindow.Render()

def ActorGlyphVisual(poly):
    """ 
    ugly class to create the color map with given opacity... to be largely ameliorated
    """

    apoly = vtkActor()
    mpoly = vtkPolyDataMapper()

    if False: # use spheres as glyph

        glyph = vtkGlyph3D()
        sphere = vtkSphereSource()
        sphere.SetRadius(1.0)
        sphere.SetThetaResolution(4)
        sphere.SetPhiResolution(4)
        if VTK_MAJOR_VERSION <= 5:
            glyph.SetSource(sphere.GetOutput())
            glyph.SetInput(poly)
        else:
            glyph.SetSourceConnection(sphere.GetOutputPort());
            glyph.SetInputData(poly)

        mpoly.SetInputConnection(glyph.GetOutputPort())

    else:

        # if True:
        if False:
            silhouette = vtkPolyDataSilhouette()
            if VTK_MAJOR_VERSION <= 5:
                silhouette.SetInput(poly)
            else:
                silhouette.SetInputData(poly)
            silhouette.SetCamera(renderer[0].GetActiveCamera())
            silhouette.SetEnableFeatureAngle(0)
            mpoly.SetInputConnection(silhouette.GetOutputPort())
            
        else:
            mpoly.ScalarVisibilityOff()

            # apoly.GetProperty().SetRepresentationToWireframe()
            # apoly.GetProperty().SetRepresentationToSurface()
            if VTK_MAJOR_VERSION <= 5:
                mpoly.SetInput(poly)
            else:
                mpoly.SetInputData(poly)
                # mpoly.SetInputConnection(silhouette.GetOutputPort())

    apoly.SetMapper(mpoly)
    apoly.GetProperty().SetColor(1,1,1)
    # apoly.GetProperty().SetOpacity(.2)

    return apoly,mpoly

def TransformPolyActor(actor,translate,rotate,scale):
    """ 
    translate-rotate-scale an actor and returns both the transformed actor and a transformed polydata
    T-R-S work in absolute terms
    """

    TN = vtkTransform()
    TN.Translate( translate )
    # TN.RotateX(rotate[0]); TN.RotateY(rotate[1]); TN.RotateZ(rotate[2])
    TN.Scale( scale )

    TN.Update()

    actor.SetUserTransform( TN )

    ptrans = vtkTransformPolyDataFilter()

    ptrans.SetTransform(TN)

    if VTK_MAJOR_VERSION<=5:
        ptrans.SetInput(actor.GetMapper().GetInput())
    else:
        ptrans.SetInputData(actor.GetMapper().GetInput())
    ptrans.Update()

    return ptrans.GetOutput()

    """

    poly = actor.GetMapper().GetInput()

    # first apply rotation abount center
    T1 = vtkTransform()
    T2 = vtkTransform()
    T3 = vtkTransform()
    T4 = vtkTransform()
    T5 = vtkTransform()

    ctr = CenterOfPoly(poly)

    T1.Translate( [-ctr[0], -ctr[1], -ctr[2]] )
    T2.RotateX(rotate[0]); T2.RotateY(rotate[1]); T2.RotateZ(rotate[2])
    T3.Scale( scale )
    T4.Translate( [+ctr[0], +ctr[1], +ctr[2]] )
    T5.Translate( translate )

    T2.Concatenate(T1)
    T3.Concatenate(T2)
    T4.Concatenate(T3)
    T5.Concatenate(T4)

    actor.SetUserTransform(T5)

    transform = vtkTransform()
    transform.SetMatrix( actor.GetMatrix() )
    ptrans = vtkTransformPolyDataFilter()
    ptrans.SetTransform(transform)
    if VTK_MAJOR_VERSION<=5:
        ptrans.SetInput(poly)
    else:
        ptrans.SetInputData(poly)
    ptrans.Update();

    if VTK_MAJOR_VERSION<=5:
        actor.GetMapper().SetInput(ptrans.GetOutput())
    else:
        # actor.GetMapper().SetInputData(ptrans.GetOutput())
        actor.GetMapper().SetInputConnection(ptrans.GetOutputPort())

    # return ptrans.GetOutput()
    return actor.GetMapper().GetInput()
    """

def TransformPoly(poly,translate,rotate,scale):
    """ 
    translate-rotate-scale a polydata
    """

    # first apply rotation abount center
    T1 = vtkTransform()
    T2 = vtkTransform()
    T3 = vtkTransform()
    T4 = vtkTransform()
    T5 = vtkTransform()

    ptrans = vtkTransformPolyDataFilter()

    if VTK_MAJOR_VERSION<=5:
        ptrans.SetInput(poly)
    else:
        ptrans.SetInputData(poly)

    ctr = CenterOfPoly(poly)

    T1.Translate( [-ctr[0], -ctr[1], -ctr[2]] )
    T2.RotateX(rotate[0]); T2.RotateY(rotate[1]); T2.RotateZ(rotate[2])
    T3.Scale( scale )
    T4.Translate( [+ctr[0], +ctr[1], +ctr[2]] )
    T5.Translate( translate )

    T2.Concatenate(T1)
    T3.Concatenate(T2)
    T4.Concatenate(T3)
    T5.Concatenate(T4)

    ptrans.SetTransform(T5)
    ptrans.Update()

    return ptrans.GetOutput()

def CenterOfPoly(poly):
    """ 
    return center of mass of polydata
    """
    centerOfMassFilter = vtkCenterOfMass()
    if VTK_MAJOR_VERSION <= 5:
        centerOfMassFilter.SetInput(poly)
    else:
        centerOfMassFilter.SetInputData(poly)
    centerOfMassFilter.SetUseScalarsAsWeights(False)
    centerOfMassFilter.Update()
    return centerOfMassFilter.GetCenter()

def DeepCopyActor(actor):
    newactor = vtkActor()
    newmapper = vtkPolyDataMapper()
    newpoly = vtkPolyData()

    newpoly.DeepCopy(actor.GetMapper().GetInput())
    if VTK_MAJOR_VERSION <= 5:
        newmapper.SetInput(newpoly)
    else:
        newmapper.SetInputData(newpoly)
    newactor.SetMapper(newmapper)

    return newactor


def StencilImg(img,imgOrigin,imgSpacing,poly,translate,rotate,scale):
    """ 
    stencil image inside polydata
    """

    pdtrans = TransformPoly(poly,translate,rotate,scale)

    polyImgSten = vtkPolyDataToImageStencil()
    polyImgSten.SetOutputOrigin(imgOrigin)
    polyImgSten.SetOutputSpacing(imgSpacing)

    if VTK_MAJOR_VERSION<=5:
        polyImgSten.SetInput( pdtrans ) 
    else:
        polyImgSten.SetInputData( pdtrans )

    polyImgSten.SetTolerance(1.0E-50)
    polyImgSten.Update()

    imsten = vtkImageStencil()

    if VTK_MAJOR_VERSION<=5:
        imsten.SetInput( img )
        imsten.SetStencil(polyImgSten.GetOutput())
    else:
        imsten.SetInputData( img )
        imsten.SetStencilConnection(polyImgSten.GetOutputPort())

    imsten.SetBackgroundValue(-100)
    imsten.Update()	

    return imsten.GetOutput()

def ReadImage(fl,downsamplelevel=1):

    # check if files exist
    if not os.path.isfile(fl):
        print 'image file :',fl,' is irregular or doesnt exist ! '
        sys.exit(1)

    # Read the file
    froot,fext = os.path.splitext(fl)

    if fext[0:4] == '.vti':
        reader = vtkXMLImageDataReader()

    elif fext[0:4] =='.mhd' or fext[0:4] =='.mha':
        reader = vtkMetaImageReader()

    else:
        print 'input img format unrecognized',fext
        sys.exit(1)

    reader.SetFileName(fl)
    reader.Update()

    if downsamplelevel>1:
        reader = downsampleImageFilter(reader,downsamplelevel)

    return reader.GetOutput()

def WriteImage(fl,img):

    # Write the file
    froot,fext = os.path.splitext(fl)

    if fext[0:4] == '.vti':
        writer = vtkXMLImageDataWriter()

    elif fext[0:4] =='.mhd' or fext[0:4] =='.mha':
        writer = vtkMetaImageWriter()

    else:
        print 'input img format unrecognized',fext
        sys.exit(1)

    writer.SetFileName(fl)

    if VTK_MAJOR_VERSION <= 5:
        writer.SetInputConnection( img.GetProducerPort() )
    else:
        writer.SetInputData( img )

    writer.Write()

def cleanUpHeart(scanImg,wlevels):

    """
    1. Threshold the input image
    2. take only the biggest part of it
    3. Transform it into a surface (remove inner structure)
    """

    scanThreshold = vtkThreshold()
    if VTK_MAJOR_VERSION <= 5:
        scanThreshold.SetInput( scanImg )
    else:
        scanThreshold.SetInputData( scanImg )
    # scanThreshold.SetInputArrayToProcess(0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes.SCALARS )
    scanThreshold.ThresholdBetween( wlevels[0], wlevels[1] )
    scanThreshold.Update()

    # breakdown into different objects
    geometryFilter = vtkGeometryFilter()
    # geometryFilter = vtkDataSetSurfaceFilter()
    if VTK_MAJOR_VERSION <= 5:
        geometryFilter.SetInput(scanThreshold.GetOutput())
    else:
        geometryFilter.SetInputData(scanThreshold.GetOutput())
    geometryFilter.Update()

    # connectivity filter to take only the biggest chunk of probe
    conFilt = vtkPolyDataConnectivityFilter()
    if vtk.VTK_MAJOR_VERSION<=5:
        conFilt.SetInput(geometryFilter.GetOutput())
    else:
        conFilt.SetInputData(geometryFilter.GetOutput())
    conFilt.SetExtractionModeToLargestRegion()
    conFilt.Update()

    print conFilt.GetOutput().GetNumberOfPoints()

    clean = vtkCleanPolyData()
    clean.SetInputData(conFilt.GetOutput())
    clean.Update()

    conFilt = vtkPolyData()
    conFilt.DeepCopy(clean.GetOutput())

    if SURFACE_PROBING:

        scanSurface = vtkDataSetSurfaceFilter()
        if VTK_MAJOR_VERSION <= 5:
            scanSurface.SetInput(conFilt)
        else:
            scanSurface.SetInputData(conFilt)
        scanSurface.Update()

        # scanSurface = fakeConvexHull(scanSurface.GetOutput())

        return scanSurface.GetOutput()

    else:

        return conFilt


if __name__ == '__main__':

    renderer= [vtkRenderer(), vtkRenderer(), vtkRenderer(), 
               vtkRenderer(), vtkRenderer(), vtkRenderer()]
    i1 = None
    i2 = None
    AScanSurface = None

    parser = argparse.ArgumentParser(description='Downsample image.')
    parser.add_argument('-i', '--input',   required=True,  help='input VTI/MHD/MHA file')
    parser.add_argument('-k', '--kprobes',   required=True,  help='probe folder')
    parser.add_argument('-v', '--view',   required=True,  help='final view [y/n]')
    parser.add_argument('-d', '--dlevel',   required=False,  help='downsampling filter level >= 1')
    parser.add_argument('-s', '--screensize',   required=False,  help='screen size (small, mid, large)')
    parser.add_argument('-p', '--postcheck_probe',   required=False,  help='postcheck_probe ')
    parser.add_argument('-t', '--postcheck_transform',   required=False,  help='postcheck_transform ')

    args = parser.parse_args()

    FILEIN = args.input
    FOLDERPROBES = args.kprobes

    if args.dlevel != None:
        DOWNSAMPLELEVEL = int(args.dlevel)
    else:
        DOWNSAMPLELEVEL = 4

    if args.screensize != None:
        SCREENSIZE = args.screensize
    else:
        SCREENSIZE = 'small'

    FINALVIEW = None
    if args.view != None:
        FINALVIEW = args.view
        if FINALVIEW != 'y' and FINALVIEW != 'n':
            print 'view can only be y or n'
            sys.exit(1)

    POSTCHECK_PROBE = None
    if args.postcheck_probe != None:
        POSTCHECK_PROBE = args.postcheck_probe

    POSTCHECK_TRANSFORM = None
    if args.postcheck_transform != None:
        POSTCHECK_TRANSFORM = args.postcheck_transform

    scanImg = ReadImage(FILEIN,DOWNSAMPLELEVEL)

    scanBounds = scanImg.GetBounds()
    scanLX = scanBounds[1]-scanBounds[0]
    scanLY = scanBounds[3]-scanBounds[2]
    scanLZ = scanBounds[5]-scanBounds[4]

    scanSurface = cleanUpHeart(scanImg,WLEVELS)


    ############################

    class Probe():
        pass

    probes = []

    # read the probe files
    f = open(os.path.join(FOLDERPROBES,'probelist.inp'),'r')

    for lines in f.readlines():

        lines = lines.split()
        line = lines[0]
        if line[0] == '#': continue

        probe = Probe()

        if POSTCHECK_PROBE != None:

            tr = open(POSTCHECK_TRANSFORM,'r')
            probe.M = vtkMatrix4x4()
            i = 0
            for l in tr.readlines():
                l = l.split()
                if i < 4:
                    probe.M.SetElement(i,0,float(l[0]))
                    probe.M.SetElement(i,1,float(l[1]))
                    probe.M.SetElement(i,2,float(l[2]))
                    probe.M.SetElement(i,3,float(l[3]))
                else:
                    probe.rr = [ float(l[0]), float(l[1]), float(l[2]) ] 
                    break
                i += 1

            probe.filename = POSTCHECK_PROBE
        else:

            # probe.filename = line[:len(line)-1]
            probe.filename = line


        if SURFACE_PROBING:
            # read the probe image and transform into a surface skin (polydata)
            # imgprobe = ReadImage(os.path.join(FOLDERPROBES,probe.filename))
            imgprobe = ReadImage(os.path.join(FOLDERPROBES,probe.filename), DOWNSAMPLELEVEL)

            threshold = vtkThreshold()
            if VTK_MAJOR_VERSION <= 5:
                threshold.SetInput( imgprobe )
            else:
                threshold.SetInputData( imgprobe )

            threshold.ThresholdBetween( PLEVELS[0], PLEVELS[1] )
            threshold.Update()

            surfaceSkinFilt = vtkDataSetSurfaceFilter()
            if VTK_MAJOR_VERSION <= 5:
                surfaceSkinFilt.SetInput(threshold.GetOutput())
            else:
                surfaceSkinFilt.SetInputData(threshold.GetOutput())
            surfaceSkinFilt.Update()

            probe.APoly, probe.MPoly = ActorGlyphVisual(surfaceSkinFilt.GetOutput())

        else:
            if True:
                imgprobe = ReadImage(os.path.join(FOLDERPROBES,probe.filename), DOWNSAMPLELEVEL)
                surfaceSkinFilt = cleanUpHeart(imgprobe,PLEVELS)
            else:
                read = vtkPolyDataReader()
                read.SetFileName( os.path.join(FOLDERPROBES,probe.filename) )
                read.Update()
                surfaceSkinFilt = read

            probe.APoly, probe.MPoly = ActorGlyphVisual(surfaceSkinFilt)

        probes.append(probe)


        if POSTCHECK_PROBE != None:
            break

    f.close()

    ###################################################
    renderer[0].SetBackground(0.0,0.0,0.0)
    renderer[1].SetBackground(0.0,0.0,0.0)
    renderer[2].SetBackground(0.0,0.0,0.0)
    renderer[3].SetBackground(0.0,0.0,0.0)
    renderer[4].SetBackground(0.0,0.0,0.0)
    renderer[5].SetBackground(0.0,0.0,0.0)

    #renderer[0].SetViewport(0.000, 0.0, 0.333, 1.0)
    #renderer[1].SetViewport(0.333, 0.0, 0.666, 1.0)
    #renderer[2].SetViewport(0.666, 0.0, 1.000, 1.0)

    renderer[0].SetViewport(0./3., 0., 1./3., .5)
    renderer[1].SetViewport(1./3., 0., 2./3., .5)
    renderer[2].SetViewport(2./3., 0., 3./3., .5)

    renderer[3].SetViewport(0./3., .5, 1./3., 1.)
    renderer[4].SetViewport(1./3., .5, 2./3., 1.)
    renderer[5].SetViewport(2./3., .5, 3./3., 1.)

    renderWindow = vtkRenderWindow()

    if SCREENSIZE == 'small':
        renderWindow.SetSize(600,200)
    elif SCREENSIZE == 'mid':
        renderWindow.SetSize(900,300)
    elif SCREENSIZE == 'large':
        renderWindow.SetSize(1500,500)

    for ren in renderer:
        renderWindow.AddRenderer( ren )

    # add volumetric images to renderes
    i1 = Imgview(scanImg, renderer[0], uint=False)
    # i2 = Imgview(scanSplatter.GetOutput(), renderer[1], uint=False)

    for i in xrange( 1, len(renderer) ):
        renderer[i].SetActiveCamera(renderer[0].GetActiveCamera())

    # renderer[0].ResetCamera()
    renderer[1].ResetCamera()

    iren = vtkRenderWindowInteractor()
    iren.SetRenderWindow(renderWindow)

    iren.AddObserver("KeyPressEvent", KeyPressEvent)

    AScanSurface = vtkActor()
    mp = vtkPolyDataMapper()

    verts1 = vtkCellArray()
    points1 = vtkPoints()
    target = vtkPolyData()
    # for i in xrange(scanSurface.GetOutput().GetNumberOfPoints()):
    # for i in xrange(0, scanSurface.GetOutput().GetNumberOfPoints(), 1):
    for i in xrange(0, scanSurface.GetNumberOfPoints(), 1):

        # pid = points1.InsertNextPoint ( scanSurface.GetOutput().GetPoints().GetPoint(i) )
        pid = points1.InsertNextPoint ( scanSurface.GetPoints().GetPoint(i) )
        verts1.InsertNextCell(1)
        verts1.InsertCellPoint(pid)

    target.SetPoints(points1)
    target.SetVerts(verts1)

    # mp.SetInputConnection(scanSurface.GetOutputPort())
    mp.SetInputData(target)
    AScanSurface.SetMapper(mp)

    # print scanSurface.GetOutput().GetNumberOfPoints()
    print 'scansurface no of pts:',scanSurface.GetNumberOfPoints()
    print '     target no of pts:',target.GetNumberOfPoints()

    renderer[1].AddActor(AScanSurface)
    renderer[1].ResetCamera()

    iren.Initialize()
    # iren.Start()

    ##############

    if POSTCHECK_PROBE:

        bestprobe = probes[0]
        bestprobe.AOptPoly.GetProperty().SetRepresentationToSurface()
        bestprobe.AOptPoly.GetProperty().SetPointSize(1)
        bestprobe.AOptPoly.GetProperty().SetColor(1,0,1)
        renderer[0].AddActor(bestprobe.AOptPoly)
        renderer[1].AddActor(bestprobe.AOptPoly)

        if POSTCHECK_TRANSFORM == None:
            print '\n\n Post-Check Transform file should be given....'
            sys.exit(1)

        TR = vtkTransform()
        TR.SetMatrix(probe.M)
        bestprobe.AOptPoly.SetUserTransform( TR )

    else:

        bestmetrics = 0
        for probe in probes:

            print 
            print '     Using probe :',probe.filename
            print 

            verts2 = vtkCellArray()
            points2 = vtkPoints()
            SOURCE = vtkPolyData()

            for i in xrange(0, probe.APoly.GetMapper().GetInput().GetNumberOfPoints(), 1):

                pid = points2.InsertNextPoint ( probe.APoly.GetMapper().GetInput().GetPoints().GetPoint(i) )
                verts2.InsertNextCell(1)
                verts2.InsertCellPoint(pid)

            SOURCE.SetPoints(points2)
            SOURCE.SetVerts(verts2)

            probe.AOptPoly, probe.MOptPoly = ActorGlyphVisual(poly = SOURCE)

            bestprobe = probe

            break

        print '\n Target heart filename :',FILEIN,'\n'

    isoTARGET = vtk.vtkMarchingCubes()
    if vtk.VTK_MAJOR_VERSION<=5:
        isoTARGET.SetInput(scanImg)
    else:
        isoTARGET.SetInputData(scanImg)
    isoTARGET.SetValue(0, MARCHINGLEVEL)
    isoTARGET.ComputeNormalsOn()
    isoTARGET.Update()

    isoPROBE = vtk.vtkMarchingCubes()
    if vtk.VTK_MAJOR_VERSION<=5:
        isoPROBE.SetInput(imgprobe)
    else:
        isoPROBE.SetInputData(imgprobe)
    isoPROBE.SetValue(0, MARCHINGLEVEL)
    isoPROBE.ComputeNormalsOn()
    isoPROBE.Update()

    TN = vtkTransform()
    TN.RotateX(10.)
    TN.RotateY(40.)
    TN.RotateZ(80.)
    TN.Update()
    TP = vtk.vtkTransformPolyDataFilter()
    TP.SetInputConnection(isoPROBE.GetOutputPort())
    TP.SetTransform(TN)
    TP.Update()

    FEAT1 = vtk.vtkFeatureEdges()
    FEAT1.SetInputData( isoTARGET.GetOutput() )
    FEAT1.BoundaryEdgesOn()
    # FEAT1.FeatureEdgesOff()
    FEAT1.SetFeatureAngle(1)
    FEAT1.ManifoldEdgesOff()
    FEAT1.ColoringOff()
    FEM1 = vtk.vtkPolyDataMapper()
    FEM1.SetInputConnection( FEAT1.GetOutputPort() )
    FEM1.SetResolveCoincidentTopologyToPolygonOffset()
    FEM1.ScalarVisibilityOff()
    FEA1 = vtk.vtkActor()
    FEA1.SetMapper(FEM1)
    FEA1.GetProperty().SetColor(1.,1.,1.)
    renderer[3].AddActor(FEA1)

    FEAT2 = vtk.vtkFeatureEdges()
    # FEAT2.SetInputData( isoPROBE.GetOutput() )
    FEAT2.SetInputData( TP.GetOutput() )
    # FEAT2.SetInputConnection( TP.GetOutputPort() )
    FEAT2.BoundaryEdgesOn()
    # FEAT2.FeatureEdgesOff()
    FEAT2.SetFeatureAngle(1)
    FEAT2.ManifoldEdgesOn()
    FEAT2.ColoringOff()
    FEM2 = vtk.vtkPolyDataMapper()
    FEM2.SetInputConnection( FEAT2.GetOutputPort() )
    FEM2.ScalarVisibilityOff()
    FEM2.SetResolveCoincidentTopologyToPolygonOffset()
    FEA2 = vtk.vtkActor()
    FEA2.SetMapper(FEM2)
    FEA2.GetProperty().SetColor(1.,0.,1.)
    renderer[4].AddActor(FEA2)

    icp = vtk.vtkIterativeClosestPointTransform()
    # icp.SetSource(isoPROBE.GetOutput())
    icp.SetSource(TP.GetOutput())
    icp.SetTarget(isoTARGET.GetOutput())
    icp.GetLandmarkTransform().SetModeToRigidBody()
    # icp.GetLandmarkTransform().SetModeToSimilarity()
    # icp.GetLandmarkTransform().SetModeToAffine()
    # icp.DebugOn()
    # icp.SetCheckMeanDistance(1)
    icp.SetMaximumNumberOfLandmarks(10000)
    icp.SetMaximumNumberOfIterations(100)
    # icp.SetMaximumMeanDistance(0.0001)
    icp.StartByMatchingCentroidsOn()
    icp.Modified()
    icp.Update()

    # icp.SetCheckMeanDistance(1)
    # icp.SetMaximumMeanDistance(0.001)
    # icp.SetMaximumNumberOfIterations(30)
    # icp.SetMaximumNumberOfLandmarks(50)

    P = vtk.vtkPolyData()
    P.DeepCopy( FEAT2.GetOutput() )
    M = vtk.vtkPolyDataMapper()
    M.SetInputData( P )
    M.ScalarVisibilityOff()
    M.SetResolveCoincidentTopologyToPolygonOffset()
    FEA2opt = vtk.vtkActor()
    FEA2opt.SetMapper(FEM2)
    FEA2opt.GetProperty().SetColor(1.,0.,1.)

    FEA2opt.SetUserTransform(icp)

    renderer[5].AddActor(FEA1)
    renderer[5].AddActor(FEA2opt)


    """
    # stencil the target image in ~20% larger space
    workactor = DeepCopyActor(bestprobe.AOptPoly)
    TransformPolyActor(workactor,[0,0,0],[0,0,0], [1.2, 1.2, 1.2] )
    stpoly = workactor.GetMapper().GetInput()

    # carve out the image within Probe_surface
    imsten = StencilImg(scanImg,scanImg.GetOrigin(),scanImg.GetSpacing(),stpoly,[0,0,0],[0,0,0],[1,1,1])

    # add image to renderer
    i3 = Imgview(imsten, renderer[2], uint=False)
    """

    if FINALVIEW == 'y':

        bestprobe.AOptPoly.GetProperty().SetColor(1,0,1)
        renderer[0].AddActor(bestprobe.AOptPoly)
        renderer[1].AddActor(bestprobe.AOptPoly)

        renderer[2].ResetCamera()

        for r in renderer: r.Render()
        renderWindow.Render()
        iren.Start()


