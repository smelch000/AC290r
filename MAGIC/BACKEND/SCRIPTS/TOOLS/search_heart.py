#!/usr/bin/env python

import math,sys,os,random,time
import argparse
from myvtk import *
import numpy as np
from downsample_image import *

# Threshold interval for detecting the heart
WLEVELS = 150, 900
# WLEVELS = 200, 900

# Threshold interval for probes
PLEVELS = 150, 900

# resolution for sampling T-R-S values
FINESSE = 10 # 20 

# number of optimization passes
NPASS = 3 

# how many neighboring voxels a point distributes
SPLATRADIUS = 0.01 

# ratio of valid probe points falling inside target image
VALIDITYPERC = 0.95 

# flag to represent the heart as :
#   True: a skin bject      (---> better if used via Delaunay3D)
#   False: a set of points falling within the Threshold window
# SURFACE_PROBING = True 
SURFACE_PROBING = False

# marking points of probes. Lowering this to 1 can degrade results ...
MASKPOINTS = 10

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

        rayCastFunction = vtkVolumeRayCastCompositeFunction()
        rayCastFunction.SetCompositeMethodToInterpolateFirst()
        self.volumeMapper = vtkVolumeRayCastMapper()
        self.volumeMapper.SetInputConnection(self.imageCast.GetOutputPort())
        # self.volumeMapper.SetInputConnection(self.v16.GetOutputPort())
        # self.volumeMapper.SetInputConnection(imagedata.GetOutputPort())
        self.volumeMapper.SetVolumeRayCastFunction(rayCastFunction)
        self.volumeMapper.SetSampleDistance(0.5)

        self.volume = vtk.vtkVolume()
        self.volume.SetMapper(self.volumeMapper)
        self.volume.SetProperty(self.volumeProperty)
        self.volume.PickableOff()	       
	        
        renderer.AddViewProp(self.volume)

def KeyPressEvent(caller,eventId):

    key = iren.GetKeySym()

    if key=='Left':
        if i1 != None:
            renderer[0].AddViewProp(i1.volume)
        for probe in probes:
            # probe.AOptPoly.GetProperty().SetRepresentationToPoints()
            probe.AOptPoly.VisibilityOff()

    elif key=='Right':
        if i1 != None:
            renderer[0].RemoveViewProp(i1.volume)
        for probe in probes:
            probe.AOptPoly.VisibilityOn()
            # probe.AOptPoly.GetProperty().SetRepresentationToWireframe()
            probe.AOptPoly.GetProperty().SetRepresentationToSurface()

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

def RotatePolyActor(actor,rotate):

    ctr = CenterOfPoly(actor.GetMapper().GetInput())

    T1 = vtkTransform()
    T1.Translate( [-ctr[0],-ctr[1],-ctr[2]] )

    T2 = vtkTransform()
    T2.RotateX(rotate[0]); T2.RotateY(rotate[1]); T2.RotateZ(rotate[2])

    T3 = vtkTransform()
    T3.Translate( ctr )

    T2.Concatenate(T1)
    T3.Concatenate(T2)
    T3.Update()

    actor.SetUserTransform( T3 )

    ptrans = vtkTransformPolyDataFilter()
    ptrans.SetTransform(T3)

    if VTK_MAJOR_VERSION<=5:
        ptrans.SetInput(actor.GetMapper().GetInput())
    else:
        ptrans.SetInputData(actor.GetMapper().GetInput())
    ptrans.Update()

    return ptrans.GetOutput()

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

def OptimizeTranslation(img,actor,window,finesse,metrics):
    """ 
    optimize over translations by accumulating over previous T-R-S
    """

    workactor = DeepCopyActor(actor)

    xarr = []; yarr = []; zarr = []
    x = window[0]; y = window[2]; z = window[4]
    dx = (window[1]-window[0])/(1.*finesse)
    dy = (window[3]-window[2])/(1.*finesse)
    dz = (window[5]-window[4])/(1.*finesse)
    while x < window[1]:
        xarr.append(x)
        x += dx
    while y < window[3]:
        yarr.append(y)
        y += dy
    while z < window[5]:
        zarr.append(z)
        z += dz

    probe = vtkProbeFilter()
    npts = workactor.GetMapper().GetInput().GetNumberOfPoints()

    mask = vtkMaskPoints()

    if VTK_MAJOR_VERSION <= 5:
        mask.SetInput( workactor.GetMapper().GetInput() )
    else:
        mask.SetInputData( workactor.GetMapper().GetInput() )

    mask.SetOnRatio(MASKPOINTS)
    mask.SetMaximumNumberOfPoints(10000)
    mask.RandomModeOff()
    mask.Update()

    npts = mask.GetOutput().GetNumberOfPoints()

    ii = -1
    for itr in xarr:
        for jtr in yarr:
            for ktr in zarr:

                ii += 1

                dr = np.random.rand(3)
                dr[0] -= .5; dr[1] -= .5; dr[2] -= .5
                dr *= dx

                translate = [itr+dr[0], jtr+dr[1], +ktr+dr[2]]

                # TransformPolyActor(actor,translate,rotateFix,scaleFix)
                # poly = TransformPolyActor(actor,translate,[0,0,0],[1,1,1])
                poly = TransformPolyActor(workactor,translate,[0,0,0],[1,1,1])

                ctr = CenterOfPoly(poly)

                mask = vtkMaskPoints()
                if VTK_MAJOR_VERSION <= 5:
                    mask.SetInput( poly )
                else:
                    mask.SetInputData( poly )
                mask.SetOnRatio(MASKPOINTS)
                mask.SetMaximumNumberOfPoints(10000)
                mask.RandomModeOff()
                mask.Update()

                # renderer[1].ResetCamera()
                # for r in renderer: r.Render()
                # renderWindow.Render()

                if ctr[0] < spltBounds[0] or ctr[0] > spltBounds[1] or \
                   ctr[1] < spltBounds[2] or ctr[1] > spltBounds[3] or \
                   ctr[2] < spltBounds[4] or ctr[2] > spltBounds[5]: 
                   continue

                if VTK_MAJOR_VERSION <= 5:
                    probe.SetSource(img.GetOutput())
                    # probe.SetInput(poly)
                    probe.SetInput(mask.GetOutput())
                else:
                    probe.SetSourceConnection(img.GetOutputPort())
                    # probe.SetInputData(poly)
                    probe.SetInputData(mask.GetOutput())

                probe.Update()

                new = HowULikeIt(probe,img)

                if new > metrics:
                    metrics = new
                    
                    T = vtkTransform()
                    T.SetMatrix( workactor.GetMatrix() )
                    actor.SetUserTransform(T)

                    renderer[1].ResetCamera()
                    for r in renderer: r.Render()
                    renderWindow.Render()

    print 'Final Transl score:',metrics

    return metrics

def OptimizeRotation(img,actor,window,finesse,metrics):
    """ 
    optimize over rotations by accumulating over previous T-R-S
    """

    workactor = DeepCopyActor(actor)

    xarr = []; yarr = []; zarr = []
    x = window[0]; y = window[2]; z = window[4]
    dx = (window[1]-window[0])/(1.*finesse)
    dy = (window[3]-window[2])/(1.*finesse)
    dz = (window[5]-window[4])/(1.*finesse)
    while x < window[1]:
        xarr.append(x)
        x += dx
    while y < window[3]:
        yarr.append(y)
        y += dy
    while z < window[5]:
        zarr.append(z)
        z += dz

    probe = vtkProbeFilter()

    npts = workactor.GetMapper().GetInput().GetNumberOfPoints()

    mask = vtkMaskPoints()
    if VTK_MAJOR_VERSION <= 5:
        mask.SetInput( workactor.GetMapper().GetInput() )
    else:
        mask.SetInputData( workactor.GetMapper().GetInput() )
    mask.SetOnRatio(MASKPOINTS)
    mask.SetMaximumNumberOfPoints(10000)
    mask.RandomModeOff()
    mask.Update()

    npts = mask.GetOutput().GetNumberOfPoints()

    ii = -1
    for itr in xarr:
        for jtr in yarr:
            for ktr in zarr:
                ii += 1

                dr = np.random.rand(3)
                dr[0] -= .5; dr[1] -= .5; dr[2] -= .5
                dr *= dx

                rotate = [itr+dr[0], jtr+dr[1], ktr+dr[2]]

                # if ii % 23 == 0:
                #     for r in renderer: r.Render()
                #     renderWindow.Render()

                poly = RotatePolyActor(workactor,rotate)

                if VTK_MAJOR_VERSION <= 5:
                    mask.SetInput( poly )
                else:
                    mask.SetInputData( poly )
                mask.SetOnRatio(MASKPOINTS)
                mask.SetMaximumNumberOfPoints(10000)
                mask.RandomModeOff()
                mask.Update()

                if VTK_MAJOR_VERSION <= 5:
                    probe.SetSource(img.GetOutput())
                    # probe.SetInput(poly)
                    probe.SetInput(mask.GetOutput())
                else:
                    probe.SetSourceConnection(img.GetOutputPort())
                    # probe.SetInputData(poly)
                    probe.SetInputData(mask.GetOutput())
                probe.Update()

                new = HowULikeIt(probe,img)

                if new > metrics:
                    # print rotate,'Rot score:', new
                    metrics = new

                    T = vtkTransform()
                    T.SetMatrix( workactor.GetMatrix() )
                    actor.SetUserTransform(T)

                    renderer[1].ResetCamera()
                    for r in renderer: r.Render()
                    renderWindow.Render()

    print 'Final Rotate score:',metrics

    return metrics

def OptimizeScaling(img,actor,window,finesse,metrics):
    """ 
    optimize over scalings by accumulating over previous T-R-S
    Warning: contrary to T-R, scaling is not a relative but absolute transformation  !!
    """

    workactor = DeepCopyActor(actor)

    xarr = []; yarr = []; zarr = []
    x = window[0]; y = window[2]; z = window[4]
    dx = (window[1]-window[0])/(1.*finesse)
    dy = (window[3]-window[2])/(1.*finesse)
    dz = (window[5]-window[4])/(1.*finesse)
    while x < window[1]:
        xarr.append(x)
        x += dx
    while y < window[3]:
        yarr.append(y)
        y += dy
    while z < window[5]:
        zarr.append(z)
        z += dz

    probe = vtkProbeFilter()

    npts = workactor.GetMapper().GetInput().GetNumberOfPoints()

    mask = vtkMaskPoints()
    if VTK_MAJOR_VERSION <= 5:
        mask.SetInput( workactor.GetMapper().GetInput() )
    else:
        mask.SetInputData( workactor.GetMapper().GetInput() )
    mask.SetOnRatio(MASKPOINTS)
    mask.SetMaximumNumberOfPoints(10000)
    mask.RandomModeOff()
    mask.Update()

    npts = mask.GetOutput().GetNumberOfPoints()

    ii = -1
    for itr in xarr:
        for jtr in [itr-2, itr, itr+2]:
        # for jtr in yarr:
            for ktr in [itr-2, itr, itr+2]:
            # for ktr in zarr:

                dr = np.random.rand(3)
                dr[0] -= .5; dr[1] -= .5; dr[2] -= .5
                dr *= dx

                # scale = [1. + itr/100., 1. + jtr/100., 1. + ktr/100.]
                scale = [1. + (itr + dr[0])/100., 1. + (jtr + dr[1])/100., 1. + (ktr+dr[2])/100.]

                poly = TransformPolyActor(workactor,[0,0,0],[0,0,0],scale)

                if VTK_MAJOR_VERSION <= 5:
                    mask.SetInput( poly )
                else:
                    mask.SetInputData( poly )
                mask.SetOnRatio(MASKPOINTS)
                mask.SetMaximumNumberOfPoints(10000)
                mask.RandomModeOff()
                mask.Update()

                if VTK_MAJOR_VERSION <= 5:
                    probe.SetSource(img.GetOutput())
                    # probe.SetInput(poly)
                    probe.SetInput(mask.GetOutput())
                else:
                    probe.SetSourceConnection(img.GetOutputPort())
                    # probe.SetInputData(poly)
                    probe.SetInputData(mask.GetOutput())
                probe.Update()

                new = HowULikeIt(probe,img)

                if new > metrics:
                    print 'Scale score:',new
                    metrics = new

                    T = vtkTransform()
                    T.SetMatrix( workactor.GetMatrix() )
                    actor.SetUserTransform(T)

                    renderer[1].ResetCamera()
                    for r in renderer: r.Render()
                    renderWindow.Render()

    print 'Final  Scale score:',metrics

    return metrics

def HowULikeIt(probe,img):

    npts = probe.GetOutput().GetNumberOfPoints()

    valid_ids = probe.GetValidPoints()
    nvalid = valid_ids.GetNumberOfTuples()

    if nvalid < VALIDITYPERC * npts or nvalid == 0: return -1

    valid_loc = 0
    new = 0
    arr = np.asarray( numpy_support.vtk_to_numpy(probe.GetOutput().GetPointData().GetScalars()), dtype=np.float)

    for i in xrange(npts):
        # print valid_ids.GetTuple1(valid_loc)
        if valid_ids.GetTuple1(valid_loc) == i:
            # print 'val:',arr[i]
            new += arr[i]
            valid_loc += 1

    new /= valid_loc

    return new

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
        conFilt.SetInputConnection(geometryFilter.GetOutputPort())

    conFilt.SetExtractionModeToLargestRegion()
    conFilt.Update()

    scanSurface = vtkDataSetSurfaceFilter()
    if VTK_MAJOR_VERSION <= 5:
        scanSurface.SetInput(conFilt.GetOutput())
    else:
        scanSurface.SetInputData(conFilt.GetOutput())
    scanSurface.Update()

    return scanSurface


if __name__ == '__main__':

    renderer= [vtkRenderer(), vtkRenderer(), vtkRenderer()]
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

    scanSplatter = vtkGaussianSplatter()
    if VTK_MAJOR_VERSION <= 5:
        scanSplatter.SetInput(scanSurface.GetOutput())
    else:
        scanSplatter.SetInputData(scanSurface.GetOutput())
 
    scanSplatter.SetModelBounds( scanImg.GetBounds() )
    scanSplatter.SetSampleDimensions( scanImg.GetDimensions() )
    scanSplatter.SetRadius( SPLATRADIUS )
    scanSplatter.SetScaleFactor(.5)
    scanSplatter.Update()

    spltBounds = scanSplatter.GetOutput().GetBounds()

    splatterLX = spltBounds[1] - spltBounds[0]
    splatterLY = spltBounds[3] - spltBounds[2]
    splatterLZ = spltBounds[5] - spltBounds[4]

    print 'Scan     Bounds  : [%.2f:%.2f] [%.2f:%.2f] [%.2f:%.2f] ' % tuple(scanBounds)
    print 'Splatter Bounds  : [%.2f:%.2f] [%.2f:%.2f] [%.2f:%.2f] ' % tuple(spltBounds)
    print 'Scan     Size    : %.2f %.2f %.2f ' % (scanLX,scanLY,scanLZ)
    print 'Splatter Size    : %.2f %.2f %.2f ' % (splatterLX,splatterLY,splatterLZ)
    print 'Scan     Dims    : ', scanImg.GetDimensions()
    print 'Splatter Dims    : ', scanSplatter.GetOutput().GetDimensions()
    print 'Scan     Spacing : %.2f %.2f %.2f ' % tuple(scanImg.GetSpacing())
    print 'Splatter Spacing : %.2f %.2f %.2f ' % tuple(scanSplatter.GetOutput().GetSpacing())

    # print 'Scan Extent : ', scanImg.GetExtent()

    """
    writer = vtkXMLImageDataWriter()
    writer.SetFileName( 'splat.vti' )
    if VTK_MAJOR_VERSION <= 5:
        writer.SetInput( scanSplatter.GetOutput() )
    else:
        writer.SetInputConnection( scanSplatter.GetOutputPort() )
    writer.Update()
    """

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
            imgprobe = ReadImage(os.path.join(FOLDERPROBES,probe.filename))

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
        else:
            read = vtkPolyDataReader()
            read.SetFileName( os.path.join(FOLDERPROBES,probe.filename) )
            read.Update()
            surfaceSkinFilt = read

        probe.AOptPoly, probe.MOptPoly = ActorGlyphVisual(poly = surfaceSkinFilt.GetOutput())

        print 'Probe : ',probe.filename,
        print ' Number of Points: %d' % surfaceSkinFilt.GetOutput().GetNumberOfPoints(),

        if POSTCHECK_PROBE == None:

            # Randomize Rotation <------------------
            rr = np.random.rand(3)
            rr[0] -= .5; rr[1] -= .5; rr[2] -= .5
            rr *= 360.

            probe.rr = rr

        print ' ==> Randomization Angles: %.3f %.3f %.3f'%tuple(probe.rr)

        # At first reset the probe at the center of image
        pctr = CenterOfPoly(probe.AOptPoly.GetMapper().GetInput())
        bctr = [spltBounds[0] + .5*splatterLX, spltBounds[2] + .5*splatterLY, spltBounds[4] + .5*splatterLZ]

        poly = TransformPolyActor(probe.AOptPoly, \
                            [bctr[0]-pctr[0],bctr[1]-pctr[1],bctr[2]-pctr[2]], \
                            [0,0,0], \
                            [1.1,1.1,1.1])
                            # [1.05,1.05,1.05]) # start with a slightly fatter probe
                            # [1,1,1])

        # print '     New CoM %8.2f %8.2f %8.2f' % CenterOfPoly(poly),\
        #         '  =   Box Ctr %8.2f %8.2f %8.2f' % tuple(bctr)

        probe.AOptPoly,probe.MOptPoly = ActorGlyphVisual(poly = poly)

        probes.append(probe)

        if POSTCHECK_PROBE != None:
            break

    f.close()

    ###################################################
    renderer[0].SetBackground(0.0,0.0,0.0)
    renderer[1].SetBackground(0.0,0.0,0.0)
    renderer[2].SetBackground(0.0,0.0,0.0)

    renderer[0].SetViewport(0.000, 0.0, 0.333, 1.0)
    renderer[1].SetViewport(0.333, 0.0, 0.666, 1.0)
    renderer[2].SetViewport(0.666, 0.0, 1.000, 1.0)

    renderWindow = vtkRenderWindow()

    if SCREENSIZE == 'small':
        renderWindow.SetSize(600,200)
    elif SCREENSIZE == 'mid':
        renderWindow.SetSize(900,300)
    elif SCREENSIZE == 'large':
        renderWindow.SetSize(1500,500)

    renderWindow.AddRenderer(renderer[0])
    renderWindow.AddRenderer(renderer[1])
    renderWindow.AddRenderer(renderer[2])

    # add volumetric images to renderes
    i1 = Imgview(scanImg, renderer[0], uint=False)
    # i2 = Imgview(scanSplatter.GetOutput(), renderer[1], uint=False)

    renderer[1].SetActiveCamera(renderer[0].GetActiveCamera())
    renderer[2].SetActiveCamera(renderer[0].GetActiveCamera())
    # renderer[0].ResetCamera()
    renderer[1].ResetCamera()

    iren = vtkRenderWindowInteractor()
    iren.SetRenderWindow(renderWindow)

    iren.AddObserver("KeyPressEvent", KeyPressEvent)

    AScanSurface = vtkActor()
    mp = vtkPolyDataMapper()
    mp.SetInputConnection(scanSurface.GetOutputPort())
    renderer[1].AddActor(AScanSurface)
    AScanSurface.SetMapper(mp)

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

        f = open('metrics.dat','a')
        m = open('transform.dat','w')
        bestmetrics = 0
        for probe in probes:

            print '     Using probe :',probe.filename
            print >> f, '# probe :',probe.filename

            probe.AOptPoly.GetProperty().SetRepresentationToSurface()
            # probe.AOptPoly.GetProperty().SetRepresentationToPoints()
            probe.AOptPoly.GetProperty().SetPointSize(1)
            probe.AOptPoly.GetProperty().SetColor(1,1,1)

            # renderer[0].AddActor(APoly)
            renderer[0].AddActor(probe.AOptPoly)
            renderer[1].AddActor(probe.AOptPoly)

            # optimization loop
            metrics = 0
            for ipass in xrange(NPASS):

                ctr = CenterOfPoly(probe.AOptPoly.GetMapper().GetInput())

                finesse = FINESSE

                WindowT = np.array([ spltBounds[0]-ctr[0], spltBounds[1]-ctr[0], \
                                    spltBounds[2]-ctr[1], spltBounds[3]-ctr[1], \
                                    spltBounds[4]-ctr[2], spltBounds[5]-ctr[2]])

                WindowR = np.array([-180., 180., -180., 180., -180., 180.])

                WindowS = np.array([ -10.,  10.,  -10.,  10.,  -10.,  10.])

                if ipass == 1 and metrics == 0:
                    pass
                    # finesse = 2 * FINESSE

                elif ipass == NPASS-1:
                    WindowT /= 3.
                    WindowR /= 3.
                    WindowS /= 3.

                metrics = OptimizeTranslation(scanSplatter,probe.AOptPoly,WindowT,finesse,metrics)

                metrics = OptimizeRotation(scanSplatter,probe.AOptPoly,WindowR,finesse,metrics)

                if ipass >= 0:
                    metrics = OptimizeScaling(scanSplatter,probe.AOptPoly,WindowS,finesse,metrics)

                for r in renderer: r.Render()
                renderWindow.Render()

                print ipass,'/',NPASS,'                                                            Optimum metrics:',metrics
                if metrics > 0:
                    print >>f, ipass,metrics
                    f.flush()

            renderer[0].RemoveActor(probe.AOptPoly)
            renderer[1].RemoveActor(probe.AOptPoly)

            if bestmetrics < metrics:
                bestmetrics = metrics
                bestprobe = probe

        print '\n Target heart filename :',FILEIN
        print 'Best fitting heart filename :',bestprobe.filename,'\n'

        print >>f, '# Best fitting heart filename :',bestprobe.filename
        print >>f, '# Best fitting heart metrics :',bestmetrics
        print >>f
        f.flush()

        for i in xrange(4): print >>m, bestprobe.AOptPoly.GetMatrix().GetElement(0,i),
        print >> m
        for i in xrange(4): print >>m, bestprobe.AOptPoly.GetMatrix().GetElement(1,i),
        print >> m
        for i in xrange(4): print >>m, bestprobe.AOptPoly.GetMatrix().GetElement(2,i),
        print >> m
        for i in xrange(4): print >>m, bestprobe.AOptPoly.GetMatrix().GetElement(3,i),
        print >> m
        print >>m, bestprobe.rr[0], bestprobe.rr[1], bestprobe.rr[2]
        m.flush()

    # stencil the target image in ~20% larger space
    workactor = DeepCopyActor(bestprobe.AOptPoly)
    TransformPolyActor(workactor,[0,0,0],[0,0,0], [1.2, 1.2, 1.2] )
    stpoly = workactor.GetMapper().GetInput()

    # carve out the image within Probe_surface
    imsten = StencilImg(scanImg,scanImg.GetOrigin(),scanImg.GetSpacing(),stpoly,[0,0,0],[0,0,0],[1,1,1])

    # add image to renderer
    i3 = Imgview(imsten, renderer[2], uint=False)

    if FINALVIEW == 'y':

        bestprobe.AOptPoly.GetProperty().SetColor(1,0,1)
        renderer[0].AddActor(bestprobe.AOptPoly)
        renderer[1].AddActor(bestprobe.AOptPoly)

        renderer[2].ResetCamera()

        for r in renderer: r.Render()
        renderWindow.Render()
        iren.Start()


