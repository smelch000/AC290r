#!/usr/bin/env python

from myvtk import *
import numpy as np
import numpy.ma as ma
import argparse,sys

# hardcoded values

HEARTWINDOWLEVELS = [100, 500]

POTATOFILE = None
# POTATOFILE = 'sphere.stl'
# POTATOFILE = 'ellipsoidSmall.stl'
# POTATOFILE = 'ellipsoid.stl'
# POTATOFILE = 'ellipsoidProlate.stl'

def TransformPoly(poly,translate,rotate,scale):

    trans = vtk.vtkTransform()
    trans.Translate(translate)
    trans.RotateX(rotate[0]); trans.RotateY(rotate[1]); trans.RotateZ(rotate[2])
    trans.Scale(scale)

    pdtrans = vtk.vtkTransformPolyDataFilter()
    pdtrans.SetTransform(trans)
    if vtk.VTK_MAJOR_VERSION<=5:
        pdtrans.SetInput(poly)
    else:
        pdtrans.SetInputData(poly)
    pdtrans.Update()

    return pdtrans

def OptimizeTranslation(img,polyelli,rotateFix,scaleFix,WindowT):

    metrics = 0
    translateOpt = [0,0,0]

    polyelli.ComputeBounds()
    bnds = polyelli.GetBounds()
    oX,oY,oZ = .5*(bnds[0] + bnds[1]), .5*(bnds[2] + bnds[3]), .5*(bnds[4] + bnds[5])

    # print '.....Elli Orig',oX,oY,oZ,'Img Ext',imgExtent,'Img Bnds',imgBounds

    # x1 = int(WindowT[0] - oX); x2 = int(WindowT[1] - oX)
    # y1 = int(WindowT[2] - oY); y2 = int(WindowT[3] - oY)
    # z1 = int(WindowT[4] - oZ); z2 = int(WindowT[5] - oZ)

    x1 = WindowT[0]; x2 = WindowT[1]
    y1 = WindowT[2]; y2 = WindowT[3]
    z1 = WindowT[4]; z2 = WindowT[5]

    dx = int( abs(x2-x1)/10. )
    dy = int( abs(y2-y1)/10. )
    dz = int( abs(z2-z1)/10. )

    for itr in xrange(x1, x2, dx):
        for jtr in xrange(y1, y2, dy):
            for ktr in xrange(z1, z2, dz):

                translate = [itr,jtr,ktr]

                imgout = StencilImg(img,imgOrigin,imgSpacing,polyelli,translate,rotateFix,scaleFix)

                arr = numpy_support.vtk_to_numpy(imgout.GetPointData().GetScalars()) 
                sub1 = ma.masked_where(arr  < HEARTWINDOWLEVELS[0], arr)
                sub2 = ma.masked_where(sub1 > HEARTWINDOWLEVELS[1], sub1)

                # new = sub2.sum()
                new = sub2.count()
                if new > metrics:
                    metrics = new
                    translateOpt = translate

    return translateOpt,metrics

def OptimizeRotation(img,polyelli,translateFix,scaleFix):

    metrics = 0
    rotateOpt = [0,0,0]
    for itr in xrange(0,90,15):
        for jtr in xrange(-90,90,30):
            for ktr in xrange(-90,90,30):

                rotate = [itr,jtr,ktr]

                imgout = StencilImg(img,imgOrigin,imgSpacing,polyelli,translateFix,rotate,scaleFix)

                arr = numpy_support.vtk_to_numpy(imgout.GetPointData().GetScalars()) 
                sub1 = ma.masked_where(arr  < HEARTWINDOWLEVELS[0], arr)
                sub2 = ma.masked_where(sub1 > HEARTWINDOWLEVELS[1], sub1)

                new = sub2.count()
                if new > metrics:
                    metrics = new
                    rotateOpt = rotate
    return rotateOpt,metrics

def OptimizeScaling(img,polyelli,translateFix,rotateFix):

    f = open('scale.dat','w')

    metrics = +1.e10
    # metrics = 0
    scaleOpt = [1,1,1]

    arr = numpy_support.vtk_to_numpy(img.GetOutput().GetPointData().GetScalars()) 
    sub0 = ma.masked_where(False, arr)
    allinbox = sub0.count()

    for itr in xrange(-20,50+1,5):
    # for itr in xrange(-20,30,10):
    # for itr in xrange(-10,20,10):

        # for jtr in xrange(-30,30+1,10): # fully isotropic search
        # for jtr in xrange(itr-20,itr+20+1,5): # quasi-isotropic search
        for jtr in xrange(itr,itr+1,1):

            # for ktr in xrange(-30,30+1,10): # fully isotropic search
            # for ktr in xrange(jtr-20,jtr+20+1,5): # quasi-isotropic search
            for ktr in xrange(itr,itr+1,1):

                scale = [1. + itr/100., 1. + jtr/100., 1. + ktr/100.]

                imgout = StencilImg(img,imgOrigin,imgSpacing,polyelli,translateFix,rotateFix,scale)

                arr = numpy_support.vtk_to_numpy(imgout.GetPointData().GetScalars()) 

                # count all pixels in ellipse
                sub0 = ma.masked_where(arr  < -2000, arr)
                allinell = sub0.count()

                # count pixels in ellipse within HU window
                sub1 = ma.masked_where(arr  < HEARTWINDOWLEVELS[0], arr)
                sub2 = ma.masked_where(sub1 > HEARTWINDOWLEVELS[1], sub1)
                wininell = sub2.count()

                pdtrans = TransformPoly(polyelli,translateOpt,rotateOpt,scale)

                props = vtk.vtkMassProperties()
                props.SetInputConnection(pdtrans.GetOutputPort())
                props.Update()
                potatoVolume = props.GetVolume()
                potatoSurface = props.GetSurfaceArea()

                # new = potatoSurface * sub2.count() / (potatoVolume - sub2.count())

                new = (allinell - wininell) * 1. / wininell - 1

                iso = (scale[0] + scale[1] + scale[2])/3.
                print >> f, iso, new, wininell, allinell

                if new > 0.5:
                    metrics = new
                    return scale,metrics

                """
                if new < metrics: # search for Minimum
                # if new > metrics: # search for Maximum
                    metrics = new
                    scaleOpt = scale
                    # print '                                     best',scale
                """
    f.close()
    # sys.exit(1)

    return scaleOpt,metrics

def StencilImg(img,imgOrigin,imgSpacing,polyelli,translate,rotate,scale):

    pdtrans = TransformPoly(polyelli,translate,rotate,scale)

    polyImgSten = vtk.vtkPolyDataToImageStencil()
    polyImgSten.SetOutputOrigin(imgOrigin)
    polyImgSten.SetOutputSpacing(imgSpacing)
    if vtk.VTK_MAJOR_VERSION<=5:
        polyImgSten.SetInput( pdtrans.GetOutput() )
    else:
        polyImgSten.SetInputData( pdtrans.GetOutput() )
    polyImgSten.SetTolerance(1.0E-50)
    polyImgSten.Update()

    imsten = vtk.vtkImageStencil()
    if vtk.VTK_MAJOR_VERSION<=5:
        imsten.SetStencil(polyImgSten.GetOutput())
        imsten.SetInput( img.GetOutput() )
    else:
        imsten.SetStencilData(polyImgSten.GetOutput())
        imsten.SetInputData( img.GetOutput() )
    imsten.SetBackgroundValue(-9999) # new
    imsten.Update()	

    return imsten.GetOutput()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Carve Potato from Image.')
    parser.add_argument('-i', '--input',   required=True,  help='input VTI file')
    parser.add_argument('-p', '--potato',   required=True,  help='probe STL file')
    args = parser.parse_args()

    FILEIN = args.input
    POTATOFILE = args.potato

    print 'reading VTI img...'
    img = vtk.vtkXMLImageDataReader() # for a vti file
    # img.SetFileName( 'P2small.vti' )
    img.SetFileName( FILEIN )
    img.Update()
    imgOrigin =  img.GetOutput().GetOrigin()
    imgSpacing = img.GetOutput().GetSpacing()
    imgExtent = img.GetOutput().GetExtent()
    imgBounds = img.GetOutput().GetBounds()

    print 'reading STL template...'
    stlcase = vtk.vtkSTLReader() 
    stlcase.SetFileName( POTATOFILE )
    stlcase.Update()
    polyelli = stlcase.GetOutput()

    translateOpt = [0,0,0]
    rotateOpt = [0,0,0]
    scaleOpt = [1,1,1]

    # for ipass in range(1):
    for ipass in range(3):

        print 'optimizing pass...',ipass

        if ipass==0:
            WindowT = [int(imgBounds[0]) + 0, int(imgBounds[1]) - 0, # tolerance from edges
                       int(imgBounds[2]) + 0, int(imgBounds[3]) - 0, 
                       int(imgBounds[4]) + 0, int(imgBounds[5]) - 0]
        elif ipass==1:
            WindowT = [-50, +50, -50, +50, -50, +50]
        elif ipass==2:
            WindowT = [-20, +20, -20, +20, -20, +20]

        translateOpt, metrics = OptimizeTranslation(img,polyelli,rotateOpt,scaleOpt,WindowT)
        print ' best translational metrics:',metrics,'       Optimum:',translateOpt,' [pixel]'

        rotateOpt, metrics = OptimizeRotation(img,polyelli,translateOpt,scaleOpt)
        print ' best rotational metrics:',metrics,'          Optimum:',rotateOpt,' [deg]'

        #if ipass>=0:
        #    scaleOpt, metrics = OptimizeScaling(img,polyelli,translateOpt,rotateOpt)
        #    print ' best scaling metrics:',int(metrics),'              Optimum:',scaleOpt,' [%]'

        pdtrans = TransformPoly(polyelli,translateOpt,rotateOpt,scaleOpt)

        polyelli = pdtrans.GetOutput()
        translateOpt = [0,0,0]
        rotateOpt = [0,0,0]
        scaleOpt = [1,1,1]


    writer = vtk.vtkSTLWriter()
    if vtk.VTK_MAJOR_VERSION<=5:
        writer.SetInput(polyelli)
    else:
        writer.SetInputData(polyelli)
    writer.SetFileName('bestEllipsoid.stl')
    writer.SetFileTypeToASCII()
    writer.Write()

    pdtrans = TransformPoly(polyelli,translateOpt,rotateOpt,[1.06,1.06,1.06])
    writer = vtk.vtkSTLWriter()
    if vtk.VTK_MAJOR_VERSION<=5:
        writer.SetInput(pdtrans.GetOutput())
    else:
        writer.SetInputData(pdtrans.GetOutputData())
    writer.SetFileName('bestEllipsoidFat.stl')
    writer.SetFileTypeToASCII()
    writer.Write()


imsten = StencilImg(img,imgOrigin,imgSpacing,polyelli,translateOpt,rotateOpt,scaleOpt)

writer = vtkXMLImageDataWriter()
writer.SetFileName( 'carved.vti' )
if vtk.VTK_MAJOR_VERSION<=5:
    writer.SetInput( imsten )
else:
    writer.SetInputData( imsten )
writer.Update()
writer.Write()


