#!/usr/bin/env python

#
# script to fill a closed, water-tight surface with mesh points:
#   check if the points fall into a single inlet or *single* outlet box
#   use a scale factor to alter mesh resolution
#
# after this phase, wrap_by_wall.py should be run to wrap with wall nodes
#

import os, sys

sys.path.append( os.path.join(os.getenv("MOEBIUS_ROOT"), "BACKEND/SCRIPTS") )

TOLERANCE=0.

from myvtk import *
import argparse
import numpy as np
from scipy import ndimage

# import TOOLS.wrap_by_wall as wrap_by_wall
import TOOLS.mesh as mesh
import TOOLS.d3q19 as d3q19

def _scaleGeom(filtgeom,SCALE):
    t = vtkTransform()
    t.Scale(SCALE, SCALE, SCALE)
    tF = vtkTransformPolyDataFilter()
    tF.SetInputConnection(filtgeom.GetOutputPort())
    tF.SetTransform(t)
    tF.Update()
    return tF

def _translateGeom(filtgeom,Tx,Ty,Tz):
    t = vtkTransform()
    t.Translate(Tx, Ty, Tz)
    tF = vtkTransformPolyDataFilter()
    tF.SetInputConnection(filtgeom.GetOutputPort())
    tF.SetTransform(t)
    tF.Update()
    return tF

def readfiles(surface, clip=False, aclip=False, inlets=False, outlets=False):

    readsurf = vtkSTLReader(); readsurf.SetFileName( surface ); readsurf.SetMerging(True); readsurf.Update()
    filtsurf = readsurf

    filtins = []
    if inlets:
        for inlet in inlets:
            print 'Inlet from file:',inlet
            filtin = vtkSTLReader(); filtin.SetFileName( inlet ); filtin.Update()
            filtins.append(filtin)

    filtouts = []
    if outlets:
        for outlet in outlets:
            print 'Outlet from file:',outlet
            filtout = vtkSTLReader(); filtout.SetFileName( outlet ); filtout.Update()
            filtouts.append(filtout)

    filtclip = None
    if clip:
        print 'Clip region from file:',clip
        filtclip = vtkSTLReader(); filtclip.SetFileName( clip ); filtclip.Update()

    filtaclip = None
    if aclip:
        print 'Aclip region from file:',aclip
        filtaclip = vtkSTLReader(); filtaclip.SetFileName( aclip ); filtaclip.Update()

    return filtsurf, filtclip, filtaclip, filtins, filtouts

def positiveoctant(filtsurf, filtclip, filtaclip, filtins, filtouts, T):

    print('positiveoctant...T:',T)

    bnd = filtsurf.GetOutput().GetBounds()
    
    for filtin in filtins:
        bndi = filtin.GetOutput().GetBounds()
        bnd = (min(bnd[0],bndi[0]),max(bnd[1],bndi[1]),min(bnd[2],bndi[2]),max(bnd[3],bndi[3]),min(bnd[4],bndi[4]),max(bnd[5],bndi[5]))
 
    for filtout in filtouts:
        bndo = filtout.GetOutput().GetBounds()
        bnd = (min(bnd[0],bndo[0]),max(bnd[1],bndo[1]),min(bnd[2],bndo[2]),max(bnd[3],bndo[3]),min(bnd[4],bndo[4]),max(bnd[5],bndo[5]))
     
    print 'Initial origin: %8.3f %8.3f %8.3f '%(bnd[0],bnd[2],bnd[4]), 'Bounds: %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f '%bnd

    # now TRANSLATING
    if not T and (bnd[0]<1 or bnd[2]<1 or bnd[4]<1): # translation to positive octant
        print "Translating to positive octant..."
        filtsurf = _translateGeom(filtsurf, -bnd[0]+10, -bnd[2]+10, -bnd[4]+10)
        for i,filtin in enumerate(filtins): #enumerate needed because the following create a new object!
            filtins[i] = _translateGeom(filtin, -bnd[0]+10, -bnd[2]+10, -bnd[4]+10)
        for i,filtout in enumerate(filtouts):
            filtouts[i] = _translateGeom(filtout, -bnd[0]+10, -bnd[2]+10, -bnd[4]+10)
        if filtclip: filtclip = _translateGeom(filtclip, -bnd[0]+10, -bnd[2]+10, -bnd[4]+10)
        if filtaclip: filtaclip = _translateGeom(filtaclip, -bnd[0]+10, -bnd[2]+10, -bnd[4]+10)
    elif T:
        print "Applying user given translation",T
        print 'T',T
        [Tx,Ty,Tz]=T
        filtsurf = _translateGeom(filtsurf, Tx, Ty, Tz)
        for i,filtin in enumerate(filtins):
            filtins[i] = _translateGeom(filtin, Tx, Ty, Tz )
        for i,filtout in enumerate(filtouts):
            filtouts[i] = _translateGeom(filtout, Tx, Ty, Tz )
        if filtclip: filtclip = _translateGeom(filtclip, Tx, Ty, Tz )
        if filtaclip: filtaclip = _translateGeom(filtaclip, Tx, Ty, Tz )

    for f in [filtsurf, filtclip, filtaclip] + filtins + filtouts:
        if f!=None: f.Update()

    bnd = filtsurf.GetOutput().GetBounds()
    for filtin in filtins:
        bndi=filtin.GetOutput().GetBounds()
        bnd=(min(bnd[0],bndi[0]),max(bnd[1],bndi[1]),min(bnd[2],bndi[2]),max(bnd[3],bndi[3]),min(bnd[4],bndi[4]),max(bnd[5],bndi[5]))
      
    for filtout in filtouts:
        bndo=filtout.GetOutput().GetBounds()
        bnd=(min(bnd[0],bndo[0]),max(bnd[1],bndo[1]),min(bnd[2],bndo[2]),max(bnd[3],bndo[3]),min(bnd[4],bndo[4]),max(bnd[5],bndo[5]))

    if bnd[0]<1 : raise Exception(">>> X-bound too low: "+str(bnd[0]))
    if bnd[2]<1 : raise Exception(">>> Y-bound too low: "+str(bnd[2]))
    if bnd[4]<1 : raise Exception(">>> Z-bound too low: "+str(bnd[4]))
    print 'Translated origin: %8.3f %8.3f %8.3f '%(bnd[0],bnd[2],bnd[4]), 'Bounds: %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f '%bnd

    for f in [filtsurf, filtclip, filtaclip] + filtins + filtouts:
        if f!=None: f.Update()

    return filtsurf, filtclip, filtaclip, filtins, filtouts

def scalesolids(filtsurf,  filtclip,filtaclip, filtins, filtouts, SCALE=1.):

    print('scalesolids...scale:',SCALE)

    # now SCALING
    filtsurf = _scaleGeom(filtsurf, SCALE)

    for i,filtin in enumerate(filtins):
       filtins[i] = _scaleGeom(filtin, SCALE)

    for i,filtout in enumerate(filtouts):
       filtouts[i] = _scaleGeom(filtout, SCALE)

    if filtclip: filtclip = _scaleGeom(filtclip, SCALE)
    if filtaclip: filtaclip = _scaleGeom(filtaclip, SCALE)

    bnd = filtsurf.GetOutput().GetBounds()
    print 'Rescaled surface origin: %8.3f %8.3f %8.3f '%(bnd[0],bnd[2],bnd[4]), 'Bounds: %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f '%bnd

    for f in [filtsurf, filtclip, filtaclip] + filtins + filtouts:
        if f!=None: f.Update()

    return filtsurf, filtclip, filtaclip, filtins, filtouts
 
def buildinside(filtsurf, filtclip, filtaclip, filtins, filtouts, wrapbywall=False, openends=False, GRIDSPACING=1):

    print('_buildinside...', 'wrapbywall', wrapbywall, 'openends:',openends, 'gridspacing:', GRIDSPACING)

    bnd = filtsurf.GetOutput().GetBounds()
    for filtin in filtins:
        bndi=filtin.GetOutput().GetBounds()
        bnd=(min(bnd[0],bndi[0]),max(bnd[1],bndi[1]),min(bnd[2],bndi[2]),max(bnd[3],bndi[3]),min(bnd[4],bndi[4]),max(bnd[5],bndi[5]))
      
    for filtout in filtouts:
        bndo=filtout.GetOutput().GetBounds()
        bnd=(min(bnd[0],bndo[0]),max(bnd[1],bndo[1]),min(bnd[2],bndo[2]),max(bnd[3],bndo[3]),min(bnd[4],bndo[4]),max(bnd[5],bndo[5]))

    print 'Origin: %8.3f %8.3f %8.3f '%(bnd[0],bnd[2],bnd[4]), 'Bounds: %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f '%bnd

    surface = filtsurf.GetOutput()

    inlets = []
    for filtin in filtins:
        inlets.append(filtin.GetOutput())
    outlets = []
    for filtout in filtouts:
        outlets.append(filtout.GetOutput())

    clip, aclip = False, False
    if filtclip: clip = filtclip.GetOutput()
    if filtaclip: aclip = filtaclip.GetOutput()

    XR = int(bnd[0])-2, int(bnd[1])+2
    YR = int(bnd[2])-2, int(bnd[3])+2
    ZR = int(bnd[4])-2, int(bnd[5])+2
    print '\nInitial discrete bounds:',XR,YR,ZR
    GRIDSPACING = int(GRIDSPACING)

    points = vtkPoints()
    if clip:
        clipbnd = clip.GetBounds()
        XR = int(max(XR[0], clipbnd[0]-2)), int(min(XR[1],clipbnd[1]+2))
        YR = int(max(YR[0], clipbnd[2]-2)), int(min(YR[1],clipbnd[3]+2))
        ZR = int(max(ZR[0], clipbnd[4]-2)), int(min(ZR[1],clipbnd[5]+2))
        print 'resetting to clip region bounds:', XR,YR,ZR
    

    ibegin = int(XR[0] / GRIDSPACING) * GRIDSPACING
    jbegin = int(YR[0] / GRIDSPACING) * GRIDSPACING
    kbegin = int(ZR[0] / GRIDSPACING) * GRIDSPACING

    iend = XR[1]+2*GRIDSPACING-(XR[1]%GRIDSPACING)
    jend = YR[1]+2*GRIDSPACING-(YR[1]%GRIDSPACING)
    kend = ZR[1]+2*GRIDSPACING-(ZR[1]%GRIDSPACING)

    print '\nDiscrete bounds:',ibegin,':',iend,'  ',jbegin,':',jend,'  ',kbegin,':',kend,'\n'

    # fill in the putative mask of points by sweeping
    print 'Initial points filling...'
    for k in xrange(kbegin, kend, GRIDSPACING):
        for j in xrange(jbegin, jend, GRIDSPACING):
            for i in xrange(ibegin, iend, GRIDSPACING):

                points.InsertNextPoint(i, j, k)

    pointsPolydata = vtkPolyData()
    pointsPolydata.SetPoints(points);


    # Test the insider points
    select_surf = vtkSelectEnclosedPoints()
    select_surf.SetTolerance(TOLERANCE)
    print '    testing tolerance:', select_surf.GetTolerance()

    select_surf.SetInputData(pointsPolydata);
    select_surf.SetSurfaceData(surface);
    select_surf.Update();

    if len(inlets)>0:
        select_inlets = []
        for inlet in inlets:
            select_inlet = vtkSelectEnclosedPoints()
            select_inlet.SetTolerance(0.)
            select_inlet.SetInputData(pointsPolydata);
            select_inlet.SetSurfaceData(inlet);
            select_inlet.Update();
            select_inlets.append(select_inlet)
    if len(outlets):
        select_outlets = []
        for outlet in outlets:
            select_outlet = vtkSelectEnclosedPoints()
            select_outlet.SetTolerance(0.)
            select_outlet.SetSurfaceData(outlet);
            select_outlet.SetInputData(pointsPolydata);
            select_outlet.Update();
            select_outlets.append(select_outlet)
    if clip: 
        select_clip = vtkSelectEnclosedPoints()
        select_clip.SetTolerance(TOLERANCE)
        select_clip.SetInputData(pointsPolydata);
        select_clip.SetSurfaceData(clip);
        select_clip.Update();
    if aclip: 
        select_aclip = vtkSelectEnclosedPoints()
        select_aclip.SetTolerance(TOLERANCE)
        select_aclip.SetInputData(pointsPolydata);
        select_aclip.SetSurfaceData(aclip);
        select_aclip.Update();
 
    xmn = 1e5; xmx =-1e5
    ymn = 1e5; ymx =-1e5
    zmn = 1e5; zmx =-1e5
    nfluid,nwall,ninlet,noutlet = 0,0,0,0
    for n in xrange( pointsPolydata.GetNumberOfPoints() ):
        if select_surf.IsInside(n):
            p = pointsPolydata.GetPoint(n)
            P = int(p[0]), int(p[1]), int(p[2])

            xmn = min(xmn,P[0]); xmx = max(xmx,P[0])
            ymn = min(ymn,P[1]); ymx = max(ymx,P[1])
            zmn = min(zmn,P[2]); zmx = max(zmx,P[2])

            if clip and not select_clip.IsInside(n): continue
            if aclip and    select_aclip.IsInside(n): continue

            found = False
            if len(inlets)>0 :
                for select_inlet in select_inlets:
                    if select_inlet.IsInside(n):
                        ninlet += 1
                        found = True
                        break

            if len(outlets) and not found:
                for select_outlet in select_outlets:
                    if select_outlet.IsInside(n):
                        noutlet += 1
                        found = True
                        break

            if not found:
                nfluid += 1

    print '.....',nfluid,0,ninlet,noutlet
    print 'Point Mins:', xmn,ymn,zmn, 'Maxs:', xmx,ymx,zmx

    nx = xmx + 2*GRIDSPACING
    ny = ymx + 2*GRIDSPACING
    nz = zmx + 2*GRIDSPACING

    msh = mesh.Mesh(nx,ny,nz,gridspacing=GRIDSPACING)

    msh.itppp_f = np.zeros(nfluid,  dtype='i8')
    msh.itppp_i = np.zeros(ninlet,  dtype='i8')
    msh.itppp_o = np.zeros(noutlet,  dtype='i8')
    msh.itppp_i_id, msh.itppp_o_id = [], []
    

    msh.nfluid,msh.nwall,msh.ninlet,msh.noutlet = 0,0,0,0

    i = 1
    inletId = {}
    for inlet in inlets:
        inletId[inlet] = i
        i += 1
    outletId = {}
    for outlet in outlets:
        outletId[outlet] = i
        i += 1

    """
    select_surf = vtkSelectEnclosedPoints()
    select_surf.SetTolerance(TOLERANCE)
    print '    testing tolerance:', select_surf.GetTolerance()

    select_surf.SetInputData(pointsPolydata);
    select_surf.SetSurfaceData(surface);
    select_surf.Update();
    """

    if wrapbywall:

        DEAD_NODE=0
        FLUID_NODE=1
        INLET_NODE=3
        OUTLET_NODE=4

        nodes = np.full([nx/GRIDSPACING+1,ny/GRIDSPACING+1,nz/GRIDSPACING+1],DEAD_NODE, np.uint)

        innerPoints=np.full([nx/GRIDSPACING+1,ny/GRIDSPACING+1,nz/GRIDSPACING+1], False,np.bool)

    for n in xrange( pointsPolydata.GetNumberOfPoints() ):
        if select_surf.IsInside(n):
            p = pointsPolydata.GetPoint(n)

            i,j,k = int(p[0]), int(p[1]), int(p[2])
            i4 = msh.i4back(i,j,k)
            #i4 =  long(long(k)*nxy2 + long(j)*nx2 + long(i))

            if wrapbywall:
                #innerPoints[i,j,k]=1
                innerPoints[i/GRIDSPACING,j/GRIDSPACING,k/GRIDSPACING]=True

            if clip and not select_clip.IsInside(n): continue
            if aclip and    select_aclip.IsInside(n): continue

            found = False
            if len(inlets)>0:
                for inlet,select_inlet in zip(inlets,select_inlets):
                    if select_inlet.IsInside(n):
                        msh.inlet.append([i,j,k])
                        msh.itppp_i[msh.ninlet] = i4
                        msh.itppp_i_id.append( inletId[inlet] )
                        msh.ninlet += 1
                        found = True
                    if wrapbywall: nodes[i/GRIDSPACING,j/GRIDSPACING,k/GRIDSPACING]=INLET_NODE
                    break

            if len(outlets)>0 and not found :
                found = False
                for outlet,select_outlet in zip(outlets,select_outlets):
                    if select_outlet.IsInside(n):
                        msh.outlet.append([i,j,k])
                        msh.itppp_o[msh.noutlet] = i4
                        msh.itppp_o_id.append( outletId[outlet] )
                        msh.noutlet += 1
                        found = True
                        if wrapbywall: nodes[i/GRIDSPACING,j/GRIDSPACING,k/GRIDSPACING]=OUTLET_NODE
                        break

            if not found:
                msh.fluid.append([i,j,k])
                msh.itppp_f[msh.nfluid] = i4
                msh.nfluid += 1
                if wrapbywall: nodes[i/GRIDSPACING,j/GRIDSPACING,k/GRIDSPACING]=FLUID_NODE

    msh.itppp_f.sort()
    msh.itppp_i.sort()
    msh.itppp_o.sort()
    msh.i_globs, msh.o_globs = [], []

    for inlet in inlets:
        msh.i_globs.append( {'id':inletId[inlet], 'bctype':'velocity', 'iodir':'0 0 +1 ', 'ioval':0.0} )
    for outlet in outlets:
        msh.o_globs.append( {'id':outletId[outlet], 'bctype':'pressure', 'iodir':'0 0 -1', 'ioval':0.0} )

    print 'Number of Fluid,Wall,Inlet,Outlet: ',msh.nfluid,msh.nwall,msh.ninlet,msh.noutlet

    if wrapbywall:
      
      
        print '\n  wall wrapping'

        # determine neighbors of fluids nodes
        kernel_mat = ndimage.generate_binary_structure(3, 2)
        assert kernel_mat.sum() == 19
        fluid_nodes_dilated = ndimage.morphology.binary_dilation(nodes==FLUID_NODE,kernel_mat)
        dead_nodes=(nodes == DEAD_NODE)
        if openends:
            wall_nodes= fluid_nodes_dilated & dead_nodes
        else:
            # determine OUTER and DEAD fluid nodes neighbors
            wall_nodes= fluid_nodes_dilated & dead_nodes & (~innerPoints)
        #wall_nodes= np.reshape(,[wall_nodes.size])
        q=0
        msh.nwall=wall_nodes.sum()
        msh.itppp_w=np.zeros(msh.nwall,np.int64)
        for (i,j,k),w in np.ndenumerate(wall_nodes):
          if w:
            msh.itppp_w[q] =msh.i4back(i*GRIDSPACING,j*GRIDSPACING,k*GRIDSPACING)
            q=q+1
        
        print 'final number of wall nodes',msh.nwall,'\n'
        msh.itppp_w.sort()

    return msh


def _viz3D():

    cubeMapper = vtkPolyDataMapper()
    cubeMapper.SetInputConnection(readsurf.GetOutputPort());
 
    cubeActor = vtkActor()
    cubeActor.SetMapper(cubeMapper);
    cubeActor.GetProperty().SetOpacity(0.5);
 
    vertexGlyphFilter = vtkVertexGlyphFilter()
    vertexGlyphFilter.AddInputData(pointsPolydata)

    vertexGlyphFilter.Update();
 
    pointsMapper = vtkPolyDataMapper()
    pointsMapper.SetInputConnection(vertexGlyphFilter.GetOutputPort());
 
    pointsActor = vtkActor()
    pointsActor.SetMapper(pointsMapper);
    pointsActor.GetProperty().SetPointSize(2);
    pointsActor.GetProperty().SetColor(0.0,0.0,1)
 
    renderer = vtkRenderer()
    renderWindow = vtkRenderWindow()
    renderWindow.AddRenderer(renderer);
    renderWindowInteractor = vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
 
    # Add the actor to the scene
    renderer.AddActor(cubeActor);
    renderer.AddActor(pointsActor)
    renderer.SetBackground(.0, 1,.0);
 
    # Render and interact
    renderWindow.Render();
    renderWindowInteractor.Start();


if __name__ == '__main__':

    print '\n\n *** ',sys.argv[0],' ***\n'
    print sys.argv[1:]
    print 

    parser = argparse.ArgumentParser()

    parser.add_argument('-s', '--surface',    required=True, help='input stl file')
    parser.add_argument('-C', '--clip',   required=False, default=False, help='clip region stl file')
    parser.add_argument('-A', '--aclip',   required=False, default=False, help='anticlip region stl file')
    parser.add_argument('-i', '--inlets',      required=False, default=False, help='inlet stl file(s)')
    parser.add_argument('-o', '--outlets',     required=False, default=False, help='outlet stl file(s)')
    parser.add_argument('-O', '--outfile',     required=True, default='bgkflag', help='outfile')
    parser.add_argument('-S', '--SCALE',      required=False, default=1.0, help='scale factor')
    parser.add_argument('-G', '--GRIDSPACING', required=False, default=1, help='grid spacing')
    parser.add_argument('-w', '--wrapbywall', required=False, action='store_true', help='wrap by wall at the end')
    parser.add_argument('-e', '--openends', required=False, action='store_true', help='wrap by wall with open ends (no clipping)')
    parser.add_argument('-V', '--viz3D', required=False, default=False, help='3d viz')

    parser.add_argument('-T', '--translate',required=False,default=False,help='if given: Tx Ty Tz')

    args = parser.parse_args()


    inlets = None
    if args.inlets:
        inlets = []
        for istr in args.inlets.split():
            print 'inlet:',istr
            inlets.append(istr)

    outlets = None
    if args.outlets:
        outlets = []
        for istr in args.outlets.split():
            print 'outlet:',istr
            outlets.append(istr)

    # build mesh inside surface
    T=args.translate
    if T != False:
      T=args.translate.split()
      T=[float(T[0]),float(T[1]),float(T[2])]


    # read all needed surfaces from file 
    filtsurf, filtclip, filtaclip, filtins, filtouts = readfiles(args.surface, args.clip, args.aclip, inlets, outlets)

    filtsurf, filtclip, filtaclip, filtins, filtouts = positiveoctant(filtsurf, filtclip, filtaclip, filtins, filtouts, T)

    filtsurf, filtclip, filtaclip, filtins, filtouts = scalesolids(filtsurf, filtclip, filtaclip, filtins, filtouts, float(args.SCALE))

    # build mesh nodes in various surfaces
    mesh = buildinside(filtsurf, filtclip, filtaclip, filtins, filtouts, args.wrapbywall, args.openends, int(args.GRIDSPACING))

    """
    i_globs, o_globs = [], []

    if inlets:
      i_globs=[]
      i_globs.append( {'id':1, 'bctype':'velocity', 'iodir':'-1  0  0', 'ioval':0.0} )
    else:
      i_globs=None
      
    if outlets:
      o_globs = []
      o_globs.append( {'id':2, 'bctype':'pressure', 'iodir':'+1  0  0', 'ioval':0.1} )
    else:
      o_globs=None

    mesh.specifyNodes(mesh.itppp_f, mesh.itppp_w, mesh.itppp_i, mesh.itppp_o, 
                      i_id = mesh.itppp_i_id, o_id = mesh.itppp_o_id, 
                      i_globs = i_globs, o_globs = o_globs)
    """
    mesh.specifyNodes(mesh.itppp_f, mesh.itppp_w, mesh.itppp_i, mesh.itppp_o, 
                      i_id = mesh.itppp_i_id, o_id = mesh.itppp_o_id)

    # write mesh
    print 'writing mesh on files',args.outfile + '*'
    mesh.writeMOEBIUSinput(args.outfile)


    # visualization
    if args.viz3D: _viz3D()


    # write file xyz
    import TOOLS.bgkflag2xyz as bgkflag2xyz
    class E(): pass
    e = E()
    e.infiles = [args.outfile]
    e.outfile = args.outfile
    e.outwall = 'wall' + str(args.GRIDSPACING)

    bgkflag2xyz.main(e)


