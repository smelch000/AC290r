#!/usr/bin/env pythonw

from myvtk import *
import sys,math,random
import argparse

random.seed(17123)

#########################

def GenerationRound(clineseed,lengths,clines,GRAPHICS,renWin,ren,OUTFILE):

    print ('Generating...')

    cline = clineseed

    # for KK in xrange(MAXADDCLINES):
    KK = -1
    while True:

        cline = cline.NewBranch( cline,lengths,clines )

        if KK  == cline.MAXADDCLINES: break
        KK += 1

        for kk in xrange(cline.nPoints, cline.MAXADDPOINTS):

            added = cline.AddPoint( cline )

            if cline.LOCALSCORING :

                if not added:
                    break

            if not cline.LOCALSCORING and kk%cline.MAXDELAYEDCUT == 0:

                scorePrevious = 0.
                for i in xrange(cline.nPoints):
                    scorePrevious += cline.scorePoints[i]
                scorePrevious /= cline.nPoints

                countBad = 0
                for i in xrange(cline.nPoints - cline.MAXDELAYEDCUT, cline.nPoints):

                    scorePoint = cline.scorePoints[i]

                    if scorePoint < (cline.LOW_CONTRAST - cline.TAG * cline.BONDLENGTH * i) or scorePoint > cline.HIGH_CONTRAST:
                        countBad += 1

                if countBad > cline.MAXDELAYEDCUT/2:
                    # refused move: remove last bit from cline - hair cut
                    _inputPoints = vtkPoints()
                    _clinePut = vtkCellArray()
                    _scorePoints = []

                    _clinePut.InsertNextCell(cline.nPoints - cline.MAXDELAYEDCUT)

                    for i in xrange( cline.nPoints - cline.MAXDELAYEDCUT):
                        _inputPoints.InsertPoint(i, cline.inputPoints.GetPoint(i))
                        _clinePut.InsertCellPoint(i)
                        _scorePoints.append( cline.scorePoints[i] )
                        print (i,cline.scorePoints[i],cline.inputPoints.GetPoint(i))

                    cline.scorePoints = _scorePoints
                    cline.inputPoints = _inputPoints
                    cline.clinePut = _clinePut
                    """
                    cline.profileData.SetLines( _clinePut )
                    cline.profileData.SetPoints( cline.inputPoints )
                    cline.TubeLines = vtkTubeFilter()
                    cline.TubeLines.SetNumberOfSides(8)
                    if VTK_MAJOR_VERSION <= 5:
                        cline.TubeLines.SetInput( cline.profileData )
                    else:
                        cline.TubeLines.SetInputData( cline.profileData )
                    cline.TubeLines.SetRadius(.2)
                    cline.TubeLines.Update()
                    """

                    break

        if KK%10 == 0: 
            cline.BendingScore(lengths,clines)
            cline.ColorLines(lengths,clines)
            cline.CleanBranches(lengths,clines)
            if GRAPHICS:
                ren.ResetCamera()
                renWin.Render()

    cline.BendingScore(lengths,clines)
    cline.ColorLines(lengths,clines)
    cline.CleanBranches(lengths,clines)
    if GRAPHICS:
        ren.ResetCamera()
        renWin.Render()

    print ('completed: path multiplicity',len(lengths),'...max lengths %.1f %.1f ... %.1f ' % tuple(lengths[:2]+[lengths[-1]]))

    cline.WriteClineFile( OUTFILE,lengths,clines )

    print ('done...press z for another generation round')

class CLine():

    def __init__( self, image, pt0, pt1, GRAPHICS, renWin=None, ren=None ):

        self.GRAPHICS = GRAPHICS

        if GRAPHICS:
            self.renWin = renWin
            self.ren = ren
        else:
            self.renWin = None
            self.ren = None

        self.image = image 

        # contrast window
        self.LOW_CONTRAST = 250
        self.HIGH_CONTRAST = 500
        self.BEST_CONTRAST = 500

        # self.TAG to grow cline below contrast window, as point insertion proceeds
        # self.TAG = 2.5
        self.TAG = 1.0

        # number of random rotations to insert new point
        self.TRYROTATE = 20
        # angular window to insert new point
        self.MAXANG = 60.
        # max number of points to add
        self.MAXADDPOINTS = 50
        # max number of clines to add
        self.MAXADDCLINES = 200
        # bond length of cline
        self.BONDLENGTH = 2.5

        # multiplicative penalty for negative bending. Sensitive parameter !!
        self.BENDINGPENALTY = 7.

        # local vs delayed scoring as points are added to cline
        self.LOCALSCORING = False
        # bond length of cline for delayed scoring
        self.MAXDELAYEDCUT = 4

        self.PROBETYPE = 'POINT'
        # self.PROBETYPE = 'SPHERE'

        #########################

        x0,y0,z0 = pt0
        x1,y1,z1 = pt1
        n = self.Norm(x1-x0,y1-y0,z1-z0)
        n = max(1e-5,n)
        x1 = x0 + self.BONDLENGTH * (x1-x0)/n
        y1 = y0 + self.BONDLENGTH * (y1-y0)/n
        z1 = z0 + self.BONDLENGTH * (z1-z0)/n
        pt1 = (x1,y1,z1)

        # print 'pt0 score:',self.ProbePoint(pt0)
        # print 'pt1 score:',self.ProbePoint(pt1)

        self.inputPoints = vtkPoints()
        self.scorePoints = []

        self.colors = vtkFloatArray()
        self.colors.SetNumberOfComponents(1);
        self.colors.SetName ("Score");

        self.inputPoints.InsertPoint(0, pt0)
        self.scorePoints.append( self.ProbePoint(pt0) )
        self.colors.InsertNextValue( (1.*self.ProbePoint(pt0) - self.LOW_CONTRAST + 100.)/(self.HIGH_CONTRAST - self.LOW_CONTRAST + 100.) )

        self.inputPoints.InsertPoint(1, pt1)
        self.scorePoints.append( self.ProbePoint(pt1) )
        self.colors.InsertNextValue( (1.*self.ProbePoint(pt1) - self.LOW_CONTRAST + 100.)/(self.HIGH_CONTRAST - self.LOW_CONTRAST + 100.) )

        self.nPoints = 2

        # create glyphs for the pivot points to make the spline more clear
        # polydata to be glyphed
        self.inputData = vtkPolyData()
        self.inputData.SetPoints( self.inputPoints )

        # Use sphere as glyph source.
        self.balls = vtkSphereSource()
        self.balls.SetRadius(.8)
        self.balls.SetPhiResolution(10)
        self.balls.SetThetaResolution(10)
        self.balls.Update()

        self.glyph = vtkGlyph3D()

        if VTK_MAJOR_VERSION <= 5:
            self.glyph.SetInputConnection(  self.inputData.GetProducerPort())
            self.glyph.SetSource(   self.balls.GetOutput())
        else:
            self.glyph.SetInputData( self.inputData)
            self.glyph.SetSourceData( self.balls.GetOutput())

        self.glyph.SetScaleFactor(1)
        # self.glyph.SetColorModeToColorByVector()
        self.glyph.SetColorModeToColorByScalar()
        # self.glyph.SetScaleModeToScaleByVector()
        self.glyph.SetScaleModeToDataScalingOff()

        lut = vtkLookupTable()
        lut.SetNumberOfTableValues(10)
        lut.Build()
        for i in xrange(10):
            lut.SetTableValue(i,i/10.,i/10.,i/10.,1)

        self.MGlyph = vtkPolyDataMapper()
        self.MGlyph.SetInputConnection( self.glyph.GetOutputPort())
        self.MGlyph.SetLookupTable(lut)

        self.AGlyph = vtkActor()
        self.AGlyph.SetMapper( self.MGlyph )
        # self.AGlyph.GetProperty().SetDiffuseColor(tomato)
        # self.AGlyph.GetProperty().SetSpecular(.3)
        # self.AGlyph.GetProperty().SetSpecularPower(30)

        # Generate the polyline for the spline.
        self.points = vtkPoints()
        self.profileData = vtkPolyData()

        # Create polyline
        self.clinePut = vtkCellArray()
        self.clinePut.InsertNextCell(self.nPoints)

        for i in xrange( self.nPoints ):
            self.clinePut.InsertCellPoint(i)
 
        self.profileData.SetPoints( self.inputPoints )
        self.profileData.SetLines( self.clinePut )

        """
        # spline = vtkKochanekSpline(); spline.SetDefaultTension(.05)
        spline = vtkCardinalSpline()
        splinef = vtkSplineFilter()
        splinef.SetInputData(profileData)
        if True:
            splinef.SetSubdivideToLength()
            splinef.SetLength(.5) # the smaller length the smoother is the spline
        else:
            splinef.SetSubdivideToSpecified()
            splinef.SetNumberOfSubdivisions(500)
        splinef.SetSpline(spline);
        splinef.Update()
        splinef.GetSpline().ClosedOn();
        splinef.Update()
        """

        # tube to the polyline
        self.TubeLines = vtkTubeFilter()
        self.TubeLines.SetNumberOfSides(8)
        if VTK_MAJOR_VERSION <= 5:
            self.TubeLines.SetInput( self.profileData )
        else:
            self.TubeLines.SetInputData( self.profileData )
        self.TubeLines.SetRadius(.2)

        """
        # tube to the spline
        TubeSpline = vtkTubeFilter()
        TubeSpline.SetNumberOfSides(8)
        if VTK_MAJOR_VERSION <= 5:
            TubeSpline.SetInput(splinef)
        else:
            TubeSpline.SetInputData(splinef.GetOutput())
        TubeSpline.SetRadius(.1)
        """

        self.MLines = vtkPolyDataMapper()
        self.MLines.SetInputConnection(   self.TubeLines.GetOutputPort())
        """
        self.MapSpline = vtkPolyDataMapper()
        self.MapSpline.SetInputConnection(TubeSpline.GetOutputPort())
        """

        self.ALines = vtkActor()
        self.ALines.SetMapper(  self.MLines)
        # self.ALines.GetProperty().SetDiffuseColor(banana)
        """
        self.ASpline = vtkActor()
        self.ASpline.SetMapper( self.MapSpline)
        self.ASpline.GetProperty().SetDiffuseColor(tomato)
        """
    def GetLength(self):

        pts = self.inputPoints
        ln = 0
        for i in xrange(1, pts.GetNumberOfPoints()):
            x0,y0,z0 = pts.GetPoint(i-1)
            x1,y1,z1 = pts.GetPoint(i)
            ln += self.Norm(x1-x0,y1-y0,z1-z0)
        return ln

    def ProbePoint(self,pp):

        xprobe = vtkProbeFilter()

        if self.PROBETYPE == 'POINT':

            pt = vtkPoints()
            pt.InsertPoint(0, pp)
            pd = vtkPolyData()
            pd.SetPoints(pt)

            if VTK_MAJOR_VERSION <= 5:
                xprobe.SetInput(pd)
                xprobe.SetSource(self.image)
            else:
                xprobe.SetInputData(pd)
                xprobe.SetSourceData(self.image)

        elif self.PROBETYPE == 'SPHERE':

            sph = vtkSphereSource()
            sph.SetThetaResolution( 5 )
            sph.SetPhiResolution( 5 )
            sph.SetCenter( pp )
            sph.SetRadius( 1.0 )

            if VTK_MAJOR_VERSION <= 5:
                xprobe.SetInputConnection(sph.GetOutputPort())
                xprobe.SetSource(self.image)
            else:
                xprobe.SetInputConnection(sph.GetOutputPort())
                xprobe.SetSourceData(self.image)

        # xprobe.PassPointArraysOn()
        xprobe.Update()

        if self.PROBETYPE == 'POINT':

            result = xprobe.GetOutput().GetPointData().GetArray('DICOMImage').GetValue(0)

        elif self.PROBETYPE == 'SPHERE':

            nv = xprobe.GetValidPoints().GetNumberOfTuples()
            if nv>0:
                result = 0
                for i in xrange(nv):
                    val = xprobe.GetOutput().GetPointData().GetArray('DICOMImage').GetValue(i)

                    #if val < LOW or val > HIGH: 
                    #    result = -10000
                    #    break

                    result += val

                result /= nv
            else:
                result = -99999

        return result # xprobe.GetOutput().GetPointData().GetArray('DICOMImage').GetValue(0)

    def ProbePolyLineOLD(self,cline):

        xprobe = vtkProbeFilter()
        if VTK_MAJOR_VERSION <= 5:
            # xprobe.SetInput(probePolyData)
            xprobe.SetInput(cline.profileData)
            xprobe.SetSource(self.image)
        else:
            # xprobe.SetInputData(probePolyData)
            xprobe.SetInputData(cline.profileData)
            xprobe.SetSourceData(self.image)
        xprobe.PassPointArraysOn()
        xprobe.Update()

        valids = xprobe.GetValidPoints()
        ntup = valids.GetNumberOfTuples()

        arr = xprobe.GetOutput().GetPointData().GetArray('DICOMImage')
        score = 0
        for pid in xrange(ntup):
            print (pid,'Probes:',arr.GetValue(pid))
            score += arr.GetValue(pid)
        return score / ntup

    def Norm(self,dx,dy,dz):
        return math.sqrt(dx**2 + dy**2 + dz**2)

    def RandomRotate(self,pt0,pt1,ANG):

        trans = vtkTransform()

        # r1,r2,r3 = vtkMath.Random(-.5, +.5),vtkMath.Random(-.5, +.5),vtkMath.Random(-.5, +.5)
        # r1,r2,r3 = vtkMath.Random(0, 1),vtkMath.Random(0, 1),vtkMath.Random(0, 1)
        r1,r2,r3 = random.random()-.5, random.random()-.5, random.random()-.5

        trans.RotateX( ANG * r1 )
        trans.RotateY( ANG * r2 )
        trans.RotateZ( ANG * r3 )

        aLine = vtkLineSource()
        aLine.SetPoint1(pt0)
        aLine.SetPoint2(pt1)

        tpd1 = vtkTransformPolyDataFilter()
        tpd1.SetInputConnection(aLine.GetOutputPort())
        tpd1.SetTransform(trans)
        tpd1.Update()

        PT0 = tpd1.GetOutput().GetPoints().GetPoint(0)
        PT1 = tpd1.GetOutput().GetPoints().GetPoint(1)

        return PT0,PT1

    def CleanTwinLines(self,lengths,clines):

        removables = []
        for i in xrange(len(lengths)-1):
            lni = lengths[i]
            cli = clines[lni]

            twinscore = 0.
            for j in xrange(i+1, len(lengths)):
                lnj = lengths[j]
                clj = clines[lnj]

                for ii in xrange(cli.nPoints-1): 

                    pi = cli.inputPoints.GetPoint(ii)
                    pj = clj.inputPoints.GetPoint(ii)
                    dx,dy,dz = pi[0]-pj[0],pi[1]-pj[1],pi[2]-pj[2]
                    twinscore += self.Norm(dx,dy,dz)

                if  twinscore < cli.nPoints * 2.:
                    removables.append(lnj)

        return set(removables)

    def CleanBranches(self,lengths,clines):

        twinremovables = self.CleanTwinLines(lengths,clines)

        # print 'removing twins:',len(twinremovables)

        mxb = 0
        for ln in lengths:
            cl = clines[ln]
            mxb = max(mxb, abs(cl.BendingScore))

        # cleanCUT = min( len(lengths), self.MAXADDCLINES / 5 )
        # cleanCUT = min( len(lengths), self.MAXADDCLINES / 2 )
        cleanCUT = min( len(lengths), self.MAXADDCLINES )

        # if len(lengths) < cleanCUT: return

        if self.GRAPHICS:
            for ln in clines:
                self.ren.RemoveActor( clines[ln].ALines )
                self.ren.RemoveActor( clines[ln].AGlyph )

            self.renWin.Render()

        clNEW = {}
        lnNEW = []
        for i in xrange( cleanCUT ) :

            ln = lengths[i]
            cl = clines[ln]

            ascore = cl.BendingScore / mxb
            if ascore < 0: continue

            if ln in twinremovables: continue

            if ln < 4 * self.BONDLENGTH: continue

            lnNEW.append( ln )
            clNEW[ln] = clines[ln]

            if self.GRAPHICS:
                self.ren.AddActor( clines[ln].ALines )
                self.ren.AddActor( clines[ln].AGlyph )
                clines[ln].AGlyph.VisibilityOff()

        if self.GRAPHICS:
            self.renWin.Render()

        lengths = lnNEW
        clines = clNEW

        if len(lengths)>=3:
            print ('CleanBranch: path multiplicity',len(lengths),'...max lengths %.1f %.1f ... %.1f ' % tuple(lengths[:2]+[lengths[-1]]))


    def BendingScore(self,lengths,clines):

        for ln in lengths:
            cl = clines[ln]

            persist = 0
            for i in xrange(cl.nPoints-1): 
            # for i in xrange(1): 

                Ai = cl.inputPoints.GetPoint(i)
                Bi = cl.inputPoints.GetPoint(i+1)

                dxi,dyi,dzi = Bi[0]-Ai[0], Bi[1]-Ai[1], Bi[2]-Ai[2]
                n = self.Norm(dxi,dyi,dzi)
                if n > 1.e-6:
                    dxi /= n
                    dyi /= n
                    dzi /= n

                for j in xrange(i+1,cl.nPoints-1): 
                    Aj = cl.inputPoints.GetPoint(j)
                    Bj = cl.inputPoints.GetPoint(j+1)

                    dxj,dyj,dzj = Bj[0]-Aj[0], Bj[1]-Aj[1], Bj[2]-Aj[2]
                    n = self.Norm(dxj,dyj,dzj)
                    if n > 1.e-6:
                        dxj /= n
                        dyj /= n
                        dzj /= n

                    cosa = (dxi*dxj + dyi*dyj + dzi*dzj)

                    # push down score of curved bits
                    if cosa < 0: cosa *= self.BENDINGPENALTY 

                    persist += cosa

            persist /= ln
            persist /= ln

            if persist <= 0:
                cl.BendingScore = -1.
            else:
                # print 'persist:',persist
                # cl.BendingScore = - ln * math.log(persist / ln)
                # cl.BendingScore = - math.log(persist)
                cl.BendingScore = - ln / math.log(persist)

    def ColorLines(self,lengths,clines):

        if len(lengths)<=1: return

        mxb = 0
        for ln in lengths:
            cl = clines[ln]
            mxb = max(mxb, abs(cl.BendingScore))

        for ln in lengths:

            cl = clines[ln]

            cl.AGlyph.VisibilityOff()
            cl.ALines.VisibilityOn()
            """
            ascore = 0.
            for i in xrange(cl.nPoints): 
                ascore += cl.scorePoints[i]/cl.nPoints
            ascore = (ascore - self.LOW_CONTRAST) / (self.HIGH_CONTRAST - self.LOW_CONTRAST)
            """
            ascore = cl.BendingScore / mxb
            # ascore = math.log(cl.BendingScore / mxb)
            # ascore /= 4.

            if ascore < 0.:
                ascore = 0.1
                cl.ALines.GetProperty().SetColor( 1., 0., 0.)

            else:
                cl.ALines.GetProperty().SetColor( ascore, ascore, ascore )
            """
            elif ascore < 0.1:
                cl.ALines.GetProperty().SetColor( 0., 1., 0.)

            elif ascore < 0.2:
                cl.ALines.GetProperty().SetColor( 0., 0., 1.)

            elif ascore < 0.3:
                cl.ALines.GetProperty().SetColor( 1., 1., 0.)

            elif ascore < 1.0:
                cl.ALines.GetProperty().SetColor( 1., 0., 1.)
            else:
                cl.ALines.GetProperty().SetColor( 0., 1., 1.)
            """

            # cl.TubeLines.SetRadius(ascore)

        if self.GRAPHICS:
            self.renWin.Render()

    def ColorBalls(self,lengths,clines):

        if len(lengths)<=1: return

        for ln in lengths:

            cl = clines[ln]

            ascore = 0.
            for i in xrange(cl.nPoints): 
                asc = cl.scorePoints[i]
                ascore += asc / cl.nPoints
                asc = (asc - self.LOW_CONTRAST) / (self.HIGH_CONTRAST - self.LOW_CONTRAST)

                #if asc <= 0:
                #    cl.AGlyph.GetProperty().SetColor( - asc/2. , 0., 0. )
                #else:
                #    cl.AGlyph.GetProperty().SetColor( 0., 2.*asc, 0. )

            cl.inputData.GetPointData().SetScalars(cl.colors)
            cl.AGlyph.VisibilityOn()
            cl.ALines.VisibilityOff()


            """
            # ascore = (ascore - self.LOW_CONTRAST) / (self.HIGH_CONTRAST - self.LOW_CONTRAST)
            ascore = (ascore - self.LOW_CONTRAST + 50.) / (self.HIGH_CONTRAST - self.LOW_CONTRAST)

            if ascore <= 0.001:
                ascore = 0.1
                cl.ALines.GetProperty().SetColor( 1., 0., 0.)

            else:
                cl.ALines.GetProperty().SetColor( ascore, ascore, ascore )
            """

            """
            elif ascore < 0.1:
                cl.ALines.GetProperty().SetColor( 0., 1., 0.)

            elif ascore < 0.2:
                cl.ALines.GetProperty().SetColor( 0., 0., 1.)

            elif ascore < 0.3:
                cl.ALines.GetProperty().SetColor( 1., 1., 0.)

            elif ascore < 1.0:
                cl.ALines.GetProperty().SetColor( 1., 0., 1.)
            else:
                cl.ALines.GetProperty().SetColor( 0., 1., 1.)
            """

            # cl.TubeLines.SetRadius(ascore)

        if self.GRAPHICS:
            self.renWin.Render()


    def NewBranch( self,cline,lengths,clines ):

        # cline ordered according to length

        ln = cline.GetLength()

        lengths.append(ln)
        lengths.sort(reverse=True)

        clines[ln] = cline

        if False:
            maxln = 0.
            for ln in lengths:
                cl = clines[ln]
                if ln > maxln:
                    maxln = ln
                    clinepick = cl
        else:
            # pick one of the first few longest clines
            icln = int(random.random() * min(10, len(clines)) )
            clinepick = clines[ lengths[icln] ]

        OPT=0
        if OPT==0:
            irnd = 0

        elif OPT==1:
            irnd = int(random.random() * (clinepick.nPoints-1))

        elif OPT==2:
            # pick a point
            while True:
                irnd = int(random.random() * clinepick.nPoints/2) # select preferentially from upstream
                # irnd = int(random.random() * 3.*clinepick.nPoints/4) # select preferentially from upstream
                if irnd < clinepick.nPoints-1:
                    break

        pt0 = clinepick.inputPoints.GetPoint( irnd )
        pt1 = clinepick.inputPoints.GetPoint( irnd+1 )

        PT0,PT1 = self.RandomRotate(pt0,pt1,self.MAXANG/10.)

        PT1 = (pt0[0] + (PT1[0]-PT0[0]), pt0[1] + (PT1[1]-PT0[1]), pt0[2] + (PT1[2]-PT0[2]))

        cline = CLine( self.image, pt0, PT1, self.GRAPHICS, self.renWin, self.ren )

        if self.GRAPHICS:
            self.ren.AddActor( cline.ALines )
            self.ren.AddActor( cline.AGlyph )
            cline.AGlyph.VisibilityOff()
            # self.renWin.Render()

        return cline

    def AddPoint(self,cline):

        pt0 = cline.inputPoints.GetPoint( cline.nPoints-2 )
        pt1 = cline.inputPoints.GetPoint( cline.nPoints-1 )

        added = False

        scorePoint = -10000
        optmetric = +10000

        # make random rotations and pick base choice
        for k in xrange(self.TRYROTATE):

            PT0,PT1 = self.RandomRotate(pt0,pt1,self.MAXANG)

            pg = (pt1[0] + (PT1[0]-PT0[0]), pt1[1] + (PT1[1]-PT0[1]), pt1[2] + (PT1[2]-PT0[2]))

            scoret = cline.ProbePoint(pg)

            metric = abs( scoret - self.BEST_CONTRAST )

            # print k,'scoret:',scoret,'metric',metric

            # if scoret > 0 and scoret > scorePoint:
            if metric < optmetric:
                optmetric = metric
                scorePoint = scoret
                pt = pg

        # print 'scorePoint:',scorePoint

        # see if the new choice is in contrast interval
        if self.LOCALSCORING:

            if scorePoint > (self.LOW_CONTRAST - self.TAG * self.BONDLENGTH * cline.nPoints) and scorePoint < self.HIGH_CONTRAST:

                cline.inputPoints.InsertPoint(cline.nPoints, pt)
                cline.scorePoints.append(scorePoint)
                cline.colors.InsertNextValue( (1.*scorePoint - self.LOW_CONTRAST + 100.)/(self.HIGH_CONTRAST - self.LOW_CONTRAST + 100.) )
                cline.nPoints += 1

                # print cline.nPoints,'Points score:',scorePoint,'optmetric:',optmetric

                cline.inputData.Modified()
                # cline.glyph.Update()

                # cline.clinePut = vtkCellArray()
                cline.clinePut.InsertNextCell(cline.nPoints)
                for i in xrange(cline.nPoints):
                    cline.clinePut.InsertCellPoint(i)

                cline.profileData.Modified()
                added = True
        else:

            cline.inputPoints.InsertPoint(cline.nPoints, pt)
            cline.scorePoints.append(scorePoint)
            cline.colors.InsertNextValue( (1.*scorePoint - self.LOW_CONTRAST + 100.)/(self.HIGH_CONTRAST - self.LOW_CONTRAST + 100.) )
            cline.nPoints += 1

            # print cline.nPoints,'Points score:',scorePoint,'optmetric:',optmetric

            cline.inputData.Modified()
            # cline.glyph.Update()

            # Create polyline
            # cline.clinePut = vtkCellArray()
            cline.clinePut.InsertNextCell(cline.nPoints)
            for i in xrange(cline.nPoints):
                cline.clinePut.InsertCellPoint(i)

            cline.profileData.Modified()
            added = True

        # self.renWin.Render()
        return added
 
    def WriteClineFile(self,OUTFILE,lengths,clines):

        print ('writing on file....',OUTFILE,'\n')

        merge = vtk.vtkAppendPolyData()
        # for cl in clines:
        for ln in lengths:
            cl = clines[ln]
            if vtk.VTK_MAJOR_VERSION <= 5:
                merge.AddInput(cl.profileData) 
            else:
                merge.AddInputData(cl.profileData) 
        merge.Update()

        pdw = vtk.vtkPolyDataWriter()
        if vtk.VTK_MAJOR_VERSION <= 5:
            # pdw.SetInput( cline.profileData ) 
            pdw.SetInput( merge.GetOutput() )
        else:
            pdw.SetInputData( merge.GetOutput() )
        pdw.SetFileName( OUTFILE )
        pdw.Write()

def InitStringy(IMGFILE):

    print ('Reading file:',IMGFILE)

    reader = vtkXMLImageDataReader()
    reader.SetFileName(IMGFILE)
    reader.Update()

    return reader.GetOutput()


def Stringy(image,P1,P2,GRAPHICS,OUTFILE):

    ren = vtkRenderer()

    if GRAPHICS:
        renWin = vtkRenderWindow()
        renWin.AddRenderer(ren)
        renWin.SetSize(500, 500)
    else:
        renWin = None

    clines = {}
    lengths = []
    icol = [0]

    # cline = CLine( image,(61.6,43.3,7.1),(63.5,41.5,7.4),GRAPHICS,self.renWin ) # LCA
    # cline = CLine( image,(55.7,77.0,12.3),(52.9,80.4,14.1),GRAPHICS,self.renWin ) # RCA
    cline = CLine( image, P1, P2, GRAPHICS, renWin, ren )

    GenerationRound(cline,lengths,clines,GRAPHICS,renWin,ren,OUTFILE)

    if GRAPHICS:

        ren.AddActor( cline.ALines )
        # ren.AddActor(cline.ASpline)
        ren.AddActor( cline.AGlyph )
        cline.AGlyph.VisibilityOff()

        iren = vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)
        iren.Initialize()

        def keyPressEvent(obj, event):

            key = iren.GetKeySym()

            if key == 'z':
                GenerationRound(cline,lengths,clines,GRAPHICS,renWin,ren,OUTFILE)

            elif key == 'x':
                icol[0] += 1
                if icol[0]%2 == 0:
                    cline.ColorLines(lengths,clines)
                else:
                    cline.ColorBalls(lengths,clines)

            elif key == 'c':
                cline.WriteClineFile( OUTFILE,lengths,clines )

        iren.AddObserver("KeyPressEvent", keyPressEvent)
        renWin.Render()
        iren.Start()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--image',   required=True,  help='input vti file')
    parser.add_argument('-p1', '--point1', required=True,  help='input 1st marker point')
    parser.add_argument('-p2', '--point2', required=True,  help='input 2nd marker point')
    parser.add_argument('-o', '--output',                  help='output polyline file')
    parser.add_argument('-g', '--graphics',                help='graphics on/off')
    args = parser.parse_args()

    IMGFILE = args.image

    P = args.point1.split(',')
    P1 = (float(P[0]),float(P[1]),float(P[2]))
    P = args.point2.split(',')
    P2 = (float(P[0]),float(P[1]),float(P[2]))

    if args.output != None:
        OUTFILE = args.output
    else:
        OUTFILE = 'polyline.vtk'

    if args.graphics != None:
        if args.graphics=='yes' or args.graphics=='y' or args.graphics=='true' or args.graphics=='on':
            GRAPHICS = True
        else:
            GRAPHICS = False
    else:
        GRAPHICS = True

    image = InitStringy(IMGFILE)

    Stringy(image,P1,P2,GRAPHICS,OUTFILE)
