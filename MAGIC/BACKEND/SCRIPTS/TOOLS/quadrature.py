#!/usr/bin/env python

from math import *

class Triangle:

    # vertex
    x1=0.; y1=0.; z1=0.
    x2=0.; y2=0.; z2=0.
    x3=0.; y3=0.; z3=0.

    # V1 = [0.,0.,0.]
    # V2 = [0.,0.,0.]
    # V3 = [0.,0.,0.]
    
    n_q = 4
    w_q = [0.281250,0.260416660,0.260416660,0.260416660]
    beta_q = [[1.0/3.0,1.0/3.0,1.0/3.0],
            [11.0/15.0,  2.0/15.0,  2.0/15.0],
            [2.0/15.0, 11.0/15.0,  2.0/15.0],
            [2.0/15.0,  2.0/15.0, 11.0/15.0]] 

    def GetQuadraturePoints(self):

        qpts = [ [0,0,0], [0,0,0], [0,0,0], [0,0,0] ]
        for iq in range(0,self.n_q):
            qpts[iq] = [self.beta_q[iq][0]*self.x1 + self.beta_q[iq][1]*self.x2 + self.beta_q[iq][2]*self.x3,
                        self.beta_q[iq][0]*self.y1 + self.beta_q[iq][1]*self.y2 + self.beta_q[iq][2]*self.y3,
                        self.beta_q[iq][0]*self.z1 + self.beta_q[iq][1]*self.z2 + self.beta_q[iq][2]*self.z3]

        return qpts


    def GetQuadrature(self,sqty=None,vqty=None):

        AA = self.x2-self.x3, self.y2-self.y3, self.z2-self.z3
        BB = self.x1-self.x3, self.y1-self.y3, self.z1-self.z3

        VX = AA[1]*BB[2] - AA[2]*BB[1]
        VY = AA[2]*BB[0] - AA[0]*BB[2]
        VZ = AA[0]*BB[1] - AA[1]*BB[0]

        area_tri = 0.5 * sqrt( VX**2 + VY**2 + VZ**2 )

        # get fluid_velocity on each quadrature point
        # u,v,w = get_fluid_velocity(self.n_q,px_q,py_q,pz_q)
        # density = get_fluid_density(self.n_q,px_q,py_q,pz_q)
        flux_tri = 0
        rho_tri = 0
        for iq in range(0,self.n_q):

            if sqty:
                rho_tri = rho_tri + self.w_q[iq] * sqty[iq]

            if vqty:
                v = vqty[iq]
                # flux_tri = flux_tri + self.w_q[iq] * (VX * vqty[0][iq] + VY * vqty[1][iq] + VZ * vqty[2][iq])
                flux_tri = flux_tri + self.w_q[iq] * (VX*v[0] + VY*v[1] + VZ*v[2])


        return area_tri,rho_tri,flux_tri

        # sum over all triangles
        # section_flux = section_flux + flux_tri
        # section_area = section_area + area_tri

