// rest populations f0 need not to be copied during streaming 
//      --> use a separate variable f0
//
// populations stored during streaming are those that are then read during the collision 
//      --> combine collide, compute hydrodynamic vars and stream into a single function
//
// save hydrodamic variables only when needed for analysis
//
// unrolling of all for loop
//      --> save overheads and use discrete speeds only for non-zero components
//
// optimize memory access pattern to improve cache utilisation

int i3D(int i, int j, int ip) {return (ip*NX*NY + NY*i + j);}

void stream_collide_save(double *f0, double *f1, double *f2, double *r, double *u, double *v, bool save)
 {
   // useful constants
  const double tauinv = 2.0/(6.0*nu+1.0); // 1/tau
  const double omtauinv = 1.0-tauinv; // 1 - 1/tau

  for(unsigned int y = 0; y < NY; ++y) {
    for(unsigned int x = 0; x < NX; ++x) {

       unsigned int xp1 = (x+1)%NX; unsigned int yp1 = (y+1)%NY; unsigned int xm1 = (NX+x-1)%NX; unsigned int ym1 = (NY+y-1)%NY;
       // direction numbering scheme 
       // 6 2 5
       // 3 0 1
       // 7 4 8
       double ft0 = f0[i3D(x,y,0)];

       // load populations from adjacent nodes
       double ft1 = f1[i3D(xm1,y, 1)]; 
       double ft2 = f1[i3D(x, ym1,2)]; 
       double ft3 = f1[i3D(xp1,y, 3)]; 
       double ft4 = f1[i3D(x, yp1,4)]; 
       double ft5 = f1[i3D(xm1,ym1,5)]; 
       double ft6 = f1[i3D(xp1,ym1,6)]; 
       double ft7 = f1[i3D(xp1,yp1,7)]; 
       double ft8 = f1[i3D(xm1,yp1,8)];

       // compute moments
       double rho = ft0+ft1+ft2+ft3+ft4+ft5+ft6+ft7+ft8; 
       double rhoinv = 1.0/rho;

       double ux = rhoinv*(ft1+ft5+ft8-(ft3+ft6+ft7)); 
       double uy = rhoinv*(ft2+ft5+ft6-(ft4+ft7+ft8));

       // only write to memory when needed
       if(save) {
              r[scalar_index(x,y)] = rho; 
              u[scalar_index(x,y)] = ux; 
              v[scalar_index(x,y)] = uy;
       }

       // now compute and relax to equilibrium

       // temporary variables
       double tw0r = tauinv*w0*rho; // w[0]*rho/tau
       double twsr = tauinv*ws*rho; // w[1-4]*rho/tau
       double twdr = tauinv*wd*rho; // w[5-8]*rho/tau

       double omusq = 1.0-1.5*(ux*ux+uy*uy); // 1- (3/2) u.u

       double tux = 3.0*ux;
       
       double tuy = 3.0*uy;
       f0[i3D(x,y,0)] = omtauinv*ft0 + tw0r*(omusq);

       double cidot3u = tux; f2[i3D(x,y,1)] = omtauinv*ft1 + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));

       f2[i3D(x,y,1)] = omtauinv*ft1
                  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
       // ... similar expressions for directions 2-4

       cidot3u = tux+tuy; f2[i3D(x,y,5)] = omtauinv*ft5
                  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
       // ... similar expressions for directions 6-8
    }
  }
}
