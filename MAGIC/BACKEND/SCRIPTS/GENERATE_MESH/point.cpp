#include "point.h"

/**
 * Give a new point (with spacing <sp>, or that given at constructor [default])
 * of the neighbor of current position (set by the set method) along the 'speed direction' <i>.
 * (under D3Q19 scheme).
 * No out of bounds check is applied
 *
 */ 
Point Point::getp(int i, unsigned int sp, bool voxel) const {
    int x1,y1,z1;
    this->getPop(i,x1,y1,z1,sp,voxel);
    if (sp==0) sp=spacing;
    return Point(x1,y1,z1,nx-1,ny-1,nz-1,sp);
}

/**
 * Give the internal 1-d node index (with spacing <sp>, or that given at constructor [default])
 * of the neighbor of node at pos (set by the set method) along the 'speed direction' <i>.
 * (under D3Q19 scheme).
 * Give -1 if outside box bounds
 *
 */ 
INT Point::getPop(int i, unsigned int sp, bool voxel) const {
  int x1,y1,z1;
  this->getPop(i,x1,y1,z1,sp,voxel);
  if (x1<1 || y1<1 || z1<1 
      ||  x1 >= nx ||  y1 >= ny ||  z1 >= nz) return (-1);
  return (INT(x1)+INT(y1)*INT(nx)+INT(z1)*INT(nx)*INT(ny));
}

/**
 * Give a new point (with spacing <sp>, or that given at constructor [default])
 * of the neighbor of current position (set by the set method) along the 'speed direction' <i>.
 * (under D3Q19 scheme). If <sp> is given, this point is forced to be ON the mesh with spacing <sp>.
 * No out of bounds check is applied.
 * If checkPeriodic is false then does NOT apply periodic condition; it is useful when checking
 * that the neighbor is outside the box (e.g. when finding outer node!)
 *
 */ 
void Point::getPop(int i, int & x1, int & y1, int & z1, unsigned int sp, bool voxel) const {
  if (voxel) i=voxeldir[i]; //i MUST BE < 8 !
  if (i==0) {x1=xx; y1=yy; z1=zz; return;}
  if (sp==0) {
    x1=xx+dx[i]*spacing;
    y1=yy+dy[i]*spacing;
    z1=zz+dz[i]*spacing;
  } else {
    x1=xx+dx[i]*sp;
    y1=yy+dy[i]*sp;
    z1=zz+dz[i]*sp;
  }
  if (checkPeriodic) {
    if (Point::px) x1= ((x1+nx-2) % (nx-1)) + 1;
    if (Point::py) y1=((y1+ny-2) % (ny-1)) + 1;
    if (Point::pz) z1= ((z1+nz-2) % (nz-1)) + 1;
  }
  // THIS IS IMPORTANT IN THE TRANSITION BETWEEN ONE RESOLUTION TO ANOTHER.
  // NOTE THAT if x1 % sp == 0 this is irrelevant
  if (sp>0) {
    if (x1>0 && (x1 & (sp-1)) != 0) {
      x1=x1 & ~(sp-1); // equiv. to x1= x1 - (x1%sp);
      if (x1<xx) x1 += sp;
    }
    if (y1>0 && (y1 & (sp-1)) != 0) {
      y1=y1 & ~(sp-1); // equiv. to x1= x1 - (x1%sp);
      if (y1<yy) y1 += sp;
    }
    if (z1>0 && (z1 & (sp-1)) != 0) {
      z1=z1 & ~(sp-1); // equiv. to x1= x1 - (x1%sp);
      if (z1<zz) z1 += sp;
    }
  }
}


/**
 * Set periodicity and compute neighbors distances
 * 
 */
void Point::setPeriodic(bool _px, bool _py, bool _pz) {
    px=_px; py=_py; pz=_pz;
    dist[0]=0.;
    for (int i=1; i<NPPEXT; i++) dist[i]=sqrt(dx[i]*dx[i]+dy[i]*dy[i]+dz[i]*dz[i]);
}
//                     0 1 2 3  4  5  6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
const int Point::dx[]={0,1,0,0,-1, 0, 0,1,1,0,-1,-1, 0, 1, 1, 0,-1,-1, 0,-1, 1,-1, 1,-1, 1,-1, 1};
const int Point::dy[]={0,0,1,0, 0,-1, 0,1,0,1, 1, 0,-1,-1, 0, 1,-1, 0,-1,-1,-1, 1, 1,-1,-1, 1, 1};
const int Point::dz[]={0,0,0,1, 0, 0,-1,0,1,1, 0, 1, 1, 0,-1,-1, 0,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1};

const int Point::voxeldir[]={0,1,2,3,7,8,9,26}; //voxel vertex directions

double Point::dist[27];

bool Point::px=false;
bool Point::py=false;
bool Point::pz=false;

