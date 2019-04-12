#ifndef POINT
#define POINT

#define NPP 19
#define NPPLESS 7
//#define NPPEXT 19
#define NPPEXT 27
#define NPPVOXEL 8

#include "genmesh.h"

#include <cmath>
class Point {

protected:
  int xx,yy,zz;
  int nx,ny,nz; //these are NOT bounds. These are the bounds + 1
  unsigned int spacing;
  bool checkPeriodic;
  
  static const int voxeldir[NPPVOXEL];
  static const int dx[NPPEXT];
  static const int dy[NPPEXT];
  static const int dz[NPPEXT];
  static double dist[NPPEXT];
public:  

  inline Point(bool _checkPeriodic=true){checkPeriodic=_checkPeriodic;}
  inline Point(int _nx, int _ny, int _nz, unsigned int sp, bool _checkPeriodic=true){
      setBounds(_nx, _ny, _nz); spacing=sp; checkPeriodic=_checkPeriodic;}
  inline Point(int _x, int _y, int _z,int _nx, int _ny, int _nz, unsigned int sp, bool _checkPeriodic=true){
      xx=_x; yy=_y; zz=_z; setBounds(_nx, _ny, _nz); spacing=sp; checkPeriodic=_checkPeriodic;}

  static bool px,py,pz; // periodicity

  inline void setBounds(int _nx, int _ny, int _nz){nx=_nx+1; ny=_ny+1; nz=_nz+1;}
  static void setPeriodic(bool _px, bool _py, bool _pz);
  inline double getSpatialDistance(int i) { return (Point::dist[i]*spacing);}
  
  inline void set(int _x, int _y, int _z){xx=_x; yy=_y; zz=_z;}
  inline void set(INT i){xx=i%INT(nx); yy=((i-INT(xx))/INT(nx))%INT(ny); zz=(i-INT(xx)-INT(yy)*INT(nx))/(INT(nx)*INT(ny));}
  inline INT get() const {return (INT(xx)+INT(yy)*INT(nx)+INT(zz)*INT(nx)*INT(ny));}
  inline int x () const {return xx;}
  inline int y () const {return yy;}
  inline int z () const {return zz;}
  inline int x (double O) const {return (int)((double)(xx)+O);}
  inline int y (double O) const {return (int)((double)(yy)+O);}
  inline int z (double O) const {return (int)((double)(zz)+O);}
  inline unsigned int & s()  {return spacing;} //read-write access to spacing
  INT getPop(int i, unsigned int sp=0, bool voxel=false) const;
  void getPop(int i, int & x1, int & y1, int & z1, unsigned int sp=0, bool voxel=false) const;
  Point getp(int i, unsigned int sp=0, bool voxel=false) const;
  //WARNING: assignment doesn't copy the box bounds and the SPACING!
  inline Point & operator=(Point const & p) {if (this !=&p) {xx=p.x(); yy=p.y(); zz=p.z();} return *this;}
  inline bool operator==(Point const & p) const {return (xx==p.x() && yy==p.y() && zz==p.z());}
  inline bool operator!=(Point const & p) const {return (xx!=p.x() || yy!=p.y() || zz!=p.z());}
  inline bool isOutOfBounds() const { return (xx<1 || yy<1 || zz<1 ||  xx >= nx ||  yy >= ny ||  zz >= nz);}
};

#endif
