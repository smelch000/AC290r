//////////////-----BOX-----//////////////
// to be used with the "overall" mesh including all highest resolution node
// (read from file).

#ifndef BOX
#define BOX

#include "mesh.h"
#include <fstream>
#include <string>
#include <vector>
#include <stack>

class Box : public Mesh {
protected:
  INT n;
  double O_x,O_y,O_z; //local origin
  Node * N; //use an array instead of a map (the map name of class Mesh is hidden)
  bool isBoundingBoxGiven;
  struct Triangle {
          double3 v1;
          double3 v2;
          double3 v3;
          double3 n;
  };

  static double dvol(Triangle const &t);
  static double dsur(Triangle const &t);  
/*
 * return true if the projection of <this> point lies within the triangle ABC and
 * the point is at a distance <toll from the triangle plane (if oriented==true, th point must be outside)
 * 
 */
  static bool isOnTriangle(Point const & p0, Triangle const &t,double toll,bool oriented);
  
  
public:
  

  //_nx=LARGEST x coord. of nodes. It can be [0,_nx] (thus the size along x is _nx+1)
  Box(int _nx, int _ny, int _nz, double ox, double oy, double oz):Mesh(_nx,_ny,_nz,0) {
    n=nx*ny*nz;
    N=new Node[n];
    isBoundingBoxGiven=true;
    O_x=ox;    
    O_y=oy;    
    O_z=oz;    
  }
  
// in the STL case if _nx>0 then nx, ny and nz, O_x, O_y, O_z, must be given and will NOT
// be calculated inside loadSTL
  Box(std::string & filename, INPUTTYPE ext, int & ierr, double scale=1.,bool voxel_based=true,
    int _nx=0, int _ny=0, int _nz=0, double ox=0., double oy=0., double oz=0.
  );
  
  virtual ~Box() { delete [] N;}
  
  //point the i-th node 
  inline Node & node(INT i) {return (this->N[i]);}

  //point the node at the position p
  inline Node & node(Point const & p) {return (this->N[INT(p.x())+INT(p.y())*INT(nx)+INT(p.z())*INT(nx)*INT(ny)]);}
  
  inline Node & node(int _x, int _y, int _z) { return (this->N[INT(_x)+INT(_y)*INT(nx)+INT(_z)*INT(nx)*INT(ny)]); }

  inline INT getSize() const {return n;}
  
  inline double getO_x() const {return O_x;}
  inline double getO_y() const {return O_y;}
  inline double getO_z() const {return O_z;}

  int count(NODETYPE t) const;

  /*
   * find the first "missed" node, i.e. the first encountered fluid node
   * whose dist was not set (i.e. that were missed by FindMinDist.) AND
   * is adjacent to either a boundary or
   * a node with dist set.
   * Return its position. Return "zero" position if it isn't found.
  */
  
  Point FindMissed();

  /*
  * find an outer node on the bounding box faces and return it in p0.
  * Return false if it couldn't.
  * WARNING: the given p0 must be one of the bounding box vertex
  * 
  */
  bool FindOuterNode(Point & p0);
  
  void FindExternal(Point const & p0, bool watchout_holes);

  void FindAdjacent(Point const & p0);
  
  void FindMinDist(Point const & p0,  Point & last, int & totnodes);

  // mark internal nodes as fluid
  // WARNING: must be called AFTER FindExternal
  // If internal == false then do the opposite, i.e. mark external nodes as fluid
  int FindInternal (Point & p0,  bool internal=true);

  //return the versor (pointing inwards) that is normal to the iso-distance surface passing on the
  // given node   
  //void getNormDirection(int x, int y, int z, double * normdir, int * searchdir);
  
  //load STL triangles from f after multiplying them by scale. Return -1 on error.
  int loadSTL(std::ifstream & f, std::vector<Triangle> *triangles, double scale, std::ostream *out=NULL,
    int _nx=0, int _ny=0, int _nz=0, double ox=0., double oy=0., double oz=0.    
  );

  /*
   *  fill the STL triangles surface with WALL nodes and assign them the various dirs. (as strings)
   * if voxel_based==true, then put walls on voxels on the triangle surface
   */
  void FillTriangles(std::vector<Triangle> *triangles, bool voxel_based); 

/*
 * set as <iotype> nodes *all* those nodes (i.e., at this stage this box fluid nodes
 * aren't found yet, there is no inside/outside distinction)
 * of <this> box that are located at the same position of any ACTIVE nodes of the given ibox
 * Moreover, if iotype=INLET_NODE, the veldir of ibox WALL nodes (because only these have it) is copied.
 * The searchdir on new iotype nodes  will be found and set later on when the fluid
 * nodes of this box will be available.
 * Return the number of iotype nodes set.
 */ 
  int IntersectIOnodes(Box & ibox,NODETYPE iotype);

/*
 * find and set all i/o nodes searchdir and check
 * and correct inlet nodes veldir verse.
 * Return false if inlet OR outlet are not in touch with fluid nodes!
 *
 */ 
  bool SetIOnodesSearchDir();

/*
 * set the maxspacing of fluid nodes of <this> box that are located at the same position
 * of any ACTIVE nodes of the given ibox.
 * Return the number of involved nodes.
 */ 
  int ForceHighResolution(Box & ibox);

/*
 * set all nodes spacing.
 * to be called AFTER FindMinDist!
 */
  void SetSpacing(int dfactor);

/*
 * return true only if at least a p0 node neighbor is of type <tipo>.
 */
  bool NodeIsAdjacentToType(Point const & p0, NODETYPE tipo, int npp=NPP);
  
  /*
   * just for diagnostics
   */
  void Dump(std::string const &filename, std::string const & id);

  
};

#endif