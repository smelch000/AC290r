//////////////////////--MESH--/////////////
#ifndef MESH
#define MESH

#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <list>
#include "genmesh.h"
#include "point.h"
#include "node.h"

class Mesh {
protected:
  //WARNING: these are the total number of nodes along dir.s NOT the bounding box.
  // i.e. nodes admitted positions are [1,nx-1]x[1,ny-1]x[1,nz-1]
  int nx,ny,nz;
  Point p;
  std::map<INT, Node> N;
  Spacing spacing;
  Mesh * meshH; //higher res. connected mesh
  Mesh * meshL; //lower res. connected mesh
  double O_x,O_y,O_z; //local origin

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
  
  /**
  * Starting from a node having a certain mindist (must be >=0) find an adjacent
  * (REALLY) dead node that is connected (via simple ray-tracing) to the box boundary,
  * so that it is certainly an outer node with a higher mindist.
  * Return its mindist and its position in p0.
  * Return 0 if it couldn't find it.
  * 
  * NOTE: this routine is based on a CONNECTIVITY approach, so (topologically)
  * unconnected regions must be explored further.
  * 
  */
  int FindOuterDeadNode(Point & p0);

  /**
  * Starting from a node at a given mindist (must be <=0) go THROUGH the nodes having
  * the same mindist until get out "on the other side" and hit a (really) DEAD node,
  * that is supposed to be INNER (so to which a lower mindist should be given).
  * Return its mindist and its position in p0.
  * Return 0 if it couldn't find it.
  * 
  * WARNING: this must be executed AFTER at least one outer layer is found and marked
  * NOTE: this routine is based on a CONNECTIVITY approach, so (topologically)
  * unconnected regions must be explored further.
  * 
  */
  int FindInnerDeadNode(Point & p0);

  /**
   * Starting from a fluid node at point <p0> of a given layer at a certain mindist,
   * set fluid nodes all those *really* dead (i.e. that do not belong
   * to ANY possibly connected mesh -- meshL or meshH )that are *connected* to p0 
   * and are on the same layer, and set the same mindist for them.
   * Return the number of added nodes or -1 on error.
   * Furthermore, it sets all nodes of the *previous* layer that are connected to p0 as
   * "coated". This is to speedup the process in AddFirstLayers routine because when
   * other starting point are considered (to handle possibly unconnected volumes)
   * such coated nodes can be a-priori skipped.
   * 
   * This is used to construct the FIRST layer starting from a node that is
   * surely outer (inner) and whose position was already found by calling
   * FindOuterDeadNode (FindInnerDeadNode)
   * and then by setting it to fluid node with the proper mindist.
   * 
  * NOTE: this routine is based on a CONNECTIVITY approach, so (topologically)
  * unconnected regions must be explored further.
   */
  int ExtendLayer(Point const & p0);
  
  /**
   * Remove from this mesh all nodes having mindist=dist and return the amount.
   *
   */ 
  int RemoveLayer(int dist);
  
  bool isSaturated(Point const & p0);

public:
/*
 * return the rank and spacing corresponding to a given mindist
 * 
 */
  static int CALCULATE_RANK(int dfactor, int dist);
  //maximum spacing
  static unsigned int MAXSPACING;


  //_nx=LARGEST x coord. of nodes. It can be [0,_nx] (thus the size along x is _nx+1)
  Mesh(int _nx, int _ny, int _nz, int rank,double ox=0., double oy=0., double oz=0.);

// in the STL case if _nx>0 then nx, ny and nz, O_x, O_y, O_z, 
// must be given and will NOT be calculated inside loadSTL
  Mesh(int rank, std::string & filename, INPUTTYPE ext, int & ierr, double scale=1.,bool voxel_based=true,
    int _nx=0, int _ny=0, int _nz=0, double ox=0., double oy=0., double oz=0.
  );
  
  inline void getBounds (int & _nx, int & _ny, int & _nz)  const {
    _nx=nx-1;_ny=ny-1;_nz=nz-1;
  }
  inline void setBounds (int _nx, int _ny, int _nz) {
    nx=_nx+1;ny=_ny+1;nz=_nz+1;
    p.setBounds(_nx,_ny,_nz);
  }
  inline void setMeshL(Mesh * m){meshL=m;}
  inline void setMeshH(Mesh * m){meshH=m;}
  
  //point the nodes 'map' 
  inline std::map<INT, Node> & node() {return (N);}
  //point the i-th node 
  inline Node & node(INT i) {return (N[i]);}

  //point the node at the position p
  inline Node & node(Point const & p) {return (
      N[INT(p.x())+INT(p.y())*INT(nx)+INT(p.z())*INT(nx)*INT(ny)]
    );
  }
  
  //point the node at the position (_x,_y,_z)
  inline Node & node(int _x, int _y, int _z) { return (N[INT(_x)+INT(_y)*INT(nx)+INT(_z)*INT(nx)*INT(ny)]); }

  inline INT getSize() const {return N.size();}

  inline Spacing getSpacing() const { return spacing;}

  inline void setRank(int r) { spacing.rank=r; spacing.s=1<<r; p.s()=spacing.s;}
  
  inline bool isDead(int _x, int _y, int _z) {
    Point p0(_x,_y,_z,nx-1,ny-1,nz-1,spacing.s);
    return (this->isDead(p0.get()));
  }

  inline bool isDead(Point const & p0) const {return (this->isDead(p0.get()));}

  inline bool isDead(INT j) const {
    std::map<INT,Node>::const_iterator it=N.find(j);
    if (it == N.end() ) return true;
    return (it->second.isDead());
  }

  /**
   * Return true is the node at the given <p0> position is dead in this mesh
   * AND in both connected meshes (meshL and meshH).
   * 
   */
  bool isReallyDead(Point const & p0) const;
  
/* 
 * replace any dead node possibly found around fluid nodes, with WALL nodes.
 * Only REALLY dead nodes are considered (i.e. those that correspond to a dead node
 * also in the connected meshes). 
 * Return the number of replacements
 */
  int RepairWalls();

/*
 * replace any dead node possibly found around fluid UNSATURED and UNCONNECTED nodes,
 * of this mesh with FLUID nodes.
 * UNCONNECTED means that there are no neighbors of lower res. mesh (meshL). 
 * Return the number of replacements (-1 in case of error).
 */
  int RepairSaturation();

  //void checkInsatureNodes();

  int count(NODETYPE t) const;
  
  void SaturateNodes(bool firsttime);
  
  void printNeighs(int x, int y, int z, std::ostream & out);


  inline double getO_x() const {return O_x;}
  inline double getO_y() const {return O_y;}
  inline double getO_z() const {return O_z;}
  
  /**
   * Starting from the nodes with given mindist <dist0> belonging to this mesh
   * with rank <rank>, find and set as fluid nodes all those really dead
   * at the 'next' layer (outer ot inner depending on the sign of the <dist0>).
   * Return the number of added nodes, in <dist0> the mindist of them and in
   * <rank> the rank of the connected mesh owning the next layer.
   * <dfactor> is needed to handle multimesh.
   * NOTE: if <dfactor> is <=0 then the rank of the next layer will be imposed
   * as equal to -<dfactor>.
   * 
   * NOTE: Connectivity problem FREE.
   * This routine is NOT based on a CONNECTIVITY approach, so (topologically)
   * unconnected regions are well treated.
   * It can be called repetitively for each layer to be added. 
   */
  int AddLayer(int & dist0, int & rank, int dfactor);

  //load STL triangles from f after multiplying them by scale. Return -1 on error.
  int loadSTL(std::ifstream & f, std::vector<Triangle> *triangles, double scale, std::ostream *out=NULL,
    int _nx=0, int _ny=0, int _nz=0, double ox=0., double oy=0., double oz=0.    
  );

  /**
   * Fill the two adjacent layers of that with the given mindist = <dist0> that
   * are proved to be outer and inner, with fluid nodes.

   * Return the total no. of nodes added and in <dist0> the mindist of the added inner (or outer)
   * node (useful to start with for adding subsequent layers) if ) if <inner>==true (or false).
   * Return -1 if any of the two layers was not fill.
   * 
   * NOTE: Connectivity problem FREE (it can be called just once for any <inner> case).

   * WARNING: to be called before any AddLayer call !
   * WARNING: TO BE CALLED with a <dist0> == 0  !!
   * 
   */
  int AddFirstLayers(int & dist0, bool inner);
  
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
  int IntersectIOnodes(Mesh & ibox,NODETYPE iotype);

/* 
 * find and set all i/o nodes searchdir and check
 * and correct inlet nodes veldir verse.
 * Return false if inlet OR outlet are not in touch with fluid nodes!
 *
 */ 
  bool SetIOnodesSearchDir();

/* 
 * 
 * Create and set the nodes of <this> mesh that are located at the same position
 * of any nodes of the given mesh <source> and of the given type <NodesToBeAdded>
 * as the same type of node.
 * If <NodesToBeAdded>==DEAD_NODE (the default) then any type is considered.
 * Coordinate ref. frames are taken into account.
 * <surce> must be of the same spacing.
 * WARNING: nodes that are outside this mesh bounding-box are NOT added.
 * WARNING: Mindist is 0 for added nodes.
 * Return the number of added nodes.
 */ 
  int AddNodesFromMesh(Mesh & source, int dist_offset=0, NODETYPE NodesToBeAdded=DEAD_NODE);

  inline bool isOnTheMesh(Point const & p0) const {return (p0.x() % spacing.s ==0 && p0.y() % spacing.s ==0
                                                    && p0.z() % spacing.s ==0);}
  
  /*
   * just for diagnostics
   */
  void Dump(std::string const &filename, std::string const & id);


  
};

#endif
