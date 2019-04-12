/////////////////--Node--////////////////
#ifndef NODE
#define NODE

#include "genmesh.h"
class Node {

protected:
  NODETYPE type;
  /*
  * dist= node mindist.
  * WARNING: dist>0 means OUTER
  *          dist<0 means INNER
  *          dist==0 means Boundary
  */
  int dist;

  //bool insature;
  //bool outer,explored;
  int indice; //used only for i/o nodes (index in the ios file)
  bool coated; //used in mesh::AddFirstLayer and mesh::ExtendLayer
  std::string veldir;
  std::string searchdir;

/*  
  //address of a node beloging to a layer at larger mindist  
  Node * nextLayerNode;
  //address of a node beloging to a layer at smaller mindist  
  Node * prevLayerNode;
*/

public:
  Node(){
    //outer=explored=false;
    dist=0;    
    indice=0;
    coated=false;
    //dummy values
    veldir="0. 0. 0.";
    searchdir="0 0 0";
    //nextLayerNode=prevLayerNode= NULL;
  }
  
  
  inline void SetFluidNode(int d=0){type = FLUID_NODE; dist=d;}

  inline bool isBoundary() const {
    return (type == WALL_NODE || type == INLET_NODE || type == OUTLET_NODE);
  }
  inline bool isDead() const {return (type == DEAD_NODE); }

  inline bool isFluid() const { return (type == FLUID_NODE); }

  //inline bool isOuter() const {return (outer);}
  //inline void setOuter() {outer=true;}
  //inline bool isExplored() const {return (explored);}
  //inline void setExplored() {explored=true;}

  inline int getDist() const {return (dist); }
  inline void setDist(int d) {dist=d;}
  
  inline void setIndex(int i) {indice=i;}
  inline int getIndex() const {return (indice);}

  inline void setVelDir(std::string const &  s) {veldir=s;}
  inline std::string getVelDir() const {return (veldir);}
  inline void setSearchDir(std::string const &  s) {searchdir=s;}
  inline std::string getSearchDir() const {return (searchdir);}
  
  inline NODETYPE getType() const {return type;}

  inline void setType(NODETYPE t) {type=t;}

  inline void setCoated(bool c) {coated=c;}
  inline bool isCoated() const {return coated;}
  
/*
  inline void setPrevLayer(Node * prev) {prevLayerNode=prev;}
  inline Node * getPrevLayer() const {return prevLayerNode;}
  inline void setNextLayer(Node * next) {nextLayerNode=next;}
  inline Node * getNextLayer() const {return nextLayerNode;}
*/  
  //inline bool isInsature() const {return insature;}
  //inline void setInsature(bool i) {insature=i;}
  
  Node & operator=(Node const & source);
  //inline void SetDeadNode() {type=DEAD_NODE; spacing.s=0;spacing.rank=-1;}
};

#endif
