#include "mesh.h"

unsigned int Mesh::MAXSPACING=128;


// in the STL case if _nx>0 then nx, ny and nz, O_x, O_y, O_z, must be given and will NOT
// be calculated inside loadSTL
// in the DAT case O_x,O_y,O_z are ignored (assumed zero)

Mesh::Mesh(int rank, std::string & filename, INPUTTYPE ext, int & ierr, double scale, bool voxel_based,
      int _nx, int _ny, int _nz, double ox, double oy, double oz

)  {

  meshL=NULL; meshH=NULL;
  
  std::ifstream f;
    
  f.open(filename.c_str(),std::ifstream::in);
  if (!f.is_open()) {
    std::cerr << "Cannot open file "<< filename << " for reading, aborting...\n";
    ierr=-1;
    return;
  }
  
//  if (filename.find(".dat") != std::string::npos) {
  if (ext == DAT) {  
    int x,y,z,ty;
    if (_nx==0) {
      //find box bounds
      nx=ny=nz=0;
      O_x=O_y=O_z=0.;
      while (!f.eof()) {
        f >> x >> y >> z >> ty;
        if (ty>1 && ty<5) {
          nx= MAX(x,nx); ny= MAX(y,ny); nz= MAX(z,nz); 
        }
      }
      //leave free space around the system so not to create false inner regions
      if (!Point::px) nx+=2;
      if (!Point::py) ny+=2;
      if (!Point::pz) nz+=2;
    } else {
      nx=_nx; ny=_ny; nz=_nz;
    }
    this->setBounds(nx,ny,nz);  //also set the point p bounds and increase nx,ny,nz    
    this->setRank(rank); //also set the point p spacing
    //--------
    
    f.clear();                 // clear fail and eof bits
    f.seekg(0, std::ios::beg); //return to the beginning
    while (!f.eof()) {
      f >> x >> y >> z >> ty;
      switch(ty) {
        case 2: this->node(x,y,z).setType(WALL_NODE);  break;
        case 3: this->node(x,y,z).setType(INLET_NODE); break;
        case 4: this->node(x,y,z).setType(OUTLET_NODE); break;
      }
    }

    f.close();
    
  } else if (ext==STL) {

    std::vector<Triangle> triangles;
    triangles.reserve(MAXNUMTRIANGLES);
    // if _nx==0 then nx, ny and nz, O_x, O_y, O_z, are calculated inside loadSTL
    ierr= this->loadSTL(f,&triangles,scale,&std::cout,
      _nx, _ny, _nz, ox, oy, oz      
    );
    f.close();
    if (ierr == -2) std::cerr << "\n\t"<<RED<<"*** Too a small Bounding Box given! Aborting... "<<NOCOL<<"\n\n";
    if (ierr <0) return;
    if (ierr>0) {
      std::cerr << "\n"<<RED<<"WARNING: "<< ierr
      << " triangles are 'unresolved'."<<NOCOL<<"\n\n";      
    }
    this->setBounds(nx,ny,nz);  //also set the point p bounds and increase nx,ny,nz
    this->setRank(rank); //also set the point p spacing
//    N=new (std::nothrow) Node[n];
    /*
    try {
      N=new  Node[n];
    } catch (std::bad_alloc&) {
      std::cerr << "\nNot enough memory to allocate main Box. " << sizeof(Node)*n/1024/1024 <<" Mb needed for "
                   <<n<<" total nodes! Aborting...\n";
      ierr=-1;
      return;
    }
    */
/*    if (N == 0) {
      std::cerr << "Not enough memory to allocate main Box. " << sizeof(Node)*n/1024/1024 <<" Mb needed! Aborting...\n";
      ierr=-1;
      return;
    }
*/
    this->FillTriangles(&triangles,voxel_based);
    
  } else {
    std::cerr << "Input file extension not handled! Aborting...\n";
    ierr=-1;
    return;
  }
  
  ierr=0;
}

  //_nx=LARGEST x coord. of nodes. It can be [0,_nx] (thus the size along x is _nx+1)
Mesh::Mesh(int _nx, int _ny, int _nz, int rank,double ox, double oy, double oz) {
    meshL=NULL; meshH=NULL;
    this->setBounds (_nx, _ny, _nz);
    this->setRank(rank);
    O_x=ox;    
    O_y=oy;    
    O_z=oz;    
}


int Mesh::count(NODETYPE t) const {
  int c;
  c=0;
  for (std::map<INT,Node>::const_iterator it=N.begin(); it!=N.end(); ++it)
    if (it->second.getType()==t) c++;
  return c;
}

bool Mesh::isReallyDead(Point const & p0) const {
    INT j;
    if (!this->isDead(p0)) return false;
    Point p1(nx-1,ny-1,nz-1,spacing.s);
    p1.set(p0.get());
    if (meshL) {
        if (meshL->isOnTheMesh(p1)) {
          if (!meshL->isDead(p1)) return false;
        } else { 
          //point p0 is not on the low res. connected mesh,
          // so see if there is an active point on any extended neighbor location
          // so to be sure that one is not within an interface between two
          // resolutions
          Point newp(nx-1,ny-1,nz-1,spacing.s);
          for (int i=1;i<NPPEXT;i++){
              j=p1.getPop(i);
              if (j<0) continue;
              if (!meshL->isDead(j)) return false;
          }
        }
    }
    if (meshH) return (meshH->isDead(p1)); 
    return true;
}


int Mesh::FindOuterDeadNode(Point & p0) {
  if (this->isDead(p0))  return 0; // node at p0 doesn't exist (dead node)
  const int dist0=this->node(p0).getDist();
  //it makes no sense to start from a negative mindist and going outward.
  if (dist0 < 0 ) return 0;
  Point padj(nx-1,ny-1,nz-1,spacing.s);  
  Point newp(nx-1,ny-1,nz-1,spacing.s,false); //DOES NOT apply periodic conditions
  Point newp2(nx-1,ny-1,nz-1,spacing.s);
  newp=p0;
  for (int i=1; i<NPP; i++) {
    padj=p0.getp(i); // save the adjacent node
    if (padj.isOutOfBounds()) continue; //try another direction
    //do ray-tracing along i direction
    while (true){
      newp2=newp.getp(i); //DOES NOT apply periodic conditions
      if (newp2.isOutOfBounds()) {p0=padj; return (dist0+spacing.s);}//got it!
      if (!this->isReallyDead(newp2)) {newp=p0; break;} //try another direction
      newp=newp2;
    }
  }
  return 0;
}


int Mesh::FindInnerDeadNode(Point & p0) {
  Point newp(nx-1,ny-1,nz-1,spacing.s);
  Point newp2(nx-1,ny-1,nz-1,spacing.s);
  if (this->isDead(p0))  return 0; // node at p0 doesn't exist (dead node)
  const int dist0=this->node(p0).getDist();
  //it makes no sense to start from a positive mindist and going inward.
  if (dist0 > 0 ) return 0;
  newp=p0;
  for (int i=1; i<NPP; i++) {
    //do ray-tracing along i direction
    while (true){
      newp2=newp.getp(i);
      if (newp2.isOutOfBounds()) {
        //try another direction
        newp=p0;
        break;
      }
      if (this->isReallyDead(newp2)) {
        //got it!
        p0=newp2;
        return (dist0-spacing.s);
      }
      if (this->isDead(newp2)) {newp=p0; break;} //got out on the wrong side, try another direction
      if (this->node(newp2).getDist() != dist0) {
        //either got out on the wrong side (>dist0),                                                                  //or there is no exit here. Try another direction.
        newp=p0;
        break;
      } 
      // go ahead only if you are walking through the same layer at mindist==dist0
      newp=newp2;
    }
  }
  return 0;
}

int Mesh::ExtendLayer(Point const & p0) {
    std::stack<Point> ptocheck;
    if (this->isDead(p0))  {
          std::cerr <<"Given node is DEAD...\n";
          return -1; // node at p0 doesn't exist! (dead node)
    }
    Point newp2(nx-1,ny-1,nz-1,spacing.s);
    Point newp3(nx-1,ny-1,nz-1,spacing.s);
    int dist0;
    int previousdist,c;

    //get layer mindist
    dist0=this->node(p0).getDist();
    //determine 'previous' mindist
    previousdist= dist0>0 ? dist0-spacing.s : dist0+spacing.s;
    c=0;
    
    ptocheck.push(p0);
    while (!ptocheck.empty())  {
        //take the first available point from the stack (and remove it)
        p=ptocheck.top();ptocheck.pop();
        //check if the neighbors of p belong to the layer
        for (int i=1; i<NPP; i++) {
            newp2=p.getp(i);
            if (newp2.isOutOfBounds()) continue;
            if (!this->isReallyDead(newp2)) continue;
            //newp2 is a really dead node. Check if it's adjacent to
            //the previous layer (i.e. if at least one of its neighbors
            // belongs to the previous layer)
            for (int j=1; j<NPP; j++) {
                newp3=newp2.getp(j);
                if (newp3.isOutOfBounds() || this->isDead(newp3)) continue;
                if (this->node(newp3).getDist() == previousdist) { 
                    //newp2 is adjacent to the previous layer so
                    // set it as belonging to this layer and put it into the stack
                    // so to check its neighbors, in turn
                    ptocheck.push(newp2);
                    this->node(newp2).SetFluidNode(dist0);
                    //set the node of previous layer as "coated" by this new layer
                    this->node(newp3).setCoated(true);
                    c++;
                    break;
                }
            }
                    
        }
    } 

    return c;
}

int Mesh::RemoveLayer(int dist) {
    std::list< std::map<INT,Node>::iterator > it_toremove;
    std::list< std::map<INT,Node>::iterator >::iterator ip;
    std::map<INT,Node>::iterator it;
    
    for(it=N.begin();  it != N.end(); it++ ){
        if ( it->second.getDist() == dist) it_toremove.push_back(it);
    }
    for(ip=it_toremove.begin();  ip != it_toremove.end(); ip++){ N.erase(*ip); }    
    return (it_toremove.size());
}

int Mesh::AddLayer(int & dist0, int & rank, int dfactor) {
    //if (dist0 == 0)  return 0; // it cannot start from a boundary layer
    int nextdist, nextrank;
    unsigned int nextspacing;
    Mesh * nextmesh;
    
    INT j;
    Point newp(nx-1,ny-1,nz-1,spacing.s);
    std::map<INT,Node> Nnew;
    std::map<INT,Node>::iterator it;
    
    if (rank != spacing.rank) {
        std::cerr << "\nError: wrong mesh rank! Aborting...\n\n";
        return (-1);
    }
    //determine next mesh rank
    
    nextdist= dist0>0 ? dist0+spacing.s : dist0-spacing.s;
    if (dfactor>0) {
      nextrank =  Mesh::CALCULATE_RANK(dfactor,nextdist);
    } else { //force current resolution
      nextrank = rank;
    }
    if (nextrank>rank && meshL) {
        //resolution has to be reduced.
        nextmesh= meshL;
        nextspacing=spacing.s << 1 ; //doubled spacing
        nextdist= dist0>0 ? dist0+nextspacing : dist0-nextspacing;
    } else {
        nextmesh=this;
        //this will cause the default behavior in getPop, i.e. it will consider
        //this mesh spacing
        nextspacing=0;
        nextrank=rank;
    }
    while (true) {
      //repeat this node "layer adding" loop when the lowest resolution
      // reached is actually empty of nodes. In that case a space empty
      // of nodes remains and so there could be saturation problems.
      // So the "layer adding" must be repeated forcing the previous (higher) resolution.
      for (it=N.begin(); it!=N.end(); ++it) {
        //if (it->second.isFluid()) {
          //only treat fluid nodes at mindist==dist0
          if (it->second.getDist() != dist0) continue;
          //find all really dead neighbors and replace it with fluid
          // nodes as belonging to the next layer
          p.set(it->first);
          for (int i=1;i<NPP;i++) {
            j=p.getPop(i,nextspacing);
            if (j<0 ) {
              //std::cerr<< "????????????????????\n";
              continue;
              
            }
            if (Nnew.count(j)>0) continue; //already considered node
            newp.set(j);
            if (this->isReallyDead(newp)) Nnew[j].SetFluidNode(nextdist);
          }
        //}
      }
      
      if (Nnew.size()==0 && nextrank>rank) {
        //Nnew.clear();
        //in this case repeat the adding layer forcing previous resolution.
        nextmesh=this; nextspacing=0;
        nextdist= dist0>0 ? dist0+spacing.s : dist0-spacing.s;
        nextrank=rank;
      } else {
        break; //done
      }
    }
            
    for (it=Nnew.begin(); it!=Nnew.end(); ++it) {
      nextmesh->node()[it->first]=it->second;
    }
    
    dist0=nextdist;
    rank=nextrank;
    
    return Nnew.size();
          
}

int Mesh::AddFirstLayers(int & startingDist, bool inner) {
  std::map<INT,Node>::iterator it;
  int addedNodesOut,addedNodesIn,nextdist,din,dout;
  Point p0(nx-1,ny-1,nz-1,spacing.s);  
  //if (inner && startingDist >0) return (-1); //error
  //if (!inner && startingDist <0) return (-1); //error
  
  addedNodesOut=0;
  int c;
  while (true){
    nextdist=0;
    c=0;
    for (it=N.begin(); it!=N.end() && nextdist ==0; ++it) {
      if (c % 100 == 0) std::cout << "outer layer: "<<c << " / "<<N.size()<<"       \r";
      c++;
      //find an active node with mindist==dist0
      if (it->second.getDist() != startingDist) continue;
      //skip nodes that are adjacent to the already extended layer
      if (it->second.isCoated()) continue;
      p0.set(it->first);
      nextdist=this->FindOuterDeadNode(p0);
    }
    if (nextdist == 0) break; //no outer dead node starting from p0 was found
    this->node(p0).SetFluidNode(nextdist);
    addedNodesOut += this->ExtendLayer(p0)+1;
    dout=nextdist;
  }
  std::cout << "outer layer: done!                                  \n";
  if (addedNodesOut ==0) return (-1);
  //reset "coated" status
  for (it=N.begin(); it!=N.end(); ++it) {
    it->second.setCoated(false);
  }
  if (!inner) {
    startingDist=dout;
    return (addedNodesOut);
  }
  addedNodesIn=0;
  while (true){
    nextdist=0;
    c=0;
    for (it=N.begin(); it!=N.end() && nextdist==0; ++it) {
      if (c % 100 == 0) std::cout << "inner layer: "<<c << " / "<<N.size()<<"       \r";
      c++;
      //find an active node with mindist==dist0
      if (it->second.getDist() != startingDist) continue;
      //skip nodes that are adjacent to the already extended layer
      if (it->second.isCoated()) continue;
      p0.set(it->first);
      nextdist=this->FindInnerDeadNode(p0);
    }
    if (nextdist == 0) break; //no outer dead node starting from p0 was found
    this->node(p0).SetFluidNode(nextdist);
    addedNodesIn += this->ExtendLayer(p0)+1;
    din=nextdist;
  }
  std::cout << "inner layer: done!                                 \n";
//   if (addedNodesIn == 0 ) {
//     startingDist=dout;
//     return (0);
//   }
  startingDist=din;
  this->RemoveLayer(dout);
  return (addedNodesIn);
}

///*
void Mesh::SaturateNodes(bool first) {
  std::map<INT,Node> Nnew;
  std::map<INT,Node>::iterator it;
  
  int i;
  INT j;
  if (spacing.s<2) return;
  if (!meshH) {
      std::cerr << "\nInternal Error! No hi-res mesh defined\n\n";
      return;
  }
  if (first) { //1st pass (accostamento - get closer to meshH nodes)
    //fluid nodes of this mesh are added at neighbors dead nodes that coincides with
    //fluid nodes of meshH (only if these have mindist >= spacing) 
    for (it=N.begin(); it!=N.end(); ++it) {
      if (it->second.isFluid()) {
        p.set(it->first);
        for (i=1;i<NPPEXT;i++) {
          j=p.getPop(i);
          if (j<0) continue;
          if (!this->isDead(j)) continue;
          if (meshH->isDead(j)) continue;
          if (!meshH->node(j).isFluid()) continue;
          //the distance from the wall must be not smaller than spacing
          if ( Nnew.count(j)==0 && abs(meshH->node(j).getDist())>= int(spacing.s) 
          ) {
            Nnew[j].SetFluidNode(meshH->node(j).getDist());
          }
        }

      }
    }
    for (it=Nnew.begin(); it!=Nnew.end(); ++it) {
      N[it->first]=it->second;
    }  
        
    
  } else { //2nd pass - SATURATION of voxels unsaturated ON THE HIGHER RES. MESH SIDE 
    Point neigh(nx-1,ny-1,nz-1,spacing.s);
    int ii;
    bool isUnsaturatedVoxel;
    Point p1(nx-1,ny-1,nz-1,spacing.s);
    for (it=N.begin(); it!=N.end(); ++it) {
      if (it->second.isFluid()) {
        p.set(it->first);
        isUnsaturatedVoxel=true;
        //first, check if it is a voxel (i.e. if it has all vertexes in this mesh)
        for (i=1;i<NPPVOXEL;i++) {
          j=p.getPop(i,0,true);
          if (j<0 || this->isDead(j)) {
              isUnsaturatedVoxel=false;
              break;
          }
        }
        if (isUnsaturatedVoxel) { 
          //second, check if it is an "unsaturated" voxel  on the HIGHER res. mesh side
          isUnsaturatedVoxel=false;
          for (i=0;i<NPPVOXEL;i++) {
            p1=p.getp(i,0,true);
            // check if at least one vertex is on the meshH side and lacks of
            // at least one neigh of its own mesh (unsaturated)
            for (ii=1;ii<NPP;ii++) {
              neigh=p1.getp(ii);
              if (neigh.isOutOfBounds()) continue;
              if (!meshH->isDead(neigh)) {
                isUnsaturatedVoxel=true;
                break;
              }
            }
            if (isUnsaturatedVoxel) {
              isUnsaturatedVoxel= !this->isSaturated(p1);
              /*
              isUnsaturatedVoxel=false;
              for (ii=1;ii<NPP;ii++) {
                neigh=p1.getp(ii);
                if (neigh.isOutOfBounds()) continue;
                if ( this->isDead(neigh)) {
                  //p e' insaturo
                  isUnsaturatedVoxel=true;
                  break;
                }
              }
              */
              if (isUnsaturatedVoxel) break;
            }
          }
          if (isUnsaturatedVoxel) { 
          // fill the unsaturated voxel
            //p.set(it->first);
            //fill vertexes
            for (i=0;i<NPPVOXEL;i++) {
              j=p.getPop(i,0,true);
              if (meshH->isDead(j) && Nnew.count(j)==0) {
                Nnew[j].SetFluidNode(this->node(j).getDist());
              }
            }
            //fill faces
            for (ii=0;ii<NPPVOXEL;ii++) {
              p1=p.getp(ii,spacing.s>>1,true);
              for (i=0;i<NPPVOXEL;i++) {
                j=p1.getPop(i,spacing.s>>1,true);
                if (meshH->isDead(j) && Nnew.count(j)==0) {
                  Nnew[j].SetFluidNode(); //SORRY, the mindist is not set...
                }
              }
            }
          }
        }
          
      }
    }
    for (it=Nnew.begin(); it!=Nnew.end(); ++it) {
      meshH->node(it->first)=it->second;
    }  
  }
  
  
}


void Mesh::printNeighs(int x, int y, int z, std::ostream & out) {
//  int x1,y1,z1;
  p.set(x,y,z);
  Point neigh(nx-1,ny-1,nz-1,spacing.s);
  if (N.count(p.get()) == 0) {
    out << x <<","<<y <<","<<z << ": empty!" << std::endl;
    return; //node doesn't exist
  }
  out<<"Spacing = "<<spacing.s<<"\n";
  for (int i=0; i<NPP; i++) {
    neigh=p.getp(i);
    if (neigh.isOutOfBounds()) {
      out << neigh.x() <<","<<neigh.y() <<","<<neigh.z() << ": OOB" << std::endl;
      continue;
    }
    if (N.count(neigh.get())>0) {
      out << neigh.x() <<","<<neigh.y() <<","<<neigh.z() << ": " <<N[neigh.get()].getType() << std::endl;
    } else {
      if (meshH) {
          if (!meshH->isDead(neigh.get()))
      out << neigh.x() <<","<<neigh.y() <<","<<neigh.z() << ": " <<meshH->node(neigh.get()).getType()<< std::endl;
      }
      if (meshL) {
          if (!meshL->isDead(neigh.get()))
      out << neigh.x() <<","<<neigh.y() <<","<<neigh.z() << ": " <<meshL->node(neigh.get()).getType()<< std::endl;
      }
      out << neigh.x() <<","<<neigh.y() <<","<<neigh.z() << ": empty!" << std::endl;
    }
  }
}

int Mesh::RepairWalls(){
  std::map<INT,Node> Nnew;
  std::map<INT,Node>::iterator it;
  int i;
  Point neigh(nx-1,ny-1,nz-1,spacing.s);
  
  for (it=N.begin(); it!=N.end(); ++it) {
    if (it->second.isFluid()) {
      p.set(it->first);
#if 0
      //first check if there is at least one boundary neigh already
      for (i=1; i<NPP;i++) {
        neigh=p.getp(i,s);
        if (neigh.isOutOfBounds()) continue;
        if (!this->isDead(neigh)) {
          if (this->node(neigh).isBoundary()) break;
        }
      }
      if (i==NPP) continue; //no boundaries, give it up (because dead node
                            // could be actually nodes of the another connected mesh
#endif
      for (i=1; i<NPP;i++) {
        neigh=p.getp(i);
        if (neigh.isOutOfBounds()) continue;
        if (!this->isReallyDead(neigh)) continue;
        if (Nnew.count(neigh.get())==0) {
          Nnew[neigh.get()].setType(WALL_NODE);
        }
      }
    }
  }
  for (it=Nnew.begin(); it!=Nnew.end(); ++it) {
    this->node(it->first)=it->second;
  }  
  
  return Nnew.size();
}

/*
Compenetrate(Mesh & mesh1) {
  std::map<INT,Node>::iterator it1,it;
  
  it1=mesh1.begin();
  it=this->N.begin();
  if (it1->second.getSpacing().s > it->second.getSpacing().s ) { //mesh1 = meshL
    
    
  } else { //mesh1 = meshH
  }  
}
*/


int Mesh::RepairSaturation(){
  std::map<INT,Node> Nnew;
  std::map<INT,Node>::iterator it;
  int i,nextdist,nneigh;
  bool closetomeshL;
  Point neigh(nx-1,ny-1,nz-1,spacing.s);
  //std::cout << ">>>>>>>" << meshL.getSpacing().s << "\n";
  //spL=meshL.getSpacing();
  if (!meshL) {
      std::cerr << "\nInternal Error! No low-res mesh defined\n\n";
      return (-1);
  }
  
  for (it=N.begin(); it!=N.end(); ++it) {
    if (it->second.isFluid()) {
      p.set(it->first);
      if (this->isSaturated(p)) continue;
      //p is an unstaure node.
      //Check if it is on the meshL side
      if (!meshL->isDead(p)) continue; //it coincides with a meshL one... nothing to do
      closetomeshL=false;
      nneigh=0;
      for (i=1; i<NPPEXT;i++) {
        neigh=p.getp(i);
        if (neigh.isOutOfBounds()) continue;
        if (!meshL->isDead(neigh)) {
          closetomeshL=true;
          //if neigh is saturated and coincides with an active node in this mesh,
          //then it can be used to saturate p 
          if (i<NPP && !this->isDead(neigh)) {
            if (meshL->isSaturated(neigh)) nneigh++;
          }
        }
      }
      if (!closetomeshL || nneigh>1) continue; //node p is far from meshL or it can be saturated... nothing to do
      
      //p is on the meshL side and cannot be saturated by means of 
      //the meshL nodes so make it satured            
      nextdist=it->second.getDist();
      nextdist=nextdist>0? nextdist+spacing.s : nextdist-spacing.s;
      for (i=1; i<NPPEXT;i++) {
        neigh=p.getp(i);
        if (neigh.isOutOfBounds()) continue;
        if (this->isDead(neigh) && Nnew.count(neigh.get())==0) {
          Nnew[neigh.get()].SetFluidNode(nextdist);
        }
      }
    }
  }
  for (it=Nnew.begin(); it!=Nnew.end(); ++it) {
    N[it->first]=it->second;
  }  
  return Nnew.size();
}

bool Mesh::isSaturated(Point const & p0){
  int i;
  Point neigh(nx-1,ny-1,nz-1,spacing.s);
  Point p1(p0.x(),p0.y(),p0.z(),nx-1,ny-1,nz-1,spacing.s);

  for (i=1; i<NPP;i++) {
    neigh=p1.getp(i);
    if (neigh.isOutOfBounds()) continue;
    if (this->isDead(neigh)) return false;
  }
  return true;
}

// now returns the number of unresolved triangles
int Mesh::loadSTL(std::ifstream & f, std::vector<Triangle> *triangles, double scale, std::ostream *out,
      int _nx, int _ny, int _nz, double ox, double oy, double oz  
) {

    std::string line;
    Triangle t;
    double volume,area;
    double3 norm;
//    double O_x,O_y,O_z;
    double xmax,ymax,zmax;
    int hires; //number of triangles having an edge smaller than resolution!

    xmax=ymax=zmax=__DBL_MIN__;
    O_x=O_y=O_z=__DBL_MAX__;
    
    area=volume=0.;
    hires=0;
    
    std::getline(f,line);/* discards the first "solid ..." line */
    while( !f.eof() && line.find("endsolid") ==std::string::npos) {

      std::getline(f,line);
      
      if (line.find("facet normal") !=std::string::npos) {
        // ignored!! because normals are often wrong! 
        //t.n.from_string(line,2);
        //if (t.n.Modulus() > 1.0E-6)
        //    t.n.Normalize(); /* safety is never enough... */

        if (f.eof()) {std::cerr << "STL file truncated, missing vertex, aborting...\n"; return (-1);}
        std::getline(f,line); //skip "outer loop"
        if (f.eof()) {std::cerr << "STL file truncated, missing vertex, aborting...\n"; return (-1);}

        std::getline(f,line);
        if (line.find("vertex") ==std::string::npos || f.eof()) {
          std::cerr <<"STL file truncated, missing vertex, aborting...\n";
          return (-1);
        }
        t.v1.from_string(line,1,true);
        if (f.eof()) {std::cerr << "STL file truncated, missing vertex, aborting...\n"; return (-1);}
        
        std::getline(f,line);
        if (line.find("vertex") ==std::string::npos || f.eof()) {
          std::cerr <<"STL file truncated, missing vertex, aborting...\n";
          return (-1);
        }
        t.v2.from_string(line,1,true);
        if (f.eof()) {std::cerr << "STL file truncated, missing vertex, aborting...\n"; return (-1);}

        std::getline(f,line);
        if (line.find("vertex") ==std::string::npos) {
          std::cerr <<"STL file truncated, missing vertex, aborting...\n";
          return (-1);
        }
        t.v3.from_string(line,1,true);
        
        if (f.eof()) {std::cerr << "STL file truncated, missing endloop, aborting...\n"; return (-1);}        
        std::getline(f,line); //skip "endloop"
        if (f.eof()) {std::cerr << "STL file truncated, missing endfacet, aborting...\n"; return (-1);}
        std::getline(f,line); //skip "endfacet"

        
        /* compute normal based on vertices order */ 
        norm=(t.v2 - t.v1).Curl(t.v3 - t.v1);
        norm.Normalize();
        t.n=norm;
        // round to 5 digits so to avoid small round-off difference
        t.n.x=floor(t.n.x*1.e5)/1.e5;
        t.n.y=floor(t.n.y*1.e5)/1.e5;
        t.n.z=floor(t.n.z*1.e5)/1.e5;


        /* take computed normal if it differs more than a thousandth from stored one;
          * s->triangles[c].n is either normalized or with length "almost" zero */
        //if (fabs(1.-t.n*norm) > 1.0E-3) t.n = norm;

        t.v1 = t.v1*scale;
        t.v2 = t.v2*scale;
        t.v3 = t.v3*scale;

        if ((t.v1-t.v2).Modulus2()<1. || (t.v1-t.v3).Modulus2()<1. || (t.v3-t.v2).Modulus2()<1.) hires++;
        
        xmax = MAX(xmax, MAX(t.v1.x, MAX(t.v2.x, t.v3.x)));
        ymax = MAX(ymax, MAX(t.v1.y, MAX(t.v2.y, t.v3.y)));
        zmax = MAX(zmax, MAX(t.v1.z, MAX(t.v2.z, t.v3.z)));
        O_x = MIN(O_x, MIN(t.v1.x, MIN(t.v2.x, t.v3.x)));
        O_y = MIN(O_y, MIN(t.v1.y, MIN(t.v2.y, t.v3.y)));
        O_z = MIN(O_z, MIN(t.v1.z, MIN(t.v2.z, t.v3.z)));
        if (out != NULL) {
          volume += dvol(t);
          area   += dsur(t);
        }
        
        triangles->push_back(t);
      }

        
            
    }
    /*
    if (triangles->size()<4 ) {
        std::cerr << "Missing triangles in STL file, aborting...\n";
        return (-1);
    }
    */
    // leave some free space on the origin side (walls are still unset) 

    Point::px ? O_x=floor(O_x-2.) : O_x=floor(O_x-3.);
    Point::py ? O_y=floor(O_y-2.) : O_y=floor(O_y-3.);
    Point::pz ? O_z=floor(O_z-2.) : O_z=floor(O_z-3.);
    nx = (int) (ceil(xmax-O_x)); ny = (int) (ceil(ymax-O_y)); nz = (int) (ceil(zmax-O_z));
    // leave some free space on the other side (walls are still unset) 
    nx+= Point::px?2:3;
    ny+= Point::py?2:3;
    nz+= Point::pz?2:3;
    if (out != NULL) {
      *out << triangles->size() << " triangles read.\n";
      *out << "Area: "<<area<<" square units (non-scaled "<< area/(scale*scale) << ")\n";
      *out << "Volume: "<<volume<<" cubic units (non-scaled "<< volume/(scale*scale*scale) << ")\n";
    }
    
    if (_nx >0) {
      ox=floor(ox);
      oy=floor(oy);
      oz=floor(oz);
      if (
        ox > O_x ||  _nx-ox < nx ||
        oy > O_y ||  _ny-oy < ny ||
        oz > O_z ||  _nz-oz < nz
      ) return (-2);//the user-given box is too small, it doesn't include the system
      O_x=ox;
      O_y=oy;
      O_z=oz;
      nx=_nx-ox; //WARNING this->setBounds MUST be called with this nx,ny,nz after exiting this routine
      ny=_ny-oy;
      nz=_nz-oz;
        
    }
    //traslate triangles and box so that they lie within the positive octant 
    for (std::vector<Triangle>::iterator it=triangles->begin(); it != triangles->end(); ++it) {
      it->v1.x -= O_x; it->v2.x -= O_x; it->v3.x -= O_x;
      it->v1.y -= O_y; it->v2.y -= O_y; it->v3.y -= O_y;
      it->v1.z -= O_z; it->v2.z -= O_z; it->v3.z -= O_z;
    }
    
    return (hires);
}

void Mesh::FillTriangles(std::vector<Triangle> *triangles, bool voxel_based) {
  int i,j,ii;
  int c=0;
  Point pv(nx-1,ny-1,nz-1,spacing.s);  
  double3 AB,AC,P;
  double ABm,ACm;
  
//  Point neigh(nx,ny,nz);
  
  for (std::vector<Triangle>::iterator it=triangles->begin(); it != triangles->end(); ++it) {  
    
    AB=it->v2-it->v1;
    AC=it->v3-it->v1;
    ACm=AC.Modulus();
    ABm=AB.Modulus();
    //scan the triangle surface on a 2-D 1-spaced grid
    for (i=0; i<=int(ACm);i++) {
      for (j=0; j<=int(ABm-double(i)*ABm/ACm);j++) {
        P=AC.Versor()*double(i)+AB.Versor()*double(j)+it->v1;
        // "round" P to a mesh point (the voxel origin)
        p.set((int)(P.x)-((int)(P.x)%spacing.s),
              (int)(P.y)-((int)(P.y)%spacing.s),
              (int)(P.z)-((int)(P.z)%spacing.s));
        if (p.isOutOfBounds()) continue;
        if ( this->isDead(p) ) {
          //put a wall node in the voxel origin
          this->node(p).setType(WALL_NODE);
          c++;
          //this node gets the triangle normal
          this->node(p).setVelDir(it->n.to_string());
        }
        if (voxel_based) {
          //put a wall on any other dead voxel vertex
          for (ii=1;ii<8;ii++) {
            pv=p.getp(ii,spacing.s,true);
            if ( this->isDead(pv) ) {
              if (pv.isOutOfBounds()) continue;
              this->node(pv).setType(WALL_NODE);
              this->node(pv).setVelDir(it->n.to_string());
              c++;
            }
          }
        }
        
      }
    }
  }
    
  std::cout << c << " WALL nodes generated from STL.\n";
  
}

int Mesh::IntersectIOnodes(Mesh & ibox,NODETYPE iotype){
  int inx,iny,inz,c=0;
  if (spacing.s != ibox.getSpacing().s) {
      return (0);
      std::cerr <<"\n"<<RED<< "WARNING: I/O box not of the right spacing."<<NOCOL<<"\n\n";
  }
  ibox.getBounds(inx,iny,inz);
  std::map<INT,Node>::iterator it;
  Point pi(inx,iny,inz,spacing.s);  
  
  //scan active nodes within the given ibox
  for (it=ibox.node().begin(); it!=ibox.node().end(); ++it) {
    //set p as the point of this box that is at the same "original" location
    //as the ibox point pi
    pi.set(it->first);
    p.set((int)((double)(pi.x())+ibox.getO_x()-O_x),
            (int)((double)(pi.y())+ibox.getO_y()-O_y),
            (int)((double)(pi.z())+ibox.getO_z()-O_z));
    if (p.isOutOfBounds()) continue;
    this->node(p).setType(iotype);
    c++;
    if (iotype==INLET_NODE && it->second.getType()==WALL_NODE)
        this->node(p).setVelDir(it->second.getVelDir());
  }
         
  return (c);
}

bool Mesh::SetIOnodesSearchDir(){
  bool iexp=false;
  bool oexp=false;
  std::map <INT,Node>::iterator it;
  //double3 searchdir,veldir;
  
  Point neigh(nx-1,ny-1,nz-1,spacing.s);
  //scan nodes within the ibox
  for (it=N.begin(); it!=N.end(); ++it) {
    p.set(it->first);
    if (it->second.getType()==INLET_NODE || it->second.getType()==OUTLET_NODE) { 
      for (int i=1; i<NPPLESS;i++) {
        neigh=p.getp(i);
        if (neigh.isOutOfBounds() || this->isDead(neigh)) continue;
        if (this->node(neigh).isFluid()) {
            iexp=iexp || it->second.getType()==INLET_NODE;
            oexp=oexp || it->second.getType()==OUTLET_NODE;
            std::stringstream iss;
            iss <<p.x()-neigh.x() <<" "<< p.y()-neigh.y()<<" "<<  p.z()-neigh.z();
            it->second.setSearchDir(iss.str());
            //searchdir.from_string(iss.str());
            //it->second.setSearchDir(searchdir.to_string());
            /*
            //check and correct the verse of veldir
            if (it->second.getType()==INLET_NODE) { 
            veldir.from_string(it->second.getVelDir());
            //this is to adjust veldir that is assumed to go from inlet nodes towards fluid nodes,
            //thus its scalar product with searchdir must be negative.
            if (veldir*searchdir>0.) {
                it->second.setVelDir((veldir*(-1.)).to_string()); //veldir e' verso fuori! Rigirala!
            }
    //                std::cout << searchdir << " | " << veldir <<std::endl; 
            }
            */
            break;
        }
      }
    }
  }
  return (iexp && oexp);
}

int Mesh::AddNodesFromMesh(Mesh & source, int dist_offset, NODETYPE NodesToBeAdded){

  /*
  if (spacing.s != source.getSpacing().s) {
      std::cerr <<"\n"<<RED<< "WARNING: source mesh not of the right spacing."<<NOCOL<<"\n\n";
      return (0);
  }
  */
  int inx,iny,inz,c=0,d;
  std::map <INT,Node>::iterator it;
  source.getBounds(inx,iny,inz);
  dist_offset=abs(dist_offset);
  Point ps(inx,iny,inz,spacing.s);  
  
  //scan active nodes within the given source
  for (it=source.node().begin(); it!=source.node().end(); ++it) {
        if (NodesToBeAdded != DEAD_NODE && it->second.getType() != NodesToBeAdded) continue;
        //set p as the point of this mesh that is at the same "original" location
        //as the source point ps
        ps.set(it->first);
        p.set((int)((double)(ps.x())+source.getO_x()-O_x),
              (int)((double)(ps.y())+source.getO_y()-O_y),
              (int)((double)(ps.z())+source.getO_z()-O_z));
        if (p.isOutOfBounds()) continue;
        this->node(p).setType(it->second.getType());
        d=it->second.getDist();
        this->node(p).setDist(d>0?d+dist_offset:d-dist_offset);
        //this->node(p).setDist(10000);
        c++;
  }
         
  return (c);
}

//---- Static methods follow


int Mesh::CALCULATE_RANK(int dfactor, int dist) {
  unsigned int s;
  //Spacing sp;
  int i;
  
//  calculate spacing
  s=(unsigned int) ((unsigned int)((abs(dist)-1)/dfactor)+1);
  s=MIN(s,Mesh::MAXSPACING);
  // the following is equivalent to 2**int(log_2(sp))
  for (i=0; s > 0; i++) {
    s = s >> 1; //one bit shift right
  }
  //sp.rank=i-1;
  //sp.s = (unsigned int)(1 << sp.rank); // rank bits shift left
  return i-1;
}

double Mesh::dvol(Triangle  const &t) {

        double    bz;
        double3   v1proj, v2proj;
        double    tri_proj_area;
        double3    vprod;

        // z coo of baricentrum //
        bz = (t.v1.z + t.v2.z + t.v3.z)/3.;

        // v2-v1 projection on plane z=0 //
        v1proj = t.v2- t.v1;
        v1proj.z = 0.;
        
        // v3-v1 projection on plane z=0 //
        v2proj = t.v3 - t.v1;
        v2proj.z = 0.;       

        tri_proj_area=v1proj.Curl(v2proj).z;

        tri_proj_area = fabs(0.5*tri_proj_area)*bz;

        return tri_proj_area * copysign(1., t.n.z);
}

double Mesh::dsur(Triangle const &t) {
        return ((t.v2 - t.v1).Curl(t.v3 - t.v1).Modulus()*0.5);
}
