#include "box.h"
// in the STL case if _nx>0 then nx, ny and nz, O_x, O_y, O_z, must be given and will NOT
// be calculated inside loadSTL
// in the DAT case O_x,O_y,O_z are ignored (assumed zero)

Box::Box(std::string & filename, INPUTTYPE ext, int & ierr, double scale, bool voxel_based,
      int _nx, int _ny, int _nz, double ox, double oy, double oz

) {

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
    p.setBounds(nx,ny,nz);
    nx++;ny++;nz++;
    //--------
    
    n=INT(nx)*INT(ny)*INT(nz);
    N=new Node[n];
    
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
    p.setBounds(nx,ny,nz);
    nx++;ny++;nz++;
    n=INT(nx)*INT(ny)*INT(nz);
//    N=new (std::nothrow) Node[n];
    try {
      N=new Node[n];
    } catch (std::bad_alloc&) {
      std::cerr << "\nNot enough memory to allocate main Box. " << sizeof(Node)*n/1024/1024 <<" Mb needed for "
                   <<n<<" total nodes! Aborting...\n";
      ierr=-1;
      return;
    }
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


// now returns the number of unresolved triangles
int Box::loadSTL(std::ifstream & f, std::vector<Triangle> *triangles, double scale, std::ostream *out,
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
        t.v1.from_string(line,1);
        if (f.eof()) {std::cerr << "STL file truncated, missing vertex, aborting...\n"; return (-1);}
        
        std::getline(f,line);
        if (line.find("vertex") ==std::string::npos || f.eof()) {
          std::cerr <<"STL file truncated, missing vertex, aborting...\n";
          return (-1);
        }
        t.v2.from_string(line,1);
        if (f.eof()) {std::cerr << "STL file truncated, missing vertex, aborting...\n"; return (-1);}

        std::getline(f,line);
        if (line.find("vertex") ==std::string::npos) {
          std::cerr <<"STL file truncated, missing vertex, aborting...\n";
          return (-1);
        }
        t.v3.from_string(line,1);
        
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
    if (triangles->size()<4 ) {
        std::cerr << "Missing triangles in STL file, aborting...\n";
        return (-1);
    }
    // leave some free space on the origin side (walls are still unset) 

    Point::px ? O_x=floor(O_x) : O_x=floor(O_x-4.);
    Point::py ? O_y=floor(O_y) : O_y=floor(O_y-4.);
    Point::pz ? O_z=floor(O_z) : O_z=floor(O_z-4.);
    nx = (int) (ceil(xmax-O_x)); ny = (int) (ceil(ymax-O_y)); nz = (int) (ceil(zmax-O_z));
    // leave some free space on the other side (walls are still unset) 
    if (!Point::px) nx+=4;
    if (!Point::py) ny+=4;
    if (!Point::pz) nz+=4;
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
        ox > O_x ||
        ox + _nx < O_x + nx-1 ||
        oy > O_y ||
        oy + _ny < O_x + ny-1 ||
        oz > O_z ||
        oz + _nz < O_z + nz-1
      ) return (-2);
      O_x=ox;
      O_y=oy;
      O_z=oz;
      nx=_nx+1;
      ny=_ny+1;
      nz=_nz+1;
        
    }
    //traslate triangles so that they lie within the positive octant 
    for (std::vector<Triangle>::iterator it=triangles->begin(); it != triangles->end(); ++it) {
      it->v1.x -= O_x; it->v2.x -= O_x; it->v3.x -= O_x;
      it->v1.y -= O_y; it->v2.y -= O_y; it->v3.y -= O_y;
      it->v1.z -= O_z; it->v2.z -= O_z; it->v3.z -= O_z;
    }
    
    return (hires);
}
#if 0
// old version
void Box::FillTriangles(std::vector<Triangle> *triangles, double reinforce_wall,bool oriented) {
  int x0;
  int y0;  
  int z0;
  int x1;
  int y1;  
  int z1;
  int c=0;
  int x,y,z;
//  Point neigh(nx,ny,nz);
  
  for (std::vector<Triangle>::iterator it=triangles->begin(); it != triangles->end(); ++it) {  
    x0=floor(MIN(it->v1.x, MIN(it->v2.x, it->v3.x)))-1;
    y0=floor(MIN(it->v1.y, MIN(it->v2.y, it->v3.y)))-1;
    z0=floor(MIN(it->v1.z, MIN(it->v2.z, it->v3.z)))-1;

    x1=ceil(MAX(it->v1.x, MAX(it->v2.x, it->v3.x)))+1;
    y1=ceil(MAX(it->v1.y, MAX(it->v2.y, it->v3.y)))+1;
    z1=ceil(MAX(it->v1.z, MAX(it->v2.z, it->v3.z)))+1;
    
    
    for (z=z0; z<=z1; z++) { 
      for (y=y0; y<=y1; y++) {
        for (x=x0; x<=x1; x++) {
          p.set(x,y,z);
          if (p.isOutOfBounds()) continue;
          // if the distance between the triangle plane and the node is less then sqrt(3)/2
          // and its projection lies on the triangle's surface, then it is a wall node!
          if ( this->node(p).isDead() && Box::isOnTriangle(p,*it,sqrt(3)*0.5*reinforce_wall,oriented)) {
            this->node(p).setType(WALL_NODE);
            c++;
//            if (tipo==OUTLET_NODE || tipo == INLET_NODE) {
            //store the veldir as equal to the triangle normal
            // this is useful if this box is an I/O box, then this veldir
            // will be copied to the main box node intersecting this triangle surface
            // It will be oriented later (in IntersectIOnodes) as towards the fluid nodes.
            //
            this->node(p).setVelDir(it->n.to_string());
            //}
          }
        }
      }
    }
        
  }
  std::cout << c << " WALL nodes generated from STL.\n";
  
}
#endif
// new version
void Box::FillTriangles(std::vector<Triangle> *triangles, bool voxel_based) {
  int i,j,ii;
  int c=0;
  Point pv(nx-1,ny-1,nz-1);  
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
        p.set(floor(P.x),floor(P.y),floor(P.z));
        if (p.isOutOfBounds()) continue;
        if ( this->node(p).isDead() ) {
          //put a wall node in the voxel origin
          this->node(p).setType(WALL_NODE);
          c++;
          //this node gets the triangle normal
          this->node(p).setVelDir(it->n.to_string());
        }
        if (voxel_based) {
          //put a wall on any other dead voxel vertex
          for (ii=1;ii<8;ii++) {
            pv=p.getp(ii,1,true);
            if ( this->node(pv).isDead() ) {
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


int Box::count(NODETYPE t) const {
  int c;
  c=0;
  for (int z=1; z<nz; z++) { 
    for (int y=1; y<ny; y++) { 
      for (int x=1; x<nx; x++) { 
        if (N[x+y*nx+z*nx*ny].getType()==t) c++;
      }
    }
  }
  return c;
}
//

// find the first "missed" node, i.e. the first encountered fluid node
// whose dist was not set (i.e. that were missed by FindMinDist.) AND
// is adjacent to either a boundary or
// a node with dist set.
// Return its position. Return "zero" position if it isn't found.
Point Box::FindMissed() {
  int i;
  Point newp(nx-1,ny-1,nz-1),p0(0,0,0,nx-1,ny-1,nz-1);
  for (int z=1; z<nz; z++) { 
    for (int y=1; y<ny; y++) { 
      for (int x=1; x<nx; x++) { 
        if (this->node(x,y,z).isFluid() && this->node(x,y,z).getDist()<0) {
          p0.set(x,y,z);
          for (i=1; i<NPPEXT; i++) {
            newp=p0.getp(i);
            if (newp.isOutOfBounds()) continue;
            if (this->node(newp).isBoundary() || (this->node(newp).isFluid() && this->node(newp).getDist()>=0)  ) //got it!
              return newp; //p0;
          }
        }
      }
    }
  }
  return p0;
}

//find all outer nodes starting from point p0 that
// MUST correspond to an outer node.
// Outer nodes will be marked as  "Outer"
// This MUST be called before FindExternal
#ifndef NORECUR
void Box::FindExternal(Point const & p0, bool watchout_holes) {
  Point newp(nx-1,ny-1,nz-1);
  this->node(p0).setOuter();
    /*
     * if there is an adjacent boundary node and p0 is "inside" it
     * i.e. if the scalar product of its position with respect to p0 and
     * its normal is positive then... don't enter!! return!!
     * Assumes outwards normals!!
     * >>>>>>>>>>>> TO BE IMPLEMENTED IN NON-RECURSIVE VER. <<<<<<<<<<<<
    */  
  if (watchout_holes) {
    for (int i=1; i<NPPLESS;i++) {
      newp=p0.getp(i);
      if (newp.isOutOfBounds()) continue;
      if ( this->node(newp).isBoundary()) {
        double3 v;
        v.from_string(this->node(newp).getVelDir());
        double3 r(newp.x()-p0.x(),newp.y()-p0.y(),newp.z()-p0.z());
        if (v*r>0) return;
      }
    }
  }
  for (int i=1; i<NPPLESS;i++) {
//  for (int i=1; i<NPPEXT;i++) {
    newp=p0.getp(i);
    if (newp.isOutOfBounds()) continue;
    if ( this->node(newp).isDead() && !this->node(newp).isOuter() ) this->FindExternal(newp,watchout_holes);
  }
}
// starting from p0, a boundary node adjacent to a fluid one,
// mark with dist=0 all boundary nodes adjacent to at least a fluid one, i.e.
// the innermost boundary layer.
// This must be called *before* FindMinDist
void Box::FindAdjacent(Point const & p0) {
  Point newp(nx-1,ny-1,nz-1);
  int i,j;
  this->node(p0).setDist(0);
  for (i=1; i<NPPEXT;i++) {
    newp=p0.getp(i);
    if (newp.isOutOfBounds()) continue;
    if (this->node(newp).isBoundary() && this->node(newp).getDist()<0 ) {
      // check if newp has at least a neigh fluid
      for (j=1; j<NPPEXT;j++) {
        if (newp.getp(j).isOutOfBounds()) continue;
        if (this->node(newp.getp(j)).isFluid()) {
          this->FindAdjacent(newp);
          break;
        }
      }
    }
  }
}

/*
 * starting from a node p0 having distance d and belonging to a layer of nodes all at dist. d,
 * set to distance d+1 all fluid neighbors in the yet unexplored adjacent layer
 */

void Box::FindMinDist(Point const & p0,  Point & last, int & totnodes) {

  Point newp(nx-1,ny-1,nz-1);
  int d;
  d=this->node(p0).getDist();
  this->node(p0).setExplored();
  for (int i=1; i<NPPEXT;i++) {
    newp=p0.getp(i);
    if (newp.isOutOfBounds() || this->node(newp).isExplored() ) continue;
//    if (newp.isOutOfBounds() || this->node(newp).isExplored() ) continue;
    if (this->node(newp).isFluid() && this->node(newp).getDist()<0 ) {
      this->node(newp).setDist(d+1); //just set distance BUT NOT SPACING! 
      last=newp;
/*
      this->node(p0).setNextLayer(&this->node(newp));
      this->node(newp).setPrevLayer(&this->node(p0));
*/
      totnodes++;
    } else if (this->node(newp).getDist() == d) FindMinDist(newp,last,totnodes);
  }
  
}


#else
//>>>>>>>>>>>>>NON RECURSIVE VERSIONS
void Box::FindExternal(Point const & p00, bool watchout_holes) {
  std::stack<const Point*> p0_s;
  std::stack<Point*> newp_s;
  std::stack<int*> i_s;
  int *i;
  Point * newp;
  Point const * p0;

  int j;
  Point newp2(nx-1,ny-1,nz-1);
  double3 v,r;
  
  p0=&p00;
  
start: 
  this->node(*p0).setOuter();

  
  if (watchout_holes) {
    /*
     * if there is an adjacent boundary node and p0 is "inside" it
     * i.e. if the scalar product of its position with respect to p0 and
     * its normal is positive then... don't enter!! return!!
     * Assumes outwards normals!!
    */  
    
    for (j=1; j<NPPLESS;j++) {
      newp2=p0->getp(j);
      if (newp2.isOutOfBounds()) continue;
      if ( this->node(newp2).isBoundary()) {
        v.from_string(this->node(newp2).getVelDir());
        r = double3(newp2.x()-p0->x(),newp2.y()-p0->y(),newp2.z()-p0->z());
        if (v*r>0) {
          if (!i_s.empty()) goto goprev; //return to the previous stack level
          return;  //there is no previous level: exit
        }
      }
    }
  }
  
  i = new int;
  newp=new Point(nx-1,ny-1,nz-1);
      
  for (*i=1; *i<NPPLESS;(*i)++) {
//  for (int i=1; i<NPPEXT;i++) {
    (*newp)=p0->getp(*i);
    if (newp->isOutOfBounds()) continue;
    if ( this->node(*newp).isDead() && !this->node(*newp).isOuter() ) {
      //stack current recursion  "level" and pass to next one
      p0_s.push(p0);
      newp_s.push(newp);
      i_s.push(i);
      p0=newp;
      goto start;
ritorno:;
    }
  }
  delete i;
  delete newp;
goprev:;
  if (!i_s.empty()) {
    i=i_s.top();i_s.pop();
    newp=newp_s.top(); newp_s.pop();
    p0=p0_s.top(); p0_s.pop();
    goto ritorno;
  }
}

void Box::FindAdjacent(Point const & p00) {
  std::stack<const Point*> p0_s;
  std::stack<Point*> newp_s;
  std::stack<int*> i_s;
  int *i;
  Point * newp;
  Point const * p0;
  int j;

  p0=&p00;

start:
  i = new int;
  newp=new Point(nx-1,ny-1,nz-1);

  this->node(*p0).setDist(0);
  for (*i=1; *i<NPPEXT;(*i)++) {
    *newp=p0->getp(*i);
    if (newp->isOutOfBounds()) continue;
    if (this->node(*newp).isBoundary() && this->node(*newp).getDist()<0 ) {
      // check if newp has at least a neigh fluid
      for (j=1; j<NPPEXT;j++) {
        if (newp->getp(j).isOutOfBounds()) continue;
        if (this->node(newp->getp(j)).isFluid()) {
          p0_s.push(p0);
          newp_s.push(newp);
          i_s.push(i);
          p0=newp;
          goto start;
        }
      }
      
    }
ritorno:  
    ;
  }
  delete i;
  delete newp;
  if (!i_s.empty()) {
    i=i_s.top();i_s.pop();
    newp=newp_s.top(); newp_s.pop();
    p0=p0_s.top(); p0_s.pop();
    goto ritorno;
  }
  
}

//starting from a node at dist=0, find all nodes' distances
void Box::FindMinDist(Point const & p00, Point & last, int & totnodes) {

  int d;
  std::stack<const Point*> p0_s;
  std::stack<Point*> newp_s;
  std::stack<int*> i_s;
  int *i;
  Point * newp;
  Point const * p0;

  p0=&p00;

start:
  i = new int;
  newp=new Point(nx-1,ny-1,nz-1);
  
  d=this->node(*p0).getDist();
  this->node(*p0).setExplored();
  for (*i=1; *i<NPPEXT;(*i)++) {
    *newp=p0->getp(*i);
    if (newp->isOutOfBounds() || this->node(*newp).isExplored() ) continue;
    if (this->node(*newp).isFluid() && this->node(*newp).getDist()<0 ) {
      this->node(*newp).setDist(d+1);
      last=*newp;
/*
      this->node(*p0).setNextLayer(&this->node(*newp));
      this->node(*newp).setPrevLayer(&this->node(*p0));
*/
      totnodes++;
    } else if (this->node(*newp).getDist() == d) {
//      FindMinDist(newp,dfactor,last,totnodes);
      p0_s.push(p0);
      newp_s.push(newp);
      i_s.push(i);
      p0=newp;
      goto start;
ritorno:
      ;
    }
  }
  delete i;
  delete newp;
  if (!i_s.empty()) {
    i=i_s.top();i_s.pop();
    newp=newp_s.top(); newp_s.pop();
    p0=p0_s.top(); p0_s.pop();
    goto ritorno;
  }
  
}



//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#endif

// mark internal nodes as fluid
// WARNING: must be called AFTER FindExternal
// If internal == false then do the opposite, i.e. mark external nodes as fluid
// Return in <p0> the first fluid node found which is adjacent to a boundary
int Box::FindInternal(Point & p0, bool internal) {
  bool found;
  found=false;
  Point newp(nx-1,ny-1,nz-1);
  int c=0;
  int i;
  for (int z=1; z<nz; z++) { 
    for (int y=1; y<ny; y++) { 
      for (int x=1; x<nx; x++) { 
        if (this->node(x,y,z).isDead() && (internal ^ this->node(x,y,z).isOuter()) ) {
          if (!found) {
            p0.set(x,y,z);
            for (i=1; i<NPPEXT; i++) {
              newp=p0.getp(i);
              if (newp.isOutOfBounds()) continue;
              if (this->node(newp).isBoundary()) {
                p0=newp; //.set(x,y,z);//got it!
                found=true;
              }
            }
          }
          this->node(x,y,z).SetFluidNode();
          c++;
        }
      }
    }
  }
  if (!found) p0.set(-1,-1,-1);
  return (c);
}

/*
 * set all nodes spacing.
 * to be called AFTER FindMinDist!
 */
void Box::SetSpacing(int dfactor) {
  for (int z=1; z<nz; z++) { 
    for (int y=1; y<ny; y++) { 
      for (int x=1; x<nx; x++) { 
        if (this->node(x,y,z).isFluid()) this->node(x,y,z).setSpacing(dfactor);
      }
    }
  }
}

/*
 * set as <iotype> nodes *all* those nodes (at this stage this box fluid nodes
 * aren't found yet, there is no inside/outside distinction)
 * of <this> box that are located at the same position of any ACTIVE nodes of the given ibox.
 * Moreover, if iotype=INLET_NODE, the veldir of ibox WALL nodes (because only these have it)
 * is copied also.
 * The searchdir on new iotype nodes  will be found and set later on when the fluid
 * nodes of this box will be available.
 * Return the number of iotype nodes set.
 */ 
int Box::IntersectIOnodes(Box & ibox,NODETYPE iotype){
  int inx,iny,inz,c=0;
  ibox.getBounds(inx,iny,inz);
  
  //scan nodes within the ibox
  for (int z=1; z<=inz; z++) { 
    for (int y=1; y<=iny; y++) { 
      for (int x=1; x<=inx; x++) {
        //set p as the point of this box that is at the same "original" location
        //as the ibox point
        p.set((int)((double)(x)+ibox.getO_x()-O_x),
              (int)((double)(y)+ibox.getO_y()-O_y),
              (int)((double)(z)+ibox.getO_z()-O_z));
        if (p.isOutOfBounds()) continue;
        if (!ibox.node(x,y,z).isDead() ) {
          this->node(p).setType(iotype);
          c++;
          if (iotype==INLET_NODE && ibox.node(x,y,z).getType()==WALL_NODE)
             this->node(p).setVelDir(ibox.node(x,y,z).getVelDir());
        } 
      }
    }
  }
         
  return (c);
}


/*
 * set the maxspacing of fluid nodes of <this> box that are located at the same position
 * of any ACTIVE nodes of the given ibox.
 * Return the number of involved nodes.
 */ 
int Box::ForceHighResolution(Box & ibox){
  int inx,iny,inz,c=0;
  ibox.getBounds(inx,iny,inz);
  
  //scan nodes within the ibox
  for (int z=1; z<=inz; z++) { 
    for (int y=1; y<=iny; y++) { 
      for (int x=1; x<=inx; x++) {
        //set p as the point of this box that is at the same "original" location
        //as the ibox point
        p.set((int)((double)(x)+ibox.getO_x()-O_x),
              (int)((double)(y)+ibox.getO_y()-O_y),
              (int)((double)(z)+ibox.getO_z()-O_z));
//        p.set(int(double(x)+ibox.getO_x()-O_x),
//              int(double(y)+ibox.getO_y()-O_y),
//              int(double(z)+ibox.getO_z()-O_z));
        if (p.isOutOfBounds()) continue;
        if (!ibox.node(x,y,z).isDead() ) {
          this->node(p).SetMaxSpacing(1);
          c++;
        } 
      }
    }
  }
         
  return (c);
}


/*
 * find and set all i/o nodes searchdir and check
 * and correct inlet nodes veldir verse.
 * Return false if inlet OR outlet are not in touch with fluid nodes!
 *
 */ 
bool Box::SetIOnodesSearchDir(){
  bool iexp=false;
  bool oexp=false;
  //double3 searchdir,veldir;
  
  Point neigh(nx-1,ny-1,nz-1);
  //scan nodes within the ibox
  for (int z=1; z<nz; z++) { 
    for (int y=1; y<ny; y++) { 
      for (int x=1; x<nx; x++) {
        p.set(x,y,z);
        if (this->node(p).getType()==INLET_NODE || this->node(p).getType()==OUTLET_NODE) { 
          for (int i=1; i<NPPLESS;i++) {
            neigh=p.getp(i);
            if (neigh.isOutOfBounds()) continue;
            if (this->node(neigh).isFluid()) {
              iexp=iexp || this->node(p).getType()==INLET_NODE;
              oexp=oexp || this->node(p).getType()==OUTLET_NODE;
              std::stringstream iss;
              iss <<p.x()-neigh.x() <<" "<< p.y()-neigh.y()<<" "<<  p.z()-neigh.z();
              this->node(p).setSearchDir(iss.str());
              //searchdir.from_string(iss.str());
              //this->node(p).setSearchDir(searchdir.to_string());
              /*
              //check and correct the verse of veldir
              if (this->node(p).getType()==INLET_NODE) { 
                veldir.from_string(this->node(p).getVelDir());
                //this is to adjust veldir that is assumed to go from inlet nodes towards fluid nodes,
                //thus its scalar product with searchdir must be negative.
                if (veldir*searchdir>0.) {
                  this->node(p).setVelDir((veldir*(-1.)).to_string()); //veldir e' verso fuori! Rigirala!
                }
//                std::cout << searchdir << " | " << veldir <<std::endl; 
              }
              */
              break;
            }
          }
        }
      }
    }
  }
  return (iexp && oexp);
}



double Box::dvol(Triangle  const &t) {

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

double Box::dsur(Triangle const &t) {
        return ((t.v2 - t.v1).Curl(t.v3 - t.v1).Modulus()*0.5);
}

/*
 * return true if the projection of <this> point lies within the triangle ABC and
 * the point is at a distance <toll from the triangle plane
 * 
 */
bool Box::isOnTriangle(Point const & p0, Triangle const &t,double toll, bool oriented) {
  double3 P(double(p0.x()),double(p0.y()),double(p0.z())),AP;
  AP=P-t.v1;
  if (oriented) {
    // put outer WALL only (the triangle normal MUST be OUTWARDS)
    if (AP*t.n <0 ||AP*t.n>toll) return false; //this point (P) is too far from the plane OR is outward
    //if (AP*t.n >0 ||AP*t.n < -toll) return false; //this point (P) is too far from the plane OR is inward
  } else {
    // put both outer and inner WALL (the triangle normal verse is irrelevant)
    if (fabs(AP*t.n) > toll) return false; //this point (P) is too far from the plane
  }
  double3 AB,AC,PA,PB,PC;
  
  P=P-(t.n*(AP*t.n)); //P is now this point's projection on the triangle plane
  AB=t.v2-t.v1; AC=t.v3-t.v1; PA=t.v1-P; PB=t.v2-P; PC=t.v3-P;

  double ar= AB.Curl(AC).Modulus();  //double area of ABC triangle
  double a= PB.Curl(PC).Modulus(); //double area of BPC triangle
  //if (a > ar) return false;
  double b= PC.Curl(PA).Modulus(); //double area of CPA triangle
  //if (b > ar ) return false;
  double g= PB.Curl(PA).Modulus(); //double area of BPA triangle
  //if (fabs(g +a +b-ar)>ar) return false; //1e-6 ) return false;
  if (fabs(g +a +b-ar)>1) return false; //1e-6 ) return false;

  return true;  
  
  
}

/*
bool Box::isOnTriangle(Point const & p0, Triangle const &t,double toll, bool oriented) {
  double3 P(double(p0.x()),double(p0.y()),double(p0.z())),AP;
  double3 BAR=(t.v1+t.v2+t.v3)/3.;
  double3 A=t.v1;
  double3 B=t.v2;
  double3 C=t.v3;
  
  double3 AB=B-A;
  double3 AC=C-A;
  double3 BC=C-B;
  
  double s=MIN(AB.Modulus(),MIN(AC.Modulus(),BC.Modulus()));
  s=(s+2.*sqrt(3.))/s;
  A=BAR+((A-BAR)*s);B=BAR+((B-BAR)*s);C=BAR+((C-BAR)*s);
  
  AP=P-A;
  if (oriented) {
    // put outer WALL only (the triangle normal MUST be OUTWARDS)
    if (AP*t.n <0 ||AP*t.n>toll) return false; //this point (P) is too far from the plane OR is outward
    //if (AP*t.n >0 ||AP*t.n < -toll) return false; //this point (P) is too far from the plane OR is inward
  } else {
    // put both outer and inner WALL (the triangle normal verse is irrelevant)
    if (fabs(AP*t.n) > toll) return false; //this point (P) is too far from the plane
  }
  double3 PA,PB,PC;
  
  P=P-t.n*(AP*t.n); //P is now this point's projection on the triangle plane
  AB=B-A; AC=C-A; PA=A-P; PB=B-P; PC=C-P;

  double ar= AB.Curl(AC).Modulus();
  double a= PB.Curl(PC).Modulus(); //double area of BPC triangle
  if (a > ar) return false;
  double b= PC.Curl(PA).Modulus(); //double area of CPA triangle
  if (b > ar ) return false;
  double g= PB.Curl(PA).Modulus(); //double area of BPA triangle
  if (fabs(g +a +b- ar)>5) return false; //1e-6 ) return false;

  return true;  
  
  
}
*/
#if 0
// this is from generate_mesh.c
bool Box::isOnTriangle(Point const & p0, Triangle const &t,double toll) {

  double3 P(double(p0.x()),double(p0.y()),double(p0.z())), U, V, nn;
  double lamda1, lamda2;
  int proj;

  V = t.v3 - t.v1;

  U = t.v2 - t.v1;

  P = P - t.v1;
  
  nn.x = fabs(t.n.x);
  nn.y = fabs(t.n.y);
  nn.z = fabs(t.n.z);
  proj = (nn.x > nn.y) ? ((nn.x > nn.z) ? 0 : 2) : ((nn.y >nn.z) ? 1 : 2);
  
  switch (proj) {
    case 0:
      V.x = V.z;
      U.x = U.z;
      P.x = P.z;
    break;
    case 1:
      V.y = V.z;
      U.y = U.z;
      P.y = P.z;
    break;
    case 2:;
  }
  lamda1 = (P.x*V.y - P.y*V.x) / (U.x*V.y - U.y*V.x);
  lamda2 = (P.x*U.y - P.y*U.x) / (V.x*U.y - V.y*U.x);

  return (lamda1 >= 0.) && (lamda2 >= 0.) && ((lamda1 + lamda2) <= 1.);
}
#endif
/*
 * find an outer node on the bounding box faces and return it in p0.
 * Return false if it couldn't.
 * WARNING: the given p0 must be one of the bounding box vertex
 * 
 */
bool Box::FindOuterNode(Point & p0) {
  Point newp(nx-1,ny-1,nz-1);
  newp=p0;
  for (int i=1; i<NPP; i++) {
    while (true){
      if (this->node(newp).isDead()) {p0=newp; return true;}//got it!
      newp=newp.getp(i);
      if (newp.isOutOfBounds()) {newp=p0; break;} //try another direction
    }
  }
  return false;
}

bool Box::NodeIsAdjacentToType(Point const & p0, NODETYPE tipo,int npp) {
  Point newp(nx-1,ny-1,nz-1);
  for (int i=1; i<npp; i++) {
    newp=p0.getp(i);
    if (newp.isOutOfBounds()) continue;
    if (this->node(newp).getType()==tipo) return true;//got it!
  }
  return false;
}


/*
 * for diagnosting and debugging purposes
 */

void Box::Dump(std::string const &filename, std::string const & id) {
  std::ofstream f;
  f.open(filename.c_str(),std::ofstream::out);
  for (int z=1; z<nz; z++) { 
    for (int y=1; y<ny; y++) { 
      for (int x=1; x<nx; x++) { 
        if ( this->node(x,y,z).getDist()>=0) { //????????
//        if (this->node(x,y,z).isFluid() && !this->node(x,y,z).isExplored()) { //????????
          f << id<<this->node(x,y,z).getDist()<<" " << x <<" " << y <<" " << z <<" " << 1 << std::endl; //??????????'
        }
      }
    }
  }
  f.close();
}
