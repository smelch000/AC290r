/*************************************************
 * 
 * genmesh: automatic multi mesh generator for MUPHY2
 *         by P. Miocchi and S. Melchionna
 *
**************************************************/ 
// g++ genmesh2.cpp readargs.cpp point.cpp node.cpp mesh.cpp
#include <iostream>
#include <string>
#include "genmesh.h"
#include "point.h"
#include "node.h"
#include "mesh.h"
#include "readargs.h"

/*
template <class T>
std::string vector2string(const T d[]) {
  std::ostringstream iss;

  iss << d[0] << " " << d[1] << " " << d[2];
  return iss.str();
};
*/

int main(int argc, char *argv[]){
  int numargs,maxrank;
  Mesh** mesh;
  int dfactor;
  int nx,ny,nz,ierr;
  int n_ibox=0,n_obox=0, n_fbox=0;
  double scale;
  int MaxRank;
  bool px,py,pz;
  std::string outfile,infile;
  int n,i;
  double velin,pressout;
  bool writeios;
  std::string ifile[MAXINLETBOXES];
  std::string ofile[MAXOUTLETBOXES];
  std::string ffile[MAXFORCINGBOXES];
  int frank[MAXFORCINGBOXES];
  double boundingbox[6];
  Mesh * iobox;
  Mesh * fbox;
  double Ox,Oy,Oz;
  
  
  std::string descr;
  std::string options[]={"i","o","f","r","x","y","z","v","p","I","O","s","H","W","X","N","F","E","b","d","t","h"};
  std::string descriptions[]={"Input File=","Output file (spacing and extension will be appended)=",
    "Minimum mesh layers=",
    "Mesh maximum rank (0 : spacing 1, 1 : spacing 2, 2 : spacing 4,...)=",
    "Periodic on x-axis","Periodic on y-axis","Periodic on z-axis",
    "Write ios file with inlet vel=",
    "Write ios file with outlet pressure=",
    "No. of inlet STL boxes read=",
    "No. of outlet STL boxes read=",
    "STL input data scale factor=",
    "STL \"holes avoiding\" mode enabled",
    "WALL \"repairing\" is OFF",
    "XYZ file output is ON",
    "Single WALL nodes instead of voxels on STL surface",
    "No. of custom resolution boxes read=",
    "External FLUID nodes",
    "Bounding box (rescaled) set by the user (Ox Oy Oz nx ny nz)=",
    "Debug probe point (x y z)=",
    "Automatic translation enabled",
    ""};
  ReadArgs inputPars(options,descriptions,21,numargs,argc, argv);
  if (numargs < 2 || numargs >21 || inputPars.isset("h")) {
    std::cout << "GENMESH - ver. "<< VERSION <<" \n\n";
    std::cout << "Usage: genmesh -i FILE.[dat|stl] -o FILENAME [options...]\n";
    std::cout << "Options:\n";
    std::cout << "\t-i <file.[dat|stl]>\n\t\tRead boundaries from <file> with a DAT (e.g. bgkflag.dat) or an STL extension.\n"
    << "\t\tIn the STL case a rescaling can be done and I/O nodes can be read separately\n\t\tfrom STL files (see below)\n\n"
    << "\t-o <outfile>\n\t\tWrite mesh(es) on file(s) named <outfile> to which spacing and extension DAT will be appended\n\n"
    << "\t[-E] Set fluid nodes externally to the boundary instead of internally.\n\n"
    << "\t[-r R -f F [-F <forcingfile1 rank1 forcingfile2 rank2...>]]\n"
    << "\t\tEnable multimesh generation. The maximum permitted mesh spacing will be =2**<R> and each mesh\n"
    << "\t\twill extend at least <F> layers (>2), with the resolution decreasing from the boundary inwards.\n"
    << "\t\tOptionally force custom resolutions (<rank1...n>) in the system regions enclosed by the boxes\n"
    << "\t\tread from the <forcingfile1...n> files (they must have the same extension as input <file>).\n\n"
    << "\t[-x] [-y] [-z] set the system as periodic along coord.s axes (default= no periodicity).\n\n"
    << "\t[-W] DON'T convert DEAD to WALL nodes around fluid nodes.\n\n"    
    << "\t[-t] Node positions in output files are automatically translated to the positive octant.\n\n"
    << "\t[-X] Save the overall mesh file <outfile>.xyz in XYZ format for visualization purposes.\n\n"    
    << "\t[-d <rank Px Py Pz> ] write in debug.xyz all first neighbors of the given (rescaled) box point.\n"    
    << "\t\tUseful for debugging purposes. Mesh with the given rank and its connected meshes are considered.\n\n"    
    << "\tThe following options are effective only with an STL input file\n\n"
    << "\t[-s S] set rescaling factor to <S>.\n\n"
    << "\t[-b <[Ox Oy Oz] nx ny nz> ] set the origin (in the STL case) and the size of the Bounding Box.\n"    
    << "\t\t In the STL case they will be multiplied by <S>. In the DAT case the origin is always assumed at (0,0,0).\n\n"    
    << "\t[-N] put single WALL nodes instead of voxels on STL surface. Walls will be of half thickness.\n\n"
    << "\t[-H] enable the \"holes avoiding\" mode to automatically detect small holes in STL surface.\n\n"
    << "\t[-v V -p P [-I <inletfile1 inletfile2 ...>] [-O <outletfile1 outletfile2 ...>]]\n"
    << "\t\tset inlet nodes vel. to <V>, outlet nodes pressure to <P> and write I/O nodes to ios file\n"
    << "\t\tI/O nodes are read from the given list of inlet and outlet STL files (no wildcards!),\n"
    << "\t\tafter rescaling them by the factor <S>.\n\n"
    <<  RED <<"\t\tWARNING: in multimesh mode, I/O nodes are supposed to belong to the finest mesh!"<<NOCOL<<"\n\n";
    return (0);
  }
  px=inputPars.isset("x"); py=inputPars.isset("y"); pz=inputPars.isset("z");
#if 0
  if (px && py && pz) {
    std::cerr << "Error! The system cannot be periodic in any direction" << std::endl; return (1);
  }
#endif
  std::cout << "GENMESH - ver. "<< VERSION <<" \n\n";

  if (inputPars.GetParameter<std::string>("i",descr,infile)) {
    std::cout << descr << infile << std::endl;
  } else {
    std::cerr << "Error! Please specify input file name (option -i).\n"; return (1);
  }
  if (inputPars.GetParameter<std::string>("o",descr,outfile)) {
    std::cout << descr << outfile << std::endl;
  } else {
    std::cerr << "Error! Please specify output base file name (option -o).\n"; return (1);
  }

  if (inputPars.isset("r") != inputPars.isset("f")) {
    std::cerr << "Error! Options -f and -r must be set together.\n"; return (1);
  }

  if (infile.find(".stl")!=std::string::npos && inputPars.isset("v") != inputPars.isset("p")) {
    std::cerr << "Error! Options -v and -p must be set together.\n"; return (1);
  }
  writeios = infile.find(".stl")!=std::string::npos && inputPars.isset("v");

  dfactor=__INT_MAX__;
  if (inputPars.GetParameter<int>("f",descr,dfactor)) {
    if (dfactor <3) {std::cerr << "Error! Mesh layers must be more than 2" << std::endl; return (1);}
    std::cout << descr << dfactor << std::endl;
  }

  MaxRank=0;
  if (inputPars.GetParameter<int>("r",descr,MaxRank)) std::cout << descr << MaxRank << std::endl;
  
  n_fbox=0;
  if (MaxRank>0 && inputPars.isset("F")) {
    std::string temp[2*MAXFORCINGBOXES];
    n_fbox=inputPars.GetParameter<std::string>("F",descr,temp,MAXFORCINGBOXES);
    if (n_fbox % 2 != 0) {
      std::cerr << "Error! Check custom resolution boxes file name and rank pairs. Aborting...\n";
      return (1);        
    }  
    for (i=0;i<n_fbox;i+=2) {
      ffile[i/2]=temp[i];
      frank[i/2]=atoi( temp[i+1].c_str());
      if (frank[i/2] >MaxRank) {
        std::cerr << "Error! custom resolution rank cannot be larger than maximum rank. Aborting...\n";
        return (1);        
      }
    }
    n_fbox /=2;
    
    for (i=0;i<n_fbox;i++) {
      if (
        (infile.find(".stl")!=std::string::npos && ffile[i].find(".stl")==std::string::npos)
        || (infile.find(".dat")!=std::string::npos && ffile[i].find(".dat")==std::string::npos)
      ) {
        std::cerr << "Error! Forcing boxes files must have the same format as input file. Aborting...\n";
        return (1);        
      }
    }
    if (n_fbox>0) {
      std::cout << descr << n_fbox << std::endl;
    }
  }


  if (writeios) {
    if ( !inputPars.isset("I") && !inputPars.isset("O") ) {
      std::cerr << "Error! At least one of the options -I and -O must be set!\n"; return (1);
    }
    if (inputPars.GetParameter<double>("v",descr,velin)) std::cout << descr << velin << std::endl;
    if (inputPars.GetParameter<double>("p",descr,pressout)) std::cout << descr << pressout << std::endl;
//    inputPars.GetParameter<int>("I",descr,n_ibox);
//    std::cout << descr << n_ibox << std::endl;
    n_ibox=inputPars.GetParameter<std::string>("I",descr,ifile,MAXINLETBOXES);
    if (n_ibox>0) {
      std::cout << descr << n_ibox << std::endl;
    }
    n_obox=inputPars.GetParameter<std::string>("O",descr,ofile,MAXOUTLETBOXES);
    if (n_obox>0) {
//    inputPars.GetParameter<int>("O",descr,n_obox);
      std::cout << descr << n_obox << std::endl;
    }
//    for (i=0; i<n_ibox; i++) ifile[i]="inlet_"+number2string(i)+"_IN.stl";
//    for (i=0; i<n_obox; i++) ofile[i]="outlet_"+number2string(i)+"_OUT.stl";
  } else if ((inputPars.isset("I") || inputPars.isset("O"))) {
    std::cerr << "\n"<<RED<<"WARNING: options -I and -O are ignored because I/O nodes generation was not enabled"<<NOCOL<<"\n\n";
  }
  scale=1.;
  if (inputPars.isset("t")) std::cout << "Automatic Box TRANSLATION enabled\n";
  if (inputPars.isset("E")) std::cout << "The fluid is EXTERNAL to the boundary surface\n";
  if (!inputPars.isset("W")) std::cout << "WALL repairing enabled\n";
  if (inputPars.isset("X")) std::cout << "XYZ file output enabled\n";
 
  Mesh::MAXSPACING=(1 << MaxRank);  //=2**MaxRank
  Point::setPeriodic(px,py,pz);

  mesh = new Mesh*[MaxRank+1];
  for (int i=0; i<= MaxRank;i++) {mesh[i]=NULL;}

  INPUTTYPE InputMode;
  

  if (infile.find(".dat")!=std::string::npos) {
    InputMode=DAT;
    std::cout << "\nloading boundary nodes..." <<std::endl;

    if (inputPars.isset("b")) {
      if (inputPars.GetParameter<double>("b",descr,boundingbox,6) != 3) {
        std::cerr <<  "\n"<<RED<<"Only bounding box x,y,z sizes must be given! Aborting..." <<NOCOL<<"\n\n";
        return (1);
      }
      std::cout << descr << 0 <<" " << 0 <<" " << 0 <<" " <<
      boundingbox[0] <<" " << boundingbox[1] <<" " << boundingbox[2] <<" " <<std::endl;
      mesh[0]= new Mesh(0,infile, InputMode, ierr, 1.,false,
                  (int)(boundingbox[0]),(int)(boundingbox[1]),(int)(boundingbox[2]),0.,0.,0.                  
                 );
    } else {
      mesh[0] = new Mesh(0,infile, InputMode, ierr, 1.);
    }
  } else if (infile.find(".stl")!=std::string::npos) {
    InputMode=STL;
    inputPars.GetParameter<double>("s",descr,scale);
    std::cout << descr << scale << std::endl;

    if (inputPars.isset("b")) {
      if (inputPars.GetParameter<double>("b",descr,boundingbox,6)<6) {
        std::cerr <<  "\n"<<RED<<"Bounding box origin coordinates and sizes must be given! Aborting..." <<NOCOL<<"\n\n";
        return (1);
      }
      std::cout << descr << (int)(boundingbox[0]*scale) <<" " << (int)(boundingbox[1]*scale) <<" " << (int)(boundingbox[2]*scale) <<" " <<
      (int)(boundingbox[3]*scale) <<" " << (int)(boundingbox[4]*scale) <<" " << (int)(boundingbox[5]*scale) <<" " <<std::endl;
    }
    
    
    if (inputPars.isset("H")) std::cout << "\"holes avoiding\" mode enabled\n";
    if (inputPars.isset("N")) std::cout << "single nodes instead of voxels generated on STL surface\n";
    std::cout << "\nloading STL boundary surface..." <<std::endl;
    if (inputPars.isset("b")) {
      mesh[0] = new Mesh(0,infile, InputMode, ierr, scale,!inputPars.isset("N"),
                  (int)(boundingbox[3]*scale),(int)(boundingbox[4]*scale),(int)(boundingbox[5]*scale),                  
                  boundingbox[0]*scale,boundingbox[1]*scale,boundingbox[2]*scale                  
                 );
    } else {
      mesh[0] =  new Mesh(0,infile, InputMode, ierr, scale,!inputPars.isset("N"));
    }
    
  } else {
    std::cerr << "Unrecognized extension for input file!\n";
    return(1);
  }
  
  if (ierr<0) return (1);
  
  if (px) std::cout <<   "Periodic on x-axis" << std::endl;
  if (py) std::cout <<   "Periodic on y-axis" << std::endl;
  if (pz) std::cout <<   "Periodic on z-axis" << std::endl;

  mesh[0]->getBounds(nx,ny,nz);
  std::cout << "Box bounds: "<<nx << ", "<<ny << ", "<<nz  << std::endl;

  if ((px || py || pz)&&(MaxRank>0)) {
    //check bounds divisibility by MaxRank
    unsigned int s=Mesh::MAXSPACING-1;
    if (((nx&s) != 0 && px)|| ((ny&s) != 0 && py)|| ((nz&s) != 0 && pz) ) {
        std::cerr << "Error! Box bound must be divisible by longest mesh spacing in the periodic direction. Aborting...\n\n";
        return (1);
    }
  }
  //std::cout << "\n>>> Box memory occupation: "<<sizeof(mesh[0])/1024/1024 << " Mbyte" << std::endl;

  if (inputPars.isset("t")) {
    //output files will be automatically translated to the positive octant
    //as required by moebius
    Ox=Oy=Oz=0.; 
  } else {
    //output files will be translated to the same reference frame of the original STL/DAT input files
    Ox=mesh[0]->getO_x();
    Oy=mesh[0]->getO_y();
    Oz=mesh[0]->getO_z();
  }

  maxrank=0; //leave this here because i/o boxes and/or forced box could be with rank>0
  //allocate and link all possible meshes
  for (i=1;i<=MaxRank;i++) {
    if (mesh[i]==NULL){
        mesh[i]=new Mesh(nx, ny, nz,i,mesh[0]->getO_x(),mesh[0]->getO_y(),mesh[0]->getO_z());
          mesh[i]->setMeshH(mesh[i-1]);
          mesh[i-1]->setMeshL(mesh[i]);
    }
  }
  
  if (writeios) {
    int i_b,d;
    int nn;
    for (i_b=0; i_b<n_ibox; i_b++) {
      
      std::cout << "\nreading "<<ifile[i_b]<<" to insert INLET nodes...\n";
      iobox = new Mesh(0,ifile[i_b], STL, ierr, scale, !inputPars.isset("N"));
      if (ierr<0) return(1);
      std::cout << "inlet "<<ifile[i_b]<<" file: filling 'coating' layer..." <<std::endl;
      d=0;
      if (iobox->AddFirstLayers(d, true)<0){
          std::cerr << "\nError! Couldn't find a coating layer. Aborting...\n\n";
          return (1);          
      }
      std::cout << "inlet "<<ifile[i_b]<<" file: setting internal nodes..." <<std::endl;
      i=0;
      nn=1;
      while (nn>0){nn=iobox->AddLayer(d,i,0);}
      if (nn<0) return(1);
      
      std::cout << "inlet "<<ifile[i_b]<<" file: inserting inlet nodes..." <<std::endl;
      mesh[0]->IntersectIOnodes(*iobox,INLET_NODE);
      delete iobox;
    }
    

    for (i_b=0; i_b<n_obox; i_b++) {
      
      std::cout << "\nreading "<<ofile[i_b]<<" to insert OUTLET nodes...\n";
      iobox = new Mesh(0,ofile[i_b], STL, ierr, scale, !inputPars.isset("N"));
      if (ierr<0) return(1);
      std::cout << "outlet "<<ofile[i_b]<<" file: filling 'coating' layers..." <<std::endl;
      d=0;
      if (iobox->AddFirstLayers(d, true)<0){
          std::cerr << "\nError! Couldn't find a coating layer. Aborting...\n\n";
          return (1);          
      }
      std::cout << "outlet "<<ofile[i_b]<<" file: setting internal nodes..." <<std::endl;
      i=0;
      nn=1;
      while (nn>0){nn=iobox->AddLayer(d,i,0);}
      if (nn<0) return(1);
      
      std::cout << "outlet "<<ofile[i_b]<<" file: inserting inlet nodes..." <<std::endl;
      mesh[0]->IntersectIOnodes(*iobox,OUTLET_NODE);
      delete iobox;
    }
    std::cout << "I/O nodes successfully inserted...\n\n";
    
    
  }
    
  int d,dstart;


  {
    int i_f,ff=0;
    int dist_offset=MAX(nx+1,MAX(ny+1,nz+1));
    int forced_rank;
    int nn;
    for (i_f=0; i_f<n_fbox; i_f++) {
      forced_rank=frank[i_f];
      std::cout << "\nreading "<<ffile[i_f]<<" to force " << " rank "<< forced_rank<< " resolution...\n";
      fbox = new Mesh(forced_rank,ffile[i_f], InputMode, ierr, scale,!inputPars.isset("N"));
      if (ierr<0) return(1);
      std::cout << "Custom resolution on "<<ffile[i_f]<<": searching for 'coating' layers..." <<std::endl;

      d=0;
      if (fbox->AddFirstLayers(d, true)<0){
          std::cerr << "\nError! Couldn't find a 'coating' layer. Aborting...\n\n";
          return (1);          
      }
      std::cout << "Custom resolution on "<<ffile[i_f]<<": setting internal nodes..." <<std::endl;
      i=forced_rank;
      do {
        nn=fbox->AddLayer(d,i,0);
        if (nn>0) maxrank = MAX(maxrank,i);
      } while (nn>0);
      if (nn<0) return(1);
      
      std::cout << ffile[i_f]<<" file: forcing resolution rank "<< forced_rank<<"..." <<std::endl;
      //this gives to custom resolution nodes a min_dist that is surely different from that that
      //will be given to normal nodes and to any other custom boxes.
      // This avoids that the custom boxes nodes will "interfere" in setting layers of normal nodes, 
      dist_offset += d+1;
      ff += mesh[forced_rank]->AddNodesFromMesh(*fbox,dist_offset,FLUID_NODE);
      
      /*
      for (i=1;i<=MaxRank;i++) {
        mesh[i]->RemoveNodesFromMesh(*fbox,FLUID_NODE);
      }
      */
      delete fbox;
      
    }
    if (ff>0) std::cout << "Resolution successfully forced on " << ff << " nodes.\n\n";

  }
  
  std::cout <<"Filling first 'coating' layers..." <<std::endl;
  d=0;
  if (mesh[0]->AddFirstLayers(d, !inputPars.isset("E"))<0){
    std::cerr << "\nError! Couldn't find coating layers. Aborting...\n\n";
    return (1);          
  }
  
  if (inputPars.isset("E")) {
    std::cout << "Setting external nodes..." <<std::endl;
  } else {
    std::cout << "Setting internal nodes..." <<std::endl;
  }
  dstart=d;

  int nl;
  i=Mesh::CALCULATE_RANK(dfactor,dstart);
  do {
    //add one new layer to the mesh internally determined according to its mindist
    //and return in <i> the rank of the mesh that could be of lower resolution with
    //respect to this mesh (of the givn rank <i>). So the added mesh MUST be already allocated!
    nl=mesh[i]->AddLayer(dstart,i,dfactor);
    if (nl<0) return (1);
    if (nl>0) {
      std::cout << "Added "<< nl <<" fluid nodes on mesh ranked " <<i<<" at mindist= " <<dstart <<std::endl;
      maxrank=MAX(i,maxrank);
    }
  }  while (nl>0);
  
  /*
  if (!box->FindOuterNode(p)) {
    p.set(1,1,1); // try the opposite vertex
    if (!box->FindOuterNode(p)) {
      std::cerr << RED<<"\nError! Couldn't find an outer node on the Bounding Box faces.\n"
      <<"Try to leave a larger empty space there! Aborting..."<<NOCOL<<"\n\n";
      return (1);
    }
  }
  //  std::cerr << RED <<"\nPeriodism not yet handled! Aborting..." << NOCOL<<"\n\n";
  //  return (1);
  
  if (n==0) {
    std::cerr << "\n"<<RED<<"WARNING: no fluid (inner) nodes found!\nMaybe a bad STL surface to WALL conversion.\nTry to re-run with a higher -w factor.\n"
    <<"If it persists... maybe a perforated surface?"<<NOCOL<<"\n\n";
    MaxRank=0; //forced to monomesh
  }
  */
  if (n_obox >0 || n_ibox >0) {
    if (!mesh[0]->SetIOnodesSearchDir()) {
      std::cerr << "\n"<<RED<<"WARNING: either INLET or OUTLET nodes aren't adjacent to fluid nodes!"<<NOCOL<<"\n\n";
    }
  }
  
  if (maxrank>0) { 
    std::cout << RED << "MULTIMESH mode"<<NOCOL <<"\n\n";
  } else {
    std::cout << RED << "\nMONOMESH mode"<<NOCOL <<"\n\n";
  }    
    
  if (!inputPars.isset("W")) {
    int wr=0;
    std::cout << "\nrepairing walls...\n";
    //for (int i=0; i<= maxrank;i++) wr += mesh[i]->RepairWalls();
    wr += mesh[0]->RepairWalls();
    std::cout << wr << " WALL nodes added\n\n";    

  }
    
#ifndef NOCOMP
  if (maxrank>0) {
    std::cout << "\ndoing compenetration..." <<std::endl;

    // 1st pass for voxel-based compenetration
    for (int i=1;i<=maxrank;i++) mesh[i]->SaturateNodes(true);  //accostamento
    // 2nd pass for complete voxel-based compenetration
    for (int i=1;i<=maxrank;i++) mesh[i]->SaturateNodes(false); //voxel filling
#if 1
    int nreptot=0;
    int nrep;
    for (int i=0;i<maxrank;i++) {
      nrep=mesh[i]->RepairSaturation();
      while(nrep >0){
        nreptot += nrep;
        std::cout << "\r"<< nreptot <<" rank "<<i<<" nodes added to complete compenetration.                  ";
        nrep=mesh[i]->RepairSaturation();
      };
      if (nrep<0) return (1);
      std::cout << "\n" ;
    }
#endif    
  }
#endif  
  //########## OUTPUT ################
  
  std::cout << "counting and writing mesh nodes...\n";
  // count meshes fluid nodes
  int nw,nmono,nf0;
  n=nw=nmono=0;
  for (int i=0; i<= maxrank;i++) {
    nf0 = mesh[i]->count(FLUID_NODE);
    nmono += nf0*pow(mesh[i]->getSpacing().s,3);
    n += nf0;
    nw += mesh[i]->count(WALL_NODE);
  }  
  int ni = mesh[0]->count(INLET_NODE);
  int no = mesh[0]->count(OUTLET_NODE);
  std::cout << std::endl << n << " fluid nodes found.\n";
  std::cout << nw << " wall nodes found.\n";
  std::cout << ni+no << " i/o nodes found.\n\n";
  
  if (maxrank>0 && n>0) std::cout << ">>> Monomesh fluid nodes / Multimesh fluid nodes= "
                                      << double(nmono)/double(n) <<"\n";
  
  // write found meshes
  std::cout << RED << "\nMax rank found = " <<maxrank << NOCOL << "\n\n"; 
  std::map<std::string,int> inlet,outlet; // the keys is the "string" equal to the node searchdir+veldir and the value the index
  std::map<int,std::string> inletdir,outletdir; // the keys is the i/o node index and the value the node normal vector (as a string)
  std::map<int,std::string> inletSdir,outletSdir; // the keys is the i/o node index and the value the node searchdir (as a string)
  std::string snorm;
  
  int indice=1;
  std::ofstream m;
  std::stringstream filename("");
  std::map<INT,Node>::iterator it;
  Point p(nx,ny,nz,1);

  for (int i=0; i<= maxrank;i++) {
//    filename <<outfile<<(1<<i)<<".dat";
    filename.str("");
    if (MaxRank>0) {
      filename <<outfile<<(1<<i)<<".dat";
    } else {
      filename <<outfile<<".dat";
    }
    m.open(filename.str().c_str(),std::ofstream::out  );
    for (it=mesh[i]->node().begin();it != mesh[i]->node().end();it++) {
      if (writeios) {
        if (it->second.getType()==INLET_NODE) {
          snorm=it->second.getSearchDir();
          //hidden nodes are NOT added to any group
          if  (snorm != "0 0 0") {
            snorm = snorm + " "+it->second.getVelDir();
            if (inlet.count(snorm)==0) {
              inlet[snorm]=indice;
              inletdir[indice]=it->second.getVelDir();
              inletSdir[indice]=it->second.getSearchDir();
              indice++;
            }
            it->second.setIndex(inlet[snorm]);
          }
        } else if (it->second.getType()==OUTLET_NODE) {
          snorm=it->second.getSearchDir();
          //hidden nodes are NOT added to any group
          if  (snorm != "0 0 0") {
            //snorm = snorm + " "+mesh[i].node(x,y,z).getVelDir();
            if (outlet.count(snorm)==0) {
              outlet[snorm]=indice;
                  // pressure nodes has dummy velocity 
              outletdir[indice]=" 0. 0. 0."; //mesh[i].node(x,y,z).getVelDir();
              outletSdir[indice]=it->second.getSearchDir();
              indice++;
            }
            it->second.setIndex(outlet[snorm]);
          }
        }
      } 
      p.set(it->first);        
      m <<" " << p.x(Ox) <<" " << p.y(Oy) <<" " << p.z(Oz) <<" " << it->second.getType() << std::endl;
    }
    m.close();
  }  

// write mesh headers
  std::cout << "writing mesh headers...\n";
  for (int i=0; i<= maxrank;i++) {
//    filename <<outfile<<(1<<i)<<".dat";
    filename.str("");
    if (MaxRank>0) {
      filename <<outfile<<(1<<i)<<".hdr";
    } else {
      filename <<outfile<<".hdr";
    }
    m.open(filename.str().c_str(),std::ofstream::out  );
    m << nx << " " << ny << " " << nz << std::endl;
    m << 5 << " " << mesh[i]->count(FLUID_NODE) <<
    " " << mesh[i]->count(WALL_NODE) <<
    " " << mesh[i]->count(INLET_NODE) <<
    " " << mesh[i]->count(OUTLET_NODE) << std::endl;
    m << mesh[i]->getSpacing().s << std::endl;    
    m.close();
  }  
  
/*
  mesh[0].printNeighs(17,6,10,std::cout);
  std::cout <<std::endl;
  mesh[1].printNeighs(8,12,6,std::cout);
  std::cout <<std::endl;
  mesh[0].printNeighs(8,12,7,std::cout);
*/  
  if (writeios && ni+no >0) {
    std::cout << "writing ios file...\n";
  
    
    filename.str("");
    if (MaxRank>0) {
      filename <<outfile<<"1.ios";
    } else {
      filename <<outfile<<".ios";
    }
    m.open(filename.str().c_str(),std::ofstream::out  );
    m << inlet.size()+outlet.size()<<"\n";
    for (std::map<std::string,int>::iterator it=inlet.begin(); it!=inlet.end(); ++it) {
      indice=it->second;
      m << indice <<" inlet flow "<< inletSdir[indice] <<" " << inletdir[indice];
      m << " "<< velin <<"\n";
    }

    for (std::map<std::string,int>::iterator it=outlet.begin(); it!=outlet.end(); ++it) {
      indice=it->second; 
      m << indice <<" outlet pressure "<< outletSdir[indice] <<" " << outletdir[indice];
      m << " "<< pressout <<"\n";
    }
    
    m << ni <<"\n";
    for (it=mesh[0]->node().begin();it != mesh[0]->node().end();it++) {
      if (it->second.getType()==INLET_NODE) {
        indice=it->second.getIndex();
        //force hidden inlet nodes to belong to the first inlet nodes group
        if (indice==0) indice=inlet.begin()->second;
        p.set(it->first);
        m<< p.x(Ox) <<" " << p.y(Oy) <<" "<< p.z(Oz) << " " << indice<<"\n"; 
      }
    }
    m << no <<"\n";
    for (it=mesh[0]->node().begin();it != mesh[0]->node().end();it++) {
      if ( it->second.getType()==OUTLET_NODE ) {
        indice=it->second.getIndex();
        //force hidden outlet nodes to belong to the first outlet nodes group
        if (indice==0) indice=outlet.begin()->second;
        p.set(it->first);
        m<< p.x(Ox) <<" " << p.y(Oy) <<" "<< p.z(Oz) << " " << indice<<"\n"; 
      }
    }
    m.close();
  }

  
// ======== WRITE XYZ FILE ============
  if (inputPars.isset("X")) {
    std::cout << "writing XYZ output file...\n";
    std::ofstream mxyz;
    //ofstream mmxyz; //???????????
    std::string nome[8]={"M1", "M2", "M4", "M8", "M16l", "M32", "M64", "M128"};
    //std::string nome[8]={"H", "Fe", "K", "Li", "Al", "He", "N", "C"};
    std::stringstream filenamexyz("");
    //stringstream filexxx(""); //?????????
    //filexxx <<outfile<<"b.xyz";  //????????
    filenamexyz <<outfile<<".xyz";
    //mmxyz.open(filexxx.str().c_str(),std::ofstream::out  ); //????????
    mxyz.open(filenamexyz.str().c_str(),std::ofstream::out  );
    mxyz << n+nw+ni+no << std::endl << std::endl;
    for (int i=0; i<= maxrank;i++) {
  //    filename <<outfile<<(1<<i)<<".dat";
      for (it=mesh[i]->node().begin();it != mesh[i]->node().end();it++) {
        switch (it->second.getType()) {
          case FLUID_NODE: mxyz << nome[i] ; break; 
          case WALL_NODE: mxyz << "W"; break; 
          case INLET_NODE: mxyz << "I"; break;
          case OUTLET_NODE: mxyz << "O"; break;
          default: mxyz << "X";
        }
        p.set(it->first);
        mxyz  <<" " << p.x(Ox) <<" " << p.y(Oy) <<" " << p.z(Oz) << "\n";
      }
    }  
    mxyz.close();
  }  
  
// ======== WRITE DEBUG FILE ============
  if (inputPars.isset("d")) {
    int probe[4];
    std::string nome[8]={"M1", "M2", "M4", "M8", "M16l", "M32", "M64", "M128"};
    std::ofstream dxyz;
    dxyz.open("debug.xyz",std::ofstream::out  );
    if (inputPars.GetParameter<int>("d",descr,probe,4) != 4) {
      std::cerr <<  "\n"<<RED<<"3 box coordinate and mesh rank must be given for debugging! Aborting..." <<NOCOL<<"\n\n";
      return (1);
    }
    dxyz <<  27*(maxrank+1)<< std::endl << std::endl;
    int i1=MAX(0,probe[0]-1);
    int i2=MIN(probe[0]+1,maxrank);
    int s1=mesh[i1]->getSpacing().s;
    int s2=mesh[i2]->getSpacing().s;
    for (int i=i1; i<= i2;i++) {
  //    filename <<outfile<<(1<<i)<<".dat";
      for (int z=probe[3]-s2; z <=probe[3]+s2; z += s1) {
        for (int y=probe[2]-s2; y <=probe[2]+s2; y += s1) {
          for (int x=probe[1]-s2; x <=probe[1]+s2; x += s1) {
            Point p(x,y,z,nx, ny, nz,mesh[i]->getSpacing().s);
            if (p.isOutOfBounds()) {
              dxyz << "X";                  
            } else {
              if (mesh[i]->isDead(p)) {
                dxyz << "X";                  
              } else {
                
                switch (mesh[i]->node(p).getType()) {
                  case FLUID_NODE: dxyz << nome[i] ; break; 
                  case WALL_NODE: dxyz << "W"; break; 
                  case INLET_NODE: dxyz << "I"; break;
                  case OUTLET_NODE: dxyz << "O"; break;
                  default: dxyz << "X";
                }
              }
            }
            dxyz  <<" " << x <<" " << y <<" " << z << "\n";
          }
        }
      }
    }  
    dxyz.close();
  }
  
//  mesh[0]->printNeighs(probe[0],probe[1],probe[2],dxyz);
//  box->Dump("dump.xyz","N"); //?????????'

  std::cout << "\nCoordinates offset: ("<<mesh[0]->getO_x()<<","<<mesh[0]->getO_y()<<","<<mesh[0]->getO_z()<<")\n\n";
  for (i=0;i<=MaxRank;i++) {delete mesh[i];}
  delete[] mesh;
  std::cout << "done!" <<std::endl;
  return(0);
  
}

