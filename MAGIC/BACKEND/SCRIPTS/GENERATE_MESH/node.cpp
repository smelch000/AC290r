/////////////////--Node--////////////////
#include "node.h"
//#include <iostream>



Node & Node::operator=(Node const & source) {
  if (&source != this) {
    type=source.getType();
    dist=source.getDist();
    //insature=source.isInsature();
    //outer=source.isOuter();
    //explored=source.isExplored();
    veldir=source.getVelDir();
    searchdir=source.getSearchDir();
  }
  return *this;
}


