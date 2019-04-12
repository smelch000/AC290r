#ifndef GENMESH
#define GENMESH
#define VERSION "2.3.3 - 09/12/18"
#define MIN(j1,j2) ( j1 < j2 ? j1 : j2 )
#define MAX(j1,j2) ( j1 > j2 ? j1 : j2 )

#include <iostream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include "double3.h"

#define MAXINLETBOXES 32
#define MAXOUTLETBOXES 32
#define MAXFORCINGBOXES 32
#define MAXNUMTRIANGLES 100000

#define INT long

#define RED "\033[1;31m"
#define NOCOL "\033[0m"

enum NODETYPE {FLUID_NODE=1, WALL_NODE, INLET_NODE, OUTLET_NODE, DEAD_NODE};
enum INPUTTYPE  {UNKNOWN, DAT, STL}; 
struct Spacing {unsigned int s; int rank;}; // rank=log_2 (s)

#endif
