#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <list>
#include <math.h>


using namespace std;

struct Point { 

    public:
        float x,y,z;

        //Point(){};
        // Point (float xx, float yy, float zz){ x = xx; y = yy; z = zz; };
        // Point& operator= (const Point &cSource);
};

bool operator< (const Point &l, const Point &r) { 
    if (l.x != r.x) return (l.x < r.x);
    else
    {
        if (l.y != r.y) return (l.y < r.y);
        else return (l.z < r.z);
    }
}


typedef map <Point, int> Point2id;
typedef map <int, Point> Id2point;

int main(int argc, char *argv[])
{

    Point2id vertex2id;
    Id2point id2vertex;

    int count=0;
    Point p1;
    vertex2id[p1] = count;
    id2vertex[count] = p1;

    return 0;
}



