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

char STLFNAME[50];
char OUTFNAME[50];
bool INTERPOLATE = false;

// we will use the Taubin lambda/mu smoothing algorithm with an optimal k_PB=0.1
float K_PB = 0.1;
float SCALE_POS = 0.3;
float SCALE_NEG = SCALE_POS/(K_PB*SCALE_POS-1.0);
int ITERATIONS = 10;

struct Point { 

    public:
        float x,y,z;

        Point(){};
        Point (float xx, float yy, float zz){ x = xx; y = yy; z = zz; };
        //Point& operator= (const Point &cSource);
};


struct Edge {

    public:
        Point p1,p2;

        Edge(){};
        Edge (const Point &pp1, const Point &pp2){p1=pp1; p2=pp2;}
        //Edge& operator= (const Edge &cSource);
};

struct Triangle {

    public:
        Edge e1,e2,e3;

        Triangle (){};
        Triangle (const Edge &ee1, const Edge &ee2, const Edge &ee3){ e1=ee1; e2=ee2; e3=ee3;}
        Triangle (const Point &pp1, const Point &pp2, const Point &pp3){ 
                Edge ee1(pp1,pp2);
                Edge ee2(pp2,pp3);
                Edge ee3(pp3,pp1);
                e1=ee1; e2=ee2; e3=ee3;
        }
};

Point const PNULL(1.,2.,3.);
Triangle const TNULL(PNULL,PNULL,PNULL);

bool operator== (const Point &l, const Point &r) { return (l.x == r.x && l.y == r.y && l.z == r.z); }
bool operator!= (const Point &l, const Point &r) { return !(l == r); }
bool operator< (const Point &l, const Point &r) { 
    if (l.x != r.x) return (l.x < r.x);
    else
    {
        if (l.y != r.y) return (l.y < r.y);
        else return (l.z < r.z);
    }
}

bool operator== (const Edge &l, const Edge &r) { return (l.p1 == r.p2); }
bool operator!= (const Edge &l, const Edge &r) { return !(l == r); }
bool operator< (const Edge &l, const Edge &r) { 
    if (l.p1 == r.p1) return (l.p2 < r.p2);
    else return (l.p1 < r.p1);
}

bool operator< (const Triangle &l, const Triangle &r) { 
    if (l.e1 != r.e1) return (l.e1 < r.e1);
    else
    {
        if (l.e2 != r.e2) return (l.e2 < r.e2);
        else return (l.e3 < r.e3);
    }
}

////////////////////
void usage(int argc, char *argv[])
{
    cout << "Usage:" << argv[0] << " -g stl_file" << endl;
	cout << "Options:" << endl;
	cout << "\t-g stl_file" << endl;
	//cout << "\t--geometry stl_file" << endl;
	cout << "\t\tSpecifies the STL file to smooth." << endl;
	cout << "\t-i" << endl;
	//cout << "\t--interpolate" << endl;
	cout << "\t\tPerforms triangle interpolation on input mesh (default: no)." << endl;
	cout << "\t-s scale_factor" << endl;
    //cout << "\t--scale scale_factor" << endl;
    cout << "\t\tSpecifies the Laplacian smoothing scale factor (default: 0.3)." << endl;
	cout << "\t-n iterations" << endl;
    //cout << "\t--numiter iterations" << endl;
    cout << "\t\tSpecifies the number of times Laplacian smoothing has to be applied (default: 5)." << endl;
	cout << "\t-o out_file" << endl;
    //cout << "\t--output out_file" << endl;
    cout << "\t\tSpecifies the output file." << endl;
    cout << "\t\tIf this option is not set the standard output is used." << endl;
}

////////////////////
void print_point(Point p) {cout << "Point: " << p.x << " " << p.y << " " << p.z << endl;}

////////////////////
/*
void print_triangle(Triangle t) {
    // cout << "Triangle: " << t.e1 << " " << t.e2 << " " << t.e3 << endl;
}
*/

// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}

////////////////////
string readline(ifstream &fstl){

     string line;

     // fstl >> line; // read one word at a time
     getline(fstl, line); // read entire line
     line = ltrim(line);
     // cout << "LINE:" << line << endl;

     return line;

}
////////////////////
void tuple10(string fields[10], string line){;

    string sub;

    istringstream iss(line);
    int i=0;
    do{
        iss >> sub;
        // if (strcmp(&sub[0],"") == 0) continue;
        if (strcmp(&sub[0],"") == 0) break;
        fields[i] = sub;
        i++;
    } while(iss);
    i--;

    //for(int ii=0; ii<i; ii++) cout << "Field: " << fields[ii] << " "; cout << "\n";

}


////////////////////
bool point_compare(Point *p1, Point *p2) {

    if(p1->x==p2->x && p1->y==p2->y && p1->z==p2->z) return 0;
    return 1;
}

////////////////////
float dot_prod(float a[],float b[]){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
////////////////////
void vec_prod(float *a, float *b, float *res){
	res[0] = a[1]*b[2]-a[2]*b[1];
    res[1] = a[2]*b[0]-a[0]*b[2];
    res[3] = a[0]*b[1]-a[1]*b[0];
}

////////////////////
bool searchmap(Point *p, map<Point, int> *mp){

    for(map<Point, int>::iterator it = mp->begin(); it != mp->end(); it++){
         if (it->first == *p) return true;
    }

    return false;
}
////////////////////
int main(int argc, char *argv[])
{


    map <Point, int> vertex2id;
    map <int, Point> id2vertex;

    for (int n = 1; n < argc; n++) {
        char *o = argv[n];
        // cout << n << ": " << argv[ n ] << " " << o << endl;
        if (strcmp(o,"-g") == 0) {
            n++;
            strcpy(STLFNAME, argv[n]);
            cout << "STL FILE: " << STLFNAME << endl;
        }
        else if (strcmp(o,"-o") == 0) {
            n++;
            strcpy(OUTFNAME, argv[n]);
            cout << "OUT FILE: " << OUTFNAME << endl;
        }
        else {
            usage(argc, &argv[0]);
            exit(EXIT_FAILURE);
        }
    }
    if (STLFNAME[0] == char(0)) {
        cout << "undefined STL input file" << endl;
        exit(1);
    }

    cout << "Reading STL file ...." << STLFNAME << endl;

    ifstream fstl(STLFNAME);
    string line;
    string sub;
    string fields[10];
    float x,y,z;

    set<Triangle> triangles;

    int dupTris=0;
    int cntlines=0;
    int count=0;

    line = readline(fstl);// solid ascii
    while (true){

        line = readline(fstl);// read facet normal

        if (cntlines%10000==0) cout << "line: " << line << endl;
        cntlines++;

        tuple10(fields, line);

        if (fields[0] == "endsolid") break;

        x = (float)atof(fields[2].c_str());
        y = (float)atof(fields[3].c_str());
        z = (float)atof(fields[4].c_str());

        float n[3] = {x,y,z};

        line = readline(fstl);// 'outer loop'

        line = readline(fstl);// 'vertex 1'
        tuple10(fields, line);
        x = (float)atof(fields[1].c_str());
        y = (float)atof(fields[2].c_str());
        z = (float)atof(fields[3].c_str());
        Point p1(x,y,z);

        if ( vertex2id.find(p1) == vertex2id.end() ) { // not found
            vertex2id[p1] = count;
            id2vertex[count] = p1;
            count = count+1;
        }

        line = readline(fstl);// 'vertex 2'
        tuple10(fields, line);
        x = (float)atof(fields[1].c_str());
        y = (float)atof(fields[2].c_str());
        z = (float)atof(fields[3].c_str());
        Point p2(x,y,z);

        if (vertex2id.find(p2) == vertex2id.end()) { // not found
            vertex2id[p2] = count;
            id2vertex[count] = p2;
            count = count+1;
        }

        line = readline(fstl);// 'vertex 3'
        tuple10(fields, line);
        x = (float)atof(fields[1].c_str());
        y = (float)atof(fields[2].c_str());
        z = (float)atof(fields[3].c_str());
        Point p3(x,y,z);

        if (vertex2id.find(p3) == vertex2id.end()) { // not found
            vertex2id[p3] = count;
            id2vertex[count] = p3;
            count = count+1;
        }

        Point p12;
        Point p23;
        Point p31;
        if(INTERPOLATE){

            // splice each triangle in 4 by using middle points in each side
            Point p12 ((p1.x+p2.x)/2, (p1.y+p2.y)/2, (p1.z+p2.z)/2);
            if (vertex2id.find(p12) == vertex2id.end()) { // (not p12 in vertex2id)
                vertex2id[p12] = count;
                id2vertex[count] = p12;
                count = count+1;
            }
	
            Point p23 ((p2.x+p3.x)/2, (p2.y+p3.y)/2, (p2.z+p3.z)/2);
            if (vertex2id.find(p23) == vertex2id.end()) { // (not p23 in vertex2id)
                vertex2id[p23] = count;
                id2vertex[count] = p23;
                count = count+1;
            }
		
            Point p31 ((p1.x+p3.x)/2, (p1.y+p3.y)/2, (p1.z+p3.z)/2);
            if (vertex2id.find(p31) == vertex2id.end()) { // (not p31 in vertex2id)
                vertex2id[p31] = count;
                id2vertex[count] = p31;
                count = count+1;
            }
        }

        // insert trinagles (p1,p2,p3) such that (p2-p1)X(p3-p1) is directed as the normal n
        float side1[3] = {p2.x-p1.x, p2.y-p1.y, p2.z-p1.z};
        float side2[3] = {p3.x-p1.x, p3.y-p1.y, p3.z-p1.z};

        float res[3];
        vec_prod(side1, side2, res);
        if (dot_prod(res, n) > 0.0){
            if(!INTERPOLATE){
                triangles.insert( Triangle (p1,p2,p3) ); 
            }
            else {
                triangles.insert( Triangle (p1,p12,p31) );
                triangles.insert( Triangle (p12,p2,p23) );
                triangles.insert( Triangle (p12,p23,p31) );
                triangles.insert( Triangle (p23,p3,p31) );
            }
        }
        else{
            if(!INTERPOLATE){
                triangles.insert( Triangle (p2,p1,p3) ); 
            }
            else {
                triangles.insert( Triangle (p12,p1,p31) );
                triangles.insert( Triangle (p2,p12,p23) );
                triangles.insert( Triangle (p23,p12,p31) );
                triangles.insert( Triangle (p3,p23,p31) );
            }
        }

        line = readline(fstl);// endloop
        line = readline(fstl);// endfacet
    }
    fstl.close();

    cout << "Using " << triangles.size() << "triangles " << count << "vertices.";
    // we no longer need vertex2id...
    //delete vertex2id;

    map <Triangle, set <Triangle> > stream;
    for(set<Triangle>::iterator t = triangles.begin(); t != triangles.end(); t++){

        list<Triangle> tmp;
        tmp.push_back(*t);
        tmp.sort();

        list<Triangle>::iterator it0,it1,it2;
        it0 = tmp.begin();
        it1=it0; it1++;
        it2=it1; it2++;

        set<Triangle> ss;
        if(stream.find(*it0) == stream.end()){
            ss.clear();
            stream[*it0] = ss;
        }
        if(stream.find(*it1) == stream.end()){
            ss.insert(*it1);
            stream[*it0] = ss;
            //degree
        }


    }



    return 0;
}



