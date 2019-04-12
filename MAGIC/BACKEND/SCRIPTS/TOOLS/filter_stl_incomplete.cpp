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
int THRESHOLD=0.;
char OUTFNAME[50];

// to use std::map the operator < must be implemented

struct Point { float r[3]; };

struct Edge { Point p1,p2; };

struct Triangle { Edge i1,i2,i3; };
// struct Triangle { int i1, i2, i3; };

struct V3 { float f1,f2,f3; };

typedef map <Point, int> Points2id;
typedef map <int, Point> Id2points;

typedef map <Edge, int> Edges2id;
typedef map <int, Edge> Id2edges;

typedef map <Triangle, int> Triangles2id;
typedef map <int, Triangle > Id2triangles;

typedef map <Triangle, vector<Triangle> > Triangle2triangles;
// typedef map <int, vector<int> > Id2ids;

// typedef map <int, vector<int> > Edges2tris;
// typedef map <int, vector<Triangle> > Edges2tris;
typedef map <Edge, vector<Triangle> > Edges2tris;

Triangle const TNULL = {-1,-1,-1};

bool operator== (const Point &l, const Point &r) { return (l.r[0] == r.r[0] &&
                                                           l.r[1] == r.r[1] &&
                                                           l.r[2] == r.r[2]); }
bool operator!= (const Point &l, const Point &r) { return !(l==r); }

bool operator== (const Edge &l, const Edge &r) { return (l.p1 == r.p2); }
bool operator!= (const Edge &l, const Edge &r) { return !(l == r); }

bool operator== (const Triangle &l, const Triangle &r) { return (l.i1 == r.i1 &&
                                                                 l.i2 == r.i2 &&
                                                                 l.i3 == r.i3); }

bool operator< (const Point &l, const Point &r) { 
    if (l.r[0] != r.r[0]) {
        return (l.r[0] < r.r[0]);
    }
    else
    {
        if (l.r[1] != r.r[1]) {
            return (l.r[1] < r.r[1]);
        }
        else
        {
            return (l.r[2] < r.r[2]);
        }
    }
    // return (l.r[0] < r.r[0] l.r[1] < r.r[1] l.r[2] < r.r[2]);
}

bool operator< (const Edge &l, const Edge &r) { 
    if (l.p1 == r.p1) {
        return (l.p2 < r.p2);
    }
    else
    {
        return (l.p1 < r.p1);
    }
    //return (l.p1 < r.p1 && l.p2 < r.p2); 
}

bool operator< (const Triangle &l, const Triangle &r) { 
    if (l.i1 != r.i1) {
        return (l.i1 < r.i1);
    }
    else
    {
        if (l.i2 != r.i2) {
            return (l.i2 < r.i2);
        }
        else
        {
            return (l.i3 < r.i3);
        }
    }
    // return (l.i1 < r.i1 && l.i2 < r.i2 && l.i3 < r.i3);
}

bool operator< (const V3 &l, const V3 &r) { 
    if (l.f1 != r.f1) {
        return (l.f1 < r.f1);
    }
    else
    {
        if (l.f2 != r.f2) {
            return (l.f2 < r.f2);
        }
        else
        {
            return (l.f3 < r.f3);
        }
    }
    // return (l.f1 < r.f1 && l.f2 < r.f2 && l.f3 < r.f3);
}

////////////////////
void usage(int argc, char *argv[])
{
    cout << "Usage:" << argv[0] << " -g stl_file [-t volume_threshold] [-o out_file]" << endl;
	cout << "Options:" << endl;
	cout << "\t-g stl_file" << endl;
	//cout << "\t--geometry stl_file" << endl;
	cout << "\t\tSpecifies the STL file to filter." << endl;
	cout << "\t-t volume_threshold" << endl;
	//cout << "\t--threshold volume_threshold" << endl;
	cout << "\t\tSpecifies the volume threshold (between 0.0 and 1.0) below which connected components will be removed (default: 0.0)." << endl;
	cout << "\t-o out_file" << endl;
    //cout << "\t--output out_file" << endl;
    cout << "\t\tSpecifies the output file." << endl;
    cout << "\t\tIf this option is not set the standard output is used." << endl;
}

////////////////////
void print_point(Point p) {cout << "Point: " << p.r[0] << " " << p.r[1] << " " << p.r[2] << endl;}

////////////////////
void print_plist(string str,Points2id points2id) {

    for(Points2id::iterator it = points2id.begin(); it != points2id.end(); it++){
        cout << str << it->first.r[0] << " " << it->first.r[1] << " " << it->first.r[2] << "  (" << it->second << endl;
    }
}

////////////////////
void print_triangle(Triangle t) {
    // cout << "Triangle: " << t.i1 << " " << t.i2 << " " << t.i3 << endl;
}

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

    if(p1->r[0]==p2->r[0] && p1->r[1]==p2->r[1] && p1->r[2]==p2->r[2]) return 0;
    return 1;
}

////////////////////
// void sort3(int i1, int i2, int i3, Triangle *tri){
void sort3(Edge i1, Edge i2, Edge i3, Triangle *tri){

    if (i1 < i2 && i1 < i3 ) {
        tri->i1 = i1;
        if (i2 < i3 ) {
            tri->i2 = i2; tri->i3 = i3;
        }
        else{
            tri->i2 = i3; tri->i3 = i2;
        }
    }
    if (i2 < i1 && i2 < i3 ) {
        tri->i1 = i2;
        if (i1 < i3 ) {
            tri->i2 = i1; tri->i3 = i3;
        }
        else{
            tri->i2 = i3; tri->i3 = i1;
        }
    }
    if (i3 < i1 && i3 < i2 ) {
        tri->i1 = i3;
        if (i1 < i2 ) {
            tri->i2 = i1; tri->i3 = i2;
        }
        else{
            tri->i2 = i2; tri->i3 = i1;
        }
    }
}

////////////////////
// bool fun_less(int k, vector<Triangle> v) { 
bool fun_less(Triangle k, vector<Triangle> v) { 
    return (v.size()<3); 
}
////////////////////
bool fun_not(Triangle k, int v) { 
    return (v==-1); 
}
////////////////////
Triangle find_first_key_1(bool (*fp)(Triangle, vector<Triangle>), Triangle2triangles dic)
{ 
    Triangle ret = TNULL;
    for (Triangle2triangles::iterator it = dic.begin(); it != dic.end(); it++){
        Triangle key = it->first;
        vector<Triangle> val = it->second;
        if (! fp(key,val)) continue;
        ret = key;
        break;
    }
    return ret;
}
////////////////////
Triangle find_first_key_2(bool (*fp)(Triangle, int), Triangles2id dic)
{ 
    Triangle ret = TNULL;
    for (Triangles2id::iterator it = dic.begin(); it != dic.end(); it++){
        Triangle key = it->first;
        int val = it->second;
        if (! fp(key,val)) continue;
        ret = key;
        break;
    }
    return ret;
}
////////////////////
int main(int argc, char *argv[])
{

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

    int p1Id, p2Id, p3Id;
    int e1Id, e2Id, e3Id;

    V3 norm;

    Points2id points2id;
    Id2points id2points;

    Edges2id edges2id;
    Id2edges id2edges;

    Id2triangles id2triangles;
    map <Triangle, int> triangles2id;

    map <int, V3> id2normals;
    Edges2tris edges2tris;

    int dupTris=0;
    int cntlines=0;

    line = readline(fstl);// solid ascii
    while (true){

        line = readline(fstl);// read facet normal

        if (cntlines%10000==0) cout << "line: " << line << endl;
        cntlines++;

        tuple10(fields, line);

        if (fields[0] == "endsolid") break;

        norm.f1 = (float)atof(fields[2].c_str());
        norm.f2 = (float)atof(fields[3].c_str());
        norm.f3 = (float)atof(fields[4].c_str());

        line = readline(fstl);// 'outer loop'
        line = readline(fstl);// 'vertex 1'

        tuple10(fields, line);
        Point p1;
        p1.r[0] = (float)atof(fields[1].c_str());
        p1.r[1] = (float)atof(fields[2].c_str());
        p1.r[2] = (float)atof(fields[3].c_str());

        line = readline(fstl);// 'vertex 2'
        tuple10(fields, line);
        Point p2;
        p2.r[0] = (float)atof(fields[1].c_str());
        p2.r[1] = (float)atof(fields[2].c_str());
        p2.r[2] = (float)atof(fields[3].c_str());

        line = readline(fstl);// 'vertex 3'
        tuple10(fields, line);
        Point p3;
        p3.r[0] = (float)atof(fields[1].c_str());
        p3.r[1] = (float)atof(fields[2].c_str());
        p3.r[2] = (float)atof(fields[3].c_str());

        // if (cntlines==50) exit(1);

        if (points2id.find(p1) == points2id.end()) {
            p1Id = points2id.size();
            points2id[p1] = p1Id;
        }
        else {
            p1Id = points2id[p1];
        }

        if (points2id.find(p2) == points2id.end()) {
            p2Id = points2id.size();
            points2id[p2] = p2Id;
        }
        else {
            p2Id = points2id[p2];
        }

        if (points2id.find(p3) == points2id.end()) {
            p3Id = points2id.size();
            points2id[p3] = p3Id;
        }
        else {
            p3Id = points2id[p3];
        }

        // ignore zero-are triangles
        if (p1Id == p2Id || p1Id == p3Id || p1Id == p3Id) {
            line = readline(fstl);// endloop
            line = readline(fstl);// endfacet
            continue;
        }

        // build the 3 edges (pA, pB) s.t. pA < pB
        Edge e1 = {min(p1Id, p2Id), max(p1Id, p2Id)};
        Edge e2 = {min(p1Id, p3Id), max(p1Id, p3Id)};
        Edge e3 = {min(p2Id, p3Id), max(p2Id, p3Id)};

        if (edges2id.find(e1) == edges2id.end()) {
            e1Id = edges2id.size();
            edges2id[e1] = e1Id;
        } else {
            e1Id = edges2id[e1];
        }

        if (edges2id.find(e2) == edges2id.end()) {
            e2Id = edges2id.size();
            edges2id[e2] = e2Id;
        } else {
            e2Id = edges2id[e2];
        }

        if (edges2id.find(e3) == edges2id.end()) {
            e3Id = edges2id.size();
            edges2id[e3] = e3Id;
        } else {
            e3Id = edges2id[e3];
        }

        // Triangle tt={1,1,1};
        Triangle tt=TNULL;
        // sort3(e1Id, e2Id, e3Id, &tt);
        sort3(e1, e2, e3, &tt);
        Triangle tri = tt;
        //print_triangle(tri);

        // discards duplicated trinagles
        if (triangles2id.find(tri) == triangles2id.end()) {

            int triId = triangles2id.size();
            triangles2id[tri] = triId;
            id2normals[triId] = norm;
	
            // if (edges2tris.find(e1Id) == edges2tris.end()) {
            if (edges2tris.find(e1) == edges2tris.end()) {
                //edges2tris[e1Id].push_back(tri);
                //edges2tris[e1Id].pop_back();
                edges2tris[e1].push_back(tri);
                edges2tris[e1].pop_back();
            }
            // edges2tris[e1Id].push_back(tri);
            edges2tris[e1].push_back(tri);
		
            // if (edges2tris.find(e2Id) == edges2tris.end()) {
            if (edges2tris.find(e2) == edges2tris.end()) {
                //edges2tris[e2Id].push_back(tri);
                //edges2tris[e2Id].pop_back();
                edges2tris[e2].push_back(tri);
                edges2tris[e2].pop_back();
            }
            // edges2tris[e2Id].push_back(tri);
            edges2tris[e2].push_back(tri);
		
            // if (edges2tris.find(e3Id) == edges2tris.end()) {
            if (edges2tris.find(e3) == edges2tris.end()) {
                //edges2tris[e3Id].push_back(tri);
                //edges2tris[e3Id].pop_back();
                edges2tris[e3].push_back(tri);
                edges2tris[e3].pop_back();
            }
            // edges2tris[e3Id].push_back(tri);
            edges2tris[e3].push_back(tri);
        }
	    else {
		    dupTris++;
        }

        line = readline(fstl);// endloop
        line = readline(fstl);// endfacet
    }
    fstl.close();

    // remove dictionaries no longer needed and invert the others
    cout << "total points:" << points2id.size() << endl;
    for (Points2id::const_iterator it = points2id.begin(); it != points2id.end(); ++it) {
        id2points[it->second] = it->first;
    }
    points2id.clear();

    cout << "total edges:" << edges2id.size() << endl;
    for (Edges2id::const_iterator it = edges2id.begin(); it != edges2id.end(); ++it) {
        id2edges[it->second] = it->first;
    }
    edges2id.clear();

    cout << "total triangles:" << triangles2id.size() << endl;
    for (Triangles2id::const_iterator it = triangles2id.begin(); it != triangles2id.end(); ++it) {
        // id2triangles[it->second].push_back(it->first);
        id2triangles[it->second] = it->first;
    }
    cout << "Removed: " << dupTris << " duplicated triangles over a total:" << triangles2id.size() << endl;
    triangles2id.clear();

    cout << "Constructing triangle adjacency list..." << endl;

    // construct the adjacency list of triangle graph
    Triangle2triangles triGraph;
    // Id2triangles triGraph;
    // Id2ids triGraph;
    Triangle tdum = {-99, -99, -99};
    int outertriangles=0;
    for (Id2triangles::const_iterator it = id2triangles.begin(); it != id2triangles.end(); ++it) {
        // push and pop an element to create an empty list
        triGraph[it->second].push_back(tdum);
        triGraph[it->second].pop_back(); 
        outertriangles++;
    }
    cout << "total triGraph:" << triGraph.size() << endl;
    cout << "outer triangles:" << outertriangles << endl;

    cout << "total edges2tris:" << edges2tris.size() << endl;

    for (Edges2tris::iterator it = edges2tris.begin(); it != edges2tris.end(); it++){

        // int eid = it->first;
        // vector <int> triIds = it->second;
        vector <Triangle> triIds = it->second;
        // cout << "triIds " << triIds[0] << " " << triIds[1] << endl;

        if (triIds.size()==1) {
            cout << "Edge belongs to a single triangle !" << endl;
        } else {
            //for(vector<int>::iterator tid1 = triIds.begin(); tid1!=triIds.end(); tid1++){
            //    for(vector<int>::iterator tid2 = triIds.begin(); tid2!=triIds.end(); tid2++){
                    // if(*tid1==*tid2) continue;
            for(vector<Triangle>::iterator tid1 = triIds.begin(); tid1!=triIds.end(); tid1++){
                for(vector<Triangle>::iterator tid2 = triIds.begin(); tid2!=triIds.end(); tid2++){
                    if((*tid1).i1==(*tid2).i1 && (*tid1).i2==(*tid2).i2 && (*tid1).i3==(*tid2).i3) continue;
                    // cout << *tid1 << " " << *tid2 << endl;
                    triGraph[*tid1].push_back(*tid2);
                };
            };
        }
    }


    // perform a Breadth-first search to iteratively remove triangles 
    // with less than 3 neighboring triangles from the adjacency list
    cout << "Removing loosely connected triangles..." << endl;
    vector <Triangle> queue;
    // list <int> queue;
    do{
        if (queue.size()==0) {
            bool (*fp)(Triangle k, vector<Triangle> v);
            fp = fun_less;
            Triangle seed = find_first_key_1(fp, triGraph);
            if (seed==TNULL) break;
            queue.clear();
            queue.push_back(seed);
        }

        // take curr triangle to remove
	    Triangle cur = *(queue.begin());
	    queue.pop_back();
	    for (vector<Triangle>::iterator neigh=triGraph[cur].begin(); neigh!=triGraph[cur].end(); neigh++)
        {
	    	// remove curr from the adjacency list of its neighbors
	    	triGraph[*neigh].erase(queue.begin());

	    	if (triGraph[*neigh].size() < 3){
	    		// if the degree of neighbor dropped below 2 AND
	    		// it's not already in the queue...
	    		if (find(queue.begin(), queue.end(), *neigh)!=queue.end()) queue.push_back(*neigh);
            }
	    // remove curr from the queue
	    triGraph.erase(*queue.begin());
        }

    } while(true);

    cout << "done" << endl;

    // perform a Breadth-first search to:
    // - find all the connected components and compute their volume
    // - associate every tringle to the connected component it belongs to
    map<int, float> component;
    Triangles2id tri2comp;

    for (Triangle2triangles::iterator tid=triGraph.begin(); tid!=triGraph.end(); tid++)
        tri2comp[tid->first] = -1;

    cout << "Analyzing connected components..." << endl;

    queue.clear();
    int compId = -1;
    do{
        if (queue.size()==0) {
            bool (*fp)(Triangle k, int v);
            fp = fun_not;
            Triangle seed = find_first_key_2(fp, tri2comp);
            if (seed == TNULL) break;

            compId+=1;
            queue.clear();
            queue.push_back(seed);
            tri2comp[seed] = compId;

            // vector<Triangle> t = id2triangles[seed];
            vector<Triangle> t = seed;
            Edge e0 = t.e1;
            Edge e1 = t.e2;
            Edge e2 = t.e3;
            Edge ee = id2edges[t0];
            // int i = id2triangles[seed][0]->first;
            // set<Triangle> s(t.begin(),t.end());
            // list<Triangle> lst;
            // for(set<Triangle>::const_iterator is = s.begin(); is!=s.end(); is++) 
            //     lst.push_back(*is);

            Edge e = id2edges[0];// + t[1] + t[2]];

            // triPointIds = list(set(id2edges[id2triangles[seed][0]] +
            //                        id2edges[id2triangles[seed][1]] +
            //                        id2edges[id2triangles[seed][2]]))
            // python set: unordered set of unique elements
            // python list: unordered set of elements
            // eg:
            // a=list(set([1,2,1])) returns a=[1,2,1]
            // a=list(set([1,2,1])) returns [1,2]

        }

    } while(true);


    for (Triangle2triangles::iterator tid = triGraph.begin(); tid != triGraph.end(); tid++){

        //vector<Triangle> t = tid->second;
        //vector<Triangle> edges = id2triangles[tid->first];

    }



    return 0;
}



