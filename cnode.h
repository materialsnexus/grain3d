#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>


#ifndef CNODE_H__
#define CNODE_H__

class CBody;
struct Vertex;

class CNode {
public:
    double x, y, z;
    double dx, dy, dz;
    double cx, cy, cz;
    int    id;
    bool   good;
    bool   erase;
    bool   needsfixing;
    double area;
    std::vector <CBody*> neighbors;
    std::vector <CNode*> corners[4]; // topological info for vertex nodes
    std::vector <std::string> labels;
    
    CNode* adj_faces[3];
    
    Vertex* vertex;
    
    
    CNode(void) :  x(0),  y(0),  z(0),
    dx(0), dy(0), dz(0),
    cx(0), cy(0), cz(0),
    id(0), erase(0),
    needsfixing(0),
    area(0) {}
    CNode operator - (CNode);
    CNode operator + (CNode);
    CNode operator *(double);
    CNode operator /(double);
    void  operator  =(CNode);
    void  operator +=(CNode);
    void  operator /=(double);
    void  operator *=(double);
    
    double norm();
    void   calc_face_area();
    double get_face_area();
    double get_face_perimeter();
    
    void   calc_curvatures();
    void   calc_motion();
    
    void   move();
    
    void   determine_node_topology();
    int    get_face_sides();
    void   remove_tri_face();
    bool   remove_face();         // maybe change this into a class function?
    void   remove_digon();
    void   remove_edge_node();
    void   remove_edge_nodes();
    void   output();
    void   center();
    void   center2(double);
    void   calc_center();
    void   move_center();
    
    double calc_change();
};


class sequence {
public:
    unsigned char number[256];
    CNode* origin;
    CNode* second;
};


class edge {
public:
    CNode* one;
    CNode* two;
    
    edge(CNode* A, CNode* B) : one(A), two(B) {}
};


class face {
public:
    CNode* one;
    CNode* two;
    CNode* thr;
    face(CNode* A, CNode* B, CNode* C) : one(A), two(B), thr(C) {}
};


struct Vertex
{
    unsigned int id;
    bool touched;
    CNode*  node;
    Vertex* origin;
    Vertex* neighbors[4];
    Vertex* originals[4];
    Vertex* turns[3];
    
    //THESE ARE INSTEAD OF THREADS, WE NEVER HAVE MORE THAN 4, MAYBE 5 OR 6
    // THREADS
    char    labels[8][16];
    int     orders[8];
    Vertex* origins[8];
    int tcount;  // THIS COUNTS HOW MANY THREADS WE'RE UP TO
};

struct triplet {
public:
    triplet(void) : v1(0), v2(0), v3(0) {}
    triplet(CNode* a, CNode* b, CNode* c) : v1(a), v2(b), v3(c) {}
    CNode *v1, *v2, *v3;
    double area();
};


#endif

