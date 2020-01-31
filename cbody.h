#include "cnode.h"
#include <vector>

#ifndef CBODY_H__
#define CBODY_H__


class CBody {
public:
    int  id;
    double volume;
    std::vector <triplet> triplets;
    
//    int faces[50];       // for p vector
//    int fcount;
    
    double center[3];
    
    CBody(void) : id(0) {}

    int    get_faces();         // counts the body's number of faces
    double get_area();          // measures the body's surface area
    double get_volume();        // measures the body's volume
    double get_pvolume();       // measures the body's projected volume, after motions
    void   calc_volume();       // measures the body's volume
    
    int    get_triplet(CNode*, CNode*);
    
    bool   is_tetrahedron();
    bool   is_football();
    void   simple_remove_face            (CNode*);
    void   replace_edges_with_node       (CNode*, CNode*, CNode*, CNode*);
    void   replace_node_on_face_with_edge(CNode*, CNode*, CNode*, CNode*);
    void   switch_edge_with_edge         (CNode*, CNode*, CNode*, CNode*);
    void   remove_tetrahedron();
    void   remove_football();
    
    void   determine_node_topologies();
    void   output();
    void   remove_body();
    void   blank();
    
    void   evolver_out();
    void   calc_center();
};



#endif

