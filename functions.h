#ifndef FUNCTIONS_H__
#define FUNCTIONS_H__

double det   (CNode*, CNode*, CNode*, CNode*);  // determinant of a matrix
double pdet   (CNode*, CNode*, CNode*, CNode*);  // determinant of a matrix
double length(CNode*, CNode*);
void   calc_and_print_stats();
int    count_bodies();
int    count_vertices();
int    count_faces();
int    count_edges();
void   calc_volumes();
void   calc_areas();
void   calc_motions();
double calc_avg_face_area();
double calc_avg_cell_area();
double calc_avg_cell_faces();
double calc_avg_face_sides();

bool   is_neighbor  (CBody*, CNode*);
bool   are_neighbors(CNode*, CNode*);
bool   are_neighbors(CBody*, CBody*);
bool   compare (CNode* one, CNode* two);
bool   bcompare (CBody* one, CBody* two);

CBody* identify_body(CNode*, CNode*);
CBody* identify_body(CNode*, CNode*, CNode*);
CNode* common_face  (CBody*, CBody*, CNode*);
CNode* common_face  (Vertex*, Vertex*, Vertex*);
bool   common_face  (Vertex*, Vertex*, Vertex*, Vertex*);
void   merge_faces  (CNode*, CNode*, CNode*, CNode*);
void   determine_neighbors(void);
void   determine_node_topologies();
void   cluster_evolver_out(int pid);

CNode* remove_edge(CNode*, CNode*);
void   remove_system_edge_nodes();
void   remove_small_bodies();
void   remove_small_faces ();
void   remove_small_edges ();
void   remove_node_from_face(CNode*, CNode*);

void   initialize_variables();
void   update_variables();
void   compactify();
void   compactify_bodies();

CNode* refine_edge(CNode*, CNode*);
void   refine_edges();

void   read_in_data(int argc, char *argv[]);
void   print_out_data(void);

void   evolver_out(int id);

int vector_search(std::vector<edge>& edges, CNode* one, CNode* two);

int g(int edge1[], int edge2[], int edge_id[], int one, int two, int size);
int f(int the_nodes[], int number, int size);

void crash();

double calc_echange(CNode*, CNode*);

void surface_evolver_system();

void fit_nodes_in_cube();

bool not_digon_neighbors(CNode* one, CNode* two);

void move_nodes();

#endif


