#include "cnode.h"
#include "cbody.h"

#include <ostream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>

#ifndef GVAR_H__
#define GVAR_H__

extern int step_count;    // step number
extern double time_scale; // size of time-step
extern double total_time; // total time of simulation

extern double avg_face_area;
extern int offcounter1;
extern int offcounter2;
extern int offcounter3;

extern int toffcounter2;
extern int toffcounter3;

extern time_t start;    // for timing purposes
extern clock_t cstart;  // for timing purposes
extern const double pi; // = 3.141592653589793;

extern unsigned int node_count; // number of nodes in system
extern unsigned int body_count; // number of bodies in system

extern double smallest_edge;
extern double smallest_face;
extern double smallest_body;

extern int nthreads; // needed for multiprocessing

extern int tetrahedra_deleted;
extern int footballs_deleted;
extern int triangles_deleted;
extern int digons_deleted;
extern int edges_deleted;
extern int edge_nodes_deleted;
extern int edge_nodes_added;

extern std::ofstream outfile1;
extern std::ofstream outfile2;

extern std::vector<CNode *> nodes;  // our set of nodes
extern std::vector<CBody *> bodies; // our set of bodies

extern unsigned int end_body;  // final number of bodies in system
extern unsigned int min_steps; // minimum number of steps taken

extern bool output_during_run;                       // Should we output during the run or after
extern unsigned int refine_edges_number;             // frequency that the simulation refines edges
extern unsigned int remove_system_edge_nodes_number; // frequency that the simulation removes system edge nodes
extern unsigned int output_summary_number;           // frequency in which to output summary statistics

extern unsigned int output_whole_system_number; // frequency in which to output whole system
extern unsigned int output_specific_neigh;      // frequency in which to output specific grain and its neighbours
extern unsigned int output_specific;            // frequency in which to output specific grain

#endif
