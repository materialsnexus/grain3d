//#define DEBUG

#include "cbody.h"
#include "cnode.h"
#include "functions.h"
#include "globals.h"

#include <algorithm>

int step_count;    // step number
double time_scale; // size of time-step
double total_time; // total time of simulation

time_t start; // for timing purposes
clock_t cstart;
const double pi = 3.141592653589793;

unsigned int node_count; // number of nodes in system
unsigned int body_count; // number of bodies in system

double smallest_edge;
double smallest_face;
double smallest_body;
double avg_face_area;

int nthreads; // needed for multiprocessing

int offcounter1;
int offcounter2;
int offcounter3;

int toffcounter2;
int toffcounter3;

int tetrahedra_deleted;
int footballs_deleted;
int triangles_deleted;
int digons_deleted;
int edges_deleted;
int edge_nodes_deleted;
int edge_nodes_added;

std::ofstream outfile1;
std::ofstream outfile2;

std::vector<CNode *> nodes;  // our set of nodes
std::vector<CBody *> bodies; // our set of bodies

unsigned int end_body;  // final number of bodies in system
unsigned int min_steps; // minimum number of steps taken

bool output_during_run;                       // Should we output during the run or after
unsigned int refine_edges_number;             // frequency that the simulation refines edges
unsigned int remove_system_edge_nodes_number; // frequency that the simulation removes system edge nodes
unsigned int output_summary_number;           // frequency in which to output summary statistics

unsigned int output_whole_system_number; // frequency in which to output whole system
unsigned int output_specific_neigh;      // frequency in which to output specific grain and its neighbours
unsigned int output_specific;            // frequency in which to output specific grain