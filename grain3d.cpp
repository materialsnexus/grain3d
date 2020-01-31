//#define DEBUG

#include "cbody.h"
#include "cnode.h"
#include "functions.h"
#include "globals.h"

#include <cstdio>

#include <ostream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <iomanip>
#include <vector>
#include <algorithm>

/////////////////////////////////////////////////////////////////
// THE MAIN PROGRAM

int main(int argc, char *argv[])
{
  read_in_data(argc, argv);
  initialize_variables();
  determine_neighbors();
  determine_node_topologies();
  fit_nodes_in_cube();

  // TIME EVOLUTION UNTIL FEW GRAINS REMAIN
  int start_body = count_bodies();
  int num_bodies = start_body;

  while (num_bodies > end_body || step_count < min_steps)
  {

    update_variables();
    if (step_count % refine_edges_number == 0)
    {
      refine_edges();
    }

    if (step_count % remove_system_edge_nodes_number == 0)
    {
      remove_system_edge_nodes();
    }

    calc_motions();
    if (count_bodies() > end_body)
    {
      remove_small_edges();
      remove_small_faces();
      remove_small_bodies();
    }
    move_nodes();

    if (output_during_run)
    {
      if (step_count % output_whole_system_number == 0)
      {
        surface_evolver_system();
      }

      if (step_count % output_specific_neigh == 0)
      {
        for (int c = 1; c <= start_body; c++)
        {
          cluster_evolver_out(c);
        }
      }

      if (step_count % output_specific == 0)
      {
        for (int c = 1; c <= start_body; c++)
        {
          evolver_out(c);
        }
      }
    }

    total_time += time_scale;
    step_count++;
    num_bodies = count_bodies();

    // output key statistics about the system
    if (step_count % output_summary_number == 0)
    {
      compactify();
      calc_volumes();
      calc_areas();
      calc_avg_face_area();
      calc_and_print_stats();
    }
  }

  if (!output_during_run)
  {
    for (int c = 1; c <= start_body; c++)
    {
      evolver_out(c);
    }
    // for (int c = 0; c < count_bodies(); c++)
    // {
    //   evolver_out(c);
    // }
  }
}
// THE MAIN PROGRAM
/////////////////////////////////////////////////////////////////

void update_variables()
{
  time_scale = 1. / (double)count_bodies() / 2500.;
  smallest_edge = 0.015 * pow(1. / (double)count_bodies(), 1. / 3.); // C=0.15 gives on average 310 triangles per shape
}

void initialize_variables()
{
#pragma omp parallel
  nthreads = omp_get_num_threads();
  omp_set_num_threads(nthreads);
  offcounter1 = 0;
  offcounter2 = 0;
  offcounter3 = 0;
  toffcounter2 = 0;
  toffcounter3 = 0;

  time_scale = 1. / (double)count_bodies() / 2500.;
  smallest_edge = 0.015 * pow(1. / (double)count_bodies(), 1. / 3.); // this gives on average 310 triangles per shape
}
