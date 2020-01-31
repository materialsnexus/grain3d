#include "cbody.h"
#include "cnode.h"
#include "functions.h"
#include "globals.h"

#include <algorithm>
#include <ostream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <iomanip>
#include <vector>

// #define DEBUG

void CBody::remove_tetrahedron(void)
{
#ifdef DEBUG
  std::cout << "Removing tetrahedron " << id << " on step " << step_count << '\t';
#endif

  std::vector<CNode *> corner_nodes;
  std::vector<CNode *> face_nodes;
  std::vector<CNode *> all_neighbors;

  for (unsigned int c = 1; c < triplets.size(); c++)
    face_nodes.push_back(triplets[c].v2);
  std::sort(face_nodes.begin(), face_nodes.end(), compare);
  face_nodes.erase(unique(face_nodes.begin(), face_nodes.end()), face_nodes.end());
  for (unsigned int d = 0; d < face_nodes.size(); d++)
    face_nodes[d]->remove_edge_nodes();

  for (unsigned int c = 1; c < triplets.size(); c++)
    corner_nodes.push_back(triplets[c].v1);
  std::sort(corner_nodes.begin(), corner_nodes.end(), compare);
  corner_nodes.erase(unique(corner_nodes.begin(), corner_nodes.end()), corner_nodes.end());

  for (unsigned int i = 0; i < 4; i++)     // four corner nodes
    for (unsigned int c = 0; c < 4; c++)   // each has four bodies
      for (unsigned int d = 0; d < 6; d++) // each has six things
        if (corner_nodes[i]->corners[c][d] != corner_nodes[0] &&
            corner_nodes[i]->corners[c][d] != corner_nodes[1] &&
            corner_nodes[i]->corners[c][d] != corner_nodes[2] &&
            corner_nodes[i]->corners[c][d] != corner_nodes[3] &&
            corner_nodes[i]->corners[c][d] != face_nodes[0] &&
            corner_nodes[i]->corners[c][d] != face_nodes[1] &&
            corner_nodes[i]->corners[c][d] != face_nodes[2] &&
            corner_nodes[i]->corners[c][d] != face_nodes[3])
          all_neighbors.push_back(corner_nodes[i]->corners[c][d]);

  std::sort(all_neighbors.begin(), all_neighbors.end(), compare);
  all_neighbors.erase(unique(all_neighbors.begin(), all_neighbors.end()), all_neighbors.end());

  CBody *adjacent_bodies[4];
  for (unsigned int c = 0; c < 4; c++)
  {
    face_nodes[c]->center2(1.0);
    if (face_nodes[c]->neighbors[0] == this)
      adjacent_bodies[c] = face_nodes[c]->neighbors[1];
    else
      adjacent_bodies[c] = face_nodes[c]->neighbors[0];
  }

  // Make new nodes, place them appropriately

  nodes.resize(nodes.size() + 1);
  nodes[node_count + 1] = new CNode;

  *nodes[node_count + 1] = *corner_nodes[0] + (*corner_nodes[1] - *corner_nodes[0]) / 4 + (*corner_nodes[2] - *corner_nodes[0]) / 4 + (*corner_nodes[3] - *corner_nodes[0]) / 4;
  node_count++;
  nodes[node_count]->id = node_count;

  for (unsigned int c = 0; c < 4; c++)
    nodes[node_count]->neighbors.push_back(adjacent_bodies[c]);

  for (unsigned int c = 0; c < 4; c++)
  {
    adjacent_bodies[c]->simple_remove_face(face_nodes[c]);
    adjacent_bodies[c]->replace_edges_with_node(face_nodes[c]->corners[0][0], face_nodes[c]->corners[0][1], face_nodes[c]->corners[0][2], nodes[node_count]);
  }

#ifdef DEBUG
  std::vector<int> numbers;
  for (unsigned int c = 0; c < 4; c++)
    numbers.push_back(adjacent_bodies[c]->id);
  numbers.push_back(id);
  std::sort(numbers.begin(), numbers.end());

  std::cout << "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t";
  for (unsigned int c = 0; c < 5; c++)
    std::cout << std::setw(3) << numbers[c] << "  ";
  std::cout << 'T' << std::endl;
#endif

  // Get rid of old stuff
  blank();

  for (unsigned int c = 0; c < 4; c++)
  {
    corner_nodes[c]->neighbors.clear();
    face_nodes[c]->neighbors.clear();
  }

  // Recompute topological information
  all_neighbors.push_back(nodes[node_count]);
  for (unsigned int c = 0; c < all_neighbors.size(); c++)
    all_neighbors[c]->determine_node_topology();
  for (unsigned int c = 0; c < all_neighbors.size(); c++)
    all_neighbors[c]->calc_motion();

  // Recompute volumes

  for (unsigned int c = 0; c < 4; c++)
    adjacent_bodies[c]->calc_volume();

  // Update system variable

  tetrahedra_deleted++;
}

void CBody::remove_football(void)
{
#ifdef DEBUG
  std::cout << "Removing football " << id << " on step " << step_count << "\t";
#endif

  std::vector<CNode *> all_neighbors;
  std::vector<CNode *> face_nodes;

#ifdef DEBUG
  std::vector<int> numbers;
  for (unsigned int c = 0; c < 4; c++)
  {
    numbers.push_back(triplets[1].v1->neighbors[c]->id);
    numbers.push_back(triplets[1].v3->neighbors[c]->id);
  }
  std::sort(numbers.begin(), numbers.end());
  numbers.erase(unique(numbers.begin(), numbers.end()), numbers.end());

  std::cout << "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t     ";
  for (unsigned int c = 0; c < 4; c++)
    std::cout << std::setw(3) << numbers[c] << "  ";
  std::cout << 'F' << std::endl;
#endif

  for (unsigned int c = 1; c < triplets.size(); c++)
    face_nodes.push_back(triplets[c].v2);
  std::sort(face_nodes.begin(), face_nodes.end(), compare);
  face_nodes.erase(unique(face_nodes.begin(), face_nodes.end()), face_nodes.end());
  for (unsigned int d = 0; d < face_nodes.size(); d++)
    face_nodes[d]->remove_edge_nodes();

  CBody *adjacent_bodies[3];
  for (unsigned int c = 0; c < 3; c++)
  {
    face_nodes[c]->center2(1.0);
    if (face_nodes[c]->neighbors[0] == this)
      adjacent_bodies[c] = face_nodes[c]->neighbors[1];
    else
      adjacent_bodies[c] = face_nodes[c]->neighbors[0];
  }

  for (unsigned int c = 0; c < 4; c++)
    for (unsigned int d = 0; d < 6; d++)
    {
      all_neighbors.push_back(triplets[1].v1->corners[c][d]);
      all_neighbors.push_back(triplets[1].v3->corners[c][d]);
    }
  std::sort(all_neighbors.begin(), all_neighbors.end(), compare);
  all_neighbors.erase(unique(all_neighbors.begin(), all_neighbors.end()), all_neighbors.end());

  // Remove faces from adjacent bodies (very messy, maybe could be reworked)
  for (unsigned int c = 0; c < 3; c++)
  {
    for (unsigned int d = 1; d < adjacent_bodies[c]->triplets.size(); d++)
    {
      if ((adjacent_bodies[c]->triplets[d].v1 == triplets[1].v1 && adjacent_bodies[c]->triplets[d].v3 == triplets[1].v3) ||
          (adjacent_bodies[c]->triplets[d].v1 == triplets[1].v3 && adjacent_bodies[c]->triplets[d].v3 == triplets[1].v1))
      {
        adjacent_bodies[c]->triplets.erase(adjacent_bodies[c]->triplets.begin() + d);
        d--;
      }
    }

    int p1 = 0, p2 = 0, p3 = 0, p4 = 0;
    for (unsigned int d = 1; d < adjacent_bodies[c]->triplets.size(); d++)
    {
      if (adjacent_bodies[c]->triplets[d].v3 == triplets[1].v1)
        p1 = d;
      if (adjacent_bodies[c]->triplets[d].v1 == triplets[1].v3)
        p2 = d;
      if (adjacent_bodies[c]->triplets[d].v3 == triplets[1].v3)
        p3 = d;
      if (adjacent_bodies[c]->triplets[d].v1 == triplets[1].v1)
        p4 = d;
    }

    adjacent_bodies[c]->triplets[p1].v3 = adjacent_bodies[c]->triplets[p2].v3;
    adjacent_bodies[c]->triplets[p3].v3 = adjacent_bodies[c]->triplets[p4].v3;

    adjacent_bodies[c]->triplets.erase(adjacent_bodies[c]->triplets.begin() + std::max(p2, p4));
    adjacent_bodies[c]->triplets.erase(adjacent_bodies[c]->triplets.begin() + std::min(p2, p4));
  }

  // Get rid of old stuff

  triplets[1].v1->neighbors.clear();
  triplets[1].v3->neighbors.clear();
  for (unsigned int c = 0; c < 3; c++)
    face_nodes[c]->neighbors.clear();
  blank();

  // Recompute topological information
  for (unsigned int c = 0; c < all_neighbors.size(); c++)
    all_neighbors[c]->determine_node_topology();
  for (unsigned int c = 0; c < all_neighbors.size(); c++)
    all_neighbors[c]->calc_motion();

  // Recompute volumes
  for (unsigned int c = 0; c < 3; c++)
    adjacent_bodies[c]->calc_volume();

  // Update system variable
  footballs_deleted++;
}

void CNode::remove_tri_face(void)
{
#ifdef DEBUG
  std::cout << "Removing triangular face " << id << '\n';
#endif
  remove_edge_nodes();

  std::vector<CNode *> up_nodes;
  std::vector<CNode *> down_nodes;
  std::vector<CNode *> all_neighbors;

  for (unsigned int c = 0; c < 4; c++)
    for (unsigned int d = 0; d < 6; d++)
    {
      if (corners[0][0]->corners[c][d] != corners[0][1] && corners[0][0]->corners[c][d] != corners[0][2] && corners[0][0]->corners[c][d] != this)
        all_neighbors.push_back(corners[0][0]->corners[c][d]);
      if (corners[0][1]->corners[c][d] != corners[0][2] && corners[0][1]->corners[c][d] != corners[0][0] && corners[0][1]->corners[c][d] != this)
        all_neighbors.push_back(corners[0][1]->corners[c][d]);
      if (corners[0][2]->corners[c][d] != corners[0][0] && corners[0][2]->corners[c][d] != corners[0][1] && corners[0][2]->corners[c][d] != this)
        all_neighbors.push_back(corners[0][2]->corners[c][d]);
    }
  std::sort(all_neighbors.begin(), all_neighbors.end(), compare);
  all_neighbors.erase(unique(all_neighbors.begin(), all_neighbors.end()), all_neighbors.end());

  double lengths = length(corners[0][0], corners[0][1]) + length(corners[0][1], corners[0][2]) + length(corners[0][2], corners[0][0]);
  lengths *= 0.1;

  CBody *bodyA = identify_body(corners[0][0], corners[0][1]);
  CBody *bodyB = identify_body(corners[0][1], corners[0][2]);
  CBody *bodyC = identify_body(corners[0][2], corners[0][0]);
  CBody *bodyU = identify_body(corners[0][0], this, corners[0][1]);
  CBody *bodyD = identify_body(corners[0][1], this, corners[0][0]);

  double iupv = bodyU->get_volume();
  double idnv = bodyD->get_volume();

  // Make new nodes, place them appropriately
  for (unsigned int c = 1; c < bodyU->triplets.size(); c++)
    if ((bodyU->triplets[c].v3 == corners[0][0] || bodyU->triplets[c].v3 == corners[0][1] || bodyU->triplets[c].v3 == corners[0][2]) && bodyU->triplets[c].v2 != this && (bodyU->triplets[c].v1 != corners[0][0] && bodyU->triplets[c].v1 != corners[0][1] && bodyU->triplets[c].v1 != corners[0][2]))
      up_nodes.push_back(bodyU->triplets[c].v1);

  for (unsigned int c = 1; c < bodyD->triplets.size(); c++)
    if ((bodyD->triplets[c].v3 == corners[0][0] || bodyD->triplets[c].v3 == corners[0][1] || bodyD->triplets[c].v3 == corners[0][2]) && bodyD->triplets[c].v2 != this && (bodyD->triplets[c].v1 != corners[0][0] && bodyD->triplets[c].v1 != corners[0][1] && bodyD->triplets[c].v1 != corners[0][2]))
      down_nodes.push_back(bodyD->triplets[c].v1);

  std::sort(up_nodes.begin(), up_nodes.end(), compare);
  std::sort(down_nodes.begin(), down_nodes.end(), compare);

  up_nodes.erase(unique(up_nodes.begin(), up_nodes.end()), up_nodes.end());
  down_nodes.erase(unique(down_nodes.begin(), down_nodes.end()), down_nodes.end());

#ifdef DEBUG
  std::cout << "Removing face between bodies " << neighbors[0]->id << " and " << neighbors[1]->id << '\n';
  std::cout << "Removing face " << std::setw(5) << id << " on step " << std::setw(5) << step_count << "; ";
  std::cout << "with nodes " << std::setw(5) << corners[0][0]->id << ", " << std::setw(5) << corners[0][1]->id << ", " << std::setw(5) << corners[0][2]->id << "; ";
#endif

  nodes.resize(nodes.size() + 2);
  nodes[node_count + 1] = new CNode();
  nodes[node_count + 2] = new CNode();

  CNode *nodeU = nodes[node_count + 1];
  CNode *nodeD = nodes[node_count + 2];

  nodeD->x = x;
  nodeD->y = y;
  nodeD->z = z;
  nodeU->x = x;
  nodeU->y = y;
  nodeU->z = z;

  for (unsigned int c = 0; c < 3; c++)
  {
    nodeU->x += lengths * (*up_nodes[c] - *this).x / (*up_nodes[c] - *this).norm();
    nodeU->y += lengths * (*up_nodes[c] - *this).y / (*up_nodes[c] - *this).norm();
    nodeU->z += lengths * (*up_nodes[c] - *this).z / (*up_nodes[c] - *this).norm();
    nodeD->x += lengths * (*down_nodes[c] - *this).x / (*down_nodes[c] - *this).norm();
    nodeD->y += lengths * (*down_nodes[c] - *this).y / (*down_nodes[c] - *this).norm();
    nodeD->z += lengths * (*down_nodes[c] - *this).z / (*down_nodes[c] - *this).norm();
  }

  center();
  *this = *corners[0][0] + (*corners[0][1] - *corners[0][0]) / 3. + (*corners[0][2] - *corners[0][0]) / 3.;

  nodeU->neighbors.push_back(bodyU);
  nodeD->neighbors.push_back(bodyD);
  nodeU->neighbors.push_back(bodyA);
  nodeU->neighbors.push_back(bodyB);
  nodeU->neighbors.push_back(bodyC);
  nodeD->neighbors.push_back(bodyA);
  nodeD->neighbors.push_back(bodyB);
  nodeD->neighbors.push_back(bodyC);

  node_count++;
  nodeU->id = node_count;
  node_count++;
  nodeD->id = node_count;

#ifdef DEBUG
  std::vector<int> numbers;
  numbers.push_back(bodyA->id);
  numbers.push_back(bodyB->id);
  numbers.push_back(bodyC->id);
  numbers.push_back(bodyU->id);
  numbers.push_back(bodyD->id);
  std::sort(numbers.begin(), numbers.end());

  std::cout << "new nodes " << std::setw(5) << nodeU->id << " and " << std::setw(5) << nodeD->id << ";\t\t\t\t\t\t\t\t\t";
  for (unsigned int c = 0; c < 5; c++)
    std::cout << std::setw(3) << numbers[c] << "  ";
  std::cout << std::endl;
#endif

  bodyU->simple_remove_face(this);
  bodyD->simple_remove_face(this);

  bodyU->replace_edges_with_node(corners[0][1], corners[0][0], corners[0][2], nodeU);
  bodyD->replace_edges_with_node(corners[0][0], corners[0][1], corners[0][2], nodeD);

  bodyA->switch_edge_with_edge(corners[0][1], corners[0][2], nodeU, nodeD);
  bodyB->switch_edge_with_edge(corners[0][2], corners[0][0], nodeU, nodeD);
  bodyC->switch_edge_with_edge(corners[0][0], corners[0][1], nodeU, nodeD);

  // Get rid of old stuff

  neighbors.clear();
  corners[0][0]->neighbors.clear();
  corners[0][1]->neighbors.clear();
  corners[0][2]->neighbors.clear();

  // Recompute topological information
  all_neighbors.push_back(nodeU);
  all_neighbors.push_back(nodeD);
  for (unsigned int c = 0; c < all_neighbors.size(); c++)
    all_neighbors[c]->determine_node_topology();

  if ((bodyU->get_volume() - iupv) / iupv < -0.2)
  {
#ifdef DEBUG
    std::cout << "Body " << bodyU->id << " changed a lot when removing triangle face, so fixing faces\n";
#endif
    std::vector<CNode *> uface_nodes;
    for (unsigned int c = 1; c < bodyU->triplets.size(); c++)
      uface_nodes.push_back(bodyU->triplets[c].v2);
    sort(uface_nodes.begin(), uface_nodes.end(), compare);
    uface_nodes.erase(unique(uface_nodes.begin(), uface_nodes.end()), uface_nodes.end());
    for (unsigned int c = 0; c < uface_nodes.size(); c++)
      uface_nodes[c]->center2(1.);
  }
#ifdef DEBUG
  else
    std::cout << "We didn't change body " << bodyU->id << " that much, just " << (bodyU->get_volume() - iupv) / iupv << '\n';
#endif

  if ((bodyD->get_volume() - idnv) / idnv < -0.2)
  {
#ifdef DEBUG
    std::cout << "Body " << bodyD->id << " changed alot when removing triangle face, so fixing faces\n";
#endif
    std::vector<CNode *> dface_nodes;
    for (unsigned int c = 1; c < bodyD->triplets.size(); c++)
      dface_nodes.push_back(bodyD->triplets[c].v2);
    std::sort(dface_nodes.begin(), dface_nodes.end(), compare);
    dface_nodes.erase(unique(dface_nodes.begin(), dface_nodes.end()), dface_nodes.end());
    for (unsigned int c = 0; c < dface_nodes.size(); c++)
      dface_nodes[c]->center2(1.);
  }
#ifdef DEBUG
  else
    std::cout << "We didn't change body " << bodyD->id << " that much, just " << (bodyD->get_volume() - idnv) / idnv << '\n';
#endif

#ifdef DEBUG
  std::cout << '\n'
            << '\n';
#endif

  std::vector<CNode *> facesn;
  for (unsigned int c = 0; c < all_neighbors.size(); c++)
    if (all_neighbors[c]->neighbors.size() == 2)
      facesn.push_back(all_neighbors[c]);

  for (unsigned int c = 0; c < all_neighbors.size(); c++)
    all_neighbors[c]->calc_motion();

  // Recompute volumes
  bodyU->calc_volume();
  bodyD->calc_volume();
  bodyA->calc_volume();
  bodyB->calc_volume();
  bodyC->calc_volume();
}

void CNode::remove_digon(void)
{
#ifdef DEBUG
  std::cout << "Removing digon " << id << " on step " << step_count << " with neighbors " << neighbors[0]->id << " and " << neighbors[1]->id << '\t';
#endif

  remove_edge_nodes(); // moved this down just recently, not sure if good idea or not

  // Identify important nodes and bodies
  CNode *nodeA = corners[0][0];
  CNode *nodeB = corners[0][1]; // this was changed from neighbors[0]->triplets[last ].v1

  std::vector<int> numbers;
  for (unsigned int c = 0; c < 4; c++)
    numbers.push_back(nodeA->neighbors[c]->id);
  std::sort(numbers.begin(), numbers.end());
  numbers.erase(unique(numbers.begin(), numbers.end()), numbers.end());

#ifdef DEBUG
  std::cout << "\t\t\t\t\t\t\t\t\t\t\t\t     ";
  for (unsigned int c = 0; c < 4; c++)
    std::cout << std::setw(3) << numbers[c] << "  ";
  std::cout << 'D' << std::endl;
#endif

  std::vector<CNode *> neighborsA;
  for (unsigned int c = 0; c < nodeA->neighbors.size(); c++)
    for (unsigned int d = 0; d < nodeA->corners[c].size(); d++)
      if (nodeA->corners[c][d] != nodeB && nodeA->corners[c][d]->neighbors.size() != 2)
        neighborsA.push_back(nodeA->corners[c][d]);

  std::vector<CNode *> neighborsB;
  for (unsigned int c = 0; c < nodeB->neighbors.size(); c++)
    for (unsigned int d = 0; d < nodeB->corners[c].size(); d++)
      if (nodeB->corners[c][d] != nodeA && nodeB->corners[c][d]->neighbors.size() != 2)
        neighborsB.push_back(nodeB->corners[c][d]);

  std::sort(neighborsA.begin(), neighborsA.end(), compare);
  std::sort(neighborsB.begin(), neighborsB.end(), compare);

  neighborsA.erase(unique(neighborsA.begin(), neighborsA.end()), neighborsA.end());
  neighborsB.erase(unique(neighborsB.begin(), neighborsB.end()), neighborsB.end());

  std::vector<CBody *> lbodies;
  for (unsigned int c = 0; c < nodeA->neighbors.size(); c++)
    if (nodeA->neighbors[c] != neighbors[0] && nodeA->neighbors[c] != neighbors[1])
      lbodies.push_back(nodeA->neighbors[c]);

  CBody *up = neighbors[0];
  CBody *down = neighbors[1];
  CBody *back = lbodies[0];
  CBody *front = lbodies[1];

  up->calc_volume();
  down->calc_volume();
  back->calc_volume();
  front->calc_volume();

  CNode *nodeH = common_face(up, front, nodeA);
  CNode *nodeI = common_face(down, front, nodeA);
  CNode *nodeJ = common_face(back, front, nodeA);
  CNode *nodeK = common_face(back, front, nodeB);
  CNode *nodeL = common_face(up, back, nodeA);
  CNode *nodeM = common_face(down, back, nodeA);

  // Make new nodes, place them appropriately

  up->simple_remove_face(this);
  down->simple_remove_face(this);

  remove_node_from_face(nodeA, nodeH);
  remove_node_from_face(nodeB, nodeH);
  remove_node_from_face(nodeA, nodeI);
  remove_node_from_face(nodeB, nodeI);
  remove_node_from_face(nodeA, nodeL);
  remove_node_from_face(nodeB, nodeL);
  remove_node_from_face(nodeA, nodeM);
  remove_node_from_face(nodeB, nodeM);

  merge_faces(nodeJ, nodeK, nodeA, nodeB);
  if (nodeK->id < nodeJ->id) // not sure we ever needed to do this stuff....
  {                          // also in the merge_faces routine....
    CNode *temp = nodeJ;
    nodeJ = nodeK;
    nodeK = temp;
  }

  // Get rid of old stuff

  neighbors.clear();

  // THESE ARE PARTS THAT NEED TO BE FIXED UP.  I HAD TROUBLE LABELING SPECIFICALLY THOSE NODES
  // THAT NEED TO BE REEVALUATED.  THAT IS, ONLY CERTAIN NODES HAVE BEEN AFFECTED, IN TERMS OF
  // TOPOLOGY, IN TERMS OF MOTIONS, IN TERMS OFF ALL OF THIS.  WE NEED TO FIGURE THAT OUT, SO
  // WE DON'T LOSE TIME LOOKING THROUGH ENTIRE BODIES OF NODES.

  up->determine_node_topologies();
  down->determine_node_topologies();
  back->determine_node_topologies();
  front->determine_node_topologies();

  // Recompute topological information
  nodeJ->center();

  for (unsigned int c = 0; c < nodeH->corners[0].size(); c++)
    nodeH->corners[0][c]->center();
  for (unsigned int c = 0; c < nodeI->corners[0].size(); c++)
    nodeI->corners[0][c]->center();
  for (unsigned int c = 0; c < nodeJ->corners[0].size(); c++)
    nodeJ->corners[0][c]->center();
  for (unsigned int c = 0; c < nodeL->corners[0].size(); c++)
    nodeL->corners[0][c]->center();
  for (unsigned int c = 0; c < nodeM->corners[0].size(); c++)
    nodeM->corners[0][c]->center();

  nodeH->center();
  nodeI->center();
  nodeJ->center();
  nodeL->center();
  nodeM->center();

#pragma omp parallel
  for (unsigned int c = omp_get_thread_num() + 1; c < nodes.size(); c += nthreads)
    nodes[c]->calc_motion();

  // Recompute volumes
  up->calc_volume();
  down->calc_volume();
  back->calc_volume();
  front->calc_volume();
}

bool CNode::remove_face()
{
  while (get_face_sides() > 3)
  {
    int index = 999999;
    double min_length = 1;

#ifdef DEBUG
    std::cout << "We have too many sides for the face to erase, so we will try erasing an edge\n";
#endif

    remove_edge_nodes();
    center2(1.0);
    center();

    // HERE WE ARE LOOKING FOR AN EDGE TO REMOVE FROM THIS FACE.
    for (unsigned int c = 0; c < corners[0].size(); c++)
    {
      double current_length = length(corners[0][c], corners[0][(c + 1) % corners[0].size()]);
      if (current_length < min_length &&
          (corners[0][c]->neighbors[0] != corners[0][(c + 1) % corners[0].size()]->neighbors[0] || // THIS CONDITION ENSURES THAT OUR EDGE THAT WE'RE CONSIDERING
           corners[0][c]->neighbors[1] != corners[0][(c + 1) % corners[0].size()]->neighbors[1] || // IS NOT PART OF A DIGON.  WE MAKE SURE THAT AT LEAST ONE OF
           corners[0][c]->neighbors[2] != corners[0][(c + 1) % corners[0].size()]->neighbors[2] || // THEIR NEIGHBORS IS DISTINCT
           corners[0][c]->neighbors[3] != corners[0][(c + 1) % corners[0].size()]->neighbors[3]))
      {
        std::vector<CNode *> top_nodes;
        for (unsigned int p = 0; p < 4; p++)
          for (unsigned int q = 1; q < 6; q += 2)
            top_nodes.push_back(corners[0][c]->corners[p][q]);
        std::sort(top_nodes.begin(), top_nodes.end(), compare);
        top_nodes.erase(unique(top_nodes.begin(), top_nodes.end()), top_nodes.end());

        std::vector<CNode *> bottom_nodes;
        for (unsigned int p = 0; p < 4; p++)
          for (unsigned int q = 1; q < 6; q += 2)
            bottom_nodes.push_back(corners[0][(c + 1) % corners[0].size()]->corners[p][q]);
        std::sort(bottom_nodes.begin(), bottom_nodes.end(), compare);
        bottom_nodes.erase(unique(bottom_nodes.begin(), bottom_nodes.end()), bottom_nodes.end());

        // THE LAST TWELVE LINES AND THIS CONDITION PREVENT US FROM REMOVING AN EDGE WHERE ONE OF THE VERTICES
        // ARE PART OF A DIGON.  THIS WILL CAUSE QUITE THE TROUBLE IN REMOVING THIS FACE.  IN THEORY I THINK THAT
        // THIS CHECK ISN'T EVEN COMPLETELY SAFE.  I THINK EVEN HAVING ONE VERTEX AS PART OF A DIGON COULD BE A
        // PROBLEM.  BUT FIXING THIS PROBLEM WOULD FORCE US TO PREVENT SOME FACES FROM BEING ERASED ALTOGETHER,
        // WHICH SOMETIMES COULD BE A PROBLEM.  MAYBE WE COULD ALLOW THIS, AND THEN HAVE A RETURN ON THIS FUNCTION,
        // WITH A ZERO FOR FAILURE TO ERASE A FACE.  IN REMOVING SMALL FACES THIS WOULD BE OK, BUT WHEN TRYING TO
        // TO REMOVE A CELL, FOR EXAMPLE, AND WE NEED TO ERASE A FACE FIRST, THEN THIS MIGHT BE A PROBLEM, BUT WE
        // COULD USE THE 0 RETURN TO LOOK FOR ANOTHER FACE TO DELETE.
        if (top_nodes.size() == 4 || bottom_nodes.size() == 4) // THIS CONDITION SEEMS TO REQUIRE THAT AT LEAST ONE OF THE
        {                                                      // TWO NODES OF THIS EDGE ARE GOOD VERTEX NODES (NOT PART
          index = c;                                           // OF A DIGON.
          min_length = current_length;
        }
        else
        {
          std::cout << "We have an edge we can't erase...\n";
          return 0;
        }
      }
    }

#ifdef DEBUG
    std::cout << "Removing edge in an attempt to remove a face " << id << "\n";
#endif

    if (index == 999999)
    {
      std::cout << "We have a problem.  We are trying to remove this face " << id << " but none of its edges are suitable.\n";
      return 0; // WE WILL NOW TEST IF WE CAN AVOID TRYING TO ERASE THIS FACE
      output();
      neighbors[0]->evolver_out();
      neighbors[1]->evolver_out();
      std::cout << "crashing in topchanges: remove_face";
      crash();
    }

#ifdef DEBUG
    std::cout << "We have decided to remove edge " << corners[0][index]->id << "-" << corners[0][(index + 1) % corners[0].size()]->id << '\n';
#endif

    remove_edge(corners[0][index], corners[0][(index + 1) % corners[0].size()]);
  }

#ifdef DEBUG
  std::cout << "K, ready to remove 2- or 3-face " << id << "\n";
#endif

  if (get_face_sides() == 2)
    remove_digon();
  if (get_face_sides() == 3)
    remove_tri_face();

  return 1;
}

void remove_node_from_face(CNode *node, CNode *face)
{
  CBody *one = face->neighbors[0];
  CBody *two = face->neighbors[1];

  int first = 0, second = 0;
  for (unsigned int c = 1; c < one->triplets.size(); c++)
  {
    if (one->triplets[c].v3 == node && one->triplets[c].v2 == face)
      first = c;
    if (one->triplets[c].v1 == node && one->triplets[c].v2 == face)
      second = c;
  }

  one->triplets[first].v3 = one->triplets[second].v3;
  one->triplets.erase(one->triplets.begin() + second);

  first = 0;
  second = 0;
  for (unsigned int c = 1; c < two->triplets.size(); c++)
  {
    if (two->triplets[c].v3 == node && two->triplets[c].v2 == face)
      first = c;
    if (two->triplets[c].v1 == node && two->triplets[c].v2 == face)
      second = c;
  }

  two->triplets[first].v3 = two->triplets[second].v3;
  two->triplets.erase(two->triplets.begin() + second);
}

void merge_faces(CNode *one, CNode *two, CNode *A, CNode *B)
{
#ifdef DEBUG
  std::cout << "Merging faces " << one->id << " and " << two->id << " on step " << step_count << std::endl;
#endif

  if (two->id < one->id)
  {
    CNode *temp = one;
    one = two;
    two = temp;
  }

  if (one->neighbors[0] != two->neighbors[0])
  {
    CBody *temp = two->neighbors[0];
    two->neighbors[0] = two->neighbors[1];
    two->neighbors[1] = temp;
  }

  CBody *nbodies[2];
  nbodies[0] = one->neighbors[0];
  nbodies[1] = one->neighbors[1];

  for (unsigned int d = 0; d < 2; d++)
  {
    int p1 = 0, p2 = 0, p3 = 0, p4 = 0;
    for (unsigned int c = 1; c < nbodies[d]->triplets.size(); c++)
    {
      if (nbodies[d]->triplets[c].v2 == one && p1 == 0)
        p1 = c;
      if (nbodies[d]->triplets[c].v2 == one)
        p2 = c;
      if (nbodies[d]->triplets[c].v2 == two && p3 == 0)
        p3 = c;
      if (nbodies[d]->triplets[c].v2 == two)
        p4 = c;
    }

    std::vector<triplet> temp(nbodies[d]->triplets.begin() + p3, nbodies[d]->triplets.begin() + p4 + 1);
    temp.insert(temp.begin(), nbodies[d]->triplets.begin() + p1, nbodies[d]->triplets.begin() + p2 + 1);

    p1 = 0, p2 = 0, p3 = 0, p4 = 0;
    for (unsigned int c = 0; c < temp.size(); c++)
    {
      if (temp[c].v3 == A)
        p1 = c;
      if (temp[c].v1 == A)
        p2 = c;
      if (temp[c].v3 == B)
        p3 = c;
      if (temp[c].v1 == B)
        p4 = c;
    }

    temp[p1].v3 = temp[p4].v3;
    temp[p2].v1 = temp[p3].v1;

    temp.erase(temp.begin() + std::max(p3, p4));
    temp.erase(temp.begin() + std::min(p3, p4));

    for (unsigned int c = 0; c < temp.size(); c++)
      if (temp[c].v2 == two)
        temp[c].v2 = one;

    CNode *tempv1;
    CNode *tempv3;
    unsigned int number = 1;
    for (unsigned int c = number; c < temp.size(); c++)
    {
      if (temp[c].v1 == temp[number - 1].v3)
      {
        if (c == number)
          number++;
        else
        {
          tempv1 = temp[c].v1;
          tempv3 = temp[c].v3;
          temp[c].v1 = temp[number].v1;
          temp[c].v3 = temp[number].v3;
          temp[number].v1 = tempv1;
          temp[number].v3 = tempv3;
          number++;
          c = number;
        }
      }
    }

    for (unsigned int k = 1; k < nbodies[d]->triplets.size(); k++)
      if (nbodies[d]->triplets[k].v2 == one || nbodies[d]->triplets[k].v2 == two)
      {
        nbodies[d]->triplets.erase(nbodies[d]->triplets.begin() + k);
        k--;
      }

    int place = 0;
    for (unsigned int k = 1; k < nbodies[d]->triplets.size(); k++)
      if (nbodies[d]->triplets[k].v2->id < one->id)
        place = k;

    nbodies[d]->triplets.insert(nbodies[d]->triplets.begin() + place + 1, temp.begin(), temp.end());
  }

  int one_id = one->id;
  *one = *A + (*B - *A) / 2;
  one->id = one_id;
  one->neighbors.push_back(nbodies[0]);
  one->neighbors.push_back(nbodies[1]);

  two->neighbors.clear();
  A->neighbors.clear();
  B->neighbors.clear();
}

void CBody::switch_edge_with_edge(CNode *e1, CNode *e2, CNode *e3, CNode *e4)
{
  CNode *top = 0;
  CNode *bottom = 0;
  for (unsigned int c = 1; c < triplets.size(); c++)
  {
    if (triplets[c].v1 == e1 && triplets[c].v3 == e2)
      top = triplets[c].v2;
    if (triplets[c].v1 == e2 && triplets[c].v3 == e1)
      bottom = triplets[c].v2;
  }

  for (unsigned int c = 1; c < triplets.size(); c++)
  {
    if (triplets[c].v3 == e1 && triplets[c].v2 == top)
      triplets[c].v3 = e3;
    if (triplets[c].v1 == e2 && triplets[c].v2 == top)
      triplets[c].v1 = e3;
    if (triplets[c].v1 == e1 && triplets[c].v3 == e2)
    {
      triplets.erase(triplets.begin() + c);
      c--;
    }

    if (triplets[c].v3 == e2 && triplets[c].v2 == bottom)
      triplets[c].v3 = e4;
    if (triplets[c].v1 == e1 && triplets[c].v2 == bottom)
      triplets[c].v1 = e4;
    if (triplets[c].v1 == e2 && triplets[c].v3 == e1)
    {
      triplets.erase(triplets.begin() + c);
      c--;
    }

    if (triplets[c].v1 == e1 && triplets[c].v2 != top && triplets[c].v2 != bottom)
      triplets[c].v1 = e3;
    if (triplets[c].v3 == e1 && triplets[c].v2 != top && triplets[c].v2 != bottom)
    {
      triplets[c].v3 = e4;
      triplets.insert(triplets.begin() + c + 1, triplet(e4, triplets[c].v2, e3));
    }

    if (triplets[c].v1 == e2 && triplets[c].v2 != top && triplets[c].v2 != bottom)
      triplets[c].v1 = e4;
    if (triplets[c].v3 == e2 && triplets[c].v2 != top && triplets[c].v2 != bottom)
    {
      triplets[c].v3 = e3;
      triplets.insert(triplets.begin() + c + 1, triplet(e3, triplets[c].v2, e4));
    }
  }
}

void CBody::replace_node_on_face_with_edge(CNode *face_node, CNode *old_node, CNode *e1, CNode *e2)
{
  for (unsigned int c = 1; c < triplets.size(); c++)
  {
    if (triplets[c].v2 == face_node && triplets[c].v1 == old_node)
      triplets[c].v1 = e2;
    if (triplets[c].v2 == face_node && triplets[c].v3 == old_node)
    {
      triplets[c].v3 = e1;
      triplets.insert(triplets.begin() + c + 1, triplet(e1, face_node, e2));
    }
  }
}

CNode *refine_edge(CNode *one, CNode *two)
{
  nodes.resize(nodes.size() + 1);
  nodes[node_count + 1] = new CNode();
  *nodes[node_count + 1] = *one + (*two - *one) / 2;
  nodes[node_count + 1]->dx = 0;
  nodes[node_count + 1]->dy = 0;
  nodes[node_count + 1]->dz = 0;

  node_count++;
  nodes[node_count]->id = node_count;

#ifdef DEBUG
  //  cout << "Adding new edge node " << node_count << " on edge " << one->id << "-" << two->id << " on step " << step_count << endl;
#endif

  std::vector<CBody *> side_bodies;
  for (unsigned int c = 0; c < one->neighbors.size(); c++)
  {
    if (is_neighbor(one->neighbors[c], two) == 1)
    {
      side_bodies.push_back(one->neighbors[c]);
    }
  }

  if (side_bodies.size() == 4)
  {
#ifdef DEBUG
    std::cout << "problem, trying to refine edge " << one->id << "-" << two->id << " of a digon on step " << step_count << "\n";
    std::cout << "size of side_bodies is " << side_bodies.size() << std::endl;
    nodes[23423423]->output();
    one->output();
    two->output();
#endif
    return 0;
  }

  for (unsigned int c = 0; c < 3; c++)
    nodes[node_count]->neighbors.push_back(side_bodies[c]);

  CNode *side_face_nodes[3];
  for (unsigned int c = 0; c < 3; c++)
  {    
    side_face_nodes[c] = common_face(side_bodies[c], side_bodies[(c + 1) % 3], one);
  }

  for (unsigned int c = 0; c < 3; c++)
  {
    // we insert triplets in appropriate places.  if we have a triplet v1 v2 v3,
    // and v1 and v3 are one and two respectively, then we want two triplets,
    // v1 v2 new, new v2 v3.  we need to do v1 v2 v3 and v3 v2 v1 for each body.
    for (unsigned int d = 1; d < side_bodies[c]->triplets.size(); d++)
    {
      if (side_bodies[c]->triplets[d].v1 == one && side_bodies[c]->triplets[d].v3 == two)
      {
        CNode *temp = side_bodies[c]->triplets[d].v3;
        side_bodies[c]->triplets[d].v3 = nodes[node_count];
        side_bodies[c]->triplets.insert(side_bodies[c]->triplets.begin() + d + 1, triplet(nodes[node_count], side_bodies[c]->triplets[d].v2, temp));
      }
      if (side_bodies[c]->triplets[d].v1 == two && side_bodies[c]->triplets[d].v3 == one)
      {
        CNode *temp = side_bodies[c]->triplets[d].v3;
        side_bodies[c]->triplets[d].v3 = nodes[node_count];
        side_bodies[c]->triplets.insert(side_bodies[c]->triplets.begin() + d + 1, triplet(nodes[node_count], side_bodies[c]->triplets[d].v2, temp));
      }
    }
  }

  for (unsigned int c = 0; c < 3; c++)
    side_face_nodes[c]->determine_node_topology();
  one->determine_node_topology();
  two->determine_node_topology();
  nodes[node_count]->determine_node_topology();

  edge_nodes_added++;
  return nodes[node_count];
}

void CNode::remove_edge_node()
{
#ifdef DEBUG
  std::cout << "About to remove edge node " << id << "   ";
  output();
#endif

  std::vector<CNode *> neighboring_nodes;              // not sure these are necessary, they may be the same as one and two in the end.
  for (unsigned int c = 0; c < corners[0].size(); c++) // maybe this is a better way of getting one and two
    if (corners[0][c]->neighbors.size() != 2)
      neighboring_nodes.push_back(corners[0][c]);

  std::vector<CNode *> next_nodes;
  for (unsigned int d = 0; d < corners[0].size(); d++)
  {
    if (corners[0][d]->neighbors.size() != 2)
      next_nodes.push_back(corners[0][d]);
  }
  std::sort(next_nodes.begin(), next_nodes.end(), compare);
  next_nodes.erase(unique(next_nodes.begin(), next_nodes.end()), next_nodes.end());

  CNode *one = 0;
  CNode *two = 0;
  CNode *sides[3];
  for (unsigned int c = 0; c < 3; c++)
    sides[c] = common_face(neighbors[c], neighbors[(c + 1) % 3], this);

  double iarea[3];
  for (unsigned int c = 0; c < 3; c++)
    iarea[c] = sides[c]->get_face_area();

  // For each body we need to get rid of two triplets.
  for (unsigned int c = 0; c < 3; c++)
  {
    for (unsigned int d = 0; d < neighbors[c]->triplets.size(); d++)
    {
      if (neighbors[c]->triplets[d].v3 == this)
        one = neighbors[c]->triplets[d].v1;
      if (neighbors[c]->triplets[d].v1 == this)
        two = neighbors[c]->triplets[d].v3;
    }

    for (unsigned int d = 1; d < neighbors[c]->triplets.size(); d++)
      if (neighbors[c]->triplets[d].v3 == this)
      {
        neighbors[c]->triplets.erase(neighbors[c]->triplets.begin() + d);
        d--;
      }

    for (unsigned int d = 1; d < neighbors[c]->triplets.size(); d++)
    {
      if (neighbors[c]->triplets[d].v3 == one && neighbors[c]->triplets[d].v1 == this)
        neighbors[c]->triplets[d].v1 = two;
      if (neighbors[c]->triplets[d].v3 == two && neighbors[c]->triplets[d].v1 == this)
        neighbors[c]->triplets[d].v1 = one;
    }
  }

  for (unsigned int c = 0; c < neighboring_nodes.size(); c++)
    neighboring_nodes[c]->determine_node_topology();
  for (unsigned int c = 0; c < 3; c++)
    sides[c]->determine_node_topology();
  for (unsigned int c = 0; c < neighboring_nodes.size(); c++)
    neighboring_nodes[c]->center();
  for (unsigned int c = 0; c < 3; c++)
    sides[c]->center();

  for (unsigned int c = 0; c < 3; c++)
  {
    if (fabs(sides[c]->get_face_area() - iarea[c]) / iarea[c] > 0.2)
    {
#ifdef DEBUG
      std::cout << "We are changing the area of face " << sides[c]->id << " by much after erasing the edge node\n";
#endif
      sides[c]->center2(1.);
      for (unsigned int c = 0; c < neighboring_nodes.size(); c++)
        neighboring_nodes[c]->center();
    }
  }

  for (unsigned int c = 0; c < neighboring_nodes.size(); c++)
    neighboring_nodes[c]->calc_motion();
  for (unsigned int c = 0; c < 3; c++)
    sides[c]->calc_motion();

  neighbors.clear();
  edge_nodes_deleted++;
}

// remove_edge_node()
///////////////////////////////////////////////////////////////

void refine_edges()
{
  double minedge = 10. * smallest_edge;
  for (unsigned int c = 1; c < nodes.size(); c++)
  {
    if (nodes[c]->neighbors.size() == 2)
    {
      for (unsigned int d = 0; d < nodes[c]->corners[0].size(); d++)
      {
        if (length(nodes[c]->corners[0][d], nodes[c]->corners[0][(d + 1) % nodes[c]->corners[0].size()]) > minedge)
        {
          refine_edge(nodes[c]->corners[0][d], nodes[c]->corners[0][(d + 1) % nodes[c]->corners[0].size()]);
        }
      }
    }
  }
}

void remove_system_edge_nodes()
{
  for (unsigned int c = 1; c < nodes.size(); c++)
  {
    if (nodes[c]->neighbors.size() == 3)
    {
      double e1x = nodes[c]->x - nodes[c]->corners[0][1]->x;
      double e1y = nodes[c]->y - nodes[c]->corners[0][1]->y;
      double e1z = nodes[c]->z - nodes[c]->corners[0][1]->z;
      double e2x = nodes[c]->x - nodes[c]->corners[0][3]->x;
      double e2y = nodes[c]->y - nodes[c]->corners[0][3]->y;
      double e2z = nodes[c]->z - nodes[c]->corners[0][3]->z;

      if (e1x > 0.5)
        e1x -= 1.0;
      if (e1x < -0.5)
        e1x += 1.0;
      if (e1y > 0.5)
        e1y -= 1.0;
      if (e1y < -0.5)
        e1y += 1.0;
      if (e1z > 0.5)
        e1z -= 1.0;
      if (e1z < -0.5)
        e1z += 1.0;

      if (e2x > 0.5)
        e2x -= 1.0;
      if (e2x < -0.5)
        e2x += 1.0;
      if (e2y > 0.5)
        e2y -= 1.0;
      if (e2y < -0.5)
        e2y += 1.0;
      if (e2z > 0.5)
        e2z -= 1.0;
      if (e2z < -0.5)
        e2z += 1.0;

      double lx = e1x - e2x;
      double ly = e1y - e2y;
      double lz = e1z - e2z;

      if (sqrt(nodes[c]->dx * nodes[c]->dx + nodes[c]->dy * nodes[c]->dy + nodes[c]->dz * nodes[c]->dz) * time_scale / sqrt(lx * lx + ly * ly + lz * lz) > 0.02)
      {
        if (nodes[c]->neighbors[0]->is_football() || nodes[c]->neighbors[1]->is_football() || nodes[c]->neighbors[2]->is_football() ||
            nodes[c]->neighbors[0]->is_tetrahedron() || nodes[c]->neighbors[1]->is_tetrahedron() || nodes[c]->neighbors[2]->is_tetrahedron())
          continue;

        nodes[c]->remove_edge_node();
      }
    }
  }
}

void remove_small_edges()
{
  int o2 = offcounter2;
  int o3 = offcounter3;

  double initial_time_scale = time_scale;
  time_scale /= 10.0;

  // I THINK I SHOULD FIX THIS STUFF UP BELOW?  I'M NOT SURE WHY, BUT BASICALLY EVERY TIME WE RUN IT, WE CHECK EVERY EDGE TWICE,
  // ONCE FOR A-B AND ONCE FOR B-A.  I ASSUME THAT'S A HUGE WASTE OF TIME, BUT I'M NOT SURE WHY I IMPLEMENTED IT THAT WAY IN THE
  // BEGINNING.  ANYWAY, WORTH LOOKING INTO.

#pragma omp parallel
  for (unsigned int c = omp_get_thread_num() + 1; c < nodes.size(); c += nthreads)
  {
    nodes[c]->needsfixing = 0;
    if (nodes[c]->neighbors.size() == 4)
      for (int d = 0; d < 4; d++)
        for (int e = 1; e < 6; e += 2)
          if (nodes[c]->corners[d][e]->neighbors.size() == 4)
            if (nodes[c]->id < nodes[c]->corners[d][e]->id)
              if (calc_echange(nodes[c], nodes[c]->corners[d][e]) < -0.02)
                nodes[c]->needsfixing = 1;
  }

  for (unsigned int c = 1; c < nodes.size(); c++)
  {
    if (nodes[c]->needsfixing == 0 || nodes[c]->neighbors.size() != 4)
      continue;

    // WE WILL CONSIDER AT MOST ONE NEIGHBORING EDGE PER NODE.  IF WE DON'T FIND ONE, WE CONTINUE TO OTHER NODES.
    CNode *other = 0;
    for (int d = 0; d < 4; d++)
      for (int e = 1; e < 6; e += 2)
        if (nodes[c]->corners[d][e]->neighbors.size() == 4)
          if (calc_echange(nodes[c], nodes[c]->corners[d][e]) < -0.02)
            other = nodes[c]->corners[d][e];

    if (other == 0)
      continue;

    std::vector<CNode *> adjacent_faces;
    std::vector<CBody *> adjacent_bodies;
    std::vector<CNode *> adjacent_nodes;

    for (unsigned int p = 0; p < 4; p++)
      for (unsigned int q = 1; q < 6; q += 2)
        adjacent_nodes.push_back(nodes[c]->corners[p][q]);
    std::sort(adjacent_nodes.begin(), adjacent_nodes.end(), compare);
    adjacent_nodes.erase(unique(adjacent_nodes.begin(), adjacent_nodes.end()), adjacent_nodes.end());
    if (adjacent_nodes.size() != 4)
      continue;

    adjacent_nodes.clear();
    for (unsigned int p = 0; p < 4; p++)
      for (unsigned int q = 1; q < 6; q += 2)
        adjacent_nodes.push_back(other->corners[p][q]);
    sort(adjacent_nodes.begin(), adjacent_nodes.end(), compare);
    adjacent_nodes.erase(unique(adjacent_nodes.begin(), adjacent_nodes.end()), adjacent_nodes.end());
    if (adjacent_nodes.size() != 4)
      continue;

    for (unsigned int f = 0; f < 4; f++)
      if (is_neighbor(nodes[c]->neighbors[f], other))
        adjacent_bodies.push_back(nodes[c]->neighbors[f]);

    // IF THE TWO END NODES SHARE 4 BODIES IN COMMON, THEY ARE PART OF A DIGON  AND WE SHOULDN'T ERASE THIS EDGE.
    // IF ANY OF THE ADJACENT BODIES ARE FOOTBALLS, WE ALSO SHOULDN'T ERASE THIS EDGE
    if (adjacent_bodies.size() == 4 || adjacent_bodies[0]->is_football() || adjacent_bodies[1]->is_football() || adjacent_bodies[2]->is_football())
      continue;

    adjacent_faces.push_back(common_face(adjacent_bodies[0], adjacent_bodies[1], nodes[c]));
    adjacent_faces.push_back(common_face(adjacent_bodies[1], adjacent_bodies[2], nodes[c]));
    adjacent_faces.push_back(common_face(adjacent_bodies[2], adjacent_bodies[0], nodes[c]));

    CNode *nodeD = nodes[c];
    CNode *nodeU = other;

    CBody *up = 0;
    CBody *down = 0;

    for (unsigned int p = 0; p < 4; p++)
      if (is_neighbor(nodeU->neighbors[p], nodeD) == 0)
        up = nodeU->neighbors[p];
    for (unsigned int p = 0; p < 4; p++)
      if (is_neighbor(nodeD->neighbors[p], nodeU) == 0)
        down = nodeD->neighbors[p];

    std::vector<CNode *> top_nodes(3);
    std::vector<CNode *> bottom_nodes(3);

    int index = 0;
    for (unsigned int p = 0; p < 4; p++)
      if (nodeU->neighbors[p] == up)
        index = p;
    for (unsigned int p = 0; p < 3; p++)
    {
      if (are_neighbors(adjacent_faces[0], nodeU->corners[index][2 * p + 1]))
        top_nodes[0] = nodeU->corners[index][2 * p + 1];
      if (are_neighbors(adjacent_faces[1], nodeU->corners[index][2 * p + 1]))
        top_nodes[1] = nodeU->corners[index][2 * p + 1];
      if (are_neighbors(adjacent_faces[2], nodeU->corners[index][2 * p + 1]))
        top_nodes[2] = nodeU->corners[index][2 * p + 1];
    }

    for (unsigned int p = 0; p < 4; p++)
      if (nodeD->neighbors[p] == down)
        index = p;
    for (unsigned int p = 0; p < 3; p++)
    {
      if (are_neighbors(adjacent_faces[0], nodeD->corners[index][2 * p + 1]))
        bottom_nodes[0] = nodeD->corners[index][2 * p + 1];
      if (are_neighbors(adjacent_faces[1], nodeD->corners[index][2 * p + 1]))
        bottom_nodes[1] = nodeD->corners[index][2 * p + 1];
      if (are_neighbors(adjacent_faces[2], nodeD->corners[index][2 * p + 1]))
        bottom_nodes[2] = nodeD->corners[index][2 * p + 1];
    }

    // CRASHED HERE!!!
    for (unsigned int q = 0; q < top_nodes.size(); q++)
      if (top_nodes[q] == 0)
      {
        adjacent_bodies[0]->output();
        adjacent_bodies[1]->output();
        adjacent_bodies[2]->output();
        up->output();
        down->output();
        adjacent_bodies[0]->evolver_out();
        adjacent_bodies[1]->evolver_out();
        adjacent_bodies[2]->evolver_out();
        std::cout << "up and down are " << up << " " << down << '\n';
        std::cout << "up and down ids are " << up->id << " " << down->id << '\n';
        up->evolver_out();
        down->evolver_out();
        std::cout << "We have top nodes " << top_nodes[0] << " " << top_nodes[1] << " " << top_nodes[2] << '\n';
        std::cout << "We have adjacent faces " << adjacent_faces[0] << " " << adjacent_faces[1] << " " << adjacent_faces[2] << '\n';
        std::cout << "We have adjacent face ids " << adjacent_faces[0]->id << " " << adjacent_faces[1]->id << " " << adjacent_faces[2]->id << '\n';
        adjacent_faces[0]->output();
        adjacent_faces[1]->output();
        adjacent_faces[2]->output();
        std::cout << "We have nodeU and nodeD " << nodeU << " " << nodeD << '\n';
        std::cout << "We have nodeU and nodeD ids " << nodeU->id << " " << nodeD->id << '\n';
        nodeU->output();
        nodeD->output();
        std::cout << "Here are the corners of nodeU: ";
        for (unsigned int y = 0; y < nodeU->neighbors.size(); y++)
        {
          for (unsigned int z = 0; z < nodeU->corners[y].size(); z++)
            std::cout << nodeU->corners[y][z]->id << '\t';
          std::cout << '\n';
        }
        for (unsigned int p = 0; p < 4; p++)
          if (nodeU->neighbors[p] == up)
            index = p;
        std::cout << "index is " << index << '\n';
        std::cout << "crashing in topchanges: remove_small_edges 1";
        crash();
      }
    std::sort(top_nodes.begin(), top_nodes.end(), compare);
    top_nodes.erase(unique(top_nodes.begin(), top_nodes.end()), top_nodes.end());

    for (unsigned int q = 0; q < bottom_nodes.size(); q++)
    {
      if (bottom_nodes[q] == 0)
      {
        std::cout << "crashing in topchanges: remove_small_edges 2";
        crash();
      }
    }
    std::sort(bottom_nodes.begin(), bottom_nodes.end(), compare);
    bottom_nodes.erase(unique(bottom_nodes.begin(), bottom_nodes.end()), bottom_nodes.end());
    if (top_nodes.size() != 3 || bottom_nodes.size() != 3)
      continue;

    double current_length = length(nodes[c], other);

    if (adjacent_faces[0]->get_face_area() > 1.5 * current_length * current_length &&
        adjacent_faces[1]->get_face_area() > 1.5 * current_length * current_length &&
        adjacent_faces[2]->get_face_area() > 1.5 * current_length * current_length)
    {
#ifdef DEBUG
      std::cout << "Ready to remove short edge " << nodes[c]->id << "-" << other->id << '\n';
      for (unsigned int p = 0; p < 3; p++)
        std::cout << adjacent_faces[p]->get_face_area() << '\t' << current_length * current_length << '\t' << current_length << '\t' << (2 * current_length * current_length - adjacent_faces[p]->get_face_area()) / adjacent_faces[p]->get_face_area() << '\n';
#endif
      remove_edge(nodes[c], other);
      edges_deleted++;
    }
  }
  time_scale = initial_time_scale;

  offcounter3 = o3;
  offcounter2 = o2;
}

void remove_small_faces()
{
  int o2 = offcounter2;
  int o3 = offcounter3;

  double initial_time_scale = time_scale;
  time_scale /= 10.0;

#pragma omp parallel
  for (unsigned int c = omp_get_thread_num() + 1; c < nodes.size(); c += nthreads)
  {
    if (nodes[c]->neighbors.size() != 2)
      continue;

    nodes[c]->erase = 0;
    if (nodes[c]->area < 0.1 * avg_face_area && (nodes[c]->calc_change() < -0.02 || nodes[c]->get_face_area() < 0.000000001 || nodes[c]->get_face_perimeter() / sqrt(nodes[c]->get_face_area()) > 16.0))
      nodes[c]->erase = 1;
    if (nodes[c]->erase == 1)
      if (nodes[c]->neighbors[0]->is_football() || nodes[c]->neighbors[1]->is_football() || 2 * nodes[c]->neighbors[0]->get_volume() < pow(nodes[c]->get_face_area(), 1.5) || 2 * nodes[c]->neighbors[1]->get_volume() < pow(nodes[c]->get_face_area(), 1.5))
        nodes[c]->erase = 0;
  }

  for (unsigned int c = 1; c < nodes.size(); c++)
  {
    if (nodes[c]->erase == 0 || nodes[c]->neighbors.size() != 2 || nodes[c]->neighbors[0]->is_football() || nodes[c]->neighbors[1]->is_football() || 2 * nodes[c]->neighbors[0]->get_volume() < pow(nodes[c]->get_face_area(), 1.5) || 2 * nodes[c]->neighbors[1]->get_volume() < pow(nodes[c]->get_face_area(), 1.5))
      continue;

#ifdef DEBUG
    std::cout << "Ready to remove small face " << nodes[c]->id << " with " << nodes[c]->get_face_sides() << " sides\n";
    //    cout << " ... quickly losing area " << current_area << '\t' << projected_area << '\t' << (projected_area-current_area)/current_area << '\n';
    //    cout << "Face " << nodes[c]->id << " is shrinking: " << (projected_area-current_area)/current_area << " and has " << nodes[c]->corners[0].size() << " and current area ";
    //    cout << "is " << current_area << " and projected area is " << projected_area << '\n';
#endif
    if (nodes[c]->get_face_sides() == 2)
      digons_deleted++;
    if (nodes[c]->get_face_sides() == 3)
      triangles_deleted++;

    CBody *b1 = nodes[c]->neighbors[0];
    CBody *b2 = nodes[c]->neighbors[1];

    double bvol1 = b1->get_volume();
    double bvol2 = b2->get_volume();

    nodes[c]->remove_face();

    if (fabs(b1->get_volume() - bvol1) / b1->get_volume() > 0.5)
    {
#ifdef DEBUG
      std::cout << "Wo, we've changed a lot for body " << b1->id << ", will try removing it.\n";
#endif
      b1->remove_body();
    }
    if (fabs(b2->get_volume() - bvol2) / b2->get_volume() > 0.5)
    {
#ifdef DEBUG
      std::cout << "Wo, we've changed a lot for body " << b2->id << ", will try removing it.\n";
#endif
      b2->remove_body();
    }
  }
  time_scale = initial_time_scale;
  offcounter3 = o3;
  offcounter2 = o2;
}

void remove_small_bodies()
{
  int o2 = offcounter2;
  int o3 = offcounter3;

  for (unsigned int c = 1; c < bodies.size(); c++)
  {
    if (bodies[c]->volume < 0.01 / body_count)
    {
      if (bodies[c]->triplets.size() < 2)
        continue;

      double current_volume = bodies[c]->get_volume();

      double initial_time_scale = time_scale;
      time_scale /= 10.0;

      double projected_volume = bodies[c]->get_pvolume();
      time_scale = initial_time_scale;

      // HERE WE NEED TO CHECK FOR A SPECIAL TYPE OF PROBLEM.  BASICALLY, WE CAN HAVE A BODY WHOSE TOPOLOGY
      // IS A FOOTBALL (THREE 2-SIDED FACES), BUT WHOSE SHAPE REALLY MORE RESEMBLES A LENS WITH AN EXTRA
      // DIGON ALONG ITS EDGE.  WE CAN'T ERASE THE DIGON BECAUSE THAT WILL CREATE A LENS AND WE DON'T ALLOW
      // THOSE BUT SINCE THE FOOTBALL IS RELATIVELY LARGE, WE WILL ALLOW THE DIGON TO COLLAPSE, ESSENTIALLY,
      // TO A POINT.  WHAT WE WILL TRY DOING HERE IS DETECTING THIS POOR SHAPE BY COMPARING THE RELATIVE
      // SIZES OF THE THREE FACES.  IF ONE IS CONSIDERABLY SMALLER THAN THE OTHER TWO (BY A FACTOR OF 10?)
      // THEN WE KNOW WE HAVE THIS PROBLEM AND DECIDE TO ERASE THE FOOTBALL

      std::vector<double> facesn(3);
      for (unsigned int p = 0; p < 3; p++)
        facesn[p] = 1.;
      if (bodies[c]->is_football())
      {
        std::vector<CNode *> facenodes;
        for (unsigned int d = 1; d < bodies[c]->triplets.size(); d++)
          facenodes.push_back(bodies[c]->triplets[d].v2);
        std::sort(facenodes.begin(), facenodes.end(), compare);
        facenodes.erase(unique(facenodes.begin(), facenodes.end()), facenodes.end());
        for (unsigned int p = 0; p < 3; p++)
          facesn[p] = facenodes[p]->get_face_area();
        std::sort(facesn.begin(), facesn.end());

#ifdef DEBUG
        std::cout << "We have football " << bodies[c]->id << " with " << facenodes.size() << " faces here: " << facenodes[0]->get_face_area() << '\t' << facenodes[1]->get_face_area() << '\t' << facenodes[2]->get_face_area() << '\n';
        for (unsigned int p = 0; p < 3; p++)
          std::cout << facesn[p] << '\t';
        std::cout << '\n';
        std::cout << facesn[1] / facesn[0] << '\t' << facesn[2] / facesn[0] << '\n';
        std::cout << '\n';
#endif
      }

      if ((projected_volume - current_volume) / current_volume < -0.04 || current_volume < 0.000000000001 || facesn[1] / facesn[0] > 100. || facesn[2] / facesn[0] > 100.)
        bodies[c]->remove_body();
    }
  }
  offcounter3 = o3;
  offcounter2 = o2;
}
