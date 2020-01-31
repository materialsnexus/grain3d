//#define DEBUG

#include "cbody.h"
#include "cnode.h"
#include "functions.h"
#include "globals.h"

#include <ostream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <iomanip>
#include <vector>
#include <algorithm>



double CBody::get_area()
{
    double surface_area = 0;
    for(unsigned int c=1; c<triplets.size(); c++)
        surface_area += triplets[c].area();
    return surface_area;
}


void CBody::calc_volume()
{
    volume = 0;
    
    if(triplets.size()<1) return;
    if(triplets.size()>1)
    {
        CNode* origin = triplets[1].v1;
        for(unsigned int c=2; c<triplets.size(); c++)
            volume += det(triplets[c].v1, triplets[c].v2, triplets[c].v3, origin)/6;
    }
}


double CBody::get_volume()
{
    if(triplets.size()<1) return 0;
    
    double volume=0;
    CNode* origin = triplets[1].v1;
    for(unsigned int c=2; c<triplets.size(); c++)
        volume += det(triplets[c].v1, triplets[c].v2, triplets[c].v3, origin)/6;
    
    return volume;
}


double CBody::get_pvolume()
{
    if(triplets.size()<1) return 0;
    
    double volume=0;
    CNode* origin = triplets[1].v1;
    for(unsigned int c=1; c<triplets.size(); c++)
        volume += pdet(triplets[c].v1, triplets[c].v2, triplets[c].v3, origin)/6;
    
    return volume;
}


void CBody::blank()
{
    id        = 0;
    volume    = 0;
    triplets.clear();
}



void CBody::output()
{
    std::cout << "Body " << id << " with " << get_faces() << " faces, " << triplets.size()-1 << " triplets on step " << step_count << ": " << std::endl;
    for(unsigned int c=1; c<triplets.size(); c++)
    {
        std::cout << triplets[c].v1->id << " ";
        std::cout << triplets[c].v2->id << " ";
        std::cout << triplets[c].v3->id << "\t";
        //    if(c%8==0) std::cout << std::endl;
        if(c<triplets.size()-1) if(triplets[c].v2!=triplets[c+1].v2) std::cout << '\n';
    }
    std::cout << "Body has volume " << '\t' << volume << '\t' << get_volume() << '\n';
    std::cout << std::endl << std::endl;
}



void CBody::simple_remove_face(CNode* face_node)
{
    for(unsigned int c=1; c<triplets.size(); c++)
        if(triplets[c].v2==face_node)
        {
            triplets.erase(triplets.begin()+c);
            c--;
        }
}


bool CBody::is_tetrahedron()
{
    int flag=1;
    int qd=0;
    for(unsigned int c=1; c<triplets.size(); c++) if(triplets[c].v1->neighbors.size()==4) qd++;
    for(unsigned int c=1; c<triplets.size(); c++) if(triplets[c].v2->get_face_sides()!=3) flag=0;
    if(qd==12 && flag==1) return 1;
    else return 0;
}


bool CBody::is_football()
{
    // Any three-faced body is a football (3 two-sided faces).  There are no other top. possibilites
    std::vector <CNode*> face_nodes;
    for(unsigned int c=1; c<triplets.size(); c++) face_nodes.push_back(triplets[c].v2);
    sort(face_nodes.begin(), face_nodes.end(), compare);
    face_nodes.erase(std::unique(face_nodes.begin(), face_nodes.end()), face_nodes.end());
    return face_nodes.size()==3;
}



int CBody::get_faces()
{
    int   faces = 0;
    CNode* last = 0;
    for(unsigned int c=1; c<triplets.size(); c++)
        if(last != triplets[c].v2) { last = triplets[c].v2; faces++;}
    return faces;
}



void CBody::determine_node_topologies()
{
    std::vector <CNode*> our_nodes;
    for(unsigned int c=1; c<triplets.size(); c++)
    {
        our_nodes.push_back(triplets[c].v1);
        our_nodes.push_back(triplets[c].v2);
    }
    sort(our_nodes.begin(), our_nodes.end(), compare);
    our_nodes.erase(std::unique(our_nodes.begin(), our_nodes.end()), our_nodes.end());
    
    for(unsigned int c=0; c<our_nodes.size(); c++)
        our_nodes[c]->determine_node_topology();
}


void CBody::replace_edges_with_node(CNode* e1, CNode* e2, CNode* e3, CNode* newnode)
{
    for(unsigned int c=1; c<triplets.size(); c++)
    {
        if((triplets[c].v1==e1 && triplets[c].v3==e2) || (triplets[c].v1==e2 && triplets[c].v3==e1)) { triplets.erase(triplets.begin()+c); c--; }
        if((triplets[c].v1==e2 && triplets[c].v3==e3) || (triplets[c].v1==e3 && triplets[c].v3==e2)) { triplets.erase(triplets.begin()+c); c--; }
        if((triplets[c].v1==e3 && triplets[c].v3==e1) || (triplets[c].v1==e1 && triplets[c].v3==e3)) { triplets.erase(triplets.begin()+c); c--; }
    }
    
    for(unsigned int c=1; c<triplets.size(); c++)
    {
        if(triplets[c].v1==e1) triplets[c].v1 = newnode;
        if(triplets[c].v3==e1) triplets[c].v3 = newnode;
        if(triplets[c].v1==e2) triplets[c].v1 = newnode;
        if(triplets[c].v3==e2) triplets[c].v3 = newnode;
        if(triplets[c].v1==e3) triplets[c].v1 = newnode;
        if(triplets[c].v3==e3) triplets[c].v3 = newnode;
    }
}




void CBody::evolver_out(void)
{
    if(get_volume()==0) return;
    
    std::vector <CNode*> the_nodes;
    
    for(unsigned int c=1; c<triplets.size(); c++)
    {
        the_nodes.push_back(triplets[c].v1);
        the_nodes.push_back(triplets[c].v2);
    }
    
    sort(the_nodes.begin(), the_nodes.end(), compare);
    the_nodes.erase(std::unique(the_nodes.begin(), the_nodes.end()), the_nodes.end());
    
    int nodeids[the_nodes.size()];
    for(unsigned int c=0; c<the_nodes.size(); c++) nodeids[c]=the_nodes[c]->id;
    
    CNode origin = *the_nodes[0];
    
    char file_name[20];
    sprintf(file_name, "%.4d-%.8d.fe", step_count, id);
    std::ofstream out_file;
    out_file.open(file_name);
    
    out_file << "vertices\n";
    for(unsigned int c=0; c<the_nodes.size(); c++)
    {
        out_file << c+1 << '\t';
        out_file << (origin + (*the_nodes[c]-origin)).x << '\t';
        out_file << (origin + (*the_nodes[c]-origin)).y << '\t';
        out_file << (origin + (*the_nodes[c]-origin)).z << '\t';
        out_file << "original " << the_nodes[c]->id << std::endl;
    }
    
    int edge_counter=0;
    int size = the_nodes.size() + triplets.size() - 2;
    int edge1[size];
    int edge2[size];
    int edge_id[size];
    
    out_file << std::endl;
    out_file << "edges\n";
    
    for(unsigned int d=1; d<triplets.size(); d++)
    {
        if(triplets[d].v1->id > triplets[d].v2->id)
        {
            edge1[edge_counter]=triplets[d].v1->id;
            edge2[edge_counter]=triplets[d].v2->id;
            edge_id[edge_counter]=edge_counter;
            out_file << edge_counter+1 << '\t';
            out_file << f(nodeids, edge1[edge_counter], the_nodes.size()) << '\t';
            out_file << f(nodeids, edge2[edge_counter], the_nodes.size()) << '\t';
            out_file << std::endl;
            edge_counter++;
        }
        if(triplets[d].v2->id > triplets[d].v3->id)
        {
            edge1[edge_counter]=triplets[d].v2->id;
            edge2[edge_counter]=triplets[d].v3->id;
            edge_id[edge_counter]=edge_counter;
            out_file << edge_counter+1 << '\t';
            out_file << f(nodeids, edge1[edge_counter], the_nodes.size()) << '\t';
            out_file << f(nodeids, edge2[edge_counter], the_nodes.size()) << '\t';
            out_file << std::endl;
            edge_counter++;
        }
        if(triplets[d].v3->id > triplets[d].v1->id)
        {
            edge1[edge_counter]=triplets[d].v3->id;
            edge2[edge_counter]=triplets[d].v1->id;
            edge_id[edge_counter]=edge_counter;
            out_file << edge_counter+1 << '\t';
            out_file << f(nodeids, edge1[edge_counter], the_nodes.size()) << '\t';
            out_file << f(nodeids, edge2[edge_counter], the_nodes.size()) << '\t';
            out_file << "color red " << std::endl;
            edge_counter++;
        }
    }
    
    out_file << std::endl;
    out_file << "faces\n";
    
    for(unsigned int d=1; d<triplets.size(); d++)
    {
        if(g(edge1, edge2, edge_id, triplets[d].v1->id, triplets[d].v2->id, size)==0||
           g(edge1, edge2, edge_id, triplets[d].v2->id, triplets[d].v3->id, size)==0||
           g(edge1, edge2, edge_id, triplets[d].v3->id, triplets[d].v1->id, size)==0)
        {
            output();
            std::cout << "crashing on body " << id << " with face ";
            std::cout << triplets[d].v1->id << ", " << triplets[d].v2->id << ", " << triplets[d].v3->id << "\n";
            out_file.close();
            crash();
        }
        out_file << d << '\t'
        << g(edge1, edge2, edge_id, triplets[d].v1->id, triplets[d].v2->id, size) << '\t'
        << g(edge1, edge2, edge_id, triplets[d].v2->id, triplets[d].v3->id, size) << '\t'
        << g(edge1, edge2, edge_id, triplets[d].v3->id, triplets[d].v1->id, size) << '\t'
        << std::endl;
    }
    
    out_file << std::endl;
    out_file << "bodies\n";
    out_file << 1 << '\t';
    
    for(unsigned int d=1; d<triplets.size(); d++)
        out_file << d << ' ';
    
    out_file << std::endl;
    out_file.close();
}




void CBody::remove_body()
{
    std::vector <CNode*> face_nodes;
    for(unsigned int c=1; c<triplets.size(); c++)
        face_nodes.push_back(triplets[c].v2);
    sort(face_nodes.begin(), face_nodes.end(), compare);
    face_nodes.erase(std::unique(face_nodes.begin(), face_nodes.end()), face_nodes.end());
    
    for(unsigned int c=0; c<face_nodes.size(); c++) face_nodes[c]->remove_edge_nodes();
    for(unsigned int c=0; c<face_nodes.size(); c++) face_nodes[c]->center2(1.);
    
    
    while(is_football()==0 && is_tetrahedron()==0)
    {
        int index=0;
        double min_area =100;
        double min_sides=100;
        
        for(unsigned int c=0; c<face_nodes.size(); c++)
        {
            double current_face_area = face_nodes[c]->get_face_area();
            if(current_face_area < min_area)
                if(face_nodes[c]->get_face_sides() <= min_sides)
                {
                    index    = c;
                    min_area = current_face_area;
                    min_sides= face_nodes[c]->get_face_sides();
                }
        }
#ifdef DEBUG
        std::cout << "We are inside remove_body(), removing a face\n";
#endif
        face_nodes[index]->remove_face();
        face_nodes.erase(face_nodes.begin() + index);
        for(unsigned int c=0; c<face_nodes.size(); c++) face_nodes[c]->center2(1.);
    }
    
    if(is_football())    remove_football();
    if(is_tetrahedron()) remove_tetrahedron();
    compactify_bodies();
}




void cluster_evolver_out(int pid)
{
    unsigned int id=0;
    for(unsigned int c=1; c<bodies.size(); c++)
        if(bodies[c]->id==pid) id=c;
    
    if(id==0) exit(-1);
    
    std::vector <CBody*> cluster;
    for(unsigned int c=1; c<bodies[id]->triplets.size(); c++) if(bodies[id]->triplets[c].v2->neighbors[0]!=bodies[id]) cluster.push_back(bodies[id]->triplets[c].v2->neighbors[0]);
    else cluster.push_back(bodies[id]->triplets[c].v2->neighbors[1]);
    
    sort(cluster.begin(), cluster.end(), bcompare);
    cluster.erase(std::unique(cluster.begin(), cluster.end()), cluster.end());
    
    std::vector <CNode*> cluster_nodes;
    for(unsigned int c=0; c<cluster.size(); c++)
        for(unsigned int d=1; d<cluster[c]->triplets.size(); d++)
        {
            cluster_nodes.push_back(cluster[c]->triplets[d].v1);
            cluster_nodes.push_back(cluster[c]->triplets[d].v2);
        }
    
    sort(cluster_nodes.begin(), cluster_nodes.end(), compare);
    cluster_nodes.erase(std::unique(cluster_nodes.begin(), cluster_nodes.end()), cluster_nodes.end());
    
    std::cout << "vertices\n";
    for(unsigned int c=0; c<cluster_nodes.size(); c++)
        std::cout << cluster_nodes[c]->id << '\t' << cluster_nodes[c]->x << '\t' << cluster_nodes[c]->y << '\t' << cluster_nodes[c]->z << '\n';
    std::cout << '\n';
    
    std::vector <edge> edges;
    for(unsigned int c=0; c<cluster_nodes.size(); c++)
        for(unsigned int d=c+1; d<cluster_nodes.size(); d++)
            if(are_neighbors(cluster_nodes[c], cluster_nodes[d]))
                edges.push_back(edge(cluster_nodes[c], cluster_nodes[d]));
    
    std::cout << "edges\n";
    for(unsigned int c=0; c<edges.size(); c++)
    {
        std::cout << c+1 << '\t' << edges[c].one->id << '\t' << edges[c].two->id << '\t';
        if(edges[c].one->neighbors.size()!=2 && edges[c].two->neighbors.size()!=2) std::cout << "color red";
        std::cout << '\n';
    }
    std::cout << '\n';
    
    int face_counter=1;
    std::cout << "faces\n";
    for(unsigned int c=0; c<cluster_nodes.size(); c++) if(cluster_nodes[c]->neighbors.size()==2)
    {
        for(unsigned int d=0; d<cluster_nodes[c]->corners[0].size(); d++)
        {
            std::cout
            << face_counter << '\t'
            << vector_search(edges, cluster_nodes[c]->corners[0][d], cluster_nodes[c]) << '\t'
            << vector_search(edges, cluster_nodes[c], cluster_nodes[c]->corners[0][(d+1)%cluster_nodes[c]->corners[0].size()]) << '\t'
            << vector_search(edges, cluster_nodes[c]->corners[0][(d+1)%cluster_nodes[c]->corners[0].size()], cluster_nodes[c]->corners[0][d]) << '\n';
            face_counter++;
        }
    }
}



int vector_search(std::vector<edge>& edges, CNode* one, CNode* two)
{
    for(unsigned int c=0; c<edges.size(); c++)
    {
        if(edges[c].one==one && edges[c].two==two) return  c+1;
        if(edges[c].one==two && edges[c].two==one) return -c-1;
    }
    
    return 0;
}



void CBody::calc_center()
{
    std::vector<CNode*> vert_nodes;
    for(unsigned int d=1; d<triplets.size(); d++) if(triplets[d].v1->neighbors.size()==4) vert_nodes.push_back(triplets[d].v1);
    sort(vert_nodes.begin(), vert_nodes.end());
    vert_nodes.erase(std::unique(vert_nodes.begin(), vert_nodes.end()), vert_nodes.end());
    
    
    double number = double(vert_nodes.size());
    
    center[0] = vert_nodes[0]->x;
    center[1] = vert_nodes[0]->y;
    center[2] = vert_nodes[0]->z;
    
    for(unsigned int d=0; d<vert_nodes.size(); d++)
    {
        double dx = vert_nodes[d]->x - vert_nodes[0]->x;
        double dy = vert_nodes[d]->y - vert_nodes[0]->y;
        double dz = vert_nodes[d]->z - vert_nodes[0]->z;
        
        if(dx <-0.5) dx += 1.0;
        if(dy <-0.5) dy += 1.0;
        if(dz <-0.5) dz += 1.0;
        if(dx > 0.5) dx -= 1.0;
        if(dy > 0.5) dy -= 1.0;
        if(dz > 0.5) dz -= 1.0;
        
        center[0] += dx/number;
        center[1] += dy/number;
        center[2] += dz/number;
    }
}











