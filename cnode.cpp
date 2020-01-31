//#define DEBUG

#include "cbody.h"
#include "cnode.h"
#include "functions.h"
#include "globals.h"

#include <algorithm>
#include <ostream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>



void CNode::move()
{
    if(neighbors.size()==0) return;   // double check whether this ever happens.  maybe yes?
    
    if(good==0)
    {
        x += cx;
        y += cy;
        z += cz;
    }
    
    else
    {
        x += dx * time_scale;
        y += dy * time_scale;
        z += dz * time_scale;
    }
}


double CNode::norm()
{
    return sqrt(x*x + y*y + z*z);
}


CNode CNode::operator* (double param) {
    CNode temp;
    temp.x = x * param;
    temp.y = y * param;
    temp.z = z * param;
    return temp;
}


CNode CNode::operator/ (double param) {
    CNode temp;
    temp.x = x / param;
    temp.y = y / param;
    temp.z = z / param;
    return temp;
}

CNode CNode::operator+ (CNode increase)
{
    CNode temp;
    temp.x = x + increase.x;
    temp.y = y + increase.y;
    temp.z = z + increase.z;
    return temp;
}

void CNode::operator+= (CNode param) {
    x+=param.x;
    y+=param.y;
    z+=param.z;
}

void CNode::operator/= (double param) {
    x/=param;
    y/=param;
    z/=param;
}

void CNode::operator*= (double param) {
    x*=param;
    y*=param;
    z*=param;
}

CNode CNode::operator- (CNode param) {
    CNode temp;
    temp.x = x - param.x;
    temp.y = y - param.y;
    temp.z = z - param.z;
    if (temp.x > 0.5) temp.x -= 1.0;
    if (temp.x <-0.5) temp.x += 1.0;
    if (temp.y > 0.5) temp.y -= 1.0;
    if (temp.y <-0.5) temp.y += 1.0;
    if (temp.z > 0.5) temp.z -= 1.0;
    if (temp.z <-0.5) temp.z += 1.0;
    return temp;
}

void CNode::operator= (CNode param) {
    id = param.id;
    x  = param.x;
    y  = param.y;
    z  = param.z;
    dx  = param.dx;
    dy  = param.dy;
    dz  = param.dz;
    
    neighbors = param.neighbors;
    for(unsigned int c=0; c<4; c++)
        corners[c] = param.corners[c];
}



void CNode::output(void)
{
    std::cout.precision(8);
    std::cout << "On step " << step_count << " outputting node " << id << " with " << corners[0].size() << " nodes and valence " << neighbors.size() << ": ";
    for(unsigned int c=0; c<neighbors.size(); c++) std::cout << neighbors[c]->id << ", ";
    std::cout << "Coordinates: ";
    std::cout.width(25);
    std::cout.precision(22);
    std::cout << std::left << x << "  ";
    std::cout.width(25);
    std::cout.precision(22);
    std::cout << std::left << y << "  ";
    std::cout.width(25);
    std::cout.precision(22);
    std::cout << std::left << z << "  ";
    std::cout.width(25);
    std::cout.precision(22);
    std::cout << std::left << dx << "  ";
    std::cout.width(25);
    std::cout.precision(22);
    std::cout << std::left << dy << "  ";
    std::cout.width(25);
    std::cout.precision(22);
    std::cout << std::left << dz << "  ";
    std::cout << '\n';
}



void CNode::remove_edge_nodes()   // HERE WE ESSENTIALLY DON'T COUNT ANY FUNNY CURVATURES THAT ARE
{                                 // TEMPORARILY CREATED AND SEEN ONLY WHEN REMOVING EDGE NODES
    int initoff1 = offcounter1;
    int initoff2 = offcounter2;
    int initoff3 = offcounter3;
    std::vector<CNode*> temp;
    for(unsigned int c=0; c<corners[0].size(); c++) if(corners[0][c]->neighbors.size()==3) temp.push_back(corners[0][c]);
    for(unsigned int c=0; c<temp.size(); c++) temp[c]->remove_edge_node();
    offcounter1=initoff1;
    offcounter2=initoff2;
    offcounter3=initoff3;
}


void CNode::center()
{
    calc_center();
    move_center();
}


void CNode::move_center()
{
    x += cx;
    y += cy;
    z += cz;
}


// This function moves this node nicely in a way that does not change the volumes of any of the adjacent bodies.
// Only possible for face and edge nodes.
void CNode::calc_center()
{
    ////////////////////////////////////////////////////////////////////
    // Centering face nodes
    if(neighbors.size()==2)
    {
        // WE PROJECT THE CENTER OF MASS OF THE FACE BOUNDARY ONTO THE PLANE-OF-ZERO-VOLUME-CHANGE.  THIS IS A GOOD
        // WAY TO 'CENTER' THE NODE.  IT IS INVARIENT UNDER REFINING THE EDGES, AND ALSO DOESN'T USE THE FACES THEMSELVES,
        // WHICH WOULD ACTUALLY MOVE WITH THIS MOTION.
        
        unsigned int size = corners[0].size();
        double edges[size][3];
        
        for(unsigned int c=0; c<size; c++)
        {
            edges[c][0] = corners[0][c]->x - x;
            edges[c][1] = corners[0][c]->y - y;
            edges[c][2] = corners[0][c]->z - z;
            
            if(edges[c][0] > 0.5) edges[c][0] -= 1.0;
            if(edges[c][0] <-0.5) edges[c][0] += 1.0;
            if(edges[c][1] > 0.5) edges[c][1] -= 1.0;
            if(edges[c][1] <-0.5) edges[c][1] += 1.0;
            if(edges[c][2] > 0.5) edges[c][2] -= 1.0;
            if(edges[c][2] <-0.5) edges[c][2] += 1.0;
        }
        
        double bounds[size];
        for(unsigned int c=0; c<size; c++)
            bounds[c]=sqrt((edges[c][0]-edges[(c+1)%size][0])*(edges[c][0]-edges[(c+1)%size][0])+
                           (edges[c][1]-edges[(c+1)%size][1])*(edges[c][1]-edges[(c+1)%size][1])+
                           (edges[c][2]-edges[(c+1)%size][2])*(edges[c][2]-edges[(c+1)%size][2]));
        
        double bx=0;
        double by=0;
        double bz=0;
        double bt=0;  // TOTAL BOUNDARY LENGTH
        
        for(unsigned int c=0; c<size; c++)
        {
            double lb = (bounds[c]+bounds[(c+1)%size])/2.;
            bx+=edges[c][0]*lb;
            by+=edges[c][1]*lb;
            bz+=edges[c][2]*lb;
            bt+=lb;
        }
        
        bx/=bt;
        by/=bt;
        bz/=bt;
        
        double dvx=0;
        double dvy=0;
        double dvz=0;
        
        if(size==0) 
        {
            std::cout << "crashing in cnode: calc_center";
            crash();
        }
        for(unsigned int c=0; c<size-1; c++)
        {
            dvx += edges[c][1]*edges[c+1][2] - edges[c][2]*edges[c+1][1];
            dvy += edges[c][2]*edges[c+1][0] - edges[c][0]*edges[c+1][2];
            dvz += edges[c][0]*edges[c+1][1] - edges[c][1]*edges[c+1][0];
        }
        
        dvx += edges[size-1][1]*edges[0][2] - edges[size-1][2]*edges[0][1];
        dvy += edges[size-1][2]*edges[0][0] - edges[size-1][0]*edges[0][2];
        dvz += edges[size-1][0]*edges[0][1] - edges[size-1][1]*edges[0][0];
        
        
        if(dvx*dvx+dvy*dvy+dvz*dvz==0)
        {
            cx=cy=cz=0.;
            double e1x, e1y, e1z;
            
            for(unsigned int c=0; c<size; c++)
            {
                e1x = corners[0][c]->x - x;
                e1y = corners[0][c]->y - y;
                e1z = corners[0][c]->z - z;
                
                if(e1x > 0.5) e1x -= 1.0;
                if(e1x <-0.5) e1x += 1.0;
                if(e1y > 0.5) e1y -= 1.0;
                if(e1y <-0.5) e1y += 1.0;
                if(e1z > 0.5) e1z -= 1.0;
                if(e1z <-0.5) e1z += 1.0;
                
                cx += e1x/size;
                cy += e1y/size;
                cz += e1z/size;
            }
        }
        
        else
        {
            double mult = (dvx*bx+dvy*by+dvz*bz)/(dvx*dvx+dvy*dvy+dvz*dvz);
            cx = bx - dvx*mult;
            cy = by - dvy*mult;
            cz = bz - dvz*mult;
        }
    }
    // Centering face nodes
    ////////////////////////////////////////////////////////////////////
    
    else if(neighbors.size()==3)
    {
        double e1x = corners[0][1]->x - x;
        double e1y = corners[0][1]->y - y;
        double e1z = corners[0][1]->z - z;
        double e2x = corners[0][3]->x - x;
        double e2y = corners[0][3]->y - y;
        double e2z = corners[0][3]->z - z;
        
        if(e1x > 0.5) e1x -= 1.0;
        if(e1x <-0.5) e1x += 1.0;
        if(e1y > 0.5) e1y -= 1.0;
        if(e1y <-0.5) e1y += 1.0;
        if(e1z > 0.5) e1z -= 1.0;
        if(e1z <-0.5) e1z += 1.0;
        if(e2x > 0.5) e2x -= 1.0;
        if(e2x <-0.5) e2x += 1.0;
        if(e2y > 0.5) e2y -= 1.0;
        if(e2y <-0.5) e2y += 1.0;
        if(e2z > 0.5) e2z -= 1.0;
        if(e2z <-0.5) e2z += 1.0;
        
        double bottom = (e2x-e1x)*(e2x-e1x)+(e2y-e1y)*(e2y-e1y)+(e2z-e1z)*(e2z-e1z);
        if(bottom==0)
        {
            cx = (e1x+e2x)/2.;
            cy = (e1y+e2y)/2.;
            cz = (e1z+e2z)/2.;
            
            return;
        }
        
        double constant = ((e1x+e2x)/2.*(e2x-e1x)+(e1y+e2y)/2.*(e2y-e1y)+(e1z+e2z)/2.*(e2z-e1z))/bottom;
        cx = constant*(e2x-e1x);
        cy = constant*(e2y-e1y);
        cz = constant*(e2z-e1z);
    }
    // Motions for edge nodes; at this point we don't have any.
    ////////////////////////////////////////////////////////////////////
    
    else  // THIS SHOULD COVER THE 4-NODES
    {
        cx = 0;
        cy = 0;
        cz = 0;
    }
}




////////////////////////////////////////////////////////////////////
// Checking area of one face
double CNode::get_face_area()
{
    if(neighbors.size()!=2 || corners[0].size()==2) return 0;
    double area = 0;
    for(unsigned int c=0; c<corners[0].size(); c++)
        area += triplet(corners[0][c],this,corners[0][(c+1)%corners[0].size()]).area();
    return area;
}
// Checking area of one faces
////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////
// Checking area of one face
void CNode::calc_face_area()
{
    area = 0;
    if(neighbors.size()!=2 || corners[0].size()==2) return;
    for(unsigned int c=0; c<corners[0].size(); c++)
        area += triplet(corners[0][c],this,corners[0][(c+1)%corners[0].size()]).area();
}
// Checking area of one faces
////////////////////////////////////////////////////////////////////



double CNode::get_face_perimeter(void)
{
    if(neighbors.size()!=2) return 0;
    double perimeter = 0;
    for(unsigned int c=0; c<corners[0].size(); c++)
        perimeter += length(corners[0][c],corners[0][(c+1)%corners[0].size()]);
    return perimeter;
}


////////////////////////////////////////////////////////////////////
// Calculating the topological information of each node
void CNode::determine_node_topology(void)
{
    sort(neighbors.begin(), neighbors.end(), bcompare);
    
    ////////////////////////////////////////////////////////////////////
    // Topology of face nodes
    if(neighbors.size()==2)
    {
        for(unsigned int n=0; n<2; n++)
        {
            corners[n].clear();
            for(unsigned int c=1; c<neighbors[n]->triplets.size(); c++)
                if(neighbors[n]->triplets[c].v2==this)
                    corners[n].push_back(neighbors[n]->triplets[c].v1);
        }
        
        return;
    }
    
    ////////////////////////////////////////////////////////////////////
    // Topology of edge nodes
    
    if(neighbors.size()==3)
    {
        CNode* temp1;
        CNode* temp2;
        CNode* temp3;
        
        for(unsigned int n=0; n<3; n++)
        {
            corners[n].clear();
            std::vector<triplet> temp;
            
            // finds the ordered, oriented edges associated with this node, for this body
            for(unsigned int c=1; c<neighbors[n]->triplets.size(); c++)
                if(neighbors[n]->triplets[c].v3==this)
                    temp.push_back(triplet(neighbors[n]->triplets[c].v2, neighbors[n]->triplets[c].v3, neighbors[n]->triplets[c].v1));
            
            for(unsigned int c=1; c<neighbors[n]->triplets.size(); c++)
                if(neighbors[n]->triplets[c].v1==this)
                    temp.push_back(triplet(neighbors[n]->triplets[c].v3, neighbors[n]->triplets[c].v1, neighbors[n]->triplets[c].v2));
            
            if(temp.size()!=4)
            {
                std::cout << "Something wrong in edge refining on step " << step_count << "\n";
                std::cout << "Attempting to compute 'topology' of node " << id << '\n';
                std::cout << "Size of temp is " << temp.size() << '\n';
                
                neighbors[n]->output();
                
                nodes[3322342342342]->output();
                exit(0);
            }
            
            for(unsigned int sort=0; sort<3; sort++) {
                for(unsigned int sort2=sort+1; sort2<4; sort2+=1){
                    if(temp[sort].v3==temp[sort2].v1) {
                        temp1=temp[sort+1].v1;
                        temp2=temp[sort+1].v2;
                        temp3=temp[sort+1].v3;
                        temp[sort+1].v1 = temp[sort2].v1;
                        temp[sort+1].v2 = temp[sort2].v2;
                        temp[sort+1].v3 = temp[sort2].v3;
                        temp[sort2].v1 = temp1;
                        temp[sort2].v2 = temp2;
                        temp[sort2].v3 = temp3;
                    }
                }
            }
            
            for(unsigned int c=0; c<4; c++) corners[n].push_back(temp[c].v1);
        }
        
        return;
    }
    
    ////////////////////////////////////////////////////////////////////
    // Topology of vertex nodes
    if(neighbors.size()==4)
    {
        CNode* temp1;
        CNode* temp2;
        CNode* temp3;
        
        for(unsigned int n=0; n<4; n++)
        {
            corners[n].clear();
            std::vector<triplet> temp;
            
            // finds the ordered, oriented edges associated with this node, for this body
            for(unsigned int c=1; c<neighbors[n]->triplets.size(); c++)
                if(neighbors[n]->triplets[c].v3==this)
                    temp.push_back(triplet(neighbors[n]->triplets[c].v2, neighbors[n]->triplets[c].v3, neighbors[n]->triplets[c].v1));
            
            for(unsigned int c=1; c<neighbors[n]->triplets.size(); c++)
                if(neighbors[n]->triplets[c].v1==this)
                    temp.push_back(triplet(neighbors[n]->triplets[c].v3, neighbors[n]->triplets[c].v1, neighbors[n]->triplets[c].v2));
            
            for(unsigned int sort=0; sort<5; sort++) {
                for(unsigned int sort2=sort+1; sort2<6; sort2+=1){
                    if(temp[sort].v3==temp[sort2].v1) {
                        temp1=temp[sort+1].v1;
                        temp2=temp[sort+1].v2;
                        temp3=temp[sort+1].v3;
                        temp[sort+1].v1 = temp[sort2].v1;
                        temp[sort+1].v2 = temp[sort2].v2;
                        temp[sort+1].v3 = temp[sort2].v3;
                        temp[sort2].v1 = temp1;
                        temp[sort2].v2 = temp2;
                        temp[sort2].v3 = temp3;
                    }
                }
            }
            
            if(temp[3].v3!=temp[4].v1)
            {
                triplet t1(temp[1].v1, temp[1].v2, temp[1].v3);
                triplet t2(temp[2].v1, temp[2].v2, temp[2].v3);
                triplet t3(temp[3].v1, temp[3].v2, temp[3].v3);
                triplet t4(temp[4].v1, temp[4].v2, temp[4].v3);
                triplet t5(temp[5].v1, temp[5].v2, temp[5].v3);
                
                temp[1].v1 = t5.v1; temp[1].v2 = t5.v2; temp[1].v3 = t5.v3;
                temp[2].v1 = t4.v1; temp[2].v2 = t4.v2; temp[2].v3 = t4.v3;
                temp[3].v1 = t1.v1; temp[3].v2 = t1.v2; temp[3].v3 = t1.v3;
                temp[4].v1 = t2.v1; temp[4].v2 = t2.v2; temp[4].v3 = t2.v3;
                temp[5].v1 = t3.v1; temp[5].v2 = t3.v2; temp[5].v3 = t3.v3;
            }
            
            for(unsigned int c=0; c<6; c++) corners[n].push_back(temp[c].v1);
        }
    }
}
// Calculating the topological information of each node
////////////////////////////////////////////////////////////////////



int CNode::get_face_sides()             // maybe change this?  oh this could be a problem
{                                       // when we're erasing stuff from a face and need
    int counter=0;
    for(unsigned int c=0; c<corners[0].size(); c++)
        if(corners[0][c]->neighbors.size()==4) counter++;
    return counter;
}


////////////////////////////////////////////////////////////////////
// Calculates motion for one node
void CNode::calc_motion()
{
    if(neighbors.size()==0) return;
    good = 1;
    
    double emin = 1.;
    
    ////////////////////////////////////////////////////////////////////
    // Motions for face nodes
    if(neighbors.size()==2)
    {
        if(corners[0].size()==2) good=0;
        
        double dvx=0;
        double dvy=0;
        double dvz=0;
        double local_curvature=0;
        
        // Calculates how much to move a face node given local curvature
        unsigned int size = corners[0].size();
        
        double edges[size][3];
        for(unsigned int c=0; c<size; c++)
        {
            edges[c][0] = corners[0][c]->x - x;
            edges[c][1] = corners[0][c]->y - y;
            edges[c][2] = corners[0][c]->z - z;
            
            if(edges[c][0] > 0.5) edges[c][0] -= 1.0;
            if(edges[c][0] <-0.5) edges[c][0] += 1.0;
            if(edges[c][1] > 0.5) edges[c][1] -= 1.0;
            if(edges[c][1] <-0.5) edges[c][1] += 1.0;
            if(edges[c][2] > 0.5) edges[c][2] -= 1.0;
            if(edges[c][2] <-0.5) edges[c][2] += 1.0;
        }
        
        // HERE WE CALCULATE ALL THE FACE NORMALS AND ALSO ALL OF THEIR NORMS (FACE AREAS)
        double normals[size][3];
        double norms  [size];
        
        for(unsigned int c=0; c<size-1; c++)
        {
            normals[c][0] = edges[c][1]*edges[c+1][2] - edges[c][2]*edges[c+1][1];
            normals[c][1] = edges[c][2]*edges[c+1][0] - edges[c][0]*edges[c+1][2];
            normals[c][2] = edges[c][0]*edges[c+1][1] - edges[c][1]*edges[c+1][0];
            norms[c] = sqrt(normals[c][0]*normals[c][0]+normals[c][1]*normals[c][1]+normals[c][2]*normals[c][2]);
            
            dvx += (edges[c][1]*edges[c+1][2] - edges[c][2]*edges[c+1][1]);
            dvy += (edges[c][2]*edges[c+1][0] - edges[c][0]*edges[c+1][2]);
            dvz += (edges[c][0]*edges[c+1][1] - edges[c][1]*edges[c+1][0]);
            
            if(norms[c]==0) good=0;
        }
        
        normals[size-1][0] = edges[size-1][1]*edges[0][2] - edges[size-1][2]*edges[0][1];
        normals[size-1][1] = edges[size-1][2]*edges[0][0] - edges[size-1][0]*edges[0][2];
        normals[size-1][2] = edges[size-1][0]*edges[0][1] - edges[size-1][1]*edges[0][0];
        norms[size-1] = sqrt(normals[size-1][0]*normals[size-1][0]+normals[size-1][1]*normals[size-1][1]+normals[size-1][2]*normals[size-1][2]);
        
        dvx += (edges[size-1][1]*edges[0][2] - edges[size-1][2]*edges[0][1]);
        dvy += (edges[size-1][2]*edges[0][0] - edges[size-1][0]*edges[0][2]);
        dvz += (edges[size-1][0]*edges[0][1] - edges[size-1][1]*edges[0][0]);
        
        if(norms[size-1]==0) good=0;
        
        dvx /= 6.;
        dvy /= 6.;
        dvz /= 6.;
        // FINISHED CALCULATING NORMALS AND FACE AREAS
        
        // NOW THAT WE HAVE THE FACE NORMALS AND EDGE LENGTHS CALCULATED, WE CAN EASILY CALCULATE THE LOCAL MEAN CURVATURE HERE
        for(unsigned int c=0; c<size-1; c++)
        {
            double inside = (normals[c][0]*normals[c+1][0]+normals[c][1]*normals[c+1][1]+normals[c][2]*normals[c+1][2])/norms[c]/norms[c+1];
            if(inside> 1.) inside = 1.000; // THIS IS TO AVOID MACHINE ERROR PROBLEMS
            if(inside<-1.) inside =-1.000; // THIS IS TO AVOID MACHINE ERROR PROBLEMS
            
            double angle = acos(inside);
            if(edges[c+1][0]*(normals[c][1]*normals[c+1][2]-normals[c][2]*normals[c+1][1]) +
               edges[c+1][1]*(normals[c][2]*normals[c+1][0]-normals[c][0]*normals[c+1][2]) +
               edges[c+1][2]*(normals[c][0]*normals[c+1][1]-normals[c][1]*normals[c+1][0]) < 0.0) angle *= -1;
            
            double elength = sqrt(edges[c+1][0]*edges[c+1][0]+edges[c+1][1]*edges[c+1][1]+edges[c+1][2]*edges[c+1][2]);
            if(elength<emin) emin = elength;
            
            local_curvature += angle*elength;
        }
        
        double inside = (normals[size-1][0]*normals[0][0]+normals[size-1][1]*normals[0][1]+normals[size-1][2]*normals[0][2])/norms[size-1]/norms[0];
        if(inside> 1.) inside = 1.000; // THIS IS TO AVOID MACHINE ERROR PROBLEMS
        if(inside<-1.) inside =-1.000; // THIS IS TO AVOID MACHINE ERROR PROBLEMS
        
        double angle = acos(inside);
        if(edges[0][0]*(normals[size-1][1]*normals[0][2]-normals[size-1][2]*normals[0][1]) +
           edges[0][1]*(normals[size-1][2]*normals[0][0]-normals[size-1][0]*normals[0][2]) +
           edges[0][2]*(normals[size-1][0]*normals[0][1]-normals[size-1][1]*normals[0][0]) < 0.0) angle *= -1;
        
        double elength = sqrt(edges[0][0]*edges[0][0]+edges[0][1]*edges[0][1]+edges[0][2]*edges[0][2]);
        if(elength<emin) emin = elength;
        
        local_curvature += angle*elength;
        local_curvature /= 2.;
        // DONE CALCULATING LOCAL MEAN CURVATURE
        
        
        if(dvx*dvx+dvy*dvy+dvz*dvz==0) good=0;
        double mult = local_curvature/(dvx*dvx+dvy*dvy+dvz*dvz);
        
        dx = mult*dvx;
        dy = mult*dvy;
        dz = mult*dvz;
        
        if(dx*dx+dy*dy+dz*dz > pow(0.1*emin/time_scale, 2) || dx!=dx) good=0;
        
        if(good==0)
        {
#ifdef DEBUG
            std::cout<<"Bad face node " << id << ", will center it\n";
#endif
            dx=dy=dz=0.;
            
            double e1x=0, e1y=0, e1z=0;
            double e2x=0, e2y=0, e2z=0;
            
            for(unsigned int c=0; c<corners[0].size(); c++)
            {
                e1x = corners[0][c]->x - x;
                e1y = corners[0][c]->y - y;
                e1z = corners[0][c]->z - z;
                
                if(e1x > 0.5) e1x -= 1.0;
                if(e1x <-0.5) e1x += 1.0;
                if(e1y > 0.5) e1y -= 1.0;
                if(e1y <-0.5) e1y += 1.0;
                if(e1z > 0.5) e1z -= 1.0;
                if(e1z <-0.5) e1z += 1.0;
                
                e2x += e1x;
                e2y += e1y;
                e2z += e1z;
            }
            
            cx = e2x/corners[0].size();
            cy = e2y/corners[0].size();
            cz = e2z/corners[0].size();
        }
        
        if(dx!=dx)
        {
            std::cout.precision(18);
            std::cout << "Good is " << good << '\n';
            std::cout << "mult is " << mult << '\n';
            std::cout << "emin is " << emin << '\n';
            std::cout << "dx is " << dx << ", " << dy << ", " << dz << '\n';
            std::cout << "dvx is " << dvx << ", " << dvy << ", " << dvz << '\n';
            std::cout << "timesc " << time_scale << '\n';
            std::cout << "localc " << local_curvature << '\n';
            for(unsigned int c=0; c<size; c++)
                std::cout << c << '\n'
                << edges[c][0] << '\t' << edges[c][1] << '\t' << edges[c][2] << '\n'
                << normals[c][0] << '\t' << normals[c][1] << '\t' << normals[c][2] << '\n'
                << norms[c] << '\n';
            
        }
    }
    // Motions for face nodes
    ////////////////////////////////////////////////////////////////////
    
    
    ////////////////////////////////////////////////////////////////////
    // Motions for edge nodes
    else if(neighbors.size()==3)
    {
        double local_curvature[4]={0,0,0};
        double dvol[4][4]={{0, 0,0,0},{0, 0,0,0},{0,0,0,0}};
        
        for(unsigned int n=0; n<3; n++)
        {
            double edges[4][4];
            for(unsigned int c=0; c<4; c++)
            {
                edges[c][0] = corners[n][c]->x - x;
                edges[c][1] = corners[n][c]->y - y;
                edges[c][2] = corners[n][c]->z - z;
                if(edges[c][0] > 0.5) edges[c][0] -= 1.0;
                if(edges[c][0] <-0.5) edges[c][0] += 1.0;
                if(edges[c][1] > 0.5) edges[c][1] -= 1.0;
                if(edges[c][1] <-0.5) edges[c][1] += 1.0;
                if(edges[c][2] > 0.5) edges[c][2] -= 1.0;
                if(edges[c][2] <-0.5) edges[c][2] += 1.0;
            }
            
            double normals[4][3];
            double norms[4];
            
            for(unsigned int c=0; c<3; c++)
            {
                normals[c][0] = edges[c][1]*edges[c+1][2] - edges[c][2]*edges[c+1][1];
                normals[c][1] = edges[c][2]*edges[c+1][0] - edges[c][0]*edges[c+1][2];
                normals[c][2] = edges[c][0]*edges[c+1][1] - edges[c][1]*edges[c+1][0];
                norms[c] = sqrt(normals[c][0]*normals[c][0]+normals[c][1]*normals[c][1]+normals[c][2]*normals[c][2]);
                
                dvol[n][0] += (edges[c][1]*edges[c+1][2] - edges[c][2]*edges[c+1][1]);
                dvol[n][1] += (edges[c][2]*edges[c+1][0] - edges[c][0]*edges[c+1][2]);
                dvol[n][2] += (edges[c][0]*edges[c+1][1] - edges[c][1]*edges[c+1][0]);
                
                if(norms[c]==0) good=0;
            }
            
            normals[3][0] = edges[3][1]*edges[0][2] - edges[3][2]*edges[0][1];
            normals[3][1] = edges[3][2]*edges[0][0] - edges[3][0]*edges[0][2];
            normals[3][2] = edges[3][0]*edges[0][1] - edges[3][1]*edges[0][0];
            norms[3] = sqrt(normals[3][0]*normals[3][0]+normals[3][1]*normals[3][1]+normals[3][2]*normals[3][2]);
            
            dvol[n][0] += (edges[3][1]*edges[0][2] - edges[3][2]*edges[0][1]);
            dvol[n][1] += (edges[3][2]*edges[0][0] - edges[3][0]*edges[0][2]);
            dvol[n][2] += (edges[3][0]*edges[0][1] - edges[3][1]*edges[0][0]);
            
            if(norms[3]==0) good=0;
            
            dvol[n][0] /= 6.;
            dvol[n][1] /= 6.;
            dvol[n][2] /= 6.;
            
            for(unsigned int c=0; c<3; c++)
            {
                double inside = (normals[c][0]*normals[c+1][0]+normals[c][1]*normals[c+1][1]+normals[c][2]*normals[c+1][2])/norms[c]/norms[c+1];
                if(inside> 1.) inside = 1.000; // THIS IS TO AVOID MACHINE ERROR PROBLEMS
                if(inside<-1.) inside =-1.000; // THIS IS TO AVOID MACHINE ERROR PROBLEMS
                
                double angle = acos(inside);
                if(edges[c+1][0]*(normals[c][1]*normals[c+1][2]-normals[c][2]*normals[c+1][1]) +
                   edges[c+1][1]*(normals[c][2]*normals[c+1][0]-normals[c][0]*normals[c+1][2]) +
                   edges[c+1][2]*(normals[c][0]*normals[c+1][1]-normals[c][1]*normals[c+1][0]) < 0.0) angle *= -1;
                
                double elength = sqrt(edges[c+1][0]*edges[c+1][0]+edges[c+1][1]*edges[c+1][1]+edges[c+1][2]*edges[c+1][2]);
                if(elength<emin) emin = elength;
                
                if (c==1) local_curvature[n] +=  angle        *elength;
                else      local_curvature[n] += (angle-pi/3.0)*elength;
            }
            
            double inside = (normals[3][0]*normals[0][0]+normals[3][1]*normals[0][1]+normals[3][2]*normals[0][2])/norms[3]/norms[0];
            if(inside> 1.) inside = 1.000; // THIS IS TO AVOID MACHINE ERROR PROBLEMS
            if(inside<-1.) inside =-1.000; // THIS IS TO AVOID MACHINE ERROR PROBLEMS
            double angle = acos(inside);
            if(edges[0][0]*(normals[3][1]*normals[0][2]-normals[3][2]*normals[0][1]) +
               edges[0][1]*(normals[3][2]*normals[0][0]-normals[3][0]*normals[0][2]) +
               edges[0][2]*(normals[3][0]*normals[0][1]-normals[3][1]*normals[0][0]) < 0.0) angle *= -1;
            double elength = sqrt(edges[0][0]*edges[0][0]+edges[0][1]*edges[0][1]+edges[0][2]*edges[0][2]);
            if(elength<emin) emin = elength;
            local_curvature[n] +=  angle*elength;
            
            local_curvature[n] /= 2.;
        }
        
        if(fabs(local_curvature[0] + local_curvature[1] + local_curvature[2])>0.000001) { good=0; offcounter2++; }
        
        double e1x = corners[0][1]->x - corners[0][3]->x;
        double e1y = corners[0][1]->y - corners[0][3]->y;
        double e1z = corners[0][1]->z - corners[0][3]->z;
        
        if(e1x > 0.5) e1x -= 1.0;
        if(e1x <-0.5) e1x += 1.0;
        if(e1y > 0.5) e1y -= 1.0;
        if(e1y <-0.5) e1y += 1.0;
        if(e1z > 0.5) e1z -= 1.0;
        if(e1z <-0.5) e1z += 1.0;
        
        double  det = dvol[0][0]*dvol[1][2]*e1y - dvol[0][0]*dvol[1][1]*e1z +
        dvol[0][1]*dvol[1][0]*e1z - dvol[0][1]*dvol[1][2]*e1x +
        dvol[0][2]*dvol[1][1]*e1x - dvol[0][2]*dvol[1][0]*e1y;
        
        dx = -((dvol[1][1]*e1z-dvol[1][2]*e1y)*local_curvature[0]
               + (dvol[0][2]*e1y-dvol[0][1]*e1z)*local_curvature[1])/det;
        dy = -((dvol[1][2]*e1x-dvol[1][0]*e1z)*local_curvature[0]
               + (dvol[0][0]*e1z-dvol[0][2]*e1x)*local_curvature[1])/det;
        dz = -((dvol[1][0]*e1y-dvol[1][1]*e1x)*local_curvature[0]
               + (dvol[0][1]*e1x-dvol[0][0]*e1y)*local_curvature[1])/det;
        
        if(dx*dx+dy*dy+dz*dz > pow(0.1*emin/time_scale, 2)) good=0;
        
        if(good==0)
        {
#ifdef DEBUG
            std::cout<<"Bad edge node " << id << ", moving on\n";
#endif
            dx=dy=dz=0.;
            
            double c1x = corners[0][1]->x - x;
            double c1y = corners[0][1]->y - y;
            double c1z = corners[0][1]->z - z;
            double c2x = corners[0][3]->x - x;
            double c2y = corners[0][3]->y - y;
            double c2z = corners[0][3]->z - z;
            
            if(c1x > 0.5) c1x -= 1.0;
            if(c1x <-0.5) c1x += 1.0;
            if(c1y > 0.5) c1y -= 1.0;
            if(c1y <-0.5) c1y += 1.0;
            if(c1z > 0.5) c1z -= 1.0;
            if(c1z <-0.5) c1z += 1.0;
            if(c2x > 0.5) c2x -= 1.0;
            if(c2x <-0.5) c2x += 1.0;
            if(c2y > 0.5) c2y -= 1.0;
            if(c2y <-0.5) c2y += 1.0;
            if(c2z > 0.5) c2z -= 1.0;
            if(c2z <-0.5) c2z += 1.0;
            
            cx = (c1x+c2x)/2.;
            cy = (c1y+c2y)/2.;
            cz = (c1z+c2z)/2.;
        }
    }
    // Motions for edge nodes; at this point we don't have any.
    ////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////
    // Motions for vertex nodes
    else if(neighbors.size()==4)
    {
        double dvol[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
        double local_curvature[4]={0,0,0,0};
        
        double edges[6][4];
        for(unsigned int n=0; n<4; n++)
        {
            for(unsigned int c=0; c<6; c++)
            {
                edges[c][0] = corners[n][c]->x - x;
                edges[c][1] = corners[n][c]->y - y;
                edges[c][2] = corners[n][c]->z - z;
                
                if(edges[c][0] > 0.5) edges[c][0] -= 1.0;
                if(edges[c][0] <-0.5) edges[c][0] += 1.0;
                if(edges[c][1] > 0.5) edges[c][1] -= 1.0;
                if(edges[c][1] <-0.5) edges[c][1] += 1.0;
                if(edges[c][2] > 0.5) edges[c][2] -= 1.0;
                if(edges[c][2] <-0.5) edges[c][2] += 1.0;
            }
            
            double normals[6][4];
            double norms[6];
            
            for(unsigned int c=0; c<5; c++)
            {
                normals[c][0] = edges[c][1]*edges[c+1][2] - edges[c][2]*edges[c+1][1];
                normals[c][1] = edges[c][2]*edges[c+1][0] - edges[c][0]*edges[c+1][2];
                normals[c][2] = edges[c][0]*edges[c+1][1] - edges[c][1]*edges[c+1][0];
                norms[c] = sqrt(normals[c][0]*normals[c][0]+normals[c][1]*normals[c][1]+normals[c][2]*normals[c][2]);
                
                dvol[n][0] += (edges[c][1]*edges[c+1][2] - edges[c][2]*edges[c+1][1]);
                dvol[n][1] += (edges[c][2]*edges[c+1][0] - edges[c][0]*edges[c+1][2]);
                dvol[n][2] += (edges[c][0]*edges[c+1][1] - edges[c][1]*edges[c+1][0]);
            }
            
            normals[5][0] = edges[5][1]*edges[0][2] - edges[5][2]*edges[0][1];
            normals[5][1] = edges[5][2]*edges[0][0] - edges[5][0]*edges[0][2];
            normals[5][2] = edges[5][0]*edges[0][1] - edges[5][1]*edges[0][0];
            norms[5] = sqrt(normals[5][0]*normals[5][0]+normals[5][1]*normals[5][1]+normals[5][2]*normals[5][2]);
            
            dvol[n][0] += (edges[5][1]*edges[0][2] - edges[5][2]*edges[0][1]);
            dvol[n][1] += (edges[5][2]*edges[0][0] - edges[5][0]*edges[0][2]);
            dvol[n][2] += (edges[5][0]*edges[0][1] - edges[5][1]*edges[0][0]);
            
            
            if(norms[0]==0 || norms[1]==0 || norms[2]==0 || norms[3]==0 || norms[4]==0 || norms[5]==0) good=0;
            
            dvol[n][0] /= 6.;
            dvol[n][1] /= 6.;
            dvol[n][2] /= 6.;
            
            
            for(unsigned int c=0; c<5; c++)
            {
                double inside = (normals[c][0]*normals[c+1][0]+normals[c][1]*normals[c+1][1]+normals[c][2]*normals[c+1][2])/norms[c]/norms[c+1];
                if(inside> 1.) inside = 1.000; // THIS IS TO AVOID MACHINE ERROR PROBLEMS
                if(inside<-1.) inside =-1.000; // THIS IS TO AVOID MACHINE ERROR PROBLEMS
                
                double angle = acos(inside);
                if(edges[c+1][0]*(normals[c][1]*normals[c+1][2]-normals[c][2]*normals[c+1][1]) +
                   edges[c+1][1]*(normals[c][2]*normals[c+1][0]-normals[c][0]*normals[c+1][2]) +
                   edges[c+1][2]*(normals[c][0]*normals[c+1][1]-normals[c][1]*normals[c+1][0]) < 0.0) angle *= -1.;
                
                double elength = sqrt(edges[c+1][0]*edges[c+1][0]+edges[c+1][1]*edges[c+1][1]+edges[c+1][2]*edges[c+1][2]);
                if(elength<emin) emin = elength;
                
                if (c==1 || c==3) local_curvature[n] +=  angle        *elength;
                else              local_curvature[n] += (angle-pi/3.0)*elength;
            }
            
            double inside = (normals[5][0]*normals[0][0]+normals[5][1]*normals[0][1]+normals[5][2]*normals[0][2])/norms[5]/norms[0];
            if(inside> 1.) inside = 1.000; // THIS IS TO AVOID MACHINE ERROR PROBLEMS
            if(inside<-1.) inside =-1.000; // THIS IS TO AVOID MACHINE ERROR PROBLEMS
            double angle = acos(inside);
            if(edges[0][0]*(normals[5][1]*normals[0][2]-normals[5][2]*normals[0][1]) +
               edges[0][1]*(normals[5][2]*normals[0][0]-normals[5][0]*normals[0][2]) +
               edges[0][2]*(normals[5][0]*normals[0][1]-normals[5][1]*normals[0][0]) < 0.0) angle *= -1.;
            double elength = sqrt(edges[0][0]*edges[0][0]+edges[0][1]*edges[0][1]+edges[0][2]*edges[0][2]);
            if(elength<emin) emin = elength;
            local_curvature[n] +=  angle        *elength;
            
            local_curvature[n] /= 2.;
        }
        
        
        if(fabs(local_curvature[0] + local_curvature[1] + local_curvature[2] + local_curvature[3])>0.000001) { good=0; offcounter3++; }
        if(local_curvature[0]!=local_curvature[0]||local_curvature[1]!=local_curvature[1]||
           local_curvature[2]!=local_curvature[2]||local_curvature[3]!=local_curvature[3]) good=0;
        
        // Now we have calculated curvatures for all 4 adjacent bodies of this vertex.
        // We have also calculated how volume changes with unit motions in directions x, y, and z;
        // We invert a matrix and solve a system of linear equations.
        double  det = dvol[0][0]*dvol[1][2]*dvol[2][1] - dvol[0][0]*dvol[1][1]*dvol[2][2] +
        dvol[0][1]*dvol[1][0]*dvol[2][2] - dvol[0][1]*dvol[1][2]*dvol[2][0] +
        dvol[0][2]*dvol[1][1]*dvol[2][0] - dvol[0][2]*dvol[1][0]*dvol[2][1];
        
        dx = -((dvol[1][1]*dvol[2][2]-dvol[1][2]*dvol[2][1])*local_curvature[0]
               + (dvol[0][2]*dvol[2][1]-dvol[0][1]*dvol[2][2])*local_curvature[1]
               + (dvol[0][1]*dvol[1][2]-dvol[0][2]*dvol[1][1])*local_curvature[2])/det;
        dy = -((dvol[1][2]*dvol[2][0]-dvol[1][0]*dvol[2][2])*local_curvature[0]
               + (dvol[0][0]*dvol[2][2]-dvol[0][2]*dvol[2][0])*local_curvature[1]
               + (dvol[0][2]*dvol[1][0]-dvol[0][0]*dvol[1][2])*local_curvature[2])/det;
        dz = -((dvol[1][0]*dvol[2][1]-dvol[1][1]*dvol[2][0])*local_curvature[0]
               + (dvol[0][1]*dvol[2][0]-dvol[0][0]*dvol[2][1])*local_curvature[1]
               + (dvol[0][0]*dvol[1][1]-dvol[0][1]*dvol[1][0])*local_curvature[2])/det;
        
        if(dx*dx+dy*dy+dz*dz > pow(0.1*emin/time_scale, 2)) good=0;
        
        if(good==0)
        {
#ifdef DEBUG
            std::cout<<"Bad vertex " << id << ", moving on\n";
#endif
            dx=dy=dz=0.;
            cx=cy=cz=0.;
            
            double minv = 1.;
            
            for(unsigned int c=0; c<4; c++)
            {
                for(unsigned int d=1; d<6; d+=2)
                {
                    double ex = corners[c][d]->x - x;
                    double ey = corners[c][d]->y - y;
                    double ez = corners[c][d]->z - z;
                    
                    if(ex > 0.5) ex -= 1.0;
                    if(ex <-0.5) ex += 1.0;
                    if(ey > 0.5) ey -= 1.0;
                    if(ey <-0.5) ey += 1.0;
                    if(ez > 0.5) ez -= 1.0;
                    if(ez <-0.5) ez+= 1.0;
                    
                    double norm = sqrt(ex*ex+ey*ey+ez*ez);
                    if(norm < minv) minv = norm;
                    
                    cx += ex/12./norm;
                    cy += ey/12./norm;
                    cz += ez/12./norm;
                }
            }
            
            cx*=minv;
            cy*=minv;
            cz*=minv;
        }
    }
    // Motions for vertex nodes
    ////////////////////////////////////////////////////////////////////
    
    // A LAST CHECK, JUST IN CASE WE STILL HAVE A MESSED UP MOTION
     if(dx!=dx || dy!=dy || dz!=dz)
     {
         output();
         for(unsigned int c=0; c<neighbors.size(); c++) neighbors[c]->evolver_out();
         std::cout << "crashing in cnode: calc_motion";
         crash();
     }
}
// All motions are now calculated...
////////////////////////////////////////////////////////////////////



// Moves the node to the center of its adjacent edges, not weighted in any way

void CNode::center2(double amount)
{
    //double amount is a variable, should be in [0,1], to tell us 'how much' to center
    // the node.  Sometimes we might want to fully center it, sometimes only
    // partly
    if(neighbors.size()==2)
    {
        double e1x=0, e1y=0, e1z=0;
        double e2x=0, e2y=0, e2z=0;
        
        for(unsigned int c=0; c<corners[0].size(); c++)
        {
            e1x = corners[0][c]->x - x;
            e1y = corners[0][c]->y - y;
            e1z = corners[0][c]->z - z;
            
            if(e1x > 0.5) e1x -= 1.0;
            if(e1x <-0.5) e1x += 1.0;
            if(e1y > 0.5) e1y -= 1.0;
            if(e1y <-0.5) e1y += 1.0;
            if(e1z > 0.5) e1z -= 1.0;
            if(e1z <-0.5) e1z += 1.0;
            
            e2x += e1x;
            e2y += e1y;
            e2z += e1z;
        }
        
        x += amount*e2x/corners[0].size();
        y += amount*e2y/corners[0].size();
        z += amount*e2z/corners[0].size();
    }
    
    else if(neighbors.size()==3)
    {
        double e1x, e1y, e1z;
        double e2x, e2y, e2z;
        
        e1x = corners[0][1]->x - x;
        e1y = corners[0][1]->y - y;
        e1z = corners[0][1]->z - z;
        e2x = corners[0][3]->x - x;
        e2y = corners[0][3]->y - y;
        e2z = corners[0][3]->z - z;
        
        if(e1x > 0.5) e1x -= 1.0;
        if(e1x <-0.5) e1x += 1.0;
        if(e1y > 0.5) e1y -= 1.0;
        if(e1y <-0.5) e1y += 1.0;
        if(e1z > 0.5) e1z -= 1.0;
        if(e1z <-0.5) e1z += 1.0;
        if(e2x > 0.5) e2x -= 1.0;
        if(e2x <-0.5) e2x += 1.0;
        if(e2y > 0.5) e2y -= 1.0;
        if(e2y <-0.5) e2y += 1.0;
        if(e2z > 0.5) e2z -= 1.0;
        if(e2z <-0.5) e2z += 1.0;
        
        x += amount*(e1x+e2x)/2.;
        y += amount*(e1y+e2y)/2.;
        z += amount*(e1z+e2z)/2.;
    }
    
    else if(neighbors.size()==4)
    {
        double e1x=0, e1y=0, e1z=0;
        double e2x=0, e2y=0, e2z=0;
        
        for(unsigned int c=0; c<4; c++)
            for(unsigned int d=0; d<6; d++)
                if(corners[c][d]->neighbors.size()!=2)
                {
                    e1x = corners[c][d]->x - x;
                    e1y = corners[c][d]->y - y;
                    e1z = corners[c][d]->z - z;
                    
                    if(e1x > 0.5) e1x -= 1.0;
                    if(e1x <-0.5) e1x += 1.0;
                    if(e1y > 0.5) e1y -= 1.0;
                    if(e1y <-0.5) e1y += 1.0;
                    if(e1z > 0.5) e1z -= 1.0;
                    if(e1z <-0.5) e1z += 1.0;
                    
                    e2x += e1x;
                    e2y += e1y;
                    e2z += e1z;
                }
        x += amount*e2x/12.;
        y += amount*e2y/12.;
        z += amount*e2z/12.;
    }
    
    center();
}




// THIS IS TO CALCULATE THE WAY IN WHICH THE AREA OF A FACE IS PROJECTED TO CHANGE
double CNode::calc_change()
{
    double initial_area=0;
    double projected_area=0;
    
    unsigned int size = corners[0].size();
    for(unsigned int c=0; c<size; c++)
    {
        unsigned int d = c+1; if(d==size) d=0;
        
        double e1x = corners[0][c]->x - x;
        double e1y = corners[0][c]->y - y;
        double e1z = corners[0][c]->z - z;
        double e2x = corners[0][d]->x - x;
        double e2y = corners[0][d]->y - y;
        double e2z = corners[0][d]->z - z;
        
        if(e1x > 0.5) e1x -= 1.0;
        if(e1x <-0.5) e1x += 1.0;
        if(e1y > 0.5) e1y -= 1.0;
        if(e1y <-0.5) e1y += 1.0;
        if(e1z > 0.5) e1z -= 1.0;
        if(e1z <-0.5) e1z += 1.0;
        if(e2x > 0.5) e2x -= 1.0;
        if(e2x <-0.5) e2x += 1.0;
        if(e2y > 0.5) e2y -= 1.0;
        if(e2y <-0.5) e2y += 1.0;
        if(e2z > 0.5) e2z -= 1.0;
        if(e2z <-0.5) e2z += 1.0;
        
        initial_area += sqrt((e1y*e2z-e1z*e2y)*(e1y*e2z-e1z*e2y)+(e1z*e2x-e1x*e2z)*(e1z*e2x-e1x*e2z)+(e1x*e2y - e1y*e2x)*(e1x*e2y - e1y*e2x))/2.0;
        
        if(corners[0][c]->good==1)
        {
            e1x += corners[0][c]->dx*time_scale;
            e1y += corners[0][c]->dy*time_scale;
            e1z += corners[0][c]->dz*time_scale;
        }
        else
        {
            e1x += corners[0][c]->cx;
            e1y += corners[0][c]->cy;
            e1z += corners[0][c]->cz;
        }
        if(corners[0][d]->good==1)
        {
            e2x += corners[0][d]->dx*time_scale;
            e2y += corners[0][d]->dy*time_scale;
            e2z += corners[0][d]->dz*time_scale;
        }
        else
        {
            e2x += corners[0][d]->cx;
            e2y += corners[0][d]->cy;
            e2z += corners[0][d]->cz;
        }
        
        projected_area += sqrt((e1y*e2z-e1z*e2y)*(e1y*e2z-e1z*e2y)+(e1z*e2x-e1x*e2z)*(e1z*e2x-e1x*e2z)+(e1x*e2y - e1y*e2x)*(e1x*e2y - e1y*e2x))/2.0;
    }
    
    return (projected_area-initial_area)/initial_area;
}



