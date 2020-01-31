//#define DEBUG

#include "cbody.h"
#include "cnode.h"
#include "functions.h"
#include "globals.h"

#include <algorithm>
#include <omp.h>

double det(CNode *v1, CNode *v2, CNode *v3, CNode *origin)
{
  double e1x = v1->x - origin->x;
  double e1y = v1->y - origin->y;
  double e1z = v1->z - origin->z;
  double e2x = v2->x - origin->x;
  double e2y = v2->y - origin->y;
  double e2z = v2->z - origin->z;
  double e3x = v3->x - origin->x;
  double e3y = v3->y - origin->y;
  double e3z = v3->z - origin->z;

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
  if (e3x > 0.5)
    e3x -= 1.0;
  if (e3x < -0.5)
    e3x += 1.0;
  if (e3y > 0.5)
    e3y -= 1.0;
  if (e3y < -0.5)
    e3y += 1.0;
  if (e3z > 0.5)
    e3z -= 1.0;
  if (e3z < -0.5)
    e3z += 1.0;

  return e1x * e2y * e3z + e2x * e3y * e1z + e3x * e1y * e2z - e1x * e3y * e2z - e2x * e1y * e3z - e3x * e2y * e1z;
}

double pdet(CNode *v1, CNode *v2, CNode *v3, CNode *origin)
{
  double e1x = v1->x - origin->x;
  double e1y = v1->y - origin->y;
  double e1z = v1->z - origin->z;
  double e2x = v2->x - origin->x;
  double e2y = v2->y - origin->y;
  double e2z = v2->z - origin->z;
  double e3x = v3->x - origin->x;
  double e3y = v3->y - origin->y;
  double e3z = v3->z - origin->z;

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
  if (e3x > 0.5)
    e3x -= 1.0;
  if (e3x < -0.5)
    e3x += 1.0;
  if (e3y > 0.5)
    e3y -= 1.0;
  if (e3y < -0.5)
    e3y += 1.0;
  if (e3z > 0.5)
    e3z -= 1.0;
  if (e3z < -0.5)
    e3z += 1.0;

  if (v1->good == 0)
  {
    e1x += v1->cx;
    e1y += v1->cy;
    e1z += v1->cz;
  }
  else
  {
    e1x += v1->dx * time_scale;
    e1y += v1->dy * time_scale;
    e1z += v1->dz * time_scale;
  }
  if (v2->good == 0)
  {
    e2x += v2->cx;
    e2y += v2->cy;
    e2z += v2->cz;
  }
  else
  {
    e2x += v2->dx * time_scale;
    e2y += v2->dy * time_scale;
    e2z += v2->dz * time_scale;
  }
  if (v3->good == 0)
  {
    e3x += v3->cx;
    e3y += v3->cy;
    e3z += v3->cz;
  }
  else
  {
    e3x += v3->dx * time_scale;
    e3y += v3->dy * time_scale;
    e3z += v3->dz * time_scale;
  }

  if (origin->good == 0)
  {
    e1x -= origin->cx;
    e2x -= origin->cx;
    e3x -= origin->cx;
    e1y -= origin->cy;
    e2y -= origin->cy;
    e3y -= origin->cy;
    e1z -= origin->cz;
    e2z -= origin->cz;
    e3z -= origin->cz;
  }
  else
  {
    e1x -= origin->dx * time_scale;
    e2x -= origin->dx * time_scale;
    e3x -= origin->dx * time_scale;
    e1y -= origin->dy * time_scale;
    e2y -= origin->dy * time_scale;
    e3y -= origin->dy * time_scale;
    e1z -= origin->dz * time_scale;
    e2z -= origin->dz * time_scale;
    e3z -= origin->dz * time_scale;
  }

  return e1x * e2y * e3z + e2x * e3y * e1z + e3x * e1y * e2z - e1x * e3y * e2z - e2x * e1y * e3z - e3x * e2y * e1z;
}

double length(CNode *one, CNode *two)
{
  double tempx = one->x - two->x;
  double tempy = one->y - two->y;
  double tempz = one->z - two->z;
  if (tempx > 0.5)
    tempx -= 1.0;
  if (tempx < -0.5)
    tempx += 1.0;
  if (tempy > 0.5)
    tempy -= 1.0;
  if (tempy < -0.5)
    tempy += 1.0;
  if (tempz > 0.5)
    tempz -= 1.0;
  if (tempz < -0.5)
    tempz += 1.0;
  return sqrt(tempx * tempx + tempy * tempy + tempz * tempz);
}

// this calculates and print out all the statistics we want
void calc_and_print_stats()
{

  /////////////////////////////////////////////////////
  // PART I PRINTS INFORMATION TO THE SCREEN

  if (step_count % 2000 == 0)
  {
    std::cout.precision(8);
    std::cout << std::endl

              << " "
              << "Steps   "
              << "Total time   "
              << "Clock  "
              << "Bodies  "

              << "<body V>   "
              << "<body A>   "
              << "<face A>  "
              << "<faces>    "
              << "<sides>  "

              << "Tetra  "
              << "Footballs  "
              << "3-gons "
              << "2-gons   "
              << "Edges     "
              << "\n";

    std::cout << "---------------------------------------------------------------------------------------------------------------------------------------------\n";
  }

  double avg_cell_area = calc_avg_cell_area();
  double avg_cell_faces = calc_avg_cell_faces();
  double avg_face_sides = calc_avg_face_sides();

  std::cout.precision(6);
  std::cout
      << std::setw(6) << step_count << "   "
      << std::fixed << total_time << "   "
      << std::setw(6) << time(0) - start << "   "
      << count_bodies() << "   "

      << 1. / count_bodies() << "   "
      << avg_cell_area << "   "
      << avg_face_area << "  "
      << avg_cell_faces << "  "
      << avg_face_sides << "  ";

  std::cout
      << std::setw(5) << tetrahedra_deleted << "   "
      << std::setw(5) << footballs_deleted << "     "
      << std::setw(5) << triangles_deleted << "  "
      << std::setw(5) << digons_deleted << "   "
      << std::setw(5) << edges_deleted << "   "
      << '\t'
      << std::setw(5) << edge_nodes_added << "   "
      << std::setw(5) << edge_nodes_deleted << "   "
      << std::setw(1) << offcounter1 << "-"
      << std::setw(1) << offcounter2 << "-"
      << std::setw(1) << offcounter3 << "   "
      << '\n';

  /*
    /////////////////////////////////////////////////////////////////////////////////
    // PART II PRINTS INFORMATION TO A DATA FILE.  WE PRINT MORE INFORMATION HERE
    //         THAN WE DO TO THE SCREEN.
    
    outfile1.open("output1.data", std::ios::app);
    outfile2.open("output2.data", std::ios::app);
    
    // FIRST WE NEED TO CALCULATE ALL INTERESTING STATISTICS
    
    outfile1.precision(15);
    outfile1
    << step_count            << '\t'
    << time_scale            << '\t'
    << total_time            << '\t'
    
    << time(0)-start         << '\t'
    << count_bodies()        << '\t'
    << 1./count_bodies()     << '\t'
    
    << avg_cell_area         << '\t'
    << avg_face_area         << '\t'
    << avg_cell_faces        << '\t'
    << avg_face_sides        << '\t'
    
    << tetrahedra_deleted    << '\t'
    << footballs_deleted     << '\t'
    << triangles_deleted     << '\t'
    << digons_deleted        << '\t'
    << edges_deleted         << '\t'
    
    << edge_nodes_added      << '\t'
    << edge_nodes_deleted    << '\t'
    << offcounter1           << '\t'
    << offcounter2           << '\t'
    << offcounter3           << '\t';
    
    //  << calc_shape_entropy() << '\t'
    
    // Let us try calculating the moments of the distribution, volume per cell
    double vmoments[5];
    double nums = double(count_bodies());
    for(int c=0; c<5; c++) vmoments[c]=0;
    for(unsigned int c=1; c<bodies.size(); c++)
    {
        for(int d=2; d<5; d++)
            vmoments[d] += pow(nums*bodies[c]->volume - 1, d);
    }
    for(int c=2; c<5; c++) vmoments[c] /= nums;
    double vstddev = sqrt(vmoments[2]);
    
    outfile1 << "VolumesMoments\t";
    outfile1 << vstddev << '\t';                // STDDEV
    outfile1 << vmoments[2] << '\t';            // VARIANCE
    outfile1 << vmoments[3] << '\t';            // 3RD CENTRAL MOMENT
    outfile1 << vmoments[4] << '\t';            // 4TH CENTRAL MOMENT
    outfile1 << vmoments[3]/vstddev/vstddev/vstddev << '\t';            // NOMALIZED 3RD CM (SKEWNESS)
    outfile1 << vmoments[4]/vstddev/vstddev/vstddev/vstddev << '\t';    // NOMALIZED 4TH CM
    
    
    int fnodes=0;
    int enodes=0;
    int vnodes=0;
    int znodes=0;
    for(unsigned int c=1; c<nodes.size(); c++)
    {
        if(nodes[c]->neighbors.size()==2)      fnodes++;
        else if(nodes[c]->neighbors.size()==3) enodes++;
        else if(nodes[c]->neighbors.size()==4) vnodes++;
        else znodes++;
    }
    outfile1 << "Nodes\t";
    outfile1 << vnodes << '\t';
    outfile1 << enodes << '\t';
    outfile1 << fnodes << '\t';
    outfile1 << znodes << '\t';
    
    // FACES PER CELL DISTRIBUTIONS
    for(unsigned int c=1; c<bodies.size(); c++) bodies[c]->fcount=0;
    int fnumbers[50];
    for(int c=1; c<50; c++)            fnumbers[c]=0;
    
    for(unsigned int c=1; c<nodes.size(); c++)
        if(nodes[c]->neighbors.size()==2)
        {
            nodes[c]->neighbors[0]->fcount++;
            nodes[c]->neighbors[1]->fcount++;
        }
    
    for(unsigned int c=1; c<bodies.size(); c++)
        fnumbers[bodies[c]->fcount]++;
    
    outfile1 << "FacesPerCell\t";
    for(int c=1; c<50; c++)            outfile1 << double(fnumbers[c])/double(bodies.size()-1) << '\t';
    
    double cell_faces_entropy=0;
    for(int c=1; c<50; c++) if(fnumbers[c]>0)
        cell_faces_entropy -= double(fnumbers[c])/double(bodies.size()-1.)*log(double(fnumbers[c])/double(bodies.size()-1.));
    outfile1 << cell_faces_entropy << '\t';
    
    
    // Let us try calculating the moments of the distribution, faces per cell
    double fmoments[5];
    for(int c=0; c<5; c++) fmoments[c]=0;
    for(unsigned int c=1; c<bodies.size(); c++)
    {
        for(int d=2; d<5; d++)
            fmoments[d] += pow(bodies[c]->fcount - avg_cell_faces, d);
    }
    for(int c=2; c<5; c++) fmoments[c] /= nums;
    double fstddev = sqrt(fmoments[2]);
    
    outfile1 << "FaceMoments\t";
    outfile1 << fstddev << '\t';
    outfile1 << fmoments[2] << '\t';
    outfile1 << fmoments[3] << '\t';
    outfile1 << fmoments[4] << '\t';
    outfile1 << fmoments[3]/fstddev/fstddev/fstddev << '\t';
    outfile1 << fmoments[4]/fstddev/fstddev/fstddev/fstddev << '\t';
    
    
    
    // VOLUME VS. FACE NUMBER
    double bvolumes[50];
    for(int c=1; c<50; c++)            bvolumes[c]=0;
    for(unsigned int c=1; c<bodies.size(); c++)
        bvolumes[bodies[c]->fcount] += bodies[c]->volume;
    outfile1 << "TotalVols\t";
    for(int c=1; c<50; c++)            outfile1 << bvolumes[c] << '\t';
    
    outfile1 << "AvgVols\t";
    for(int c=1; c<50; c++)  if(fnumbers[c]!=0)  outfile1 << bodies.size()*bvolumes[c]/double(fnumbers[c]) << '\t';
    else outfile1 << 0 << '\t';
    
    
    // EDGES PER FACE DISTRIBUTIONS
    int total_faces=0;
    int enumbers[50];
    double fareas[50];
    
    for(int c=1; c<50; c++)
    {
        enumbers[c]=0;
        fareas[c]=0;
    }
    
    for(unsigned int c=1; c<nodes.size(); c++)  if(nodes[c]->neighbors.size()==2)
    {
        total_faces++;
        int faces = nodes[c]->get_face_sides();
        enumbers[faces]++;
        fareas[faces] += nodes[c]->area;
    }
    
    outfile1 << "EdgesPerFace\t";
    for(int c=1; c<50; c++)            outfile1 << double(enumbers[c])/double(total_faces) << '\t';
    
    double face_edges_entropy=0;
    for(int c=1; c<50; c++) if(enumbers[c]>0)
        face_edges_entropy -= double(enumbers[c])/double(total_faces)*log(double(enumbers[c])/double(total_faces));
    outfile1 << face_edges_entropy << '\t';
    
    
    
    
    // Let us try calculating the moments of the distribution, edges per face
    double emoments[5];
    double ftotal=0;
    for(int c=0; c<5; c++) emoments[c]=0;
    for(unsigned int c=1; c<nodes.size(); c++) if(nodes[c]->neighbors.size()==2)
    {
        for(int d=2; d<5; d++)
            emoments[d] += pow(nodes[c]->get_face_sides() - avg_face_sides, d);
        ftotal += 1.0;
    }
    for(int c=2; c<5; c++) emoments[c] /= ftotal;
    double estddev = sqrt(emoments[2]);
    
    outfile1 << "EdgeMoments\t";
    outfile1 << estddev << '\t';
    outfile1 << emoments[2] << '\t';
    outfile1 << emoments[3] << '\t';
    outfile1 << emoments[4] << '\t';
    outfile1 << emoments[3]/estddev/estddev/estddev << '\t';
    outfile1 << emoments[4]/estddev/estddev/estddev/estddev << '\t';
    
    
    
    // AREA VS. EDGE NUMBER
    outfile1 << "TotalAreas\t";
    for(int c=1; c<50; c++)             outfile1 << fareas[c] << '\t';
    outfile1 << "AvgAreas\t";
    for(int c=1; c<50; c++)
        if(enumbers[c]!=0)                 outfile1 << fareas[c]/enumbers[c] << '\t';
        else                               outfile1 << 0.0 << '\t';
    
    
    // VOLUME DISTRIBUTIONS
    int volumes[100];
    int bodyc = bodies.size()-1;
    for(int c=0; c<100; c++)            volumes[c]=0;
    for(unsigned int c=1; c<bodies.size(); c++)
        volumes[int(bodies[c]->volume * bodyc * 10.)]++;
    outfile1 << "Volumes\t";
    for(int c=0; c<100; c++)            outfile1 << volumes[c]/double(bodies.size()-1) << '\t';
    
    
    // NOW WE ARE LOOKING AT SIXTUPLETS OF FACES AROUND EVERY VERTEX
    // Sum of face edges over vertices
    int vtotals[100];
    for(int c=0; c<100; c++) vtotals[c]=0;
    for(unsigned int c=1; c<nodes.size(); c++)
        if(nodes[c]->neighbors.size()==4)
            vtotals[common_face(nodes[c]->neighbors[0], nodes[c]->neighbors[1], nodes[c])->get_face_sides()
                    + common_face(nodes[c]->neighbors[0], nodes[c]->neighbors[2], nodes[c])->get_face_sides()
                    + common_face(nodes[c]->neighbors[0], nodes[c]->neighbors[3], nodes[c])->get_face_sides()
                    + common_face(nodes[c]->neighbors[1], nodes[c]->neighbors[2], nodes[c])->get_face_sides()
                    + common_face(nodes[c]->neighbors[1], nodes[c]->neighbors[3], nodes[c])->get_face_sides()
                    + common_face(nodes[c]->neighbors[2], nodes[c]->neighbors[3], nodes[c])->get_face_sides()]++;
    outfile1 << "vtotals\t";
    for(int c=20; c<50; c++)
        outfile1 << vtotals[c]/double(vnodes) << '\t';
    
    
    // Sum of cell faces over vertices
    int ftotals[100];
    for(int c=0; c<100; c++) ftotals[c]=0;
    for(unsigned int c=1; c<nodes.size(); c++)
        if(nodes[c]->neighbors.size()==4)
            ftotals[nodes[c]->neighbors[0]->fcount + nodes[c]->neighbors[1]->fcount + nodes[c]->neighbors[2]->fcount + nodes[c]->neighbors[3]->fcount]++;
    outfile1 << "ftotals\t";
    for(int c=20; c<100; c++)
        outfile1 << ftotals[c]/double(vnodes) << '\t';
    
    
    // WEAIRE-ABOAV 3D - FIRST NEIGHBORS
    int neighboringfaces[50];
    int neighboringfacesc[50];
    for(int c=1; c<50; c++) neighboringfaces[c] = neighboringfacesc[c] = 0;
    for(unsigned int c=1; c<nodes.size(); c++)
        if(nodes[c]->neighbors.size()==2)
        {
            neighboringfaces[nodes[c]->neighbors[0]->fcount] += nodes[c]->neighbors[1]->fcount;
            neighboringfaces[nodes[c]->neighbors[1]->fcount] += nodes[c]->neighbors[0]->fcount;
            neighboringfacesc[nodes[c]->neighbors[0]->fcount]++;
            neighboringfacesc[nodes[c]->neighbors[1]->fcount]++;
        }
    
    outfile1 << "WA3D\t";
    for(int c=1; c<50; c++)
        if(neighboringfacesc[c]!=0)
            outfile1 << double(neighboringfaces[c])/double(neighboringfacesc[c]) << '\t';
        else outfile1 << 0 << '\t';
    
    
    // WEAIRE-ABOAV 3D - SECOND NEIGHBORS
    for(int c=1; c<50; c++) neighboringfaces[c] = neighboringfacesc[c] = 0;
    for(unsigned int c=1; c<nodes.size(); c++)
    {
        if(nodes[c]->neighbors.size()==4)
        {
            std::vector <CNode*> nnodes;
            for(int d=0; d<4; d++)
                for(int e=1; e<6; e+=2)
                    nnodes.push_back(nodes[c]->corners[d][e]);
            sort(nnodes.begin(), nnodes.end(), compare);
            nnodes.erase(std::unique(nnodes.begin(), nnodes.end()), nnodes.end());
            
            for(unsigned int d=0; d<nnodes.size(); d++)
            {
                if(not_digon_neighbors(nnodes[d], nodes[c])==0) continue;   // HERE WE ARE SKIPPING A DIGON
                CNode* last = nodes[c];
                CNode* next = nnodes[d];
                
                while(next->neighbors.size()!=4)
                {
                    if(next->corners[0][1] == last)
                    {
                        last = next;
                        next = next->corners[0][3];
                    }
                    else
                    {
                        last = next;
                        next = next->corners[0][1];
                    }
                }
                
                // WE WILL HAVE A CRASH PROBLEM IF WE ARE ON THE EDGE OF A DIGON,
                // BECAUSE identify_body WILL CRASH.  HMM, WE NEED TO WRITE SOMETHING
                // TO CHECK IF ALL NEIGHBORS ARE THE SAME.  IF YES, THEN SKIP THIS.
                CBody* one=identify_body(nnodes[d], nodes[c]);
                CBody* two=identify_body(last, next);
                
                if(nodes[c]->id < next->id)
                {
                    neighboringfaces[one->fcount] += two->fcount;
                    neighboringfaces[two->fcount] += one->fcount;
                    neighboringfacesc[one->fcount]++;
                    neighboringfacesc[two->fcount]++;
                }
            }
        }
    }
    
    outfile1 << "WA3D2\t";
    for(int c=1; c<50; c++)
        if(neighboringfacesc[c]!=0)
            outfile1 << double(neighboringfaces[c])/double(neighboringfacesc[c]) << '\t';
        else outfile1 << 0 << '\t';
    
    
    // WEAIRE-ABOAV 3D - THIRD NEIGHBORS
    for(int c=1; c<50; c++) neighboringfaces[c] = neighboringfacesc[c] = 0;
    for(unsigned int c=1; c<nodes.size(); c++)
    {
        if(nodes[c]->neighbors.size()==4)
        {
            std::vector <CNode*> nnodes;
            for(int d=0; d<4; d++)
                for(int e=1; e<6; e+=2)
                    nnodes.push_back(nodes[c]->corners[d][e]);
            sort(nnodes.begin(), nnodes.end(), compare);
            nnodes.erase(std::unique(nnodes.begin(), nnodes.end()), nnodes.end());
            
            for(int d=0; d<3; d++)
            {
                for(int e=d+1; e<4; e++)
                {
                    
                    // FIRST WE GET THE D NODE, THEN WE'LL GET THE E NODE
                    CNode* last = nodes[c];
                    CNode* last2 =nodes[c];
                    
                    CNode* next = nnodes[d];
                    CNode* next2= nnodes[e];
                    
                    while(next->neighbors.size()!=4)
                    {
                        if(next->corners[0][1] == last)
                        {
                            last = next;
                            next = next->corners[0][3];
                        }
                        else
                        {
                            last = next;
                            next = next->corners[0][1];
                        }
                    }
                    
                    while(next2->neighbors.size()!=4)
                    {
                        if(next2->corners[0][1] == last2)
                        {
                            last2 = next2;
                            next2 = next2->corners[0][3];
                        }
                        else
                        {
                            last2 = next2;
                            next2 = next2->corners[0][1];
                        }
                    }
                    
                    if(not_digon_neighbors(last, next)==0||not_digon_neighbors(last2, next2)==0) continue;   // HERE WE ARE SKIPPING A DIGON
                    
                    CBody* one = identify_body(last, next);
                    CBody* two = identify_body(last2, next2);
                    
                    
                    if(next->id < next2->id && one->id != two->id)
                    {
                        neighboringfaces[one->fcount] += two->fcount;
                        neighboringfaces[two->fcount] += one->fcount;
                        neighboringfacesc[one->fcount]++;
                        neighboringfacesc[two->fcount]++;
                    }
                }
            }
        }
    }
    
    outfile1 << "WA3D3\t";
    for(int c=1; c<50; c++)
        if(neighboringfacesc[c]!=0)
            outfile1 << double(neighboringfaces[c])/double(neighboringfacesc[c]) << '\t';
        else outfile1 << 0 << '\t';
    
    
    
    if(step_count%1000==0)
    {
        int vtable2[32][32];    // THIS STORES NUMBERS OF PAIRS OF FACES
        for(int c=0; c<32; c++) for(int d=0; d<32; d++) vtable2[c][d]=0;
        
        int vtable[7][50000];   // THIS SHOULD STORE ALL SETS OF SIXTUPLETS
        int vhighest=0;
        
        for(int c=0; c<7; c++) for(int d=0; d<10000; d++) vtable[c][d]=0;
        for(unsigned int c=1; c<nodes.size(); c++)
        {
            if(nodes[c]->neighbors.size()==4)
            {
                int f1 = common_face(nodes[c]->neighbors[0], nodes[c]->neighbors[1], nodes[c])->get_face_sides();
                int f2 = common_face(nodes[c]->neighbors[2], nodes[c]->neighbors[3], nodes[c])->get_face_sides();
                int f3 = common_face(nodes[c]->neighbors[0], nodes[c]->neighbors[2], nodes[c])->get_face_sides();
                int f4 = common_face(nodes[c]->neighbors[1], nodes[c]->neighbors[3], nodes[c])->get_face_sides();
                int f5 = common_face(nodes[c]->neighbors[0], nodes[c]->neighbors[3], nodes[c])->get_face_sides();
                int f6 = common_face(nodes[c]->neighbors[1], nodes[c]->neighbors[2], nodes[c])->get_face_sides();
                
                int temp1;
                int temp2;
                if(f2<f1) { temp1=f1; f1=f2; f2=temp1; }
                if(f4<f3) { temp1=f3; f3=f4; f4=temp1; }
                if(f6<f5) { temp1=f5; f5=f6; f6=temp1; }
                
                if(f3<f1 || (f1==f3 && f4<f2)) { temp1=f1; temp2=f2; f1=f3; f2=f4; f3=temp1; f4=temp2; }
                if(f5<f1 || (f1==f5 && f6<f2)) { temp1=f1; temp2=f2; f1=f5; f2=f6; f5=temp1; f6=temp2; }
                if(f5<f3 || (f3==f5 && f6<f4)) { temp1=f3; temp2=f4; f1=f5; f2=f6; f5=temp1; f6=temp2; }
                // HERE WE PUT TOGETHER THE SIXTUPLET, THEN ADD IT TO SOME TABLE...
                vtable2[f1][f2]++;
                vtable2[f3][f4]++;
                vtable2[f5][f6]++;
                
                // FIRST WE CHECK IF IT'S ALREADY IN THE TABLE
                bool intable=0;
                for(int d=0; d<vhighest; d++)
                {
                    if(f1==vtable[1][d]&&
                       f2==vtable[2][d]&&
                       f3==vtable[3][d]&&
                       f4==vtable[4][d]&&
                       f5==vtable[5][d]&&
                       f6==vtable[6][d])
                    {
                        intable=1;
                        vtable[0][d]++;
                    }
                }
                if(intable==0)
                {
                    vtable[1][vhighest]=f1;
                    vtable[2][vhighest]=f2;
                    vtable[3][vhighest]=f3;
                    vtable[4][vhighest]=f4;
                    vtable[5][vhighest]=f5;
                    vtable[6][vhighest]=f6;
                    vtable[0][vhighest]++;
                    vhighest++;
                }
            }
        }
        
        for(int d=0; d<vhighest-1; d++)
            for(int e=d+1; e<vhighest; e++)
            {
                if(vtable[0][d]<vtable[0][e])
                {
                    int t0 = vtable[0][d]; vtable[0][d]=vtable[0][e]; vtable[0][e]=t0;
                    int t1 = vtable[1][d]; vtable[1][d]=vtable[1][e]; vtable[1][e]=t1;
                    int t2 = vtable[2][d]; vtable[2][d]=vtable[2][e]; vtable[2][e]=t2;
                    int t3 = vtable[3][d]; vtable[3][d]=vtable[3][e]; vtable[3][e]=t3;
                    int t4 = vtable[4][d]; vtable[4][d]=vtable[4][e]; vtable[4][e]=t4;
                    int t5 = vtable[5][d]; vtable[5][d]=vtable[5][e]; vtable[5][e]=t5;
                    int t6 = vtable[6][d]; vtable[6][d]=vtable[6][e]; vtable[6][e]=t6;
                }
            }
        
        outfile2 << step_count << " steps, " << total_time << " total time, " << vnodes << " vnodes, and " << vhighest << " types\n";
        outfile2.precision(4);
        for(int d=0; d<vhighest; d++)
            outfile2 << d << '\t' << std::setw(10) << vtable[0][d] << '\t' << '\t' << vtable[1][d] <<  '\t' << vtable[2][d] <<  '\t' << vtable[3][d] <<
            '\t' << vtable[4][d] <<  '\t' << vtable[5][d] <<  '\t' << vtable[6][d] << '\n';
        
        outfile2 << '\n';
        
        outfile2.precision(4);
        for(int c=2; c<15; c++)
        {
            for(int d=2; d<15; d++)
                outfile2 << std::setw(10) << vtable2[c][d]/double(vnodes*3.0) << "  ";
            outfile2 << '\n';
        }
        outfile2 << '\n'<<'\n';
    }
    
    outfile1 << '\n';                          // END OF LINE
    
    outfile1.close();
    outfile2.close();

    */

  // THESE ARE STATS WE WANT TO RESTART EVERY TIME STATS ARE PRINTED OUT
  toffcounter2 += offcounter2;
  toffcounter3 += offcounter3;
  offcounter1 = 0;
  offcounter2 = 0;
  offcounter3 = 0;
  edge_nodes_added = 0;
  edge_nodes_deleted = 0;
}

double calc_avg_face_area(void)
{
  double total = 0;
  int count = 0;

  for (unsigned int c = 1; c < nodes.size(); c++)
    if (nodes[c]->neighbors.size() == 2)
    {
      total += nodes[c]->area;
      count++;
    }

  avg_face_area = total / double(count);
  return avg_face_area;
}

double calc_avg_cell_area()
{
  double total = 0;
  for (unsigned int c = 1; c < nodes.size(); c++)
    total += nodes[c]->area;
  return 2. * double(total) / double(bodies.size() - 1);
}

double calc_avg_cell_faces()
{
  int total = 0;
  for (unsigned int c = 1; c < nodes.size(); c++)
    if (nodes[c]->neighbors.size() == 2)
      total++;
  return 2. * double(total) / double(bodies.size() - 1);
}

double calc_avg_face_sides()
{
  double total = 0;
  double count = 0;

  for (unsigned int c = 1; c < nodes.size(); c++)
    if (nodes[c]->neighbors.size() == 2)
    {
      total += nodes[c]->get_face_sides();
      count++;
    }

  return total / count;
}

////////////////////////////////////////////////////////////////////
// Read in data
void read_in_data(int argc, char *argv[])
{
  // if not enough arguments exit
  if (argc < 3)
  {
    std::cout << "usage: " << argv[0] << " <vertex_filename>"
              << " <param_filename>\n";
    exit(0);
  }

  std::ifstream paramFile(argv[2]);
  if (paramFile.is_open())
  {
    std::string line;
    while (getline(paramFile, line))
    {
      line.erase(std::remove_if(line.begin(), line.end(), isspace),
                 line.end());
      if (line[0] == '#' || line.empty())
        continue;
      auto delimiterPos = line.find("=");
      auto name = line.substr(0, delimiterPos);
      auto value = line.substr(delimiterPos + 1);
      if (name == "end_body")
      {
        end_body = std::stoi(value);
      }
      else if (name == "refine_edges_number")
      {
        refine_edges_number = std::stoi(value);
      }
      else if (name == "remove_system_edge_nodes_number")
      {
        remove_system_edge_nodes_number = std::stoi(value);
      }
      else if (name == "output_summary_number")
      {
        output_summary_number = std::stoi(value);
      }
      else if (name == "output_whole_system_number")
      {
        output_whole_system_number = std::stoi(value);
      }
      else if (name == "output_specific_neigh")
      {
        output_specific_neigh = std::stoi(value);
      }
      else if (name == "output_specific")
      {
        output_specific = std::stoi(value);
      }
      else if (name == "min_steps")
      {
        min_steps = std::stoi(value);
      }
      else if (name == "output_during_run")
      {
        if (value == "true")
        {
          output_during_run = true;
        }
        else
        {
          output_during_run = false;
        }
      }
    }
  }
  else
  {
    std::cerr << "Couldn't open param file: " << argv[2] << "\n";
  }

  int clocks;
  std::ifstream data_file(argv[1]);
  if (!data_file.is_open())
  {
    std::cout << "Could not open data file" << argv[1] << "\n";
    exit(0);
  }

  data_file >> step_count; //step_count++;
  data_file >> total_time;
  data_file >> clocks;
  start = time(NULL) - clocks;
  cstart = clock();
  data_file >> triangles_deleted;
  data_file >> digons_deleted;
  data_file >> edges_deleted;
  data_file >> tetrahedra_deleted;
  data_file >> footballs_deleted;
  data_file >> edge_nodes_deleted;
  data_file >> edge_nodes_added;

  data_file >> node_count;
  nodes.resize(node_count + 1);
  for (int c = 1; c <= node_count; c++)
  {
    nodes[c] = new CNode;
    data_file >> nodes[c]->x;
    data_file >> nodes[c]->y;
    data_file >> nodes[c]->z;
    nodes[c]->id = c;
  }

  data_file >> body_count;
  bodies.resize(body_count + 1);

  for (int c = 1; c <= body_count; c++)
  {
    bodies[c] = new CBody;
    int body_nodes;
    data_file >> bodies[c]->id;
    data_file >> body_nodes;
    bodies[c]->triplets.resize(body_nodes + 1);
    for (int d = 1; d <= body_nodes; d++)
    {
      int temp;
      data_file >> temp;
      bodies[c]->triplets[d].v1 = nodes[temp];
      data_file >> temp;
      bodies[c]->triplets[d].v2 = nodes[temp];
      data_file >> temp;
      bodies[c]->triplets[d].v3 = nodes[temp];
    }
  }
}
// Read in data
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
// Print out the data
void print_out_data(void)
{
  char file_name[30];
  sprintf(file_name, "%.8d-%.8d.data", step_count, body_count);
  std::ofstream outfile;
  outfile.open(file_name);

  outfile.precision(20);

  outfile << step_count << std::endl;
  outfile << total_time << std::endl;
  outfile << time(0) - start << std::endl;
  outfile << triangles_deleted << std::endl;
  outfile << digons_deleted << std::endl;
  outfile << edges_deleted << std::endl;
  outfile << tetrahedra_deleted << std::endl;
  outfile << footballs_deleted << std::endl;
  outfile << edge_nodes_deleted << std::endl;
  outfile << edge_nodes_added << std::endl;

  outfile << node_count << std::endl;

  outfile.fill('0');
  outfile << std::fixed;

  for (int c = 1; c <= node_count; c++)
    outfile << std::scientific << nodes[c]->x << "\t" << nodes[c]->y << "\t" << nodes[c]->z << "\n";

  outfile << std::endl
          << body_count << std::endl;
  for (int c = 1; c <= body_count; c++)
  {
    outfile << bodies[c]->id << " ";

    if (bodies[c]->triplets.size() != 0)
      outfile << bodies[c]->triplets.size() - 1;
    else
      outfile << 0;

    outfile << std::endl;
    for (unsigned int d = 1; d < bodies[c]->triplets.size(); d++)
      outfile << bodies[c]->triplets[d].v1->id << "\t" << bodies[c]->triplets[d].v2->id << "\t" << bodies[c]->triplets[d].v3->id << "\n";

    outfile << std::endl;
  }
}
// Print out data
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
// Determines each vertex's valence and neighbors
void determine_neighbors(void)
{
  for (unsigned int c = 1; c < nodes.size(); c++)
  {
    nodes[c]->neighbors.clear();
  }

  for (int c = 1; c <= body_count; c++)
  {
    for (unsigned int d = 1; d < bodies[c]->triplets.size(); d++)
    {
      bodies[c]->triplets[d].v1->neighbors.push_back(bodies[c]);
      bodies[c]->triplets[d].v2->neighbors.push_back(bodies[c]);
    }
  }

  for (unsigned int c = 1; c < nodes.size(); c++)
  {
    sort(nodes[c]->neighbors.begin(), nodes[c]->neighbors.end(), bcompare);
  }

  for (unsigned int c = 1; c < nodes.size(); c++)
  {
    nodes[c]->neighbors.erase(unique(nodes[c]->neighbors.begin(), nodes[c]->neighbors.end()), nodes[c]->neighbors.end());
  }
}
// Determines each vertex's valence and neighbors
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
// Given three nodes identifies the body on which this triplet lies.
// The middle node, two, is a face node, so it has two neighbors,
// and the body in question will be one of those two.
CBody *identify_body(CNode *one, CNode *two, CNode *three)
{
  CBody *body = two->neighbors[0];
  for (unsigned int c = 1; c < body->triplets.size(); c++)
    if (one == body->triplets[c].v1 &&
        two == body->triplets[c].v2 &&
        three == body->triplets[c].v3)
      return body;
  return two->neighbors[1];
}
////////////////////////////////////////////////////////////////////

CBody *identify_body(CNode *one, CNode *two)
{
  for (int c = 0; c < 4; c++)
    if (is_neighbor(two->neighbors[c], one) == 0)
      return two->neighbors[c];

  std::cout << "Error in identify body on step " << step_count << ", no body matches criteria\n";
  for (int c = 0; c < 4; c++)
  {
    one->neighbors[c]->evolver_out();
    two->neighbors[c]->evolver_out();
  }
  one->output();
  two->output();

  std::cout << "crashing in functions: identify_body";
  crash();
  return 0;
}

bool is_neighbor(CBody *body, CNode *node)
{
  for (unsigned int c = 0; c < node->neighbors.size(); c++)
    if (node->neighbors[c] == body)
      return 1;
  return 0;
}

int count_bodies()
{
  int body_count = 0;
  for (unsigned int c = 1; c < bodies.size(); c++)
    if (bodies[c]->triplets.size() > 1)
      body_count++;
  return body_count;
}

int count_vertices()
{
  int vertex_count = 0;
  for (unsigned int c = 1; c < nodes.size(); c++)
    if (nodes[c]->neighbors.size() == 4)
      vertex_count++;
  return vertex_count;
}

int count_faces()
{
  int face_count = 0;
  for (unsigned int c = 1; c < nodes.size(); c++)
    if (nodes[c]->neighbors.size() == 2)
      face_count++;
  return face_count;
}

int count_edges()
{
  return count_faces() + count_vertices() - count_bodies();
}

// Removes deleted nodes in system
void compactify()
{
  for (unsigned int c = 1; c < nodes.size(); c++)
    if (nodes[c]->neighbors.size() == 0)
    {
      delete nodes[c];
      nodes[c] = 0;
    }

  nodes.erase(remove(nodes.begin() + 1, nodes.end(), (CNode *)0), nodes.end());
  for (unsigned int c = 1; c < nodes.size(); c++)
    nodes[c]->id = c;

  node_count = nodes.size() - 1;
}

void compactify_bodies() // removes leftover bodies from the system
{
  for (unsigned int c = 1; c < bodies.size(); c++)
    if (bodies[c]->id == 0)
    {
      delete bodies[c];
      bodies[c] = 0;
    }

  bodies.erase(remove(bodies.begin() + 1, bodies.end(), (CBody *)0), bodies.end());
  body_count = bodies.size() - 1;
}

bool are_neighbors(CNode *one, CNode *two)
{
  for (unsigned int c = 0; c < one->neighbors.size(); c++)
    for (unsigned int d = 0; d < one->corners[c].size(); d++)
      if (one->corners[c][d] == two)
        return 1;
  return 0;
}

CNode *common_face(CBody *one, CBody *two, CNode *node)
{
  for (unsigned int c = 0; c < node->neighbors.size(); c++)
    for (unsigned int d = 0; d < node->corners[c].size(); d += 2)
      if ((node->corners[c][d]->neighbors[0] == one && node->corners[c][d]->neighbors[1] == two) ||
          (node->corners[c][d]->neighbors[0] == two && node->corners[c][d]->neighbors[1] == one))
        return node->corners[c][d];

  std::cout << "We have a major common_face problem, trying to find a face between body " << one->id << " and " << two->id << " that shares node " << node->id << std::endl;
  one->output();
  two->output();
  node->output();
  one->evolver_out();
  two->evolver_out();
  std::cout << "crashing in functions: common_face";
  crash();
  return 0;
}

int g(int edge1[], int edge2[], int edge_id[], int one, int two, int size)
{
  for (int c = 0; c < size; c++)
    if (edge1[c] == one && edge2[c] == two)
      return edge_id[c] + 1;
    else if (edge1[c] == two && edge2[c] == one)
      return -edge_id[c] - 1;
  return 0;
}

int f(int the_nodes[], int number, int size)
{
  for (int c = 0; c < size; c++)
    if (the_nodes[c] == number)
      return c + 1;
  return 0;
}

void evolver_out(int id)
{
  CBody *body = 0;
  for (unsigned int c = 1; c < bodies.size(); c++)
  {
    if (bodies[c]->id == id)
    {
      body = bodies[c];
      break;
    }
  }
  if (body == 0)
    return;

  if (body->get_volume() == 0)
    return;

  std::vector<CNode *> the_nodes;

  for (unsigned int c = 1; c < body->triplets.size(); c++)
  {
    the_nodes.push_back(body->triplets[c].v1);
    the_nodes.push_back(body->triplets[c].v2);
  }

  sort(the_nodes.begin(), the_nodes.end(), compare);
  the_nodes.erase(std::unique(the_nodes.begin(), the_nodes.end()), the_nodes.end());

  int nodeids[the_nodes.size()];
  for (unsigned int c = 0; c < the_nodes.size(); c++)
    nodeids[c] = the_nodes[c]->id;

  CNode origin = *the_nodes[0];

  char file_name[20];
  sprintf(file_name, "%.4d-%.8d.fe", step_count, id);
  std::ofstream out_file;
  out_file.open(file_name);

  out_file << "vertices\n";
  for (unsigned int c = 0; c < the_nodes.size(); c++)
  {
    out_file << c + 1 << '\t';
    out_file << (origin + (*the_nodes[c] - origin)).x << '\t';
    out_file << (origin + (*the_nodes[c] - origin)).y << '\t';
    out_file << (origin + (*the_nodes[c] - origin)).z << '\t';
    out_file << "original " << the_nodes[c]->id << std::endl;
  }

  int edge_counter = 0;
  int size = the_nodes.size() + body->triplets.size() - 2;
  int edge1[size];
  int edge2[size];
  int edge_id[size];

  out_file << std::endl;
  out_file << "edges\n";

  for (unsigned int d = 1; d < body->triplets.size(); d++)
  {
    if (body->triplets[d].v1->id > body->triplets[d].v2->id)
    {
      edge1[edge_counter] = body->triplets[d].v1->id;
      edge2[edge_counter] = body->triplets[d].v2->id;
      edge_id[edge_counter] = edge_counter;
      out_file << edge_counter + 1 << '\t';
      out_file << f(nodeids, edge1[edge_counter], the_nodes.size()) << '\t';
      out_file << f(nodeids, edge2[edge_counter], the_nodes.size()) << '\t';
      out_file << std::endl;
      edge_counter++;
    }
    if (body->triplets[d].v2->id > body->triplets[d].v3->id)
    {
      edge1[edge_counter] = body->triplets[d].v2->id;
      edge2[edge_counter] = body->triplets[d].v3->id;
      edge_id[edge_counter] = edge_counter;
      out_file << edge_counter + 1 << '\t';
      out_file << f(nodeids, edge1[edge_counter], the_nodes.size()) << '\t';
      out_file << f(nodeids, edge2[edge_counter], the_nodes.size()) << '\t';
      out_file << std::endl;
      edge_counter++;
    }
    if (body->triplets[d].v3->id > body->triplets[d].v1->id)
    {
      edge1[edge_counter] = body->triplets[d].v3->id;
      edge2[edge_counter] = body->triplets[d].v1->id;
      edge_id[edge_counter] = edge_counter;
      out_file << edge_counter + 1 << '\t';
      out_file << f(nodeids, edge1[edge_counter], the_nodes.size()) << '\t';
      out_file << f(nodeids, edge2[edge_counter], the_nodes.size()) << '\t';
      out_file << "color red " << std::endl;
      edge_counter++;
    }
  }

  out_file << std::endl;
  out_file << "faces\n";

  for (unsigned int d = 1; d < body->triplets.size(); d++)
  {
    out_file << d << '\t'
             << g(edge1, edge2, edge_id, body->triplets[d].v1->id, body->triplets[d].v2->id, size) << '\t'
             << g(edge1, edge2, edge_id, body->triplets[d].v2->id, body->triplets[d].v3->id, size) << '\t'
             << g(edge1, edge2, edge_id, body->triplets[d].v3->id, body->triplets[d].v1->id, size) << '\t'
             << std::endl;
  }

  out_file << std::endl;
  out_file << "bodies\n";
  out_file << 1 << '\t';

  for (unsigned int d = 1; d < body->triplets.size(); d++)
    out_file << d << ' ';

  out_file << std::endl;
  out_file.close();
}

bool compare(CNode *one, CNode *two)
{
  if (one || two)
    return one->id < two->id;
  else
    std::cout << "crashing in functions: compare, not one or two";
  crash();
  return 0;
}

bool bcompare(CBody *one, CBody *two)
{
  return one->id < two->id;
}

void crash()
{
  std::cout << "crashing\n";
  nodes[23423423423423]->output();
  exit(0);
}

double calc_echange(CNode *one, CNode *two)
{
  double e1x = one->x - two->x;
  double e1y = one->y - two->y;
  double e1z = one->z - two->z;

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

  double initial_length = sqrt(e1x * e1x + e1y * e1y + e1z * e1z);

  if (one->good == 1)
  {
    e1x += one->dx * time_scale;
    e1y += one->dy * time_scale;
    e1z += one->dz * time_scale;
  }
  else
  {
    e1x += one->cx;
    e1y += one->cy;
    e1z += one->cz;
  }

  if (two->good == 1)
  {
    e1x -= two->dx * time_scale;
    e1y -= two->dy * time_scale;
    e1z -= two->dz * time_scale;
  }
  else
  {
    e1x -= two->cx;
    e1y -= two->cy;
    e1z -= two->cz;
  }

  if (sqrt(e1x * e1x + e1y * e1y + e1z * e1z) == initial_length || initial_length < 0.00001)
    return -1;
  else
    return (sqrt(e1x * e1x + e1y * e1y + e1z * e1z) - initial_length) / initial_length;
}

bool not_digon_neighbors(CNode *one, CNode *two)
{
  if (one->neighbors.size() != two->neighbors.size())
    return 1;
  for (unsigned int c = 0; c < one->neighbors.size(); c++)
    if (one->neighbors[c] != two->neighbors[c])
      return 1;
  return 0;
}

void surface_evolver_system()
{
  std::vector<edge> edges;
  std::vector<face> faces;

  for (unsigned int c = 1; c < nodes.size(); c++)
  {
    std::vector<CNode *> neighbors;
    for (unsigned int d = 0; d < nodes[c]->neighbors.size(); d++)
      for (unsigned int e = 0; e < nodes[c]->corners[d].size(); e++)
        if (nodes[c]->id < nodes[c]->corners[d][e]->id)
          neighbors.push_back(nodes[c]->corners[d][e]);
    sort(neighbors.begin(), neighbors.end(), compare);
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());

    for (unsigned int d = 0; d < neighbors.size(); d++)
      edges.push_back(edge(nodes[c], neighbors[d]));
  }

  char file_name[30];
  sprintf(file_name, "%.6d-%.6d.fe", step_count, count_bodies());
  std::ofstream out_file;
  out_file.open(file_name);

  out_file << "define body attribute xseed real[3]\n";
  out_file << "torus_filled\n";
  out_file << "periods\n";
  out_file << "1 0 0\n";
  out_file << "0 1 0\n";
  out_file << "0 0 1\n";
  out_file << "\n";
  out_file << "vertices\n";

  for (unsigned int c = 1; c < nodes.size(); c++)
    out_file << nodes[c]->id << '\t' << nodes[c]->x << '\t' << nodes[c]->y << '\t' << nodes[c]->z << '\n';

  out_file << "\n";
  out_file << "edges\n";

  for (unsigned int c = 0; c < edges.size(); c++)
  {
    out_file << c + 1 << '\t' << edges[c].one->id << '\t' << edges[c].two->id << '\t';

    if (edges[c].one->x - edges[c].two->x > 0.5)
      out_file << "+" << '\t';
    else if (edges[c].one->x - edges[c].two->x < -0.5)
      out_file << "-" << '\t';
    else
      out_file << "*" << '\t';

    if (edges[c].one->y - edges[c].two->y > 0.5)
      out_file << "+" << '\t';
    else if (edges[c].one->y - edges[c].two->y < -0.5)
      out_file << "-" << '\t';
    else
      out_file << "*" << '\t';

    if (edges[c].one->z - edges[c].two->z > 0.5)
      out_file << "+" << '\n';
    else if (edges[c].one->z - edges[c].two->z < -0.5)
      out_file << "-" << '\n';
    else
      out_file << "*" << '\n';
  }

  int face_counter = 1;
  out_file << "\n";
  out_file << "faces\n";

  int sedges = edges.size() / 1000;
  int eindex[sedges];
  for (int c = 0; c < sedges; c++)
    eindex[c] = edges[c * 1000].one->id;

  for (unsigned int c = 1; c < nodes.size(); c++)
  {
    if (nodes[c]->neighbors.size() != 2)
      continue;

    for (unsigned int d = 0; d < nodes[c]->corners[0].size(); d++)
    {
      CNode *corners[3];

      corners[0] = nodes[c]->corners[0][d];
      corners[1] = nodes[c];
      corners[2] = nodes[c]->corners[0][(d + 1) % nodes[c]->corners[0].size()];

      std::vector<int> cindices;
      for (int e = 0; e < 3; e++)
        cindices.push_back(corners[e]->id);
      sort(cindices.begin(), cindices.end());

      int minindex = 0;
      int midindex = 0;

      for (int c = 0; c < sedges; c++)
      {
        if (eindex[c] < cindices[0])
          minindex++;
        if (eindex[c] < cindices[1])
          midindex++;
      }

      int start = (minindex - 1) * 1000;
      if (start < 0)
        start = 0;
      int end = (midindex + 1) * 1000;
      if (end > int(edges.size()))
        end = edges.size();

      int indices[3] = {0, 0, 0};
      for (int e = start; e < end; e++)
      {
        if (corners[0] == edges[e].one && corners[1] == edges[e].two)
          indices[0] = e + 1;
        if (corners[0] == edges[e].two && corners[1] == edges[e].one)
          indices[0] = -e - 1;
        if (corners[1] == edges[e].one && corners[2] == edges[e].two)
          indices[1] = e + 1;
        if (corners[1] == edges[e].two && corners[2] == edges[e].one)
          indices[1] = -e - 1;
        if (corners[2] == edges[e].one && corners[0] == edges[e].two)
          indices[2] = e + 1;
        if (corners[2] == edges[e].two && corners[0] == edges[e].one)
          indices[2] = -e - 1;
      }

      faces.push_back(face(corners[0], corners[1], corners[2]));
      out_file << face_counter << '\t' << indices[0] << '\t' << indices[1] << '\t' << indices[2] << '\n';
      face_counter++;
    }
  }

  out_file << "\n";
  out_file << "bodies\n";

  for (unsigned int c = 1; c < bodies.size(); c++)
  {
    out_file << c << '\t';

    for (unsigned int e = 0; e < faces.size(); e++)
    {
      if (faces[e].two->neighbors[0] == bodies[c])
        out_file << e + 1 << " ";
      if (faces[e].two->neighbors[1] == bodies[c])
        out_file << -(signed int)e - 1 << " ";
    }

    out_file << '\n';
  }

  out_file.close();
}

void fit_nodes_in_cube()
{
  for (unsigned int c = 1; c < nodes.size(); c++)
  {
    if (nodes[c]->x < 0.)
      nodes[c]->x += 1.0;
    if (nodes[c]->x >= 1.)
      nodes[c]->x -= 1.0;
    if (nodes[c]->y < 0.)
      nodes[c]->y += 1.0;
    if (nodes[c]->y >= 1.)
      nodes[c]->y -= 1.0;
    if (nodes[c]->z < 0.)
      nodes[c]->z += 1.0;
    if (nodes[c]->z >= 1.)
      nodes[c]->z -= 1.0;
  }
}

void calc_volumes()
{
#pragma omp parallel
  for (unsigned int c = omp_get_thread_num() + 1; c < bodies.size(); c += nthreads)
    bodies[c]->calc_volume();
}

void calc_areas()
{
#pragma omp parallel
  for (unsigned int c = omp_get_thread_num() + 1; c < nodes.size(); c += nthreads)
    nodes[c]->calc_face_area();
}

void determine_node_topologies()
{
  for (unsigned int c = 1; c < nodes.size(); c++)
    nodes[c]->determine_node_topology(); // CAN BE PARALLELIZED I THINK
}

void calc_motions()
{
#pragma omp parallel
  for (unsigned int c = omp_get_thread_num() + 1; c < nodes.size(); c += nthreads)
    nodes[c]->calc_motion();
}

void move_nodes()
{
#pragma omp parallel
  for (unsigned int c = omp_get_thread_num() + 1; c < nodes.size(); c += nthreads)
    nodes[c]->move();
#pragma omp parallel
  for (unsigned int c = omp_get_thread_num() + 1; c < nodes.size(); c += nthreads)
    nodes[c]->calc_center();
#pragma omp parallel
  for (unsigned int c = omp_get_thread_num() + 1; c < nodes.size(); c += nthreads)
    nodes[c]->move_center();
}

double triplet::area()
{
  double e1x = v1->x - v2->x;
  double e1y = v1->y - v2->y;
  double e1z = v1->z - v2->z;
  double e2x = v3->x - v2->x;
  double e2y = v3->y - v2->y;
  double e2z = v3->z - v2->z;

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

  return sqrt((e1y * e2z - e1z * e2y) * (e1y * e2z - e1z * e2y) + (e1z * e2x - e1x * e2z) * (e1z * e2x - e1x * e2z) + (e1x * e2y - e1y * e2x) * (e1x * e2y - e1y * e2x)) / 2.0;
}

// THIS IS USED ONLY TO REMOVE AN EDGE BETWEEN TWO QUAD-NODES.  NEITHER NODE IS AN EDGE NODE.
CNode *remove_edge(CNode *nodeU, CNode *nodeD)
{
#ifdef DEBUG
  std::cout << "Removing edge " << std::setw(5) << nodeU->id << " - " << std::setw(5) << nodeD->id << " on step " << step_count << "; ";
#endif

  // Identify important nodes and bodies
  CBody *up = identify_body(nodeD, nodeU);   // check
  CBody *down = identify_body(nodeU, nodeD); // check

  std::vector<CNode *> all_neighbors;
  std::vector<CNode *> top_nodes(3);
  std::vector<CNode *> top_face_nodes(3);
  std::vector<CNode *> bottom_nodes(3);
  std::vector<CNode *> bottom_face_nodes(3);
  std::vector<CNode *> side_face_nodes(3);
  std::vector<CBody *> side_bodies(3);

  int index = 0;
  for (unsigned int c = 0; c < 4; c++)
    if (nodeU->neighbors[c] == up)
      index = c;
  for (unsigned int c = 0; c < 3; c++)
  {
    if (nodeU->corners[index][2 * c]->neighbors[0] == up)
      side_bodies[c] = nodeU->corners[index][2 * c]->neighbors[1];
    else
      side_bodies[c] = nodeU->corners[index][2 * c]->neighbors[0];
  }

  for (unsigned int c = 0; c < 3; c++)
    side_face_nodes[c] = common_face(side_bodies[c], side_bodies[(c + 1) % 3], nodeU);
  double fareas[3];
  for (unsigned int c = 0; c < 3; c++)
    fareas[c] = side_face_nodes[c]->get_face_area();
  double perims[3];
  for (unsigned int c = 0; c < 3; c++)
    perims[c] = side_face_nodes[c]->get_face_perimeter();

  for (unsigned int c = 0; c < 3; c++)
  {
    if (are_neighbors(side_face_nodes[0], nodeU->corners[index][2 * c + 1]))
      top_nodes[0] = nodeU->corners[index][2 * c + 1];
    if (are_neighbors(side_face_nodes[1], nodeU->corners[index][2 * c + 1]))
      top_nodes[1] = nodeU->corners[index][2 * c + 1];
    if (are_neighbors(side_face_nodes[2], nodeU->corners[index][2 * c + 1]))
      top_nodes[2] = nodeU->corners[index][2 * c + 1];

    if (is_neighbor(side_bodies[0], nodeU->corners[index][2 * c]))
      top_face_nodes[0] = nodeU->corners[index][2 * c];
    if (is_neighbor(side_bodies[1], nodeU->corners[index][2 * c]))
      top_face_nodes[1] = nodeU->corners[index][2 * c];
    if (is_neighbor(side_bodies[2], nodeU->corners[index][2 * c]))
      top_face_nodes[2] = nodeU->corners[index][2 * c];
  }

  for (unsigned int c = 0; c < 4; c++)
    if (nodeD->neighbors[c] == down)
      index = c;
  for (unsigned int c = 0; c < 3; c++)
  {
    if (is_neighbor(side_bodies[0], nodeD->corners[index][2 * c]))
      bottom_face_nodes[0] = nodeD->corners[index][2 * c];
    if (is_neighbor(side_bodies[1], nodeD->corners[index][2 * c]))
      bottom_face_nodes[1] = nodeD->corners[index][2 * c];
    if (is_neighbor(side_bodies[2], nodeD->corners[index][2 * c]))
      bottom_face_nodes[2] = nodeD->corners[index][2 * c];

    if (are_neighbors(side_face_nodes[0], nodeD->corners[index][2 * c + 1]))
      bottom_nodes[0] = nodeD->corners[index][2 * c + 1];
    if (are_neighbors(side_face_nodes[1], nodeD->corners[index][2 * c + 1]))
      bottom_nodes[1] = nodeD->corners[index][2 * c + 1];
    if (are_neighbors(side_face_nodes[2], nodeD->corners[index][2 * c + 1]))
      bottom_nodes[2] = nodeD->corners[index][2 * c + 1];
  }

  double tfareas[3];
  for (unsigned int c = 0; c < 3; c++)
    tfareas[c] = top_face_nodes[c]->get_face_area();
  double tperims[3];
  for (unsigned int c = 0; c < 3; c++)
    tperims[c] = top_face_nodes[c]->get_face_perimeter();

  double bfareas[3];
  for (unsigned int c = 0; c < 3; c++)
    bfareas[c] = bottom_face_nodes[c]->get_face_area();
  double bperims[3];
  for (unsigned int c = 0; c < 3; c++)
    bperims[c] = bottom_face_nodes[c]->get_face_perimeter();

  // Make new nodes, place them appropriately
  nodes.resize(nodes.size() + 4);
  for (unsigned int c = 1; c < 5; c++)
    nodes[node_count + c] = new CNode();

  CNode *nodeA = nodes[node_count + 1]; // the three side nodes
  CNode *nodeB = nodes[node_count + 2];
  CNode *nodeC = nodes[node_count + 3];
  CNode *nodeF = nodes[node_count + 4]; // the face node
  *nodeF = *nodeU + (*nodeD - *nodeU) / 2;

  *nodeA = (*top_nodes[0] - *nodeU) / 2. + (*bottom_nodes[0] - *nodeD) / 2.;
  *nodeB = (*top_nodes[1] - *nodeU) / 2. + (*bottom_nodes[1] - *nodeD) / 2.;
  *nodeC = (*top_nodes[2] - *nodeU) / 2. + (*bottom_nodes[2] - *nodeD) / 2.;

  // These lines compute the area of a triangle whose three 'inside edges' are nodeA, nodeB, and nodeC.  We want
  // to 'normalize' the area of the new triangle formed so that its area is 4*smallest_face.  This should prevent
  // it from being erased on the next step and should also hopefully be small enough

  double edge_length = length(nodeU, nodeD);

  *nodeA *= (edge_length / nodeA->norm());
  *nodeB *= (edge_length / nodeB->norm());
  *nodeC *= (edge_length / nodeC->norm());

#ifdef DEBUG
  vector<int> numbers;
  numbers.push_back(up->id);
  numbers.push_back(down->id);
  numbers.push_back(side_bodies[0]->id);
  numbers.push_back(side_bodies[1]->id);
  numbers.push_back(side_bodies[2]->id);
  sort(numbers.begin(), numbers.end());
#endif

  *nodeA += *nodeF;
  *nodeB += *nodeF;
  *nodeC += *nodeF;
  *nodeF = *nodeA + (*nodeB - *nodeA) / 3. + (*nodeC - *nodeA) / 3.;

  nodeF->neighbors.push_back(up);
  nodeF->neighbors.push_back(down);

  nodeA->neighbors.push_back(up);
  nodeB->neighbors.push_back(up);
  nodeC->neighbors.push_back(up);
  nodeA->neighbors.push_back(down);
  nodeB->neighbors.push_back(down);
  nodeC->neighbors.push_back(down);

  nodeA->neighbors.push_back(side_bodies[0]);
  nodeA->neighbors.push_back(side_bodies[1]);
  nodeB->neighbors.push_back(side_bodies[1]);
  nodeB->neighbors.push_back(side_bodies[2]);
  nodeC->neighbors.push_back(side_bodies[2]);
  nodeC->neighbors.push_back(side_bodies[0]);

  node_count++;
  nodeA->id = node_count;
  node_count++;
  nodeB->id = node_count;
  node_count++;
  nodeC->id = node_count;
  node_count++;
  nodeF->id = node_count;

#ifdef DEBUG
  std::cout << "new nodes " << nodeA->id << ", " << nodeB->id << ", " << nodeC->id << ", and face node " << nodeF->id << "; ";

  std::cout << "\t\t\t\t\t\t\t\t\t";
  for (unsigned int c = 0; c < 5; c++)
    std::cout << std::setw(3) << numbers[c] << "  ";
  std::cout << 'E' << endl;
#endif

  side_bodies[0]->switch_edge_with_edge(nodeU, nodeD, nodeA, nodeC);
  side_bodies[1]->switch_edge_with_edge(nodeU, nodeD, nodeB, nodeA);
  side_bodies[2]->switch_edge_with_edge(nodeU, nodeD, nodeC, nodeB);

  up->replace_node_on_face_with_edge(top_face_nodes[0], nodeU, nodeA, nodeC);
  up->replace_node_on_face_with_edge(top_face_nodes[1], nodeU, nodeB, nodeA);
  up->replace_node_on_face_with_edge(top_face_nodes[2], nodeU, nodeC, nodeB);

  down->replace_node_on_face_with_edge(bottom_face_nodes[0], nodeD, nodeC, nodeA);
  down->replace_node_on_face_with_edge(bottom_face_nodes[1], nodeD, nodeA, nodeB);
  down->replace_node_on_face_with_edge(bottom_face_nodes[2], nodeD, nodeB, nodeC);

  up->triplets.push_back(triplet(nodeA, nodeF, nodeB));
  up->triplets.push_back(triplet(nodeB, nodeF, nodeC));
  up->triplets.push_back(triplet(nodeC, nodeF, nodeA));
  down->triplets.push_back(triplet(nodeB, nodeF, nodeA));
  down->triplets.push_back(triplet(nodeA, nodeF, nodeC));
  down->triplets.push_back(triplet(nodeC, nodeF, nodeB));

  // Get rid of old stuff
  nodeU->neighbors.clear();
  nodeD->neighbors.clear();

  // Recompute topological information
  for (unsigned int c = 0; c < 3; c++)
  {
    top_nodes[c]->determine_node_topology();
    top_face_nodes[c]->determine_node_topology();
    bottom_nodes[c]->determine_node_topology();
    bottom_face_nodes[c]->determine_node_topology();
    side_face_nodes[c]->determine_node_topology();
  }

  nodeA->determine_node_topology();
  nodeB->determine_node_topology();
  nodeC->determine_node_topology();
  nodeF->determine_node_topology();

  // CUT THIS?
  for (unsigned int c = 0; c < 50; c++)
  {
    nodeA->center2(0.001);
    nodeB->center2(0.001);
    nodeC->center2(0.001);
    nodeF->center2(0.100);
  }

  double initial_scale = time_scale;
  time_scale /= 100.;

  for (unsigned int c = 0; c < 3; c++)
    side_face_nodes[c]->center();

  for (unsigned int c = 0; c < 3; c++)
  {
    if (fabs(side_face_nodes[c]->get_face_area() - fareas[c]) / fareas[c] > 0.2 || fabs(side_face_nodes[c]->get_face_perimeter() - perims[c]) / perims[c] > 0.2)
    {
#ifdef DEBUG
      std::cout << "We have changed the area of face " << side_face_nodes[c]->id << " much while removing an edge.  Was " << side_face_nodes[c]->area << " or "
                << fareas[c] << " and now is " << side_face_nodes[c]->get_face_area() << "\n";
#endif
      side_face_nodes[c]->center2(1.);
    }
#ifdef DEBUG
    //    else cout<<"Changed face " << side_face_nodes[c]->id << " area by " << (side_face_nodes[c]->get_face_area()-fareas[c])/fareas[c] << " and perimeter by " << (side_face_nodes[c]->get_face_perimeter()-perims[c])/perims[c] << '\t' << side_face_nodes[c]->get_face_perimeter() << '\n';
#endif
  }

  for (unsigned int c = 0; c < 3; c++)
  {
    if (fabs(top_face_nodes[c]->get_face_area() - tfareas[c]) / tfareas[c] > 0.2 || fabs(top_face_nodes[c]->get_face_perimeter() - tperims[c]) / tperims[c] > 0.2)
    {
#ifdef DEBUG
      std::cout << "We have changed the area of top face " << top_face_nodes[c]->id << " much while removing an edge.  Was " << top_face_nodes[c]->area << " or "
                << tfareas[c] << " and now is " << top_face_nodes[c]->get_face_area() << "\n";
#endif
      top_face_nodes[c]->center2(1.);
    }
  }

  for (unsigned int c = 0; c < 3; c++)
  {
    if (fabs(bottom_face_nodes[c]->get_face_area() - bfareas[c]) / bfareas[c] > 0.2 || fabs(bottom_face_nodes[c]->get_face_perimeter() - bperims[c]) / bperims[c] > 0.2)
    {
#ifdef DEBUG
      std::cout << "We have changed the area of bottom face " << bottom_face_nodes[c]->id << " much while removing an edge.  Was " << bottom_face_nodes[c]->area << " or "
                << bfareas[c] << " and now is " << bottom_face_nodes[c]->get_face_area() << "\n";
#endif
      bottom_face_nodes[c]->center2(1.);
    }
  }

  nodeF->center2(1.0);
  time_scale = initial_scale;

  // Recompute volumes
  up->calc_volume();
  down->calc_volume();
  side_bodies[0]->calc_volume();
  side_bodies[1]->calc_volume();
  side_bodies[2]->calc_volume();

  for (unsigned int c = 0; c < 3; c++)
    if (side_face_nodes[c]->corners[0].size() == 3)
      side_face_nodes[c]->center2(1.);
  nodeF->center2(1.0);

  for (unsigned int c = 0; c < 4; c++)
    for (unsigned int d = 0; d < 6; d++)
    {
      all_neighbors.push_back(nodeA->corners[c][d]);
      all_neighbors.push_back(nodeB->corners[c][d]);
      all_neighbors.push_back(nodeC->corners[c][d]);
    }
  sort(all_neighbors.begin(), all_neighbors.end(), compare);
  all_neighbors.erase(std::unique(all_neighbors.begin(), all_neighbors.end()), all_neighbors.end());
  for (unsigned int c = 0; c < all_neighbors.size(); c++)
    all_neighbors[c]->calc_motion();

  nodeF->area = nodeF->get_face_area();

#ifdef DEBUG
  std::cout << "New face has area " << nodeF->area << endl;
  std::cout << "Done removing edge " << std::setw(5) << nodeU->id << " - " << std::setw(5) << nodeD->id << " on step " << step_count << "\n\n";
  std::cout << '\n'
            << '\n';
#endif
  return nodeF;
}
