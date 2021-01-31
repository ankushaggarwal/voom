#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <unistd.h>
#include <tvmet/Vector.h>
#include "Node.h"
#include "SCElastic.h"
#include "EvansElastic.h"
#include "LoopShellBody.h"
#include "C0MembraneBody.h"
#include "ShapeTri3.h"
#include "Model.h"
#include "NodeBin.h"

using namespace std;
using namespace voom;

//! This function creates a shell model with LoopShellBody  for bending and  C0MembraneBodies for stretching. Then it reads the configurations and calculates strain invariants and curvatures

//void readconfig_vtk(char* filename, float* x, float* y, float* z, int npts);
void readconfig_vtk(char* filename, float* x, float* y, float* z, int npts, const std::vector<std::vector<int> > & nodeid, const std::vector<std::vector<double> > & weights);

void printParaview(const Model::BodyContainer bdc, char* filename, const std::vector<double>& curv, const std::vector<double>& I1, const std::vector<double>& I2);

void calc_projection(const char* coarseMeshFileName,const std::vector< DeformationNode<3>* > defNodes,std::vector<std::vector<int> > & nodeid,std::vector<std::vector<double> > & weights);

int main(int argc, char* argv[])
{
  if(argc < 2){
    std::cout<<"Usage calc_strains inputfile targetfile_prefix inputCoarseFile [-n start_config -N end_config]." << std::endl;
    return(0);
  }
  
  bool verbose=true;
  for(int i=0; i<argc; i++) {
    std::cout << std::setw(8) << i << "\t"
	      << argv[i] << std::endl;
  }
  // create input stream
  std::ifstream ifs;
  std::string inputFileName = argv[1];
  std::string targetFileNamePrefix = argv[2];
  std::string coarseMeshFileName = argv[3];
  ifs.open( inputFileName.c_str(), std::ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << inputFileName << endl;
    exit(0);
  }

  bool badCommandLine=false;
  int StartConfig=0, EndConfig=0;
  
  for(char option; (option=getopt(argc,argv,"n:N:")) != EOF; ) 
    switch (option) {
    default :
      badCommandLine=true;
      break;
    case 'n' :
      StartConfig = std::atoi(optarg);
      std::cout << "Starting No. of configurations: " << StartConfig << std::endl;
      break;
    case 'N' :
      EndConfig = std::atoi(optarg);
      std::cout << "Ending No. of configurations: " << EndConfig << std::endl;
      break;
    }

  if( badCommandLine || argc < 2 ) {
    std::cout<<"Usage calc_strains inputfile targetfile_prefix inputCoarseFile [-n start_config -N end_config]." << std::endl;
    return(0);
  }
  
  // create vector of nodes
  int dof=0;
  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;
  char key;  ifs>>key;
  
  // Input .vtk file containing nodes and connectivities
  std::string token;
  ifs >> token; 
  while( token != "POINTS" ) ifs >> token;
  int npts=0;
  ifs >> npts; 
  defNodes.reserve(npts);
  ifs >> token;// skip number type
  std::cout<<"npts = "<<npts<<std::endl;
  // read in points
  for(int i=0; i<npts; i++) {
    int id=i;
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  assert(nodes.size()!=0);
  std::cout << "Number of nodes: " <<nodes.size() << std::endl;
  
  // read in triangle connectivities
  while( token != "POLYGONS" ) ifs >> token;
  std::vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int,3> c;
  int ntri=0; ifs >> ntri;
  connectivities.reserve(ntri);
  cout << "Number of triangles: " << ntri << endl;
  int tmp; ifs>>tmp; std::cout<<" tmp = "<<tmp<<std::endl;
  for (int i = 0; i<ntri; i++){
    int tmp=0;
    ifs >> tmp;
    if(tmp!=3) std::cout<<"Some mistake reading the elements connectivity from file. Check again."<<std::endl;
    for(int a=0; a<3; a++) ifs >> c[a];
    connectivities.push_back(c);
  }

  // material constants
  // bending
  double Y=1.;
  double KC=1.;
  double KG=1.;
  double nu=1./3.;
  double C0=0.;
  // stretching
  double mu=1.;
  double kS=1.;

  // create bending body
  int quadOrder=1;
  
  typedef SCElastic BendingMaterial;
  BendingMaterial bending(KC,KG,C0);

  typedef LoopShellBody<BendingMaterial> LSB;
  LSB * bending_body = new LSB(bending, connectivities, nodes, quadOrder);
  bending_body->setOutput(paraview);

  // create stretching body which doesn't have any bending energy
  std::vector< std::vector<int> > s_connectivities;
  for(int i=0; i<connectivities.size(); i++) {
    std::vector<int> c(3);
    for(int j=0; j<3; j++) c[j]=connectivities[i](j);
    s_connectivities.push_back(c);
  }
  typedef EvansElastic StretchingMaterial;
  StretchingMaterial stretch(0.,0.,0.,mu,kS);

  typedef C0MembraneBody<TriangleQuadrature,StretchingMaterial,ShapeTri3> membrane;
  membrane * stretch_body = new membrane(stretch,s_connectivities,nodes,quadOrder);
  stretch_body->setOutput(paraview);

  std::cout << std::endl;
  std::cout << "Bending body energy = " << bending_body->energy() << std::endl;
  std::cout << "Stretching body energy = " << stretch_body->energy() << std::endl;
  std::cout << std::endl;

  // Create Model
  Model::BodyContainer bdc;
  bdc.push_back(bending_body);
  bdc.push_back(stretch_body);
  
  //Calculate the barycentric coordinates of these nodes w.r.t. a coarser mesh
  std::vector<int> initial_id(3,0);
  std::vector<std::vector<int> > nodeid(npts,initial_id);
  std::vector<double> initial_weight(3,0.);
  std::vector<std::vector<double> > weights(npts,initial_weight);
  calc_projection(coarseMeshFileName.c_str(),defNodes,nodeid,weights);
  /*//write the areas of elements needed in modulus calculation
  std::ofstream outfile;
  outfile.open("area.dat");
  typedef membrane::MembraneElementContainer elements;
  const elements temp_elements=stretch_body->elements();
  for(int i=0;i<temp_elements.size();i++){
    double area=temp_elements[i]->area();
    outfile << area << std::endl;
  }
  outfile.close();*/

/*  //read the trajectory and output curvatures and invariants
  std::ofstream file1, file2, file3;
  file1.open("inv1.dat",std::ios::binary);
  file2.open("inv2.dat",std::ios::binary);
  file3.open("curv.dat",std::ios::binary);

//  char tab='\t', newline='\n';
  int col=ntri+1, row=EndConfig-StartConfig+1;

  file1.write((char*)&row,sizeof(int)); file1.write((char*)&col,sizeof(int));
  file2.write((char*)&row,sizeof(int)); file2.write((char*)&col,sizeof(int));
  file3.write((char*)&row,sizeof(int)); file3.write((char*)&col,sizeof(int));
*/
  //read in the configuration files in a loop
  int di;
  if(StartConfig<=EndConfig) di=1;
  else di=-1;
  for(int i=StartConfig;i<=EndConfig;i+=di){
    char inpfile[50];
    std::cout<<"Reading "<<i<<"th configuration."<<std::endl;
    float x[npts], y[npts], z[npts];
    //sprintf(inpfile,"../surf-traj/%d.dat",i);
    sprintf(inpfile,"%s%02d_surface.vtk",targetFileNamePrefix.c_str(),i);
    readconfig_vtk(inpfile,x,y,z,npts,nodeid,weights);
    for(int inode=0;inode<npts;inode++){
      DeformationNode<3>::Point newx;
      newx(0)=x[inode]; newx(1)=y[inode]; newx(2)=z[inode];
      defNodes[inode]->setPoint(newx);
    }
    /*std::ofstream dispfile;
    dispfile.open("disp.dat");
    for(int inode=0;inode<npts;inode++){
      Vector3D nodalDisp;
      nodalDisp = defNodes[inode]->point() - defNodes[inode]->position();
      dispfile << nodalDisp(0) << '\t' << nodalDisp(1) << '\t' << nodalDisp(2) << std::endl;
    }
    dispfile.close();*/
    bending_body->compute(true,false,false); 
    std::vector<double> I1(ntri,0.), I2(ntri,0.), curv(ntri,0.);
    std::cout<<"Calculating the invariants"<<std::endl;
    bending_body->cal_curv(curv);
    stretch_body->cal_inv(I1,I2);
    //Printing the output files
    char outfile[50];
    //sprintf(outfile,"curv-%d",i);
    //bending_body->printParaview(outfile);
    sprintf(outfile,"strains-%02d.vtk",i);
    printParaview(bdc, outfile, curv, I1, I2);
    //sprintf(outfile,"inv-%d",i);
    //stretch_body->printParaviewInv(outfile);
  }
  std::cout<<"Everything done!"<<std::endl;
  return 0;
}

void readconfig_vtk(char* filename, float* x, float* y, float* z, int npts, const std::vector<std::vector<int> > & nodeid, const std::vector<std::vector<double> > & weights){
  ifstream inpfile;
  inpfile.open(filename);
  if (!inpfile) {
    cout << "Cannot open input file: " << filename << endl;
    exit(0);
  }
  int id;
  std::string token;
  inpfile >> token; 
  while( token != "POINTS" ) inpfile >> token;
  int coarse_npts=0;
  inpfile >> coarse_npts;
  std::vector<double> xtemp, ytemp, ztemp; 
  xtemp.resize(coarse_npts,0.);
  ytemp.resize(coarse_npts,0.);
  ztemp.resize(coarse_npts,0.);
  inpfile >> token;// skip number type
  // read in points
  for(int i=0; i<coarse_npts; i++){
    inpfile >>  xtemp[i] >> ytemp[i] >> ztemp[i];
    //std::cout<<xtemp[i]<<"\t"<<ytemp[i]<<"\t"<<ztemp[i]<<std::endl;
  }
  inpfile.close();
  //std::cout<<"Calculated points"<<std::endl;
  for(int i=0; i<npts; i++){
    x[i]=y[i]=z[i]=0.;
    for( int j=0; j<3; j++){
      x[i] += xtemp[nodeid[i][j]]*weights[i][j];
      y[i] += ytemp[nodeid[i][j]]*weights[i][j];
      z[i] += ztemp[nodeid[i][j]]*weights[i][j];
      //std::cout<<nodeid[i][j]<<"\t"<<xtemp[nodeid[i][j]]<<"\t"<<ytemp[nodeid[i][j]]<<"\t"<<ztemp[nodeid[i][j]]<<"\t"<<weights[i][j]<<std::endl;
    } 
    //std::cout<<weights[i][0] <<"\t"<<weights[i][1]<<"\t"<<weights[i][2]<<std::endl; 
    //std::cout<<x[i]<<"\t"<<y[i]<<"\t"<<z[i]<<std::endl;
    //std::cout<<"-------------------------------"<<std::endl;
  }
}

void printParaview(const Model::BodyContainer bdc, char* filename, const std::vector<double>& curv, const std::vector<double>& I1, const std::vector<double>& I2){
  ofstream ofs;
  ofs.open(filename);
  if (!ofs){
    std::cout << "Cannot open output file: " << filename << endl;
    exit(0);
  }

  int nbd=bdc.size();
  if(nbd<0) return;
  int npts=bdc[0]->nodes().size();
  for(int i=0;i<nbd;i++)
    if(bdc[i]->nodes().size() != npts)
      std::cout<<"print_body Error: Bodies have different node numbers"<<std::endl;

  //    Node Section
  ofs << "# vtk DataFile Version 2.0\n"
      << "Test example" << std::endl
      << "ASCII" << std::endl
      << "DATASET POLYDATA" << std::endl
      << "POINTS  " << npts << "  double" << std::endl;

  // output nodal postions
  Body::NodeContainer nodes=bdc[0]->nodes();
  for(int i=0;i<npts;i++){
      NodeBase* n=nodes[i];
      ofs << std::setprecision(16) 
          << n->getPoint(0) << " "
          << n->getPoint(1) << " "
          << n->getPoint(2) << std::endl;
    }
  
  //output connectivities
  Body::ElementContainer BendingElements = bdc[0]->elements();
  // Element Section
  int nel = BendingElements.size();
  ofs << "POLYGONS  " << nel << "  "
      << 4*nel << endl;
  for(int e = 0; e < nel; e++)
  {
    if((bdc[0]->active(e)) ) 
    ofs << 3 << "  "
        << setw(10) << BendingElements[e]->baseNodes()[0] -> id()
        << setw(10) << BendingElements[e]->baseNodes()[1] -> id()
        << setw(10) << BendingElements[e]->baseNodes()[2] -> id()
        << endl;
  }
 
  // Cell section
  ofs << "CELL_DATA    " << nel << endl;
  // Bending energy
  ofs << "SCALARS    MeanCurvature    double    1" << endl;
  ofs << "LOOKUP_TABLE default" << endl;
  for(int e = 0; e < nel; e++) {
      ofs << curv[e] << endl;
  }
  ofs << std::endl;
  
  ofs << "SCALARS    ShearInvariant    double    1" << endl;
  ofs << "LOOKUP_TABLE default" << endl;
  for(int e = 0; e < nel; e++) {
      ofs << I1[e] << endl;
  }
  ofs << std::endl;

  ofs << "SCALARS    AreaInvariant    double    1" << endl;
  ofs << "LOOKUP_TABLE default" << endl;
  for(int e = 0; e < nel; e++) {
      ofs << I2[e] << endl;
  }
  ofs << std::endl;
  //ofs << endl << "POINT_DATA " << npts << endl
  //    << "VECTORS displacements double" << endl;
  ofs.close();
}

void calc_projection(const char* coarseMeshFileName,const std::vector< DeformationNode<3>* > fineNodes,std::vector<std::vector<int> > & nodeid,std::vector<std::vector<double> > & weights){

  std::ifstream ifs;
  ifs.open( coarseMeshFileName, std::ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << coarseMeshFileName << endl;
    exit(0);
  }
  typedef DeformationNode<3> DefNode;
  // create vector of nodes
  int dof=0;
  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;
  char key;  ifs>>key;
  
  // Input .vtk file containing nodes and connectivities
  std::string token;
  ifs >> token; 
  while( token != "POINTS" ) ifs >> token;
  int npts=0;
  ifs >> npts; 
  defNodes.reserve(npts);
  ifs >> token;// skip number type
  std::cout<<"calc_prjection::npts = "<<npts<<std::endl;
  // read in points
  for(int i=0; i<npts; i++) {
    int id=i;
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  
  // read in triangle connectivities
  while( token != "polygons" ) ifs >> token;
  std::vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int,3> c;
  int ntri=0; ifs >> ntri;
  connectivities.reserve(ntri);
  cout << "calc_projection::Number of triangles: " << ntri << endl;
  int tmp; ifs>>tmp; 
  for (int i = 0; i<ntri; i++){
    int tmp=0;
    ifs >> tmp;
    if(tmp!=3) std::cout<<"calc_projection::Some mistake reading the elements connectivity from file. Check again."<<std::endl;
    for(int a=0; a<3; a++) ifs >> c[a];
    connectivities.push_back(c);
  }

    //! Centroid positions of coarse nodes
    std::vector<std::vector<double> > CoarseCentroids;
    std::vector<double> x(3,0.0);
    CoarseCentroids.resize(connectivities.size(),x);

    //update the neighbor search data structure
    for(int i=0;i<connectivities.size();i++){
      CoarseCentroids[i][0] = (defNodes[connectivities[i][0]]->getPoint(0) + defNodes[connectivities[i][1]]->getPoint(0) + defNodes[connectivities[i][2]]->getPoint(0) )/3.;

      CoarseCentroids[i][1] = (defNodes[connectivities[i][0]]->getPoint(1) + defNodes[connectivities[i][1]]->getPoint(1) + defNodes[connectivities[i][2]]->getPoint(1) )/3.;

      CoarseCentroids[i][2] = (defNodes[connectivities[i][0]]->getPoint(2) + defNodes[connectivities[i][1]]->getPoint(2) + defNodes[connectivities[i][2]]->getPoint(2) )/3.;
      }
      double spacing=3.;
      NodeBin grid(spacing,spacing,spacing,CoarseCentroids);

      //Now start a check for each fine node for nearest element in coarse mesh 
      for(int i=0; i<fineNodes.size();i++){
        const Vector3D & xstar_pos = fineNodes[i]->point();
        std::vector<double> temp(3,0.0); temp[0]=xstar_pos(0); temp[1]=xstar_pos(1); temp[2]=xstar_pos(2);
        //for each slave node, find the "close" neighbors from the search data structure
        std::vector<int> list;
        list=grid.GetList(temp);
        //find the closest of these neighbors
        double mindist = std::numeric_limits<double>::max();
        int minid=-1;
        double minN1=0.,minN2=0.,minN3=0.;
        if(list.size()==0) std::cout << "calc_projection::No neighbour was found::Problem" << std::endl;
        for(int j=0;j<list.size();j++){
          DefNode* x1=defNodes[connectivities[list[j]][0]]; const Vector3D & x1_pos = x1->point();
          DefNode* x2=defNodes[connectivities[list[j]][1]]; const Vector3D & x2_pos = x2->point();
          DefNode* x3=defNodes[connectivities[list[j]][2]]; const Vector3D & x3_pos = x3->point();
          //projection direction of node to the triangle is the normal of the triangle
          Vector3D dir = tvmet::cross(x1_pos-x2_pos,x1_pos-x3_pos);
          dir = dir/tvmet::norm2(dir);
          double dist=dot((xstar_pos-x1_pos),dir);
          //closest point
          Vector3D xhat; xhat = xstar_pos - dist*dir;
          Vector3D totalarea = tvmet::cross(x2_pos-x1_pos,x3_pos-x1_pos);
          double totalarea_normsq = tvmet::dot(totalarea,totalarea);
          double N1,N2,N3;
          N1 = tvmet::dot(tvmet::cross(x2_pos-xhat,x3_pos-xhat),totalarea)/totalarea_normsq;
          N2 = tvmet::dot(tvmet::cross(x3_pos-xhat,x1_pos-xhat),totalarea)/totalarea_normsq;
          N3 = tvmet::dot(tvmet::cross(x1_pos-xhat,x2_pos-xhat),totalarea)/totalarea_normsq;
          if(N1>-0.1 && N1<1.1 && N2>-0.1 && N2<1.1 && N3>-0.1 && N3<1.1){ //checking the barycentric coordinates
            if(mindist>std::abs(dist)) {mindist=std::abs(dist); minid=list[j];
                                        minN1=N1; minN2=N2; minN3=N3;}
          }
        }
        if(minid==-1) std::cout<<"calc_projection::No projection was found::Problem"<<std::endl;
        nodeid[i][0]=connectivities[minid][0];
        nodeid[i][1]=connectivities[minid][1];
        nodeid[i][2]=connectivities[minid][2];
        weights[i][0]=minN1;
        weights[i][1]=minN2;
        weights[i][2]=minN3;
      } 
}
