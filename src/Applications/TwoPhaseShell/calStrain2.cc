#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <unistd.h>
#include "Node.h"
#include "SCElastic.h"
#include "EvansElastic.h"
#include "LoopShellBody.h"
#include "C0MembraneBody.h"
#include "ShapeTri3.h"
#include "Model.h"

using namespace std;
#include <tvmet/Vector.h>
using namespace voom;

//! This function creates a shell model with LoopShellBody  for bending and  C0MembraneBodies for stretching. Then it reads the configurations and calculates invariants and curvatures. The difference from calStrain is that this version transfers the data from cells to nodes before writing the binary files

void readconfig_vtk(char* filename, float* x, float* y, float* z, int npts, double* thick=NULL);
void cell_to_nodes(std::vector<double> &data, std::vector<double> &data2, std::vector<std::vector<int> > &conn);

int main(int argc, char* argv[])
{
  if(argc < 2){
    std::cout<<"Usage homog2D modelName [-n start_config -N end_config]." << std::endl;
    return(0);
  }
  
  bool verbose=true;
  for(int i=0; i<argc; i++) {
    std::cout << std::setw(8) << i << "\t"
	      << argv[i] << std::endl;
  }

  bool badCommandLine=false;
  int StartConfig=0, EndConfig=0, dConfig=1;
  
  for(char option; (option=getopt(argc,argv,"n:N:d:")) != EOF; ) 
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
    case 'd' :
      dConfig = std::atoi(optarg);
      std::cout << "Increment of configurations: " << dConfig << std::endl;
      break;
    }

  if( badCommandLine || argc < 2 ) {
    std::cout<<"Usage homog2D modelName tableFile [-n start_config -N end_config]." << std::endl;
    return(0);
  }

  // create input stream
  std::ifstream ifs;
  std::string inputFileName = argv[optind];
  ifs.open( inputFileName.c_str(), std::ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << inputFileName << endl;
    exit(0);
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
  while( token != "polygons" ) ifs >> token;
  std::vector< tvmet::Vector<int,3> > connectivities;
  std::vector< std::vector<int> > reverse_connectivities(npts,std::vector<int>(0));
  tvmet::Vector<int,3> c;
  int ntri=0; ifs >> ntri;
  connectivities.reserve(ntri);
  cout << "Number of triangles: " << ntri << endl;
  int tmp; ifs>>tmp; std::cout<<" tmp = "<<tmp<<std::endl;
  for (int i = 0; i<ntri; i++){
    int tmp=0;
    ifs >> tmp;
    if(tmp!=3) std::cout<<"Some mistake reading the elements connectivity from file. Check again."<<std::endl;
    for(int a=0; a<3; a++) 
    {
        ifs >> c[a];
        reverse_connectivities[c[a]].push_back(i);
    }
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

  //read the trajectory and output curvatures and invariants
  std::ofstream file1, file2, file3, file4;
  file1.open("inv12.dat",std::ios::binary);
  file2.open("inv22.dat",std::ios::binary);
  file3.open("curv2.dat",std::ios::binary);
  file4.open("thick2.dat",std::ios::binary);

//  char tab='\t', newline='\n';
  int col=npts+1, row=(EndConfig-StartConfig)/dConfig+1;

  file1.write((char*)&row,sizeof(int)); file1.write((char*)&col,sizeof(int));
  file2.write((char*)&row,sizeof(int)); file2.write((char*)&col,sizeof(int));
  file3.write((char*)&row,sizeof(int)); file3.write((char*)&col,sizeof(int));
  file4.write((char*)&row,sizeof(int)); file4.write((char*)&col,sizeof(int));

  //read in the configuration files in a loop
  for(int i=StartConfig;i<=EndConfig;i+=dConfig){
    char inpfile[100];
    std::cout<<"Reading "<<i<<"th configuration."<<std::endl;
    float x[npts], y[npts], z[npts];
    double thick[npts];

    sprintf(inpfile,"%s%d%s",inputFileName.substr(0,16).c_str(),i,".vtk");
    readconfig_vtk(inpfile,x,y,z,npts,thick);
    std::cout<<inpfile<<" File has been read"<<std::endl;

    for(int inode=0;inode<npts;inode++){
      DeformationNode<3>::Point newx;
      newx(0)=x[inode]; newx(1)=y[inode]; newx(2)=z[inode];
      defNodes[inode]->setPoint(newx);
    }
    
    std::vector<double> I1(ntri,0.), I2(ntri,0.), curv(ntri,0.);
    std::vector<double> I12(npts,0.), I22(npts,0.), curv2(npts,0.);
    std::cout<<"Calculating the invariants"<<std::endl;
    bending_body->compute(true,false,false);
    bending_body->cal_curv(curv);
    stretch_body->cal_inv(I1,I2);
    cell_to_nodes(I1,I12,reverse_connectivities);
    cell_to_nodes(I2,I22,reverse_connectivities);
    cell_to_nodes(curv,curv2,reverse_connectivities);
    //bd.cal_stretch_ratios(I1,I2,I3);
    double di=double(i);
    //std::cout<<"double of i is "<<di<<std::endl;
    file1.write((char*)&di,sizeof(double));
    file2.write((char*)&di,sizeof(double));
    file3.write((char*)&di,sizeof(double));
    file4.write((char*)&di,sizeof(double));
    
    for(int inode=0;inode<npts;inode++){
      file1.write((char*)&I12[inode],sizeof(double)); 
      file2.write((char*)&I22[inode],sizeof(double));
      file3.write((char*)&curv2[inode],sizeof(double));
      file4.write((char*)&thick[inode],sizeof(double));
    }
    sprintf(inpfile,"strain%02d.vtk",i);
    ofstream ofs(inpfile); 
    
    //output in vtk format
    ofs << endl << "CELL_DATA " << ntri << endl;
    ofs << "SCALARS I1 float" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for(int iele=0;iele<ntri;iele++){
      ofs << I1[iele] << endl;
    }
    ofs << "SCALARS I2 float" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for(int iele=0;iele<ntri;iele++){
      ofs << I2[iele] << endl;
    }
    ofs << "SCALARS curv float" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for(int iele=0;iele<ntri;iele++){
      ofs << curv[iele] << endl;
    }
    
    ofs.close();

  }
  std::cout<<"Everything done!"<<std::endl;
  return 0;
}

void readconfig_vtk(char* filename, float* x, float* y, float* z, int npts, double* thick){
  ifstream inpfile;
  inpfile.open(filename);
  if (!inpfile) {
    cout << "Cannot open input file: " << filename << endl;
    exit(0);
  }
  int id;

  // Input .dat file containing nodes and weights
  std::string token;
  inpfile >> token;
  while( token != "POINTS" ) inpfile >> token;
  inpfile >> id;
  if (id != npts) {
    cout << "Error: the number of nodes is not the same" << endl;
    exit(0);
  }


  inpfile >> token;

  for(int i=0; i<npts; i++){
    inpfile >> x[i] >> y[i] >> z[i];
  }
  
  if(!thick){
      inpfile.close();
      return; //if thickness array is NULL pointer
  }

  while( token != "Thickness" ) inpfile >> token;
  inpfile >> token;
  inpfile >> id;
  if (id != npts) {
    cout << "Error: the number of nodes is not the same" << endl;
    exit(0);
  }
  inpfile >> token;

  for(int i=0; i<npts; i++){
    inpfile >> thick[i];
  }
  inpfile.close();
}

void cell_to_nodes(std::vector<double> &data, std::vector<double> &data2, std::vector<std::vector<int> > &conn){
    for(int i=0;i<conn.size();i++){
        data2[i]=0;
        for(int j=0;j<conn[i].size();j++)
            data2[i] += data[conn[i][j]];
        data2[i] /= conn[i].size();
    }
    return;
}
