//#include <string>
// This file reads in point clouds at different time points 
// in .vtk format
// then calculates the deformation gradient at each node (w.r.t. the first time point)
// then calculates the three invariants at each node
// and writes them to the file invariants1.dat *2.dat and *3.dat
// in the format : at each configuration a row of invariants at all the nodes
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace std;
#include "MFShape.h"
#include "NodeBase.h"
#include "Node.h"
#include "MFBody.h"
#include "Model.h"
#include "CompNeoHookean.h"
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>
using namespace voom;

void readconfig_vtk(char* filename, float* x, float* y, float* z, int npts);

int main(int argc, char* argv[])
{
  if( argc < 2 ) {
    cout << "Usage: homog inputFileName [-n start_config -N end_config -s support_size - k kernel_option -a average_number]." << endl;
    return(0);
  }
  
  bool verbose=true;
  
  for(int i=0; i<argc; i++) {
    std::cout << std::setw(8) << i << "\t"
	      << argv[i] << std::endl;
  }
  //Input file name and other parameters from run command
  bool datfile = true;
  bool pdbfile = false;
  bool badCommandLine=false;
  double support_size, dx, dy, dz;
  int StartConfig=0, EndConfig=0, dConfig=1;
  bool cubic_kernel=false, modified_kernel=false, singular_kernel=false;
  
  for(char option; (option=getopt(argc,argv,"n:N:d:s:k:x:y:z:")) != EOF; ) 
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
    case 's' :
      support_size = std::atof(optarg);
      std::cout << "Global support size: " << support_size <<std::endl;
      break;
    case 'x' :
      dx = std::atof(optarg);
      std::cout << "Node spacing in x: " << dx <<std::endl;
      break;
    case 'y' :
      dy = std::atof(optarg);
      std::cout << "Node spacing in y: " << dy <<std::endl;
      break;
    case 'z' :
      dz = std::atof(optarg);
      std::cout << "Node spacing in z: " << dz <<std::endl;
      break;
    case 'k' :
      int kern=std::atoi(optarg);
      if( kern==1){
	cubic_kernel=true;
	std::cout<<"Main:: Cubic kernel will be used"<<std::endl;
      }
      else if(kern==2){
        modified_kernel=true;
        std::cout<<"Main:: Modified kernel will be used"<<std::endl;
      }
      else if(kern==3){
        singular_kernel=true;
        std::cout<<"Main:: Singular kernel will be used"<<std::endl;
      }
      else{
        std::cout<<"Only three kernel options"<<std::endl;
        exit(0);
      }
      break;
    }
  if( badCommandLine || argc < 2 ) {
    cout << "Usage: homog inputFileName [-n start_configs -N end_config -s support_size - k kernel_option -a average_number]." << endl;
    return(0);
  }

  double maxd, mind;
  maxd = dx > dy ? (dx>dz?dx:dz): (dy>dz? dy:dz);
  mind = dx < dy ? (dx<dz?dx:dz): (dy<dz? dy:dz);
  std::cout << "Max and min are " << maxd << " " << mind << std::endl;

  string objFileName = argv[optind];
  string modelName = argv[optind];
  string inputFileName = modelName;
  // create input stream
  ifstream ifs;
  ifs.open( inputFileName.c_str(), ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << inputFileName << endl;
    exit(0);
  }
  // create vector of nodes
  int dof=0;
  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;
  vector<double> node_volume,supp_size, supp_size_hat;
  double supphatfactor=0.9;
  if(modified_kernel) std::cout<<"A factor of "<<supphatfactor<<" is used for scaling the support size for primitive kernel"<<std::endl;
  char key;  ifs>>key;
  
  // Input .dat file containing nodes and weights
  std::string token;
  ifs >> token;
  while( token != "POINTS" ) ifs >> token;
  int npts=0;
  ifs >> npts; 
  defNodes.reserve(npts);
  ifs >> token;
  
  // read in points
  for(int i=0; i<npts; i++) {
    int id=i;
    int temp;
    double weight=dx*dy*dz;
    double supphat=mind;
    DeformationNode<3>::Point x;
    //ifs >> temp >> x(0) >> x(1) >> x(2) >> weight>> supphat;
    ifs >> x(0) >> x(1) >> x(2); // >> weight>> supphat;
    
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
    node_volume.push_back(weight);
    supp_size.push_back(8.*maxd);
    supp_size_hat.push_back(supphat*supphatfactor);
  }
  cout << "Number of nodes: " <<nodes.size() << endl;
  ifs.close();
  // create material, values don't matter here
  double rho = 1.0;
  double E = 200.0;
  double nu = 0.4;
  double viscosity=1.0e-2;
  typedef CompNeoHookean MaterialType;
  MaterialType protein( rho, E, nu );
  
  // create Body
  typedef MFBody<MaterialType,MFShape> body;
  body bd;
  
  if(cubic_kernel)
    bd.initializeBody(protein,nodes,node_volume,supp_size);
  if(modified_kernel)
    bd.initializeBody(protein,nodes,node_volume,supp_size,supp_size_hat);
  if(singular_kernel)
    bd.initializeBody(protein,nodes,node_volume,supp_size,1.5);
  
  std::cout << "Body has been created and the shapefunctions have been calculated"<<std::endl;
  bd.setOutput(paraview);
  bd.compute(true,true,false);
  
  std::cout << std::endl;
  std::cout << "Body energy = " << bd.energy() << std::endl;
  std::cout << "Body volume = " << bd.volume() << std::endl;
  std::cout << std::endl;
  
  //open the output files
  ofstream file1, file2, file3;
  file1.open("inv1.dat",ios::binary);
  file2.open("inv2.dat",ios::binary);
  file3.open("inv3.dat",ios::binary);
  char tab='\t', newline='\n';
  int col=npts+1, row=(EndConfig-StartConfig)/dConfig+1;
  file1.write((char*)&row,sizeof(int)); file1.write((char*)&col,sizeof(int));
  file2.write((char*)&row,sizeof(int)); file2.write((char*)&col,sizeof(int));
  file3.write((char*)&row,sizeof(int)); file3.write((char*)&col,sizeof(int));
  float x[npts], y[npts], z[npts];

  //read in the configuration files in a loop
  for(int i=StartConfig;i<=EndConfig;i+=dConfig){
    char inpfile[100];
    std::cout<<"Reading "<<i<<"th configuration."<<std::endl;

    //sprintf(inpfile,"../bav04-point-clouds/seg09_to_%02d_bav04_root_for_strain.vtk",i);
    sprintf(inpfile,"%s%02d%s",inputFileName.substr(0,9).c_str(),i,inputFileName.substr(11).c_str());
    std::cout << inpfile << std::endl;
    readconfig_vtk(inpfile,x,y,z,npts);
    
    std::cout<<"File has been read"<<std::endl;
    for(int inode=0;inode<npts;inode++){
      DeformationNode<3>::Point newx;
      newx(0)=x[inode]; newx(1)=y[inode]; newx(2)=z[inode];
      defNodes[inode]->setPoint(newx);
    }
    std::vector<double> I1(npts,0.), I2(npts,0.), I3(npts,0.);
    std::cout<<"Calculating the invariants"<<std::endl;
    bd.cal_invariants(I1,I2,I3);
    //bd.cal_stretch_ratios(I1,I2,I3);
    double di=double(i);
    std::cout<<"double of i is "<<di<<std::endl;
    file1.write((char*)&di,sizeof(double));
    file2.write((char*)&di,sizeof(double));
    file3.write((char*)&di,sizeof(double));
    
    //file1 << i; file2 << i; file3 << i;
    for(int inode=0;inode<npts;inode++){
      //std::cout<<defNodes[inode]->getPoint(0)<<std::endl;
      //	file1 << "\t" << I1[inode];
      //	file2 << "\t" << I2[inode];
      //	file3 << "\t" << I3[inode];
      //file1.write(&tab,sizeof(char));
      //file2.write(&tab,sizeof(char));
      //file3.write(&tab,sizeof(char));
      file1.write((char*)&I1[inode],sizeof(double)); 
      file2.write((char*)&I2[inode],sizeof(double));
      file3.write((char*)&I3[inode],sizeof(double));
    }
    //file1 << std::endl; file2 << std::endl; file3 << std::endl;
    //file1.write(&newline,sizeof(char));
    //file2.write(&newline,sizeof(char));
    //file3.write(&newline,sizeof(char));

    sprintf(inpfile,"strain%02d.vtk",i);
    ofstream ofs(inpfile); 
    
    //output in vtk format
    ofs << endl << "POINT_DATA " << npts << endl;
    ofs << "SCALARS I1 float" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for(int inode=0;inode<npts;inode++){
      ofs << I1[inode] << endl;
    }
    ofs << "SCALARS I2 float" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for(int inode=0;inode<npts;inode++){
      ofs << I2[inode] << endl;
    }
    ofs << "SCALARS J float" << endl;
    ofs << "LOOKUP_TABLE default" << endl;
    for(int inode=0;inode<npts;inode++){
      ofs << I3[inode] << endl;
    }
    
    ofs.close();

  }
  std::cout<<"Everything done!"<<std::endl;
  return 0;
}

void readconfig_vtk(char* filename, float* x, float* y, float* z, int npts){
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
  inpfile.close();
}

