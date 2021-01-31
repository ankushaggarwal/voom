//#include <string>
// This file reads in the two conformation states
// either in .dat or .pdf formats
// then calculates the deformation gradient at each node (w.r.t. the first state)
// and prints it's inverse to Fhatinv.dat
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

void readconfig_pdb(char* filename, float* x, float* y, float* z, int npts);
void readconfig_dat(char* filename, float* x, float* y, float* z, int npts);
int main(int argc, char* argv[])
{
  if( argc < 2 ) {
    cout << "Usage: calF inputFileName [-s support_size - k kernel_option]." << endl;
    return(0);
  }
  
  bool verbose=true;
  
  for(int i=0; i<argc; i++) {
    std::cout << std::setw(8) << i << "\t"
	      << argv[i] << std::endl;
  }
  //Input file name and other parameters from run command
  string objFileName = argv[1];
  ifstream ifs;
  string modelName = argv[1];
  string inputFileName = modelName;
  bool datfile = true;
  bool pdbfile = false;
  bool badCommandLine=false;
  double support_size;
  bool cubic_kernel=false, modified_kernel=false, singular_kernel=false;
  
  for(char option; (option=getopt(argc,argv,"n:N:s:k:")) != EOF; ) 
    switch (option) {
    default :
      badCommandLine=true;
      break;
    case 's' :
      support_size = std::atof(optarg);
      std::cout << "Global support size: " << support_size <<std::endl;
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
    cout << "Usage: homog inputFileName [-n start_configs -N end_config -s support_size - k kernel_option]." << endl;
    return(0);
  }
  
  // create input stream
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
  
  // read in points
  for(int i=0; i<npts; i++) {
    int id=i;
    int temp;
    double weight;
    double supphat;
    DeformationNode<3>::Point x;
    ifs >> temp >> x(0) >> x(1) >> x(2) >> weight>>supphat;
    
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
    node_volume.push_back(weight);
    supp_size.push_back(support_size);
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
  ofstream file1, file2, file3, fileF;
  file1.open("inv1.dat");
  file2.open("inv2.dat");
  file3.open("inv3.dat");
  fileF.open("Fhatinv.dat");
  //fileF << npts <<std::endl;
  //read in the configuration files in a loop
  for(int i=0;i<=0;i++){
    char inpfile[20];
    std::cout<<"Reading "<<i<<"th configuration."<<std::endl;
    float x[npts], y[npts], z[npts];
    if(datfile){
      //std::cout<<"using .dat files"<<std::endl;
      sprintf(inpfile,"%d.dat",i);
      readconfig_dat(inpfile,x,y,z,npts);}
    
    if(pdbfile){
      //std::cout<<"using .pdb files"<<std::endl;
      sprintf(inpfile,"../calpha-traj/%d.pdb",i);
      readconfig_pdb(inpfile,x,y,z,npts);} 
    std::cout<<"File has been read"<<std::endl;
    for(int inode=0;inode<npts;inode++){
      DeformationNode<3>::Point newx;
      newx(0)=x[inode]; newx(1)=y[inode]; newx(2)=z[inode];
      defNodes[inode]->setPoint(newx);
    }
    //std::vector<double> I1(npts,0.), I2(npts,0.), I3(npts,0.);
    std::cout<<"Calculating the invariants"<<std::endl;
    Tensor3D zero(0.);
    std::vector<Tensor3D> Fhat(npts,zero);
    std::vector<double> I1(npts,0.), I2(npts,0.), I3(npts,0.);
    bd.cal_defGrad(Fhat);
    bd.cal_invariants(I1,I2,I3);
    //bd.cal_stretch_ratios(I1,I2,I3);
    //calculate the inverse of F 
    //file1 << i; file2 << i; file3 << i;
    ofstream jac;
    jac.open("jacobian.dat");
    for(int inode=0;inode<npts;inode++){
      Tensor3D inv1;
      double det=determinant(Fhat[inode]);
      jac << det << std::endl;
      if(det<0.) std::cout<<"Determinant < 0 : "<<det<<std::endl;
      invert(Fhat[inode],inv1);
      fileF << "\t" << inv1(0,0) << "\t" << inv1(0,1) << "\t" << inv1(0,2) << "\t" << inv1(1,0) << "\t" << inv1(1,1) << "\t" << inv1(1,2) << "\t" << inv1(2,0) <<"\t" << inv1(2,1) << "\t" << inv1(2,2) << std::endl;
      //std::cout<<defNodes[inode]->getPoint(0)<<std::endl;
      	file1 << "\n" << I1[inode];
      	file2 << "\n" << I2[inode];
      	file3 << "\n" << det;
    }
    //file1 << std::endl; file2 << std::endl; file3 << std::endl;
    jac.close();
    // print out a paraview file
    char name[20];
    sprintf(name,"deformed%d",i);
    bd.print(name);
  }
  fileF.close();
  std::cout<<"Everything done!"<<std::endl;
  return 0;
}

void readconfig_dat(char* filename, float* x, float* y, float* z, int npts){
  ifstream inpfile;
  inpfile.open(filename);
  int id;
  
  for(int i=0; i<npts; i++){
    inpfile >> id >> x[i] >> y[i] >> z[i];
  }
  inpfile.close();
}

void readconfig_pdb(char* filename, float* x, float* y, float* z, int npts){
  ifstream inpfile;
  inpfile.open(filename);
  int i=0;
  string line;
  
  while (! inpfile.eof() )
    {
      getline (inpfile,line);
      if (line.length()>0 && line.substr(0,4)=="ATOM"){
	x[i] = atof(line.substr(30,8).c_str()); 
	y[i] = atof(line.substr(38,8).c_str()); 
	z[i] = atof(line.substr(46,8).c_str()); 
	i++;
      }
    }
	inpfile.close();
}
