//#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include "Contact.h"
#include "ViscousRegularizer.h"

#if defined(_OPENMP)
#include <omp.h>
#endif


// using namespace tvmet;
 using namespace std;
#include "MFShape.h"
//#include "MFShape.cc"
#include "NodeBase.h"
#include "Node.h"
#include "MFBody.h"
#include "Model.h"
#include "CompNeoHookean.h"
#include "MooneyRivlin.h"
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>
#include "Solver.h"
#include "Lbfgsb.h"
 using namespace voom;

int main(int argc, char* argv[])
{
  if( argc < 2 ) {
      std::cout << "Usage: normal-mode modelName [-n nNormalModes -e eigenvectors/no-eigenvectors -s support size]." << std::endl;
      return 0;
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

  bool badCommandLine=false;
  double support_size;
  int nModes;
  bool eigenvector=false;
  bool cubic_kernel=false, modified_kernel=false, singular_kernel=false;

  for(char option; (option=getopt(argc,argv,"n:e:s:k:")) != EOF; ) 
    switch (option) {
    case 'n' :
      nModes = std::atoi(optarg);
	std::cout << "Number of normal modes to be computed: "<<nModes<< std::endl;
      break;
    case 'e' :
      if( std::string(optarg) == std::string("yes") ) {
	eigenvector = true;
	std::cout << "Eigenvectors will also be computed." << std::cout;
      }
      else if( std::string(optarg) != std::string("no") ) {
	std::cout << "Yes or no for eigenvectors.  Exiting."
		  << std::endl;
	return 0;
      }	
      break;
    case 'k' :{
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
      }}
      break;
    default :
      badCommandLine=true;
      break;
    case 's' :
      support_size = std::atof(optarg);
      std::cout << "Global support size: " << support_size << std::endl;
      break;
    }

  if( badCommandLine || argc < 2 ) {
      std::cout << "Usage: capsid InputFileName [-c -H finalHeight -b boundaryConditionFlag -R restartFilename -r Eratio -s SupportSize]." << std::endl;
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
  std::vector<double> node_volume, supp_size, supp_size_hat;
  char key;  ifs>>key;
  double supphatfactor=0.9;
  if(modified_kernel) std::cout<<"A factor of "<<supphatfactor<<" is used for scaling the support size for primitive kernel"<<std::endl;
  
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
    ifs >> temp >> x(0) >> x(1) >> x(2) >> weight >> supphat;
    //x /= 10.0;
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
    //node_volume.push_back(weight/1000);
    node_volume.push_back(weight);
    supp_size.push_back(support_size);
    supp_size_hat.push_back(supphat*supphatfactor);
  }
  cout << "Number of nodes: " <<nodes.size() << endl;

  // create material
  double rho = 1.0;
  double E = 200.0;
  double nu = 0.4;
  double viscosity=1.0e-2;

  ifstream inp("parameters.inp");
  inp >> E >> nu >> viscosity;
  if(verbose) 
    std::cout << "Input rho: " << rho << std::endl
	      << "Input E: " << E << std::endl
	      << "Input nu: " << nu << std::endl
	      << "Input viscosity: " << viscosity << std::endl;

  typedef CompNeoHookean MaterialType;
  MaterialType protein( rho, E, nu );
  //typedef MooneyRivlin MaterialType;
  //MaterialType protein( E, 0, 8./3.*E );
  //MaterialType protein( 0.1471, 0.1226, 0.7010 );
  if(verbose) std::cout << "Material has been created." << std::endl;

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
  
  ifs.close();

  std::cout << std::endl;
  std::cout << "Body energy = " << bd.energy() << std::endl;
  std::cout << "Body volume = " << bd.volume() << std::endl;
  std::cout << std::endl;

  //set the node wise modulus
  /*std::vector<double> C10(npts,0.),C01(npts,0.),D1(npts,0.);
  ifs.open("first.mat");
  for(int i=0;i<npts;i++){
    ifs >> C10[i];
    if(C10[i]<0) C10[i]=0.;}
  ifs.close(); 
  ifs.open("second.mat");
  for(int i=0;i<npts;i++){
    ifs >> C01[i]; 
    if(C01[i]<0) C10[i]=0.;}
  ifs.close(); 
  ifs.open("bulk.mat");
  for(int i=0;i<npts;i++){
    ifs >> D1[i]; 
    if(D1[i]<0) D1[i]=0.;}
  ifs.close(); 
  bd.setModulus(C10,C01,D1);*/

  //density calculated from atomic mass
  rho=17298.712055/bd.volume()*1.66e-27;
  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(&bd);
  Model model(bdc, nodes);
  //node_volume.clear(); node_volume.resize(npts,bd.volume()/npts);
  //model.checkConsistency(true,false);

  //std::cout << "Calculating rank now." << std::endl;
  //model.checkRank(model.dof()-6,true);
  //std::cout << "Rank calculated." << std::endl;
  //return 0;
  std::cout << "Calculating normal modes now." << std::endl;
//  model.normalModes_all(eigenvector,nModes,node_volume,rho,inputFileName);
  rho=1.;
  node_volume.clear(); node_volume.resize(npts,1.);
  model.normalModes(eigenvector,nModes,node_volume,rho,inputFileName);
  //model.normalModes_all(eigenvector,nModes,node_volume,rho,inputFileName);
  return 0;
}
