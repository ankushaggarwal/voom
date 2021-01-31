#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>
#include "MFShape.h"
#include "NodeBase.h"
#include "Node.h"
#include "MFBody.h"
#include "Model.h"
#include "CompNeoHookean.h"
#include "MooneyRivlin.h"
#include "Solver.h"
#include "Lbfgsb.h"
#include "Contact.h"
#include "ViscousRegularizer.h"
#include "RigidHemisphereAL.h"
#include "RigidPlateAL.h"
#include "RadialPressure.h"
#include "FixCenter.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace voom;
using namespace std;

int main(int argc, char* argv[])
{
  if( argc < 2 ) {
      cout << "Usage: pressurize_load InputFileName [-s support_size - k kernel_option]." << endl;
      return(0);
  }

  for(int i=0; i<argc; i++) {
    std::cout << std::setw(8) << i << "\t"
	      << argv[i] << std::endl;
  }

//Input file name and other parameters from run command
  string objFileName = argv[1];
  ifstream ifs;
  string modelName = argv[1];
  string inputFileName = modelName;
  double inside_radius=6.;
  bool badCommandLine=false;
  double support_size;
  bool cubic_kernel=false, modified_kernel=false, singular_kernel=false;

  for(char option; (option=getopt(argc,argv,"s:k:")) != EOF; ) 
    switch (option) {
    case 's' :
      support_size = std::atof(optarg);
      std::cout << "Global support size: " << support_size << std::endl;
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
    }

  if( badCommandLine || argc < 2 ) {
      cout << "Usage: pressurize_load InputFileName [-s support_size - k kernel_option]." << endl;
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
  std::vector< DeformationNode<3>* > defNodes,inside_nodes;
  std::vector<double> node_volume, supp_size, supp_size_hat, node_area;
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
    ifs >> temp >> x(0) >> x(1) >> x(2) >> weight>>supphat;
 
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
    if(sqrt(x(0)*x(0)+x(1)*x(1)+x(2)*x(2))<inside_radius*1.01)
      inside_nodes.push_back(n);
    node_volume.push_back(weight);
    double radius=pow(weight*3/4./M_PI,1./3.);
    node_area.push_back(M_PI*radius*radius);
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
  std::cout << "Input rho: " << rho << std::endl
	      << "Input E: " << E << std::endl
	      << "Input nu: " << nu << std::endl
	      << "Input viscosity: " << viscosity << std::endl;

  //typedef CompNeoHookean MaterialType;
  //MaterialType protein( rho, E, nu );
  typedef MooneyRivlin MaterialType;
  MaterialType protein( E, 0, 100.*E );
  std::cout << "Material has been created." << std::endl;

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

  //add the pressure load element
  RadialPressure inside_pressure(inside_nodes,node_area);
  bd.addElement(&inside_pressure);

  //add a penalty function to prevent rigid body motion (of the center)
  FixCenter fixcenter(defNodes,10000.);
  bd.addElement(&fixcenter);

  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(&bd);
  Model model(bdc, nodes);

  //model.checkConsistency(true,false);
  //model.checkRank(model.dof()-6,true);
  //return 0;

  //Input lbfgsb solver parameters
  int m=5;
  double factr=1.0e+7;
  double pgtol=1.0e-5;
  int iprint = 0;
  ifstream lbfgsbinp("lbfgsb.inp");
  lbfgsbinp >> iprint >> factr >> pgtol >> m ;
	 
  std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl;
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint );

  ofstream PvsR;
  PvsR.open("PvsR.dat");

  // %%%%%%%%%%%%%%%%%%%%%%
  // Begin loading loop
  // %%%%%%%%%%%%%%%%%%%%%%
  double increment=0.005;
  int Iter=0;
  for(double factor=0.;factor<1.5;factor+=increment){
  inside_pressure.setP(E*factor);
  solver.solve(&model);
  if(Iter%4==0){
  char name[20];
  sprintf(name,"disp%d",Iter/4);
  model.print(name);}
  Iter++;
  double r=inside_pressure.MeanRadius();
  PvsR<<factor<<"\t"<<r<<std::endl;
  }
  return 0;
}
