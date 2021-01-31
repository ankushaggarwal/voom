// -*- C++ -*-
//----------------------------------------------------------------------
//
//                Melissa M. Gibbons & William S. Klug
//                University of California Los Angeles
//                   (C) 2006-2008 All Rights Reserved
//
//----------------------------------------------------------------------

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <tvmet/Vector.h>
#include "Node.h"
#include "CompNeoHookean.h"
#include "Capsid3DBody.h"
#include "TetQuadrature.h"
#include "ShapeTet4CP.h"
#include "Model.h"
#include "Solver.h"
#include "Lbfgsb.h"
#include "SimulatedAnnealing.h"
//#include "Rigid.h" 
#include "Contact.h"
#include "ViscousRegularizer.h"

using namespace tvmet;
using namespace std;
using namespace voom;

int main(int argc, char* argv[])
{
  if( argc < 2 ) {
      std::cout << "Usage: normal-mode modelName [-n nNormalModes -e eigenvectors/no-eigenvectors]." << std::endl;
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
  string inputFileName = modelName + ".vtk";

  bool badCommandLine=false;
  int nModes;
  bool eigenvector=false;

  for(char option; (option=getopt(argc,argv,"n:e:")) != EOF; ) 
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
    default :
      badCommandLine=true;
      break;
    }

  if( badCommandLine || argc < 2 ) {
      cout << "Usage: capsid InputFileName [-c -H finalHeight -b boundaryConditionFlag -R restartFilename]." << endl;
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
  char key;  ifs>>key;
  
  // Input .vtk file containing nodes and connectivities
  std::string token;
  ifs >> token; 
  while( token != "POINTS" ) ifs >> token;
  int npts=0;
  ifs >> npts; 
  defNodes.reserve(npts);
  ifs >> token;// skip number type

  // read in points
  for(int i=0; i<npts; i++) {
    int id=i;
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    //x /= 10.0;
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  cout << "Number of nodes: " <<nodes.size() << endl;

  // read in tetrahedral connectivities
  while( token != "CELLS" ) ifs >> token;
  std::vector< std::vector<int> > connectivities;
  std::vector<int> c(4);
  int ntet=0; ifs >> ntet;
  connectivities.reserve(ntet);
  cout << "Number of tetrahedrons: " << ntet << endl;
  int ntmp=0;
  ifs >> ntmp;
  if(ntmp != 5*ntet) {
    cout << "Non-tetrahedral elements?" << endl;
    return 0;
  }
  for (int i = 0; i<ntet; i++){
    int tmp=0;
    ifs >> tmp;
    ifs >> tmp; c[0]=tmp;
    ifs >> tmp; c[1]=tmp;
    ifs >> tmp; c[2]=tmp;
    ifs >> tmp; c[3]=tmp;
    connectivities.push_back(c);
  }
  if(verbose) cout << connectivities.size() << endl;

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
  if(verbose) std::cout << "Material has been created." << std::endl;

  // create Body
  unsigned int quadOrder = 1;
  
  typedef Capsid3DBody<TetQuadrature,MaterialType,ShapeTet4> Capsid;
  Capsid bd(protein, connectivities, nodes, quadOrder);

  bd.setOutput(paraview);
  bd.compute(true,true,false);

  ifs.close();

  std::cout << std::endl;
  std::cout << "Body energy = " << bd.energy() << std::endl;
  std::cout << "Body volume = " << bd.volume() << std::endl;
  std::cout << std::endl;

  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(&bd);
  Model model(bdc, nodes);

  //model.checkConsistency(true,false);
  //model.checkRank(model.dof()-6,true);
  //return 0;

  std::cout << "Calculating normal modes now." << std::endl;
//  model.normalModes_all(eigenvector,nModes,node_volume,rho,inputFileName);
  rho=1.;
  std::vector<double> node_volume; node_volume.resize(npts,1.);
  model.normalModes_all(eigenvector,nModes,node_volume,rho,inputFileName);
  return 0;
}

