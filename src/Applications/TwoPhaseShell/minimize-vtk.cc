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
#include "SkewHexonBody.h"
#include "ShapeTri3.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "CGfast.h"

#include <vtkCell.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>

using namespace voom;

//! This function creates a shell model with LoopShellBody  for bending and 2 C0MembraneBodies for stretching, one of pentamers and one of hexamers. Hexmaers' modified form of stretching energy which has two-minima is used to model phase change of HK97 from symmetric state to sheared state

int main(int argc, char* argv[])
{
  if(argc < 2){
    std::cout<<"Usage minimize modelName [-c c1 -d c2]." << std::endl;
    return(0);
  }
  
  bool verbose=true;
  for(int i=0; i<argc; i++) {
    std::cout << std::setw(8) << i << "\t"
	      << argv[i] << std::endl;
  }
  const double PI=3.14159265;
  double c1=0.,c2=0.;
  const double Rcapsid = 60;
  int option;
  bool badCommandLine=false;
  std::string modelName = argv[1];

  while( (option=getopt(argc,argv,"c:d:")) != -1 ) {
    std::cout << "option = " << char(option) << std::endl;
    switch (option) {
    case 'c':
      c1 = std::atof(optarg);
      std::cout << "c1: " << c1 << std::endl;
      break;
    case 'd':
      c2 = std::atof(optarg);
      std::cout << "c2: " << c2 << std::endl;
      break;
    case '?':
      std::cout << "Option " << char(optopt) << " is not supported." 
		<< std::endl;
      return 1;
    default :
      badCommandLine=true;
      break;
    }
  }
  std::string inputFileName = modelName + ".vtk";
  vtkDataSetReader * reader = vtkDataSetReader::New();
  reader->SetFileName( inputFileName.c_str() );
  // send through normals filter to ensure that triangle orientations
  // are consistent
  vtkPolyDataNormals * normals = vtkPolyDataNormals::New();
  normals->SetInput( reader->GetOutput() );
  normals->ConsistencyOn();
  normals->SplittingOff();
  normals->AutoOrientNormalsOn();
  vtkPolyData * mesh = normals->GetOutput();
  mesh->Update();
  std::cout << "mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
	    << std::endl;

  // create vector of nodes
  int dof=0;
  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;
  double Ravg = 0;
  
  // read in points
  for(int a=0; a<mesh->GetNumberOfPoints(); a++) {
    int id=a;
    DeformationNode<3>::Point x;
    mesh->GetPoint(a, &(x[0]));
    Ravg += tvmet::norm2(x);
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  assert(nodes.size()!=0);
  Ravg /= nodes.size();
  cout << "Number of nodes: " <<nodes.size() << endl
       << "Ravg = " << Ravg << endl;

  // read in triangle connectivities
  std::vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int, 3> c;
  int ntri=mesh->GetNumberOfCells();
  connectivities.reserve(ntri);
  if(verbose) cout << "Number of triangles: " <<ntri << endl;
  for (int i = 0; i<ntri; i++){
    assert(mesh->GetCell(i)->GetNumberOfPoints() == 3);
    for(int a=0; a<3; a++) c[a] = mesh->GetCell(i)->GetPointId(a);
    connectivities.push_back(c);
  }

//std::cout<<reader->GetOutput()->GetCellData()->GetValue(0)<<std::endl;
 /* // read in the shear directions
vtkCellData *cd = reader->GetCellData("shear_direction").GetValue(0)

  while( token != "shear_direction") ifs >> token;
  std::vector<double> shear_angle;
  std::vector< std::vector<int> > pent_connectivities, hex_connectivities;
  ifs >> token;
  int junk=0; ifs >> junk;
  ifs >> token >> token;
  for (int i=0; i<ntri; i++){
    double shear=0.;
    ifs >> shear;
    if(shear<0) //pentamer
    {
      std::vector<int> c(3);
      for(int j=0;j<3;j++) c[j]=connectivities[i](j);
      pent_connectivities.push_back(c);
    }
    else  //hexamer
    {
      std::vector<int> c(3);
      for(int j=0;j<3;j++) c[j]=connectivities[i](j);
      hex_connectivities.push_back(c);
      shear_angle.push_back(shear*PI/180.);
    }
  }
  std::cout << "Number of penton elements = "<<pent_connectivities.size() << std::endl
            << "Number of hexon elements = "<<hex_connectivities.size() << std::endl;
*/
  // rescale size to make the capsid of radius in nanometers
    for(int i=0; i<defNodes.size(); i++) {
      DeformationNode<3>::Point x;
      x = defNodes[i]->point();
      x *= Rcapsid/Ravg;
      defNodes[i]->setPoint(x);
      defNodes[i]->setPosition(x);
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
  double eta=0.25;
  if(verbose) 
    std::cout << " Y: " << Y << std::endl
	      << " nu: " << nu << std::endl
	      << " KC: " << KC << std::endl
	      << " KG: " << KG << std::endl
	      << " C0: " << C0 << std::endl
              << " mu: " << mu << std::endl
              << " kS: " << kS << std::endl
              << " eta:" << eta<<std::endl;

  // create bending body
  int quadOrder=1;
  
  typedef SCElastic BendingMaterial;
  BendingMaterial bending(KC,KG,C0);

  typedef LoopShellBody<BendingMaterial> LSB;
  LSB * bending_body = new LSB(bending, connectivities, nodes, quadOrder);
  bending_body->setOutput(paraview);

  // create stretching body which doesn't have any bending energy
  // pentamers
/*  typedef EvansElastic PentStretchingMaterial;
  PentStretchingMaterial stretch(0.,0.,0.,mu,kS);

  typedef C0MembraneBody<TriangleQuadrature,PentStretchingMaterial,ShapeTri3> pent_membrane;
  pent_membrane * pentamer_body = new pent_membrane(stretch,pent_connectivities,nodes,quadOrder);
  pentamer_body->setOutput(paraview);

  // hexamers
  typedef SkewHexonBody<TriangleQuadrature,ShapeTri3> hex_membrane;
  hex_membrane * hex_body = new hex_membrane(hex_connectivities,nodes,shear_angle,eta,mu,kS,c1,c2,quadOrder);

  // create Model
  Model::BodyContainer bdc;
  //bdc.push_back(bending_body);
  bdc.push_back(pentamer_body);
  bdc.push_back(hex_body);
  Model model(bdc,nodes);

  // set the solver parameters
  int m=5;
  double factr=1.0e1;
  double pgtol=1.5e-5;
  int iprint  =0;
  int maxIter=500;
  std::ifstream lbfgsbinp("lbfgsb.inp");
  lbfgsbinp >> iprint >> factr >> pgtol >> m;
  if(verbose) 
    std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl;
  
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxIter );
  // solve to minimize the energy
  solver.solve(&model);*/
}
