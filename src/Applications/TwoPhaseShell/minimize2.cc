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

using namespace voom;

//! This function creates a shell model with LoopShellBody  for bending and 2 C0MembraneBodies for stretching, one of pentamers and one of hexamers. Hexmaers' modified form of stretching energy which has two-minima is used to model phase change of HK97 from symmetric state to sheared state

const double EDGELEN=2.0;

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
  const double Rcapsid = 30.;
  int option;
  bool badCommandLine=false;
  std::string modelName = argv[1];
  std::ifstream ifs;

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
  // create input stream
  ifs.open( inputFileName.c_str(), std::ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << inputFileName << endl;
    exit(0);
  }

  // create vector of nodes
  int dof=0;
  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;
  double Ravg = 0;
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
    Ravg += tvmet::norm2(x);
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  assert(nodes.size()!=0);
  Ravg /= nodes.size();
  std::cout << "Number of nodes: " <<nodes.size() << std::endl
            << "Ravg = " << Ravg << std::endl;
  
  // set reference edgelength 
  C0Membrane::setRefLength(2.0);
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

  // read in the shear directions
  while( token != "shear_direction") ifs >> token;
  std::vector<double> shear_angle;
  std::vector< std::vector<int> > pent_connectivities, hex_connectivities;
  ifs >> token;
  //int junk=0; ifs >> junk;
  //ifs >> token >> token;
  //ifs >> token ;
  std::cout<<"------------------------------------------------------------------------------------"<<std::endl
           <<"REMEBER TO CHANGE IT IN THE CODE FOR USUAL DATA FORMAT"<<std::endl
           <<"------------------------------------------------------------------------------------"<<std::endl;
  std::cout<<"Calculating and printing the shear angles"<<std::endl;
  for (int i=0; i<ntri; i++){
    double shear=0.;
    double TOL=2e-1;
    tvmet::Vector<double,3> shear_dir;
    ifs >> shear_dir(0) >> shear_dir(1) >> shear_dir(2);
    //std::cout<<shear_dir<<std::endl;
    double shear_dir_norm = tvmet::norm2(shear_dir); 
    if(shear_dir_norm==0) //pentamer
    {
      std::vector<int> c(3);
      for(int j=0;j<3;j++) c[j]=connectivities[i](j);
      pent_connectivities.push_back(c);
    }
    else  //hexamer
    {
      shear_dir = shear_dir/shear_dir_norm;
      std::vector<int> c(3);
      for(int j=0;j<3;j++) c[j]=connectivities[i](j);
      hex_connectivities.push_back(c);
      //calculate the shear angle
      tvmet::Vector<double,3> e31(defNodes[c[0]]->point()-defNodes[c[2]]->point()), 
                              e32(defNodes[c[1]]->point()-defNodes[c[2]]->point()),
                              e12(defNodes[c[1]]->point()-defNodes[c[0]]->point());
      //normalize the edge vectors
      e31=e31/tvmet::norm2(e31);
      e32=e32/tvmet::norm2(e32);
      e12=e12/tvmet::norm2(e12);
      if(tvmet::norm2(tvmet::cross(shear_dir,e31))<TOL)
        if(tvmet::dot(shear_dir,e31)>0) shear=0.;
        else shear=180.;
      else if(tvmet::norm2(tvmet::cross(shear_dir,e32))<TOL)
        if(tvmet::dot(shear_dir,e32)>0) shear=60.;
        else shear=240.;
      else if(tvmet::norm2(tvmet::cross(shear_dir,e12))<TOL)
        if(tvmet::dot(shear_dir,e12)>0) shear=120.;
        else shear=300.;
      else {std::cout<<"Problem. The shear direction is not parallel to any triangle edge"<<std::endl; exit(1);}
      shear_angle.push_back(shear*PI/180.); 
      //std::cout<<shear<<std::endl;
    }
  }
  std::cout << "Number of penton elements = "<<pent_connectivities.size() << std::endl
            << "Number of hexon elements = "<<hex_connectivities.size() << std::endl;

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
  std::ifstream inp("parameters.inp");
  inp >> Y >> nu >> KC >> KG >> C0 >> mu >> kS >> eta;
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
  typedef EvansElastic PentStretchingMaterial;
  PentStretchingMaterial stretch(0.,0.,0.,mu,kS);

  typedef C0MembraneBody<TriangleQuadrature,PentStretchingMaterial,ShapeTri3> pent_membrane;
  pent_membrane * pentamer_body = new pent_membrane(stretch,pent_connectivities,nodes,quadOrder);
  pentamer_body->setOutput(paraview);

  // hexamers
  typedef SkewHexonBody<TriangleQuadrature,ShapeTri3> hex_membrane;
  // correct the shear modulus
  mu = mu/(2*(eta*eta+c2)+2*c1);
  hex_membrane * hex_body = new hex_membrane(hex_connectivities,nodes,shear_angle,eta,mu,kS,c1,c2,quadOrder);

  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(bending_body);
  bdc.push_back(pentamer_body);
  bdc.push_back(hex_body);
  Model model(bdc,nodes);

  // consistency check
  //model.checkConsistency(true,false);
  //model.checkRank(model.dof()-6,true);
  //return 0;
  // set the solver parameters
  int m=5;
  double factr=1.0e1;
  double pgtol=1.5e-5;
  int iprint  =0;
  int maxIter=5000;
  std::ifstream lbfgsbinp("lbfgsb.inp");
  lbfgsbinp >> iprint >> factr >> pgtol >> m;
  if(verbose) 
    std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl;
  
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxIter );
  // solve to minimize the energy
  bdc[0]->printParaview("initial");
  solver.solve(&model);
  bdc[0]->printParaview("minimized");
  //model.checkRank(model.dof()-6,true);

  //calculate the center
  tvmet::Vector<double,3> center(0,0,0);
  for(int i=0; i<defNodes.size(); i++) {
    center += defNodes[i]->point();
  }
  center /= nodes.size();
  std::cout << "Center = " << center << std::endl;

  //calculate the average radius and aspherity
  Ravg = 0;
  for(int i=0; i<defNodes.size(); i++) {
    Ravg += tvmet::norm2(defNodes[i]->point());
  }
  Ravg /= nodes.size();
  std::cout << "Ravg = " << Ravg << std::endl;
  std::cout << "Edge Length = " << EDGELEN << std::endl;

  double deltaR=0.;
  for(int i=0; i<defNodes.size(); i++) {
      double Rtemp = tvmet::norm2(defNodes[i]->point());
      deltaR += (Rtemp-Ravg)*(Rtemp-Ravg);
    }
  double aspherity = sqrt(deltaR/nodes.size())/Ravg;
  std::cout<<"Aspherity = "<<aspherity<<std::endl;
}
