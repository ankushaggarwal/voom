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
#include "Lbfgsb.h"
#include "CGfast.h"

using namespace std;
using namespace voom;

//! This function creates a shell model with LoopShellBody  for bending and  C0MembraneBodies for stretching. Then it calculates energy for each SH applying a small deformation

void SHdisp(int l,int m,double x,double y,double z,std::vector<double>& dx);
int main(int argc, char* argv[])
{
  if(argc < 2){
    std::cout<<"Usage SHenergy modelName [-k KC -m MU -f OUTFILENAME]." << std::endl;
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
  ifs.open( inputFileName.c_str(), std::ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << inputFileName << endl;
    exit(0);
  }
  bool badCommandLine=false;
  double inp_KC=0.,inp_mu=0.;
  std:string outfilename="default_result.dat";

  for(char option; (option=getopt(argc,argv,"k:m:f:")) != EOF; ) 
    switch (option) {
    default:
      badCommandLine=true;
      break;
    case 'k' :
      inp_KC = std::atof(optarg);
      std::cout << "KC: " << inp_KC << std::endl;
      break;
    case 'm' :
      inp_mu = std::atof(optarg);
      std::cout << "mu: " << inp_mu << std::endl;
      break;
    case 'f' :
      outfilename = optarg;
      std::cout << "Result will be output in "<<outfilename << std::endl;
    }
  if( badCommandLine || argc < 2 ) {
    std::cout<<"Usage SHenergy modelName [-k KC -m MU -f OUTFILENAME]." << std::endl;
    return(0);
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
  // material constants and their default values
  // bending
  double KC=600;
  double KG=-1.;
  double C0=0;
  // stretching
  double mu=500;
  double kS=2.*mu;

  // now input user-defined material constants
  //std::ifstream inp("parameters.inp");
  KC = inp_KC;
  mu = inp_mu;
  KG=-KC; //KG value does not change the result, but adds stability
  kS=2*mu; //using nu=1/3, we get this result
  // print parameters out to screen
  if(verbose) 
    std::cout << " KC: " << KC << std::endl
	      << " KG: " << KG << std::endl
	      << " C0: " << C0 << std::endl
              << " mu: " << mu << std::endl
              << " kS: " << kS << std::endl;

  // create bending body
  int quadOrder=1;
  
  typedef EvansElastic BendingMaterial;
  BendingMaterial bending(KC,KG,C0,mu,kS);

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

  bending_body->compute(true,false,false);
  stretch_body->compute(true,false,false);
  double initBE=bending_body->energy(),initSE=stretch_body->energy();
  std::cout << std::endl;
  std::cout << "Bending body energy = " << bending_body->energy() << std::endl;
  std::cout << "Stretching body energy = " << stretch_body->energy() << std::endl;
  std::cout << std::endl;

  // create Model
  //std::string name="resultSH";
  Model::BodyContainer bdc;
  bdc.push_back(bending_body);
  //bdc.push_back(stretch_body);
  Model model(bdc,nodes);
  
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

  //print the total energies of three bodies
  std::cout << std::endl;
  std::cout << "Bending body energy = " << bending_body->energy() << std::endl;
  std::cout << "Stretching body energy = " << stretch_body->energy() << std::endl;
  std::cout << std::endl;


  //print the displacements
  std::ofstream ofs("disp.dat");
  if(!ofs)  std::cout<<"Cannot open output file for displacements."<<std::endl;
  ofs << "POINT_DATA " << defNodes.size() << std::endl
      << "VECTORS displacements double" << std::endl;
  /*     << "LOOKUP_TABLE displacements" << endl; */
  for(int i=0; i<defNodes.size(); i++) {
    Vector3D nodalDisp;
    nodalDisp = defNodes[i]->point() - defNodes[i]->position();
    ofs << std::setprecision(16) 
        << nodalDisp(0)
        << '\t' <<nodalDisp(1)
        << '\t' <<nodalDisp(2) << std::endl;
  }
  ofs.close();
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

  double deltaR=0.;
  for(int i=0; i<defNodes.size(); i++) {
      double Rtemp = tvmet::norm2(defNodes[i]->point());
      deltaR += (Rtemp-Ravg)*(Rtemp-Ravg);
    }
  double aspherity = sqrt(deltaR/nodes.size())/Ravg;
  std::cout<<"Aspherity = "<<aspherity<<std::endl;
  
  // Calculate spherical harmonics
  model.normalModes_SH(true,400,outfilename );
  //model.normalModes(true,10,name);
  //return 0;

  //apply deformation of spherical harmonics and calculate the energy
/*  double mag=Ravg*1e-2;
  for(int l=0;l<5;l++){
    //analytical energy
    double ana_energy=1./4.*((8.*kS)+KC/Ravg/Ravg*l*(l-1.)*(l+1.)*(l+2.))*mag*mag;
    std::cout<<std::endl<<"l= "<<l<<" and Analytical energy="<<ana_energy<<std::endl;
  }*/
  std::cout<<"Everything done!"<<std::endl;
  return 0;
}

