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
#include "WarpingPenalty.h"
#include "LoopShellBody.h"
#include "C0MembraneBody.h"
#include "SkewHexonBody.h"
#include "NormalPressure.h"
#include "ShapeTri3.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "CGfast.h"
#include "DeformableContact.h"
//#include "print-body.cc"

using namespace voom;

void print_body(Model::BodyContainer bdc);

/*! This program creates a shell model with LoopShellBody for bending,
   C0MembraneBodies for stretching and one NormalPressure for simulating the normal pressure
 */
int main(int argc, char* argv[])
{
  if(argc < 2){
    std::cout<<"Usage pressure modelName" << std::endl;
    return(0);
  }
 
  // get command-line arguments
  bool verbose=true;
  for(int i=0; i<argc; i++) {
    std::cout << std::setw(8) << i << "\t"
	      << argv[i] << std::endl;
  }
  std::string modelName = argv[1];
  std::ifstream ifs;
  int option;
  bool badCommandLine=false;
  bool restart=false;
  std::string restartFileName;
  while( (option=getopt(argc,argv,"r:")) != -1 ) {
    std::cout << "option = " << char(option) << std::endl;
    switch (option) {
      // -r restartfilename 
      //
      // whether restart or not with displacements 
      // added to the current position of nodes
    case 'r':
      restart=true;
      restartFileName = optarg;
      std::cout << "Restarting from file: " << restartFileName << std::endl;
      break;
      // 
      // No other command-line options.
      //
    case '?':
      std::cout << "Option " << char(optopt) << " is not supported." 
		<< std::endl;
      return 1;
    default :
      badCommandLine=true;
      break;
    }
  }

  //
  // Read in the mesh, in (ascii) legacy vtk format.
  //
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

  //different format connectivity for stretching body
  std::vector< std::vector<int> > s_connectivities;
  for(int i=0; i<connectivities.size(); i++) {
    std::vector<int> c(3);
    for(int j=0; j<3; j++) c[j]=connectivities[i](j);
    s_connectivities.push_back(c);
  }

  //read the two contact leaflets
  std::vector< std::vector<int> > large_leaflet_connect;
  while(token!="large_leaflet") ifs >> token;
  int datasize;
  ifs >> token >> datasize >> token;
  if(datasize!=ntri) std::cerr <<"Error while reading large leaflet component, data size not right" << std::endl;
  for(int i=0;i<ntri;i++){
    int flag;
    ifs >> flag;
    if(flag==1)large_leaflet_connect.push_back(s_connectivities[i]);
  }
  //small leaflet
  std::vector<int> small_leaflet;
  while(token!="small_leaflet") ifs >> token;
  ifs >> token >> datasize >> token;
  if(datasize!=npts) std::cerr <<"Error while reading large leaflet component, data size not right" << std::endl;
  for(int i=0;i<npts;i++){
    int flag;
    ifs >> flag;
    if(flag==1)small_leaflet.push_back(i);
  }

  //create a container with "variable" nodes (without the annulus)
  std::vector< NodeBase* > vary_nodes;
  std::vector<int> annulus;
  annulus.resize(npts,0);
  while(token!="annulus") ifs >> token;
  ifs >> token >> datasize >> token;
  if(datasize!=npts) std::cerr <<"Error while reading annulus component, data size not right" << std::endl;
  int nAnnulus=0;
  for(int i=0;i<npts;i++){
    int flag;
    ifs >> flag;
    annulus[i]=flag;
    if(flag==0) vary_nodes.push_back(defNodes[i]);
    else nAnnulus++;
  }
  std::cout<<nAnnulus<<" nodes are restrained at the annulus"<<std::endl;

  ifs.close();

  // Read in the starting displacement guess when restarting the simulation
  if(restart){
    std::ifstream guess(restartFileName.c_str());
    if (!guess) {
      cout << "Cannot open restart file: " << restartFileName << endl;
      return 0;
    }
    guess >> token >> token >> token >> token >> token;
    for(int i=0; i<npts; i++) {
      int id=i;
      double dx=0.,dy=0.,dz=0.;
      guess >> dx >> dy >> dz;
      defNodes[i]->addPoint(0,dx);
      defNodes[i]->addPoint(1,dy);
      defNodes[i]->addPoint(2,dz);
    }
    guess.close();
  }
  // default material constants

  // bending
  double Y=1.;
  double KC=0.1;
  double KG=-0.1;
  double nu=1./3.;
  double C0=0.;
  // stretching
  double mu=1.;
  double kS=2.;

  // now input (if file is available) user-defined material constants
  std::ifstream inp("parameters.inp");
  inp >> Y >> nu >> KC >> KG >> C0 >> mu >> kS;
  
  // print parameters out to screen
  if(verbose) 
    std::cout << " Y: " << Y << std::endl
	      << " nu: " << nu << std::endl
	      << " KC: " << KC << std::endl
	      << " KG: " << KG << std::endl
	      << " C0: " << C0 << std::endl
              << " mu: " << mu << std::endl
              << " kS: " << kS << std::endl;


  // create bending body; SCElastic only has curvature terms, no in-plane strain terms.
  int quadOrder=1;
  
  // bending material object, to be copied when body generates new elements
  typedef SCElastic BendingMaterial;
  BendingMaterial bending(KC,KG,C0);

  // Bending body will generate a mesh of Loop Subdivision shell
  // elements from the connectivity array and the container of node pointers.
  typedef LoopShellBody<BendingMaterial> LSB;
  LSB * bending_body = new LSB(bending, connectivities, nodes, quadOrder);

  // Use the body's output routines to generate vtk files with results
  bending_body->setOutput(paraview);
  bending_body->SetRefMeanCurvature();

  // Stretching body using EvansElastic Material, which separates area strain from pure
  // shear strain (at constant area).  

  typedef EvansElastic StretchingMaterial;
  StretchingMaterial stretch(0.,0.,0.,mu,kS);

  // Use generic C0MembraneBody for stretching of pents, with the EvansElastic material.
  typedef C0MembraneBody<TriangleQuadrature,StretchingMaterial,ShapeTri3> membrane;
  membrane * stretch_body = new membrane(stretch,s_connectivities,nodes,quadOrder);
  stretch_body->setOutput(paraview);

  //create the normal pressue element
  std::cout<<"Adding normal pressure element"<<std::endl;
  NormalPressure pressure_elem(defNodes,s_connectivities);
  bending_body->pushBack(&pressure_elem);
  std::cout<<"Added normal pressue element"<<std::endl;
  
  //create the deformable contact element
  //std::cout<<"Adding contact element"<<std::endl;
  DeformableContact contact(defNodes,large_leaflet_connect,small_leaflet,3,0.8,0.001);
  std::cout<<"Added contact element"<<std::endl;
  std::cout<<"Initial energy of the contact = " << contact.energy() << std::endl;
  bending_body->pushBack(&contact);
  //contact.printParaview("contact");

  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(bending_body);
  bdc.push_back(stretch_body);
  Model model(bdc,nodes);

  // consistency check
  //model.checkConsistency(true,false);
  //model.checkRank(model.dof()-6,true);
  //return 0;
  // set the solver parameters
  int m=5;
  double factr=1.0e7;
  double pgtol=1.5e-5;
  int iprint  =1;
  int maxIter=5000;
  std::ifstream lbfgsbinp("lbfgsb.inp");
  lbfgsbinp >> iprint >> factr >> pgtol >> m;
  if(verbose) 
    std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl;
  
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxIter );

  // set up bounds for solver
  blitz::Array<int,1> nbd(3*nodes.size());
  blitz::Array<double,1> lo(3*nodes.size());
  blitz::Array<double,1> hi(3*nodes.size());
  nbd = 2; //both upper and lower bound will be assigned
  lo = 0.0;
  hi = 0.0;
  double delta=5;

  for(int a=0; a<defNodes.size(); a++) {
    Vector3D pos;
    pos = defNodes[a]->point();
    if(annulus[a]==1){
      for(int i=0;i<3;i++){
        hi(3*a+i)=pos(i);
        lo(3*a+i)=pos(i);
      }
    }
    else{
      for(int i=0;i<3;i++){
        hi(3*a+i)=pos(i)+delta;
        lo(3*a+i)=pos(i)-delta;
      }
    }
  }

  solver.setBounds(nbd,lo,hi);
//////////////////////////////////////////////////////////

  // solve to minimize the energy
  //bdc[0]->printParaview("initial");
  double final_pressure=0.5,deltaP=5e-4;
  for(double pressure=0.0265+deltaP;std::abs(pressure)<=std::abs(final_pressure);pressure += deltaP){
    std::cout<<"Pressure = "<<pressure<<std::endl;
    pressure_elem.setPressure(pressure);
    solver.solve(&model);
  }
  //print the displacements
  std::ofstream ofs("disp.dat");
  if(!ofs)  std::cout<<"Cannot open output file for displacements."<<std::endl;
  ofs << "POINT_DATA " << defNodes.size() << std::endl
      << "VECTORS displacements double" << std::endl;
  //     << "LOOKUP_TABLE displacements" << endl; 
  for(int i=0; i<defNodes.size(); i++) {
    Vector3D nodalDisp;
    nodalDisp = defNodes[i]->point() - defNodes[i]->position();
    ofs << std::setprecision(16) 
        << nodalDisp(0) 
        << '\t' <<nodalDisp(1) 
        << '\t' <<nodalDisp(2) 
        << std::endl;
  }
  ofs.close();
  contact.printParaview("contact");
  //bdc[1]->printParaview("minimized");
  //print the postprocessing strains and stress
/*  hex_body->printParaviewEigVec("Hexbody");
  pentamer_body->printParaviewEigVec("Pentbody");

  hex_body->printParaview2("Hexbody2");
  pentamer_body->printParaview2("Pentbody2");

  //model.checkRank(model.dof()-6,true);

  //print the total energies of three bodies
  std::cout <<"Bending energy = "<<bending_body->totalStrainEnergy()<<std::endl;
  std::cout <<"Pentamer Stretching energy = "<<pentamer_body->totalStrainEnergy()<<std::endl;
  std::cout <<"Hexamer Stretching energy = "<<hex_body->totalStrainEnergy()<<std::endl;
*/
  return 0;
}
