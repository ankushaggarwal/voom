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
//#include "print-body.cc"

using namespace voom;

void print_body(Model::BodyContainer bdc);
//const double EDGELEN=2.0;

/*! This program creates a shell model with LoopShellBody for bending
  and 2 C0MembraneBodies for stretching, one of pentamers and one of
  hexamers. Hexmaers' modified form of stretching energy which has
  two-minima is used to model phase change of HK97 from symmetric
  state to sheared state
 */
int main(int argc, char* argv[])
{
  if(argc < 2){
    std::cout<<"Usage minimize modelName [-c c1 -d c2]." << std::endl;
    return(0);
  }
 
  // get command-line arguments
  bool verbose=true;
  for(int i=0; i<argc; i++) {
    std::cout << std::setw(8) << i << "\t"
	      << argv[i] << std::endl;
  }
  const double PI=3.14159265;
  double c1=0.,c2=0.,in_eta=0.;
  const double Rcapsid = 30.;
  int option;
  bool badCommandLine=false;
  std::string modelName = argv[1];
  std::ifstream ifs;

  while( (option=getopt(argc,argv,"c:d:n:")) != -1 ) {
    std::cout << "option = " << char(option) << std::endl;
    switch (option) {
      // -c c1 
      //
      // gives the depth of the first minimum of a double-well energy,
      // if such an energy is used. If single-well energy is used,
      // then this is irrelevant.
    case 'c':
      c1 = std::atof(optarg);
      std::cout << "c1: " << c1 << std::endl;
      break;
      // -d c2 
      //
      // gives the depth of the second minimum of a double-well energy,
      // if such an energy is used. If single-well energy is used,
      // then this is irrelevant.
    case 'd':
      c2 = std::atof(optarg);
      std::cout << "c2: " << c2 << std::endl;
      break;
      // -n eta
      //
      // The prescribed pre-shear value.
    case 'n':
      in_eta=std::atof(optarg);
      std::cout << "input eta: " << in_eta << std::endl;
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
  double av_edgelen=0.;
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
    //calculate the three edge lengths
    const DeformationNode<3>::Point& x1=defNodes[c[0]]->point(), x2=defNodes[c[1]]->point(), x3=defNodes[c[2]]->point();
    av_edgelen += (tvmet::norm2(x1-x2) + tvmet::norm2(x1-x3) + tvmet::norm2(x3-x2) )/3.;
  }
  av_edgelen /= ntri;

  // Read in the starting displacement guess, required in some cases where solution is sensitive
  std::ifstream guess("initial-guess.dat");
  if(guess){
    std::cout << "Reading an initial guess of displacement."<<std::endl;
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
  double KC=1.;
  double KG=-1.;
  double C0=0.;
  // stretching
  double mu=100.;
  double kS=200.;

  // now input (if file is available) user-defined material constants
  std::ifstream inp("parameters.inp");
  if(inp)
  inp >> KC >> KG >> C0 >> mu >> kS;
  
  // print parameters out to screen
  if(verbose) 
    std::cout << " KC: " << KC << std::endl
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
  //bending_body->SetRefMeanCurvature();

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

  // Define the edge length of the canonical 
  const double EDGELEN=av_edgelen;
  bending_body->SetRefConfiguration(EDGELEN);
  stretch_body->SetRefConfiguration(EDGELEN);

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

  //print the total energies of three bodies
  std::cout <<"Bending energy = "<<bending_body->totalStrainEnergy()<<std::endl;
  std::cout <<"Stretching energy = "<<stretch_body->totalStrainEnergy()<<std::endl;


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
  return 0;
}

void print_body(Model::BodyContainer bdc){
  std::string fileName = "name.vtk";
  std::ofstream ofs(fileName.c_str());
  if (!ofs) {
    std::cout << "can not open output ("
              << fileName
	      << ") file." << std::endl;
    exit(0);
  }

  int nbd=bdc.size();
  if(nbd<0) return;
  int npts=bdc[0]->nodes().size();
  for(int i=0;i<nbd;i++)
    if(bdc[i]->nodes().size() != npts)
      std::cout<<"print_body Error: Bodies have different node numbers"<<std::endl;

  //    Node Section
  ofs << "# vtk DataFile Version 2.0\n"
      << "Test example" << std::endl
      << "ASCII" << std::endl
      << "DATASET POLYDATA" << std::endl
      << "POINTS  " << npts << "  double" << std::endl;
    
  // output nodal postions
  /*Body::NodeContainer nodes=bdc[0]->nodes();
  for(int i=0;i<nbd;i++){
      //DeformationNode<3>*
      NodeBase n=nodes[i];
      const Vector3D & nodalPos =  (*n)->position();
      ofs << std::setprecision(16) 
	  << nodalPos(0) << "  "
	  << nodalPos(1) << "  "
	  << nodalPos(2) << std::endl;
    }
  */
  return;
}
