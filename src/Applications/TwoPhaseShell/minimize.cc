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

  //
  // Read in the shear directions.  These are vectors in the current
  // configuration, pointing from vertex 1 to vertex 4 of each hexon,
  // as defined in Joe Rudnick's mathematica notebook. This vector
  // should point roughly along one of the edges of each triangular
  // element.
  //
  // Since we will define a "canonical" equilateral triangle reference
  // configuration for each element, we need to find which edge the
  // vector is closest to pointing along, and then choose that edge to
  // define the "shear angle" defining the reference pre-shear.
  //
  while( token != "shear_direction") ifs >> token; 
  std::vector<double> shear_angle;
  std::vector< std::vector<int> > pent_connectivities, hex_connectivities;
  ifs >> token;
  //
  // Since we are not using VTK libraries to do IO, sometimes need to
  // hack things to skip irrelevant tokens when format changes with a new mesh.
  //
  //int junk=0; ifs >> junk;
  //ifs >> token >> token;
  //ifs >> token ;
  std::cout<<"------------------------------------------------------"<<std::endl
           <<"REMEBER TO CHANGE IT IN THE CODE FOR USUAL DATA FORMAT"<<std::endl
           <<"------------------------------------------------------"<<std::endl;
  std::cout<<"Calculating and printing the shear angles"<<std::endl;
  for (int i=0; i<ntri; i++){
    double shear=0.;
    double TOL=2.5e-1;
    tvmet::Vector<double,3> shear_dir;
    ifs >> shear_dir(0) >> shear_dir(1) >> shear_dir(2);
    //std::cout<<shear_dir<<std::endl;
    double shear_dir_norm = tvmet::norm2(shear_dir); 
    //
    // pentamers don't have shear directions
    //
    if(shear_dir_norm==0) 
    {
      std::vector<int> c(3);
      for(int j=0;j<3;j++) c[j]=connectivities[i](j);
      pent_connectivities.push_back(c);
    }
    // 
    // hexamers do
    //
    else 
    {
      shear_dir = shear_dir/shear_dir_norm;
      std::vector<int> c(3);
      for(int j=0;j<3;j++) c[j]=connectivities[i](j);
      hex_connectivities.push_back(c);
      //
      // Calculate the shear angle
      //
      
      // edge vectors in current config.
      tvmet::Vector<double,3> e31(defNodes[c[0]]->point()-defNodes[c[2]]->point()), 
                              e32(defNodes[c[1]]->point()-defNodes[c[2]]->point()),
                              e12(defNodes[c[1]]->point()-defNodes[c[0]]->point());

      //normalize the edge vectors
      e31=e31/tvmet::norm2(e31);
      e32=e32/tvmet::norm2(e32);
      e12=e12/tvmet::norm2(e12);

      // Which edge is clossest to the shear vector?  Take cross
      // product.  If zero (or < TOL) then two vectors are parallel or
      // anti-parallel.  
      //
      // Compute angles in degrees.

      if(tvmet::norm2(tvmet::cross(shear_dir,e31))<TOL) {

        if(tvmet::dot(shear_dir,e31)>0) shear=0.; // parallel with edge 31
        else shear=180.; // anti-parallel

      } else if(tvmet::norm2(tvmet::cross(shear_dir,e32))<TOL){

        if(tvmet::dot(shear_dir,e32)>0) shear=60.; // parallel with 32
        else shear=240.; // anti-parallel

      } else if(tvmet::norm2(tvmet::cross(shear_dir,e12))<TOL) {

        if(tvmet::dot(shear_dir,e12)>0) shear=120.; // parallel with 12
        else shear=300.; // anti-parallel

      } 
      // If TOL is too small, it bombs.
      else {std::cout<<"Problem. The shear direction is not parallel to any triangle edge"<<std::endl; exit(1);}

      // Convert from degrees to radians.
      shear_angle.push_back(shear*PI/180.); 
      //std::cout<<shear<<std::endl;
    }
  }

  std::cout << "Number of penton elements = "<<pent_connectivities.size() << std::endl
            << "Number of hexon elements = "<<hex_connectivities.size() << std::endl;

  // Read in the starting displacement guess, required in some cases where solution is sensitive
  std::ifstream guess("initial-guess.dat");
  if(guess){
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
  double KC=1.;
  double KG=-1.;
  double nu=1./3.;
  double C0=0.;
  // stretching
  double mu=1.;
  double kS=1.;
  double eta=0.25;

  // now input (if file is available) user-defined material constants
  std::ifstream inp("parameters.inp");
  inp >> Y >> nu >> KC >> KG >> C0 >> mu >> kS;
  
  // assign pre-shear
  eta = in_eta;

  // print parameters out to screen
  if(verbose) 
    std::cout << " Y: " << Y << std::endl
	      << " nu: " << nu << std::endl
	      << " KC: " << KC << std::endl
	      << " KG: " << KG << std::endl
	      << " C0: " << C0 << std::endl
              << " mu: " << mu << std::endl
              << " kS: " << kS << std::endl
              << " eta:" << eta<<std::endl;


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

  // Create two stretching bodies, one for pents and one for hex's,
  // which don't have any bending energy (bending moduli set to zero)
  //
  // Use EvansElastic Material, which separates area strain from pure
  // shear strain (at constant area).  

  //
  // pentamers
  //

  typedef EvansElastic PentStretchingMaterial;
  PentStretchingMaterial stretch(0.,0.,0.,mu,kS);

  // Use generic C0MembraneBody for stretching of pents, with the EvansElastic material.
  typedef C0MembraneBody<TriangleQuadrature,PentStretchingMaterial,ShapeTri3> pent_membrane;
  pent_membrane * pentamer_body = new pent_membrane(stretch,pent_connectivities,nodes,quadOrder);
  pentamer_body->setOutput(paraview);


  //
  // hexamers
  //
  
  // Use specialezed SkewHexonBody (modified version of
  // C0MembraneBody), which has material hard-coded as
  // ModifiedEvansElastic, which takes care of the multiplicative F=AG
  // "Eigen-strain" decomposition.

  typedef SkewHexonBody<TriangleQuadrature,ShapeTri3> hex_membrane;

  // Correct the shear modulus if using double-well potential.
  //
  // mu = mu/(2*(eta*eta+c2)+2*c1);
  //
  // Not necessary if using single-well.
  hex_membrane * hex_body 
    = new hex_membrane(hex_connectivities,nodes,shear_angle,eta,mu,kS,c1,c2,quadOrder);

  // Define the edge length of the canonical 
  const double EDGELEN=av_edgelen;
  bending_body->SetRefConfiguration(EDGELEN);
  pentamer_body->SetRefConfiguration(EDGELEN);
  hex_body->SetRefConfiguration(EDGELEN);

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
  //print the postprocessing strains and stress
  hex_body->printParaviewEigVec("Hexbody");
  pentamer_body->printParaviewEigVec("Pentbody");

  hex_body->printParaview2("Hexbody2");
  pentamer_body->printParaview2("Pentbody2");

  //model.checkRank(model.dof()-6,true);

  //print the total energies of three bodies
  std::cout <<"Bending energy = "<<bending_body->totalStrainEnergy()<<std::endl;
  std::cout <<"Pentamer Stretching energy = "<<pentamer_body->totalStrainEnergy()<<std::endl;
  std::cout <<"Hexamer Stretching energy = "<<hex_body->totalStrainEnergy()<<std::endl;


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
  std::cout << "Edge Length = " << EDGELEN << std::endl;

  double deltaR=0.;
  for(int i=0; i<defNodes.size(); i++) {
      double Rtemp = tvmet::norm2(defNodes[i]->point());
      deltaR += (Rtemp-Ravg)*(Rtemp-Ravg);
    }
  double aspherity = sqrt(deltaR/nodes.size())/Ravg;
  std::cout<<"Aspherity = "<<aspherity<<std::endl;
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
