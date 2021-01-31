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
#include "Contact.h"
#include "ViscousRegularizer.h"
#include "RigidHemisphereAL.h"
#include "RigidPlateAL.h"
//#include "print-body.cc"

using namespace voom;

void print_body(Model::BodyContainer bdc);
double findmin(std::vector< DeformationNode<3>* > &defNodes, int dim, double rad=0.);
double findmax(std::vector< DeformationNode<3>* > &defNodes, int dim, double rad=0.);
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
  double Hfinal=0;
  bool compressFlag=false;
  int option;
  bool badCommandLine=false;
  std::string modelName = argv[1];
  std::ifstream ifs;

  while( (option=getopt(argc,argv,"c:d:n:H:")) != -1 ) {
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
    case 'H' :
      Hfinal = std::atof(optarg);
      std::cout << "Hfinal: " << Hfinal << std::endl;
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
  double KC=1.;
  double KG=-1.;
  double C0=0.;
  // stretching
  double mu=1.;
  double kS=1.;
  double eta=0.25;
  double viscosity=0; //1.0e-9;

  // now input (if file is available) user-defined material constants
  std::ifstream inp("parameters.inp");
  if(inp) inp >> KC >> KG >> C0 >> mu >> kS;
  
  // assign pre-shear
  eta = in_eta;

  // print parameters out to screen
  if(verbose) 
    std::cout << " KC: " << KC << std::endl
	      << " KG: " << KG << std::endl
	      << " C0: " << C0 << std::endl
              << " mu: " << mu << std::endl
              << " kS: " << kS << std::endl
              << " viscosity: " << viscosity << std::endl
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
  //const double EDGELEN=2.0;
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
  int maxIter=20000;
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

  compressFlag=true;
  bool extrapolate=false;
  if( !compressFlag ) {
    std::cout << "All done.  Bye now." << std::endl;
    return 0;
  }
  // Start indenting
  std::cout << "Compressing capsid." << std::endl;

  // create indentor and plate
  double k1 = 1.0e4;
  double afmR = 10.0;
  double afmX = 0.0;
  double afmY = 0.0;
  double friction_coeff = 0.0;

  //inp >> k1 >> afmR >> afmX >> afmY >> friction_coeff;
  if(verbose) 
    std::cout << "Input k1: " << k1 << std::endl
	      << "Input afmR: " << afmR << std::endl
	      << "Center position of afm tip (x): " << afmX << std::endl
	      << "Center position of afm tip (y): " << afmY << std::endl
	      << "Friction Coefficient: " << friction_coeff << std::endl;

  //stiffness of indentor
  double k_AL=k1;
  // find top and bottom of capsid
  double Zmin1=findmin(defNodes,1);
  double Zmax1=findmax(defNodes,1,afmR);
  std::cout << "Zmin1 = " << Zmin1 << std::endl
	    << "Zmax1 = " << Zmax1 << std::endl;

  tvmet::Vector<double,3> xc; xc = afmX, afmY, Zmax1; 

  RigidHemisphereAL * afm 
    = new RigidHemisphereAL( defNodes, k_AL, afmR, xc, friction_coeff );
  afm->updateContact();
  
  bending_body->pushBack( afm );
  std::cout << "Added afm to body." << std::endl;

  bool up=true;
  RigidPlateAL * glass
    = new RigidPlateAL( defNodes, k_AL, Zmin1, up, friction_coeff );
  bending_body->pushBack( glass );
  std::cout << "Added glass substrate to the body." << std::endl;

  //add some viscosity for regularization
  ViscousRegularizer vr(bending_body->nodes(), viscosity);
  bending_body->pushBack( &vr ); 

  int steps = 10;
  double vrTol = pgtol;
  inp >> steps;
  //inp >> vrTol;
  if(verbose) 
    std::cout << "Input steps: "<<steps<<std::endl
	      << "Input vrTol: "<<vrTol<<std::endl;


  const double originalHeight=Zmax1-Zmin1;
  double Zbegin=Zmin1;
  double Zend=Zmax1-Hfinal;
  double dZ = (Zend-Zbegin)/steps;
  double vrEnergy = vr.energy();
  double bdEnergy = bending_body->energy()+hex_body->totalStrainEnergy()+pentamer_body->totalStrainEnergy();
  

  std::ofstream FvsZ;
  if(verbose) FvsZ.open("FvsZ.dat");
  int Iter = 0;
  
  std::ofstream restart;
  restart.open("all-input-parameters-for-restart.dat");
  restart 
           << std::setw( 24 ) << std::setprecision(16) << Hfinal
           << std::setw( 24 ) << std::setprecision(16) << eta
           << std::setw( 24 ) << std::setprecision(16) << av_edgelen
           << std::setw( 24 ) << std::setprecision(16) << KC
           << std::setw( 24 ) << std::setprecision(16) << KG
           << std::setw( 24 ) << std::setprecision(16) << C0
           << std::setw( 24 ) << std::setprecision(16) << mu
           << std::setw( 24 ) << std::setprecision(16) << kS
           << std::setw( 24 ) << std::setprecision(16) << viscosity
           << std::setw( 24 ) << std::setprecision(16) << quadOrder
           << std::setw( 24 ) << std::setprecision(16) << m
           << std::setw( 24 ) << std::setprecision(16) << factr
           << std::setw( 24 ) << std::setprecision(16) << pgtol
           << std::setw( 24 ) << std::setprecision(16) << vrTol
           << std::setw( 24 ) << std::setprecision(16) << iprint
           << std::setw( 24 ) << std::setprecision(16) << maxIter
           << std::setw( 24 ) << std::setprecision(16) << k_AL
           << std::setw( 24 ) << std::setprecision(16) << afmR
           << std::setw( 24 ) << std::setprecision(16) << afmX
           << std::setw( 24 ) << std::setprecision(16) << afmY
           << std::setw( 24 ) << std::setprecision(16) << friction_coeff
           << std::setw( 24 ) << std::setprecision(16) << Zmax1
           << std::setw( 24 ) << std::setprecision(16) << Zmin1
           << std::setw( 24 ) << std::setprecision(16) << dZ
           << std::setw( 24 ) << std::setprecision(16) << steps;
  restart.close();
//  added for extrapolation of initial guess
  blitz::Array<double,1> x_prev(model.dof());
  blitz::Array<double,1> u_prev(model.dof());
  model.getField(solver);
  for(int i=0; i<model.dof(); i++ ) x_prev(i) = solver.field(i);
  u_prev = 0.0;

  // %%%%%%%%%%%%%%%%%%%%%%
  // Begin indentation loop
  // %%%%%%%%%%%%%%%%%%%%%%
  //ofs.open("restart1.dat",std::ios::binary);
  ofs.open("restart1.dat");
  if(!ofs)  std::cout<<"Cannot open output file for displacements."<<std::endl;
  ofs << 0 << std::endl;
  for(int i=0; i<defNodes.size(); i++) {
    Vector3D nodalDisp;
    nodalDisp = defNodes[i]->point() - defNodes[i]->position();
    ofs << std::setprecision(16) 
        << nodalDisp(0)
        << '\t' <<nodalDisp(1)
        << '\t' <<nodalDisp(2) << std::endl;
    //ofs.write((char*)&nodalDisp(0),sizeof(double));
    //ofs.write((char*)&nodalDisp(1),sizeof(double));
    //ofs.write((char*)&nodalDisp(2),sizeof(double));
  }
  ofs.close();

  for(double Z = Zbegin; Z<Zend+0.5*dZ; Z+=dZ) {
    
    // move glass up by dZ
    glass->setZ(Z);

    // add scaled version of previous displacement as an initial guess
    if(extrapolate){
    for(int i=0; i<model.dof(); i++ ) {
      if (!afm->active_node(i/3)){
          if(!glass->active_node(i/3)){
             solver.field(i) = x_prev(i) + 0.99*dZ*u_prev(i);
          }
          else if(i%3==2)
             solver.field(i) = x_prev(i) + 0.99*dZ*u_prev(i);
     }
     else solver.field(i) = x_prev(i);
    }
    model.putField(solver);
    }

    std::cout << std::endl
	      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	      << std::endl
	      << " Z = " << Z 
	      << " zeta = " << Z-Zbegin 
	      << " Iter = " << Iter 
	      << std::endl;

    int viterMax = 5;
    for(int viter=0; viter<viterMax; viter++) {

      if(verbose) std::cout << std::endl 
			    << "VISCOUS ITERATION: " << viter 
			    << "\t viscosity = " << vr.viscosity()
			    << std::endl
			    << std::endl;
//       vr.step();
      glass->updateContact();
      afm->updateContact();
//??      model.computeAndAssemble(solver,false,true,false);
      solver.solve( &model );
      vrEnergy = vr.energy();
      bdEnergy = bending_body->energy()+hex_body->totalStrainEnergy()+pentamer_body->totalStrainEnergy();

      if(verbose) {
	std::cout << "ENERGY:" << std::endl
		  << "viscous energy = " << vrEnergy << std::endl
		  << "   body energy = " << bdEnergy << std::endl
		  << "  total energy = " << solver.function() << std::endl
		  << std::endl;
      }

      // update contact
      if(verbose) {
	std::cout << "CONTACT:" << std::endl;
	std::cout << "    top active = " << afm->active() << std::endl
		  << " bottom active = " << glass->active() << std::endl
		  << "     top force = " << afm->FZ() << std::endl 
		  << "  bottom force = " << glass->FZ() << std::endl
 		  << "  top penetration  = " << afm->penetration() << std::endl
 		  << "bottom penetration = " << glass->penetration() << std::endl;
      }
	
      // step forward in "time", relaxing viscous energy & forces 
      vr.step();

      if(vrEnergy < vrTol*bdEnergy && solver.projectedGradientNorm()<=pgtol && std::max( afm->penetration(), glass->penetration() ) < 1.0e-2*std::abs(dZ) ) { 
	// viscous energy is small enough; exit
	break;
      }
    }

    double height = Zmax1-Z;
    // add up forces on top and bottom
    double Ztop = Zmax1;
    double Zbot = Zmax1 - height;       

    if(verbose) {
      std::cout << "Contact and Viscous energy converged." << std::endl
		<< "        height = " << height << std::endl
		<< "   indentation = " << originalHeight-height << std::endl
		<< std::endl;
      }

    if(extrapolate){
    // Compute displacement
    for(int i=0; i<model.dof(); i++ ) {
      u_prev(i) = ( solver.field(i) - x_prev(i) )/dZ;
    }

    // Save current (successful) state as previous
    for(int i=0; i<model.dof(); i++ ) x_prev(i) = solver.field(i);
    }

    // This is where Paraview file (if desired) is saved
    char name[20]; 
    if(Iter%4==0){
      sprintf(name,"disp%d",Iter/4);
      model.print(name);
    }
Iter++;


    // Output the resulting energy and forces
    if(verbose) {
      FvsZ << std::setw( 24 ) << std::setprecision(16) 
           << originalHeight-height
           << std::setw( 24 ) << std::setprecision(16)
	   << glass->FZ()
	   << std::setw( 24 ) << std::setprecision(16)
           << afm->FZ()
	   << std::setw( 24 ) << std::setprecision(16)
           << solver.function()
	   << std::endl;
	}

    //print the displacements for restart
    std::ofstream ofs;
    if(Iter%2==0) ofs.open("restart.dat",std::ios::binary);
    //if(Iter%2==0) ofs.open("restart.dat");
    else ofs.open("restart1.dat",std::ios::binary);
    //else ofs.open("restart1.dat");
    if(!ofs)  std::cout<<"Cannot open output file for displacements."<<std::endl;
    ofs << Iter << std::endl;
    for(int i=0; i<defNodes.size(); i++) {
      Vector3D nodalDisp;
      nodalDisp = defNodes[i]->point() - defNodes[i]->position();
      double dx=nodalDisp(0), dy=nodalDisp(1), dz=nodalDisp(2);
      //ofs << std::setprecision(16) 
      //    << nodalDisp(0)
      //    << '\t' <<nodalDisp(1)
      //    << '\t' <<nodalDisp(2) << std::endl;
      ofs.write((char*)&dx,sizeof(double));
      ofs.write((char*)&dy,sizeof(double));
      ofs.write((char*)&dz,sizeof(double));
    }
    ofs.close();
  }
  FvsZ.close();

  std::cout << "Indentation Complete." << std::endl;

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
double findmin(std::vector< DeformationNode<3>* > &defNodes, int dim, double rad){

  if(rad==0){ 
  double Zmin=std::numeric_limits<double>::max();
  for (int a=0; a<defNodes.size(); a++){
    double Z = defNodes[a]->getPoint(2);
    Zmin = std::min(Zmin, Z);
  }

  std::cout << "Zmin = " << Zmin << std::endl;
    return Zmin;
   }
  //if(rad!=0)
  double Zmin1=std::numeric_limits<double>::max();

  for (int a=0; a<defNodes.size(); a++){
    double Z = defNodes[a]->getPoint(2);
    double r=sqrt(defNodes[a]->getPoint(0)*defNodes[a]->getPoint(0) + defNodes[a]->getPoint(1)*defNodes[a]->getPoint(1));
    ////this probably needs correction
    Zmin1 = std::min(Zmin1, sqrt(rad*rad-r*r)+Z);
  }

  std::cout << "Zmin1 = " << Zmin1 << std::endl;
  return Zmin1;
}

double findmax(std::vector< DeformationNode<3>* > &defNodes, int dim, double rad){

  if(rad==0){ 
  double Zmax=-std::numeric_limits<double>::max();
  for (int a=0; a<defNodes.size(); a++){
    double Z = defNodes[a]->getPoint(2);
    Zmax = std::max(Zmax, Z);
  }

  std::cout << "Zmax = " << Zmax << std::endl;
    return Zmax;
   }

  double Zmax1=-std::numeric_limits<double>::max();  
  //if(rad!=0)
  for (int a=0; a<defNodes.size(); a++){
    double Z = defNodes[a]->getPoint(2);
    double r=sqrt(defNodes[a]->getPoint(0)*defNodes[a]->getPoint(0) + defNodes[a]->getPoint(1)*defNodes[a]->getPoint(1));
    Zmax1 = std::max(Zmax1, sqrt(rad*rad-r*r)+Z);
  }

  std::cout << "Zmax1 = " << Zmax1 << std::endl;
  return Zmax1;
}
