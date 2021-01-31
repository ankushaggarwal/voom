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
#include "Contact.h"
#include "ViscousRegularizer.h"
#include "RigidHemisphereAL.h"
#include "RigidPlateAL.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace voom;

double findmin(std::vector< DeformationNode<3>* > &defNodes, int dim, double rad=0.);
double findmax(std::vector< DeformationNode<3>* > &defNodes, int dim, double rad=0.);

int main(int argc, char* argv[])
{
  if( argc < 2 ) {
      std::cout << "Usage: indent InputFileName [-c compress/no -H final_height]." << std::endl;
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

  bool compressFlag=false;
  double Hfinal=0;
  bool badCommandLine=false;
  int option;
  while( (option=getopt(argc,argv,"c:H:")) != -1 ) {
    std::cout << "option = " << char(option) << std::endl;
    switch (option) {
    case 'c' :
      if( std::string(optarg) == std::string("yes") ) {
	compressFlag = true;
	std::cout << "Compress capsid." << std::endl;
      }
      break;
    case 'H' :
      Hfinal = std::atof(optarg);
      std::cout << "Hfinal: " << Hfinal << std::endl;
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
  //read in the material properties -- mu, ks and kappa
  std::vector<double> mu_vector(ntri,0.),ks_vector(ntri,0.),kappa_vector(ntri,0.);
  while( token != "mu" ) ifs >> token;
  ifs >> token >> token >> token ;
  for(int i=0;i<ntri;i++)
    ifs >> mu_vector[i];
  while( token != "ks" ) ifs >> token;
  ifs >> token >> token >> token ;
  for(int i=0;i<ntri;i++)
    ifs >> ks_vector[i];
  while( token != "kappa" ) ifs >> token;
  ifs >> token >> token >> token ;
  for(int i=0;i<ntri;i++)
    ifs >> kappa_vector[i];
  ifs.close();

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
  double viscosity=1.0e-2;

  // now input (if file is available) user-defined material constants
  std::ifstream inp("parameters.inp");
  inp >> Y >> nu >> KC >> KG >> C0 >> mu >> kS >> viscosity;
  
  // print parameters out to screen
  if(verbose) 
    std::cout << " Y: " << Y << std::endl
	      << " nu: " << nu << std::endl
	      << " KC: " << KC << std::endl
	      << " KG: " << KG << std::endl
	      << " C0: " << C0 << std::endl
              << " mu: " << mu << std::endl
              << " kS: " << kS << std::endl
              << " viscosity: " << viscosity << std::endl;

  // create bending body
  int quadOrder=1;
  
  typedef SCElastic BendingMaterial;
  BendingMaterial bending(KC,KG,C0);

  typedef LoopShellBody<BendingMaterial> LSB;
  LSB * bending_body = new LSB(bending, connectivities, nodes, quadOrder);
  bending_body->setOutput(paraview);
  bending_body->SetRefMeanCurvature();
  bending_body->SetBendingModulus(kappa_vector);

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
  stretch_body->SetElasticModulus(mu_vector,ks_vector);

  bending_body->compute(true,false,false);
  stretch_body->compute(true,false,false);
  std::cout << std::endl;
  std::cout << "Bending body energy = " << bending_body->energy() << std::endl;
  std::cout << "Stretching body energy = " << stretch_body->energy() << std::endl;
  std::cout << std::endl;
  bool extrapolate = false;
  if(extrapolate) std::cout<<"This uses extrapolation"<<std::endl;

  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(bending_body);
  bdc.push_back(stretch_body);
  Model model(bdc,nodes);
  
  //Input lbfgsb solver parameters
  int m=5;
  double factr=1.0e+7;
  double pgtol=1.0e-5;
  int iprint = 0;
  std::ifstream lbfgsbinp("lbfgsb.inp");
  lbfgsbinp >> iprint >> factr >> pgtol >> m ;
  if(verbose) 
    std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl;
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint );//(true);
  model.print("initial");

  if( !compressFlag ) {
    std::cout << "All done.  Bye now." << std::endl;
    return 0;
  }
  
  std::cout << "Compressing capsid." << std::endl;

  // create indentor and plate
  double k1 = 1.0e4;
  double afmR = 10.0;
  double afmX = 0.0;
  double afmY = 0.0;
  double friction_coeff = 0.0;

  inp >> k1 >> afmR >> afmX >> afmY >> friction_coeff;
  if(verbose) 
    std::cout << "Input k1: " << k1 << std::endl
	      << "Input afmR: " << afmR << std::endl
	      << "Center position of afm tip (x): " << afmX << std::endl
	      << "Center position of afm tip (y): " << afmY << std::endl
	      << "Friction Coefficient: " << friction_coeff << std::endl;

  //stiffness of indentor
  double k_AL=k1;
  // find top and bottom of capsid
  double Zmin = findmin(defNodes,1);
  double Zmax = findmax(defNodes,1);
  std::cout << "Zmin = " << Zmin << std::endl
	    << "Zmax = " << Zmax << std::endl;
  // this is written for correcting the initial flat region
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
    = new RigidPlateAL( defNodes, k_AL, Zmin, up, friction_coeff );
  bending_body->pushBack( glass );
  std::cout << "Added glass substrate to the body." << std::endl;

  //add some viscosity for regularization
  ViscousRegularizer vr(bending_body->nodes(), viscosity);
  bending_body->pushBack( &vr ); 

  int steps = 10;
  double vrTol = pgtol;
  inp >> steps;
  inp >> vrTol;
  if(verbose) 
    std::cout << "Input steps: "<<steps<<std::endl
	      << "Input vrTol: "<<vrTol<<std::endl;


  const double originalHeight=Zmax-Zmin;
  double Zbegin=Zmin;
  double Zend=Zmax-Hfinal;
  double dZ = (Zend-Zbegin)/steps;
  double vrEnergy = vr.energy();
  double bdEnergy = bending_body->energy()+stretch_body->energy();
  

  std::ofstream FvsZ;
  if(verbose) FvsZ.open("FvsZ.dat");
  int Iter = 0;

//  added for extrapolation of initial guess
  blitz::Array<double,1> x_prev(model.dof());
  blitz::Array<double,1> u_prev(model.dof());
  model.getField(solver);
  for(int i=0; i<model.dof(); i++ ) x_prev(i) = solver.field(i);
  u_prev = 0.0;

  // %%%%%%%%%%%%%%%%%%%%%%
  // Begin indentation loop
  // %%%%%%%%%%%%%%%%%%%%%%

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
	      << std::endl;

    int viterMax = 10;
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
      bdEnergy = bending_body->energy()+stretch_body->energy();

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

    double height = Zmax-Z;
    // add up forces on top and bottom
    double Ztop = Zmax;
    double Zbot = Zmax - height;       

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

  }
  FvsZ.close();

  std::cout << "Indentation Complete." << std::endl;

  return 0;
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
