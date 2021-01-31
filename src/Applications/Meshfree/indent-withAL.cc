//#include <string>
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

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace voom;
using namespace std;

double findmin(std::vector< DeformationNode<3>* > &defNodes, int dim, double rad=0.);
double findmax(std::vector< DeformationNode<3>* > &defNodes, int dim, double rad=0.);


int main(int argc, char* argv[])
{
  if( argc < 2 ) {
      cout << "Usage: indent_withAL InputFileName [-c compress/no -H final_height -s support_size - k kernel_option]." << endl;
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
  string inputFileName = modelName;

  bool compressFlag=false;
  double Hfinal=0;
  bool badCommandLine=false;
  double support_size;
  bool cubic_kernel=false, modified_kernel=false, singular_kernel=false;

  for(char option; (option=getopt(argc,argv,"c:H:b:R:r:s:k:")) != EOF; ) 
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
      cout << "Usage: indent_withAL InputFileName [-c compress/no -H finalHeight -s SupportSize -k kernel_option]." << endl;
      return(0);
  }
  if( compressFlag == true && Hfinal <= 0.0 ) {
    std::cout << "Please specify a positive value of Hfinal." << std::endl;
    return 0;
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
  std::vector<double> node_volume, supp_size, supp_size_hat, node_radius;
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
 
   //if(change of units from angstrom to nanometer)
   //double fac=10.;
   //x /= fac;
   //std::cout<<"Input dimensions changed by a factor of "<<fac<<std::endl;

    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
    //node_volume.push_back(weight/fac/fac/fac);
    node_volume.push_back(weight);
    node_radius.push_back(pow(weight*3/4./3.14,1./3.));
    supp_size.push_back(support_size);
    supp_size_hat.push_back(supphat*supphatfactor);
  }
  cout << "Number of nodes: " <<nodes.size() << endl;
  //if(!SNNIcontact)
  node_radius.clear(); node_radius.resize(npts,0.);

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

  //typedef CompNeoHookean MaterialType;
  //MaterialType protein( rho, E, nu );
  typedef MooneyRivlin MaterialType;
  MaterialType protein( E, E, 8./3.*E );
  if(verbose) std::cout << "Material has been created." << std::endl;

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
  
  bool extrapolate = false;
  if(extrapolate) std::cout<<"This uses extrapolation"<<std::endl;
  
  ifs.close();

  std::cout << std::endl;
  std::cout << "Body energy = " << bd.energy() << std::endl;
  std::cout << "Body volume = " << bd.volume() << std::endl;
  std::cout << std::endl;

  //set the node wise modulus
  std::vector<double> C10(npts,0.),C01(npts,0.),D1(npts,0.);
  ifs.open("first.mat");
  for(int i=0;i<npts;i++){
    ifs >> C10[i];
    C10[i] *= 41.4;
    if(C10[i]<0) C10[i]=0.;}
  ifs.close(); 
  ifs.open("second.mat");
  for(int i=0;i<npts;i++){
    ifs >> C01[i]; 
    C01[i] *= 41.4;
    if(C01[i]<0) C01[i]=0.;}
  ifs.close(); 
  ifs.open("bulk.mat");
  for(int i=0;i<npts;i++){
    ifs >> D1[i]; 
    D1[i] *= 41.4;
    if(D1[i]<0) D1[i]=0.;}
  ifs.close(); 
  bd.setModulus(C10,C01,D1);

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
  
  bd.addElement( afm );
  std::cout << "Added afm to body." << std::endl;

  bool up=true;
  RigidPlateAL * glass
    = new RigidPlateAL( defNodes, k_AL, Zmin, up, friction_coeff );
  bd.addElement( glass );
  std::cout << "Added glass substrate to the body." << std::endl;

  //add some viscosity for regularization
  ViscousRegularizer vr(bd.nodes(), viscosity);
  bd.addElement( &vr ); 

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
  double bdEnergy = bd.energy();
  

  ofstream FvsZ;
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
      bdEnergy = bd.energy();

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
