//#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include "ContactSNNI.h"
#include "ViscousRegularizer.h"

#if defined(_OPENMP)
#include <omp.h>
#endif


// using namespace tvmet;
 using namespace std;
#include "MFShape.h"
//#include "MFShape.cc"
#include "NodeBase.h"
#include "Node.h"
#include "MFBody.h"
#include "Model.h"
#include "CompNeoHookean.h"
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>
#include "Solver.h"
#include "Lbfgsb.h"
 using namespace voom;

int main(int argc, char* argv[])
{
  if( argc < 2 ) {
      cout << "Usage: indent modelName [-c compress/no -H final_height -b boundary_conditions -r Eratio -s support_size]." << endl;
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
  double Eratio = 1;
  bool restart=false;
  string restartFileName;
  bool boundaryConditionFlag = false;
  bool badCommandLine=false;
  double support_size;

  for(char option; (option=getopt(argc,argv,"c:H:b:R:r:s:")) != EOF; ) 
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
    case 'b' :
      if( std::string(optarg) == std::string("yes") ) {
	boundaryConditionFlag = true;
      }
      else if( std::string(optarg) != std::string("no") ) {
	std::cout << "Yes or no for boundary conditions.  Exiting."
		  << std::endl;
	return 0;
      }	
      std::cout << "boundaryConditionFlag: " << boundaryConditionFlag
		<< std::endl;
      break;
    case 'R' :
      restart=true;
      restartFileName = optarg;
      std::cout << "Restarting from file: " << restartFileName
		<< std::endl;
      break;
    default :
      badCommandLine=true;
      break;
    case 'r' :
      Eratio = std::atof(optarg);
      std::cout << "Ratio of E: " << Eratio << std::endl;
      break;
    case 's' :
      support_size = std::atof(optarg);
      std::cout << "Global support size: " << support_size << std::endl;
      break;
    }

  if( badCommandLine || argc < 2 ) {
      cout << "Usage: capsid InputFileName [-c -H finalHeight -b boundaryConditionFlag -R restartFilename -r Eratio -s SupportSize]." << endl;
      return(0);
  }
  if( Hfinal <= 0.0 ) {
    std::cout << "Please specify a positive value of Hfinal." << std::endl;
    return 0;
  }

  
  // create input stream
  ifs.open( inputFileName.c_str(), ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << inputFileName << endl;
    exit(0);
  }

  //
  // create vector of nodes
  int dof=0;

  // mesh free parameters
//  double support_size = 1.4;
//  std::cout<<"Global support size of "<<support_size<<" is used."<<std::endl;

  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;
  std::vector<double> node_volume, supp_size, node_radius;
  char key;  ifs>>key;
  
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
    DeformationNode<3>::Point x;
    ifs >> temp >> x(0) >> x(1) >> x(2) >> weight;
    //x /= 10.0;
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
    //node_volume.push_back(weight/1000);
    node_volume.push_back(weight);
    node_radius.push_back(pow(weight*3/4./3.14,1./3.));
    supp_size.push_back(support_size);
  }
  cout << "Number of nodes: " <<nodes.size() << endl;


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
  typedef MFBody<MaterialType,MFShape> body;
  body bd(protein,nodes,node_volume,supp_size);
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

  if( !compressFlag ) {
    std::cout << "All done.  Bye now." << std::endl;
    return 0;
  }
  
  std::cout << "Compressing capsid." << std::endl;
  // find top and bottom of capsid
  double Zmin=std::numeric_limits<double>::max();
  double Zmax=-std::numeric_limits<double>::max();
  double Zavg = 0.0;
  for (int a=0; a<defNodes.size(); a++){
    double Z = defNodes[a]->getPoint(2);
    Zmin = std::min(Zmin, Z);
    Zmax = std::max(Zmax, Z);
    Zavg += Z;
  }
  Zavg /= defNodes.size();

  std::cout << "Zmin = " << Zmin << std::endl
	    << "Zmax = " << Zmax << std::endl
	    << "Zavg = " << Zavg << std::endl;
    

  // create indentor and plate
  double k1 = 1.0e4;
  double afmR = 10.0;
  double afmX = 0.0;
  double afmY = 0.0;
  bool frictionalContact = false;

  inp >> k1 >> afmR >> afmX >> afmY >> frictionalContact;
  if(verbose) 
    std::cout << "Input k1: " << k1 << std::endl
	      << "Input afmR: " << afmR << std::endl
	      << "Center position of afm tip (x): " << afmX << std::endl
	      << "Center position of afm tip (y): " << afmY << std::endl
	      << "Frictional Contact: " << frictionalContact << std::endl;

  // this is written for correcting the initial flat region
  double Zmin1=std::numeric_limits<double>::max();
  double Zmax1=-std::numeric_limits<double>::max();  
  double Zavg1 = 0.0;

  for (int a=0; a<defNodes.size(); a++){
    double Z = defNodes[a]->getPoint(2);
    double r=sqrt(defNodes[a]->getPoint(0)*defNodes[a]->getPoint(0) + defNodes[a]->getPoint(1)*defNodes[a]->getPoint(1));
    Zmin1 = std::min(Zmin1, Z);
    Zmax1 = std::max(Zmax1, sqrt(afmR*afmR-r*r)+Z);
    Zavg1 += Z; 
  }
  Zavg1 /= defNodes.size();

  std::cout << "Zmin1 = " << Zmin1 << std::endl
	    << "Zmax1 = " << Zmax1 << std::endl;

	    tvmet::Vector<double,3> xc; xc = afmX, afmY, Zmax1; 
  //tvmet::Vector<double,3> xc; xc = afmX, afmY, Zmax+afmR; 


  double friction=20.0;
//  if(afmR>0) {
  RigidHemisphereContactSNNI * afm 
    = new RigidHemisphereContactSNNI( defNodes,node_radius, afmR, xc, friction, pgtol );
// }
//  else {
//  RigidPlateContactSNNI * afm
//    = new RigidPlateContactSNNI( defNodes,node_radius, Zmax, false, friction, pgtol);
//  }

  bd.pushBackContact( afm );

  bool up=true;
  RigidPlateContactSNNI * glass
    = new RigidPlateContactSNNI( defNodes,node_radius, Zmin, up, friction, pgtol );
  bd.pushBackContact( glass );

  //add some viscosity for regularization
  ViscousRegularizer vr(bd.nodes(), viscosity);
  bd.pushBack( &vr ); 

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
  if( restart ) inp >> dZ;
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
		  << "  bottom force = " << glass->FZ() << std::endl;
      }
	
      // step forward in "time", relaxing viscous energy & forces 
      vr.step();

      if(vrEnergy < vrTol*bdEnergy && solver.projectedGradientNorm()<=pgtol ) { 
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

