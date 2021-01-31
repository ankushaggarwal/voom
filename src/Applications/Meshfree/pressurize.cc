#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include "ViscousRegularizer.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace std;
#include "MFShape.h"
#include "NodeBase.h"
#include "Node.h"
#include "MFBodyDoubleWell.h"
#include "MFBody.h"
#include "Model.h"
#include "CompNeoHookean.h"
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>
#include "Solver.h"
#include "Lbfgsb.h"
#include "RadialForce.h"
using namespace voom;

void readconfig_dat(char* filename, float* x, float* y, float* z, int npts);

int main(int argc, char* argv[])
{
  if( argc < 2 ) {
      cout << "Usage: pressurize InputFileName [-R final_radius -s support_size - k kernel_option]." << endl;
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

  double Rfinal=0;
  bool badCommandLine=false;
  double support_size;
  bool cubic_kernel=false, modified_kernel=false, singular_kernel=false;

  for(char option; (option=getopt(argc,argv,"c:H:b:R:r:s:k:")) != EOF; ) 
    switch (option) {
    case 'R' :
      Rfinal = std::atof(optarg);
      std::cout << "Rfinal: " << Rfinal << std::endl;
      break;
    default :
      badCommandLine=true;
      break;
    case 's' :
      support_size = std::atof(optarg);
      std::cout << "Global support size: " << support_size << std::endl;
      break;
    case 'k' :
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
      }
      break;
    }

  if( badCommandLine || argc < 2 ) {
      cout << "Usage: capsid InputFileName [-c compress/no -H finalHeight -b boundaryConditionFlag -R restartFilename -r Eratio -s SupportSize -k kernel_option]." << endl;
      return(0);
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
 
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
    node_volume.push_back(weight);
    node_radius.push_back(pow(weight*3/4./3.14,1./3.));
    supp_size.push_back(support_size);
    supp_size_hat.push_back(supphat*supphatfactor);
  }
  cout << "Number of nodes: " <<nodes.size() << endl;

  // create material
  double rho = 1.0;
  double E = 2.0;
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
  typedef MFBodyDoubleWell<MaterialType,MFShape> body;
  //typedef MFBody<MaterialType,MFShape> body;
  body bd;

  // vector of deformation gradients
  std::vector<Tensor3D> Fhat;
  Fhat.reserve(npts);
  ifstream Fhatfile;
  Fhatfile.open("Fhatinv.dat");
  if(!Fhatfile) std::cout << "Error: Fhatinv.dat could not be opened"<<std::endl;
  for(int i=0;i<npts;i++){
    double f11,f12,f13,f21,f22,f23,f31,f32,f33;
    Fhatfile >> f11 >> f12 >> f13 >> f21 >> f22 >> f23 >> f31 >> f32 >> f33;
    Fhat[i]=f11,f12,f13,f21,f22,f23,f31,f32,f33;}
  //    protein.setDeformationGradient(Fhat[i]);
  //    protein.updateState(true,false,false);
  //    const double EDnorm=protein.energyDensity();
  //    if(determinant(Fhat[i])<0.) std::cout<<"very low EDnorm = "<<EDnorm<<std::endl
    //                           <<" i = " <<i <<std::endl
    //                           <<" Fhatinv = "<< Fhat[i] << std::endl
    //                           <<" determinant = "<<determinant(Fhat[i])<<std::endl;}
 //   Fhat[i]=1.1,0.,0.,0.,1.1,0.,0.,0.,1.1;}

  if(cubic_kernel)
  bd.initializeBody(protein,nodes,node_volume,supp_size,Fhat);
  if(modified_kernel)
  bd.initializeBody(protein,nodes,node_volume,supp_size,supp_size_hat,Fhat);
  if(singular_kernel)
  bd.initializeBody(protein,nodes,node_volume,supp_size,1.5,Fhat);

/*  if(cubic_kernel)
  bd.initializeBody(protein,nodes,node_volume,supp_size);
  if(modified_kernel)
  bd.initializeBody(protein,nodes,node_volume,supp_size,supp_size_hat);
  if(singular_kernel)
  bd.initializeBody(protein,nodes,node_volume,supp_size,1.5);*/
  std::cout << "Body has been created and the shapefunctions have been calculated"<<std::endl;
  bd.setOutput(paraview);
  bd.compute(true,true,false);
  
  ifs.close();

  std::cout << std::endl;
  std::cout << "Body energy = " << bd.energy() << std::endl;
  std::cout << "Body volume = " << bd.volume() << std::endl;
  std::cout << std::endl;
  //return 0;
  //displace the nodes
  //for(int i=0;i<npts;i++){
  //  DeformationNode<3>::Point x;
  //  x(0)=1.1*defNodes[i]->getPosition(0);
  //  x(1)=1.1*defNodes[i]->getPosition(1);
  //  x(2)=1.1*defNodes[i]->getPosition(2);
  //  defNodes[i]->setPoint(x);}

  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(&bd);
  Model model(bdc, nodes);

  //model.checkConsistency(true,false);
  //model.checkRank(model.dof()-6,true);
  //return 0;

  //Input lbfgsb solver parameters
  int m=5;
  double factr=1.0e+1;
  double pgtol=1.0e-5;
  int iprint = 0;
  ifstream lbfgsbinp("lbfgsb.inp");
  lbfgsbinp >> iprint >> factr >> pgtol >> m ;
  if(verbose) 
    std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl;
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint );


////////////////////////////////////////////////////// starting from EI-II

  float x[npts], y[npts], z[npts];
  //std::cout<<"using .dat files"<<std::endl;
  char inpfile[20];
  sprintf(inpfile,"%d.dat",0);
  readconfig_dat(inpfile,x,y,z,npts);
  Rfinal=0.;
  for(int inode=0;inode<npts;inode++){
    DeformationNode<3>::Point newx;
    newx(0)=x[inode]; newx(1)=y[inode]; newx(2)=z[inode];
    Rfinal += sqrt(x[inode]*x[inode] + y[inode]*y[inode] + z[inode]*z[inode]);
    //defNodes[inode]->setPoint(newx);
  }
  Rfinal /= npts;
/////////////////////////////////////////////////////////////  
  //add the Pressurize element
  double k1 = 5.0e5;
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
  RadialForce pressure(bd.nodes(),k1);
  bd.addElement(&pressure);
  //bd.compute(true,false,false);
  //std::cout<<"Mean Radius = " << pressure.MeanRadius() << " and Energy = "<<bd.energy()-pressure.energy()<<std::endl;
  //return 0;

  int steps = 1;
  inp >> steps;
  double Rbegin=pressure.MeanRadius();
  double dRho=1./steps;
  double dR=(Rfinal-Rbegin)/steps;
  if(verbose) 
    std::cout << "Input steps: "<<steps   <<std::endl
              << "Rbegin: "     << Rbegin <<std::endl
              << "dR: "         << dR     <<std::endl;
  
  double bdEnergy = bd.energy();
  ofstream EvsR;
  if(verbose) EvsR.open("EvsR.dat");
  // %%%%%%%%%%%%%%%%%%%%%%
  // Begin pressurization loop
  // %%%%%%%%%%%%%%%%%%%%%%
  int Iter = 0;
  double R = Rbegin;
  for(int istep=0; istep<=steps; istep++, R+=dR) {
    pressure.setR(R);
    //set the starting guess as the interpolation between the two states
    for(int inode=0;inode<npts;inode++){
      DeformationNode<3>::Point newx;
      double Rho=istep*dRho;
      newx(0)=(1-Rho)*defNodes[inode]->getPosition(0) + Rho*x[inode];
      newx(1)=(1-Rho)*defNodes[inode]->getPosition(1) + Rho*y[inode];
      newx(2)=(1-Rho)*defNodes[inode]->getPosition(2) + Rho*z[inode];
      defNodes[inode]->setPoint(newx);
    }
    //solver.solve( &model );
    bd.compute(true,true,false);
    double bdEnergyi = bd.energy()-pressure.energy();
    //solver.solve( &model );
    bdEnergy=bd.energy();
    // Calculate the RMSD from this structure to the other structure
    double rmsd=0.;
    for(int inode=0;inode<npts;inode++){
      double dx=defNodes[inode]->getPoint(0)-x[inode];
      double dy=defNodes[inode]->getPoint(1)-y[inode];
      double dz=defNodes[inode]->getPoint(2)-z[inode];
      //DeformationNode<3>::Point newx=defNodes[inode]->point();
//      rmsd += (x[inode]-newx(0))*(x[inode]-newx(0)) + (y[inode]-newx(1))*(y[inode]-newx(1)) + (z[inode]-newx(2))*(z[inode]-newx(2));    
      rmsd += dx*dx + dy*dy + dz*dz;
    }
    rmsd /= npts;
    rmsd = sqrt(rmsd);

    if(verbose) {
	std::cout << "Iteration # " << Iter << std::endl
                  << "R = " << R << std::endl
		  << "   body energy = " << bdEnergy - pressure.energy()<< std::endl
                  << " rmsd = " << rmsd << std::endl
                  << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
      }
    // This is where Paraview file (if desired) is saved
    char name[20]; 
    //if(Iter%4==0){
      sprintf(name,"disp%d",Iter);
      model.print(name);
    //}
    Iter++;
    
   
    // Output the resulting energy and forces
    if(verbose) {
      EvsR << std::setw( 24 ) << std::setprecision(16) 
           << R
	   << std::setw( 24 ) << std::setprecision(16)
	   << pressure.MeanRadius()
           << std::setw( 24 ) << std::setprecision(16)
           << bdEnergyi
           << std::setw( 24 ) << std::setprecision(16)
           << bdEnergy
           << std::setw( 24 ) << std::setprecision(16)
           << bdEnergy-pressure.energy()
           << std::setw( 24 ) << std::setprecision(16)
           << rmsd
	   << std::endl;
	}
  }
  std::cout << "Pressurization Complete." << std::endl;
  return 0;
}

void readconfig_dat(char* filename, float* x, float* y, float* z, int npts){
  ifstream inpfile;
  inpfile.open(filename);
  int id;
  
  for(int i=0; i<npts; i++){
    inpfile >> id >> x[i] >> y[i] >> z[i];
  }
  inpfile.close();
}
