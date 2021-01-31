//#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace std;
#include "MFShape.h"
#include "NodeBase.h"
#include "Node.h"
#include "MFBody.h"
#include "Model.h"
#include "CompNeoHookean.h"
#include "MooneyRivlin.h"
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>
using namespace voom;

void SHdisp(int l,int m,double x,double y,double z,std::vector<double>& dx);
int main(int argc, char* argv[])
{
  if( argc < 2 ) {
      cout << "Usage: SHenergy InputFileName [-s support_size - k kernel_option]." << endl;
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

  bool badCommandLine=false;
  double support_size;
  bool cubic_kernel=false, modified_kernel=false, singular_kernel=false;

  for(char option; (option=getopt(argc,argv,"s:k:")) != EOF; ) 
    switch (option) {
    case 's' :
      support_size = std::atof(optarg);
      std::cout << "Global support size: " << support_size << std::endl;
      break;
    default:
      badCommandLine=true;
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
      cout << "Usage: SHenergy InputFileName [-s support_size - k kernel_option]." << endl;
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
  double Ravg=0.;
  defNodes.reserve(npts);

  // read in points
  for(int i=0; i<npts; i++) {
    int id=i;
    int temp;
    double weight;
    double supphat;
    DeformationNode<3>::Point x;
    ifs >> temp >> x(0) >> x(1) >> x(2) >> weight>>supphat;
    Ravg += tvmet::norm2(x);
 
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
  assert(nodes.size()!=0);
  Ravg /= nodes.size();
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

  //typedef CompNeoHookean MaterialType;
  //MaterialType protein( rho, E, nu );
  typedef MooneyRivlin MaterialType;
  MaterialType protein( E, 0, 8./3.*E );
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
  
  ifs.close();

  std::cout << std::endl;
  std::cout << "Body energy = " << bd.energy() << std::endl;
  std::cout << "Body volume = " << bd.volume() << std::endl;
  std::cout << std::endl;

  //set the node wise modulus
  std::vector<double> C10(npts,0.),C01(npts,0.),D1(npts,0.);
  ifs.open("first.mat");
  for(int i=0;i<npts;i++) ifs >> C10[i];
  ifs.close(); 
  ifs.open("second.mat");
  for(int i=0;i<npts;i++) ifs >> C01[i]; 
  ifs.close(); 
  ifs.open("bulk.mat");
  for(int i=0;i<npts;i++) ifs >> D1[i]; 
  ifs.close(); 
  bd.setModulus(C10,C01,D1);

  //apply deformation of spherical harmonics and calculate the energy
  double mag=1e-1;
  for(int l=0;l<5;l++){
    double lenergy=0.;
    for(int m=0;m<2*l+1;m++){
      for(int inode=0;inode<npts;inode++){
        tvmet::Vector<double,3> newx(defNodes[inode]->getPosition(0),defNodes[inode]->getPosition(1),defNodes[inode]->getPosition(2));
        std::vector<double> dx(3,0.);
        SHdisp(l,m,newx(0),newx(1),newx(2),dx);
        newx(0) += mag*dx[0]/sqrt(2.*l+1.); newx(1) += mag*dx[1]/sqrt(2.*l+1.); newx(2) += mag*dx[2]/sqrt(2.*l+1.);
        defNodes[inode]->setPoint(newx);
      }
      bd.compute(true,false,false);
      double energy=bd.energy();
      lenergy += energy;
      std::cout<<"l="<<l<<", m="<<m<<": Energy="<<energy<<std::endl;
      //print the paraview files
      //char name[20];
      //sprintf(name,"l%dm%d",l,m);
      //bd.printParaview(name);
    }
    //analytical energy
    //double ana_energy=1./4.*((8.*kS)+KC/Ravg/Ravg*l*(l-1.)*(l+1.)*(l+2.))*mag*mag;
    std::cout<<std::endl<<"l="<<l<<": Energy="<<lenergy<<std::endl;
    //                    <<"Analytical energy="<<ana_energy<<std::endl;
  }

  std::cout<<"Everything done!"<<std::endl;
  return 0;
}

void SHdisp(int l,int m,double x,double y,double z,std::vector<double>& dx){
double ur;
double r=sqrt(x*x+y*y+z*z);
if(l==0)
  ur=1./2./sqrt(M_PI);

if(l==1){
  if(m==0)
    ur=sqrt(3./4./M_PI)*x/r;
  if(m==1)
    ur=sqrt(3./4./M_PI)*y/r;
  if(m==2)
    ur=sqrt(3./4./M_PI)*z/r;
}

if(l==2){
  if(m==0)
    ur=1./4.*sqrt(5./M_PI)*(-x*x-y*y+2*z*z)/r/r;
  if(m==1)
    ur=1./2.*sqrt(15./M_PI)*y*z/r/r;
  if(m==2)
    ur=1./2.*sqrt(15./M_PI)*z*x/r/r;
  if(m==3)
    ur=1./2.*sqrt(15./M_PI)*x*y/r/r;
  if(m==4)
    ur=1./4.*sqrt(15./M_PI)*(x*x-y*y)/r/r;
}

if(l==3){
  if(m==0)
    ur=1./4.*sqrt(7./M_PI)*z*(2*z*z-3*x*x-3*y*y)/r/r/r;
  if(m==1)
    ur=1./4.*sqrt(35./2./M_PI)*y*(3*x*x-y*y)/r/r/r;
  if(m==2)
    ur=1./4.*sqrt(35./2./M_PI)*x*(x*x-3*y*y)/r/r/r;
  if(m==3)
    ur=1./4.*sqrt(105./M_PI)*z*(x*x-y*y)/r/r/r;
  if(m==4)
    ur=1./2.*sqrt(105./M_PI)*x*y*z/r/r/r;
  if(m==5)
    ur=1./4.*sqrt(21./2./M_PI)*y*(4*z*z-x*x-y*y)/r/r/r;
  if(m==6)
    ur=1./4.*sqrt(21./2./M_PI)*x*(4*z*z-x*x-y*y)/r/r/r;
}

if(l==4){
  if(m==0)
    ur=3./16.*sqrt(1./M_PI)*(35.*z*z*z*z - 30.*z*z*r*r + 3.*r*r*r*r)/(r*r*r*r);
  if(m==1)
    ur=3./4.*sqrt(5./2./M_PI)*x*z*(7*z*z-3*r*r)/(r*r*r*r);
  if(m==2)
    ur=3./4.*sqrt(5./2./M_PI)*y*z*(7*z*z-3*r*r)/(r*r*r*r);
  if(m==3)
    ur=3./8.*sqrt(5./M_PI)*(x*x-y*y)*(7*z*z-r*r)/(r*r*r*r);
  if(m==4)
    ur=3./4.*sqrt(5./M_PI)*x*y*(7*z*z-r*r)/(r*r*r*r);
  if(m==5)
    ur=3./4.*sqrt(35./2./M_PI)*x*z*(x*x-3*y*y)/(r*r*r*r);
  if(m==6)
    ur=3./4.*sqrt(35./2./M_PI)*x*z*(x*x-3*y*y)/(r*r*r*r);
  if(m==7)
    ur=3./4.*sqrt(35./2./M_PI)*y*z*(3*x*x-y*y)/(r*r*r*r);
  if(m==8)
    ur=3./16.*sqrt(35./M_PI)*(x*x*x*x-6*x*x*y*y+y*y*y*y)/(r*r*r*r);
  if(m==9)
    ur=3./4.*sqrt(35./M_PI)*x*y*(x*x-y*y)/(r*r*r*r);
}

if(l>4) std::cout<<"Not coded yet"<<std::endl;


dx[0]=x/r*ur;
dx[1]=y/r*ur;
dx[2]=z/r*ur;
}
