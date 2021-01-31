//#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include "Contact.h"
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
typedef tvmet::Vector<double,3> CoordinateArray;
bool verbose = true;

//geometry
float Lx = 2.0, Ly = 2.0, Lz = 2.0;

//mesh free points
int nx = 3, ny=3, nz = 3;
int npts = nx*ny*nz;

//distance between mesh points
float dx = Lx/(nx-1.), dy = Ly/(ny-1.), dz = Lz/(nz-1.0);

// mesh free parameters
double support_size = 3*dx;

//array of mesh points
CoordinateArray x[npts];

for(int k=0;k<nz;k++){
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      int temp = k*nx*ny+i+j*nx;
      x[temp] = i*dx, j*dy, k*dz;
    }
  }
}
//weights/voumes associated with nodes
vector<double> node_volume;
vector<bool> node_boundary(npts,false);
for(int i=0;i<npts;i++)  node_volume.push_back(dx*dy*dz);

for(int i=0;i<npts;i++){
  if(i%nx==0 || i%nx==nx-1) {node_volume[i] = node_volume[i]/2; node_boundary[i]=true;}
  if((i/nx)%ny==0 || (i/nx)%ny==ny-1) {node_volume[i] = node_volume[i]/2; node_boundary[i]=true;}
  if((i/nx/ny)%nz==0 || (i/nx/ny)%nx==nz-1) {node_volume[i] = node_volume[i]/2; node_boundary[i]=true;}
  }

// create vector of nodes
int dof=0;
std::vector< NodeBase* > int_nodes;
std::vector< DeformationNode<3>* > int_defNodes;
std::vector< NodeBase* > nodes;
std::vector< DeformationNode<3>* > defNodes;
defNodes.reserve(npts);
for(int i=0; i<npts; i++) {
  int id=i;
  DeformationNode<3>::Point temp;
  temp(0) = x[i](0);
  temp(1) = x[i](1);
  temp(2) = x[i](2);
  
  NodeBase::DofIndexMap idx(3);
  if (node_boundary[i]==false)
    {for(int j=0; j<3; j++) idx[j]=dof++;}
  else
    {for(int j=0; j<3; j++) idx[j]=0;}
  DeformationNode<3>* n = new DeformationNode<3>(id,idx,temp);
  nodes.push_back( n );
  defNodes.push_back( n );
  if (node_boundary[i]==false){
     int_nodes.push_back( n );
     int_defNodes.push_back( n );
     }  
}
 cout << "Number of nodes: " <<nodes.size() << endl;
 cout << "Number of internal nodes: " <<int_nodes.size() << endl;

// deform the initial shape, set position of the boundary nodes with a random constant deformation gradient F
    srand(1234);
    double F[3][3];
    /*F[0][0] = 0.5*rand()/RAND_MAX+1.;
    F[0][1] = 0.5*rand()/RAND_MAX+1.;
    F[0][2] = 0.5*rand()/RAND_MAX+1.;
    F[1][0] = 0.5*rand()/RAND_MAX+1.;
    F[1][1] = 0.5*rand()/RAND_MAX+1.;
    F[1][2] = 0.5*rand()/RAND_MAX+1.;
    F[2][0] = 0.5*rand()/RAND_MAX+1.;
    F[2][1] = 0.5*rand()/RAND_MAX+1.;
    F[2][2] = 0.5*rand()/RAND_MAX+1.;*/
    /*F[0][0] = 2.1;
    F[0][1] = 1;
    F[0][2] = 1.;
    F[1][0] = 1.;
    F[1][1] = 2.1;
    F[1][2] = 1.;
    F[2][0] = 1.;
    F[2][1] = 1.;
    F[2][2] = 2.1;*/
    F[0][0] = 2.1156;
    F[0][1] = 1.1084;
    F[0][2] = 1.22378;
    F[1][0] = 1.24631;
    F[1][1] = 1.28468;
    F[1][2] = 1.23689;
    F[2][0] = 1.23746;
    F[2][1] = 1.01683;
    F[2][2] = 1.47702;
    std::cout<<"F is"<<std::endl<<F[0][0]<<"    "<<F[0][1]<<"    "<<F[0][2]<<std::endl<<F[1][0]<<"    "<<F[1][1]<<"    "<<F[1][2]<<std::endl<<F[2][0]<<"    "<<F[2][1]<<"    "<<F[2][2]<<std::endl;
    double JF = F[0][0]*(F[1][1]*F[2][2]-F[1][2]*F[2][1]) - F[0][1]*(F[1][0]*F[2][2]-F[2][0]*F[1][2]) + F[0][2]*(F[1][0]*F[2][1]-F[2][0]*F[1][1]);
    std::cout<<"And the determinant of F is "<<JF<<std::endl;
    
    for(int i=0; i<defNodes.size(); i++) {
    if(node_boundary[i]==true){
      DeformationNode<3>::Point x;
      x(0)=F[0][0]*defNodes[i]->getPosition(0) + F[0][1]*defNodes[i]->getPosition(1) + F[0][2]*defNodes[i]->getPosition(2);
      x(1)=F[1][0]*defNodes[i]->getPosition(0) + F[1][1]*defNodes[i]->getPosition(1) + F[1][2]*defNodes[i]->getPosition(2);
      x(2)=F[2][0]*defNodes[i]->getPosition(0) + F[2][1]*defNodes[i]->getPosition(1) + F[2][2]*defNodes[i]->getPosition(2);
      //x += defNodes[i]->position();
      defNodes[i]->setPoint(x);
      }
    else {
      DeformationNode<3>::Point x;
      x(0)=F[0][0]*defNodes[i]->getPosition(0) + F[0][1]*defNodes[i]->getPosition(1) + F[0][2]*defNodes[i]->getPosition(2);
      x(1)=F[1][0]*defNodes[i]->getPosition(0) + F[1][1]*defNodes[i]->getPosition(1) + F[1][2]*defNodes[i]->getPosition(2);
      x(2)=F[2][0]*defNodes[i]->getPosition(0) + F[2][1]*defNodes[i]->getPosition(1) + F[2][2]*defNodes[i]->getPosition(2);
      //x += defNodes[i]->position();
      std::cout<<"-----------------------------------"<<std::endl<<"Final position of center should be"<<std::endl<<x(0)<<"     "<<x(1)<<"     "<<x(2)<<std::endl;
      }
    }
  //create material
  double rho = 1.0;
  double E = 200.0;
  double nu = 0.4;
  double viscosity=1.0e-2;
  
  typedef CompNeoHookean MaterialType;
  MaterialType protein( rho, E, nu );
  if(verbose) std::cout << "Material has been created." << std::endl;
  
  //create body
  typedef MFBody<MaterialType,MFShape> body;
  body bd(protein,nodes,node_volume,support_size);
  bd.setOutput(paraview);
  bd.compute(true,true,false);
  
  std::cout<<"Volume of the body is "<<bd.volume()<<std::endl;
  std::cout<<"Strain energy of the body is "<<bd.totalStrainEnergy()<<std::endl;
  
  float xfor=0., yfor=0., zfor=0.;
  for (int i=0;i<npts;i++){
  if(node_boundary[i]==false){
     std::cout<<"Force at node "<<i<<" is: "<<defNodes[i]->getForce(0)<<" "<<defNodes[i]->getForce(1)<<" "<<defNodes[i]->getForce(2)<<std::endl;
     xfor += defNodes[i]->getForce(0);
     yfor += defNodes[i]->getForce(1);
     zfor += defNodes[i]->getForce(2);
     }
     }
     std::cout<<"Total forces are "<<xfor<<"  "<<yfor<<"  "<<zfor<<std::endl;
  
  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(&bd);
  Model model(bdc, int_nodes);
  std::cout<<"Checking the consistency"<<std::endl;
  //model.checkConsistency(true,false);
  //model.checkRank(model.dof()-6,true);
  //return 0;  
  model.print("initial");
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


  //add some viscosity for regularization
  ViscousRegularizer vr(bd.nodes(), viscosity);
  bd.pushBack( &vr ); 

  int steps = 100;
  double vrTol = pgtol;

  if(verbose) 
    std::cout << "Input steps: "<<steps<<std::endl
	      << "Input vrTol: "<<vrTol<<std::endl;


  double vrEnergy = vr.energy();
  double bdEnergy = bd.energy();

  solver.solve( &model );
  vrEnergy = vr.energy();
  bdEnergy = bd.energy();

  // This is where Paraview file (if desired) is saved
  model.print("deformed");
  
    for(int i=0; i<defNodes.size(); i++) {
    if(node_boundary[i]==false){
      DeformationNode<3>::Point x;
      x=defNodes[i]->point();
      std::cout<<"------------------------"<<std::endl<<"and the calculated position is "<<std::endl<<x(0)<<"     "<<x(1)<<"     "<<x(2)<<std::endl;
      }
    }
  return 0;
}
