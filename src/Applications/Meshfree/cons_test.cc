#include <iostream>
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
using namespace voom;

int main(int argc, char* argv[])
{
typedef tvmet::Vector<double,3> CoordinateArray;

//geometry
float Lx = 1.0, Ly = 1.0, Lz = 1.0;

//mesh free points
int nx=10, ny=10, nz=10;
int npts = nx*ny*nz;

//distance between mesh points
float dx = Lx/(nx-1.0), dy = Ly/(ny-1.0), dz = Lz/(nz-1.0);

// mesh free parameters
double global_supp_size = 3*dx;

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
vector<double> support_size(npts,global_supp_size);
for(int i=0;i<npts;i++)  node_volume.push_back(dx*dy*dz);

for(int i=0;i<npts;i++){
  if(i%nx==0 || i%nx==nx-1) {node_volume[i] = node_volume[i]/2; node_boundary[i]=true;}
  if((i/nx)%ny==0 || (i/nx)%ny==ny-1) {node_volume[i] = node_volume[i]/2; node_boundary[i]=true;}
  if((i/nx/ny)%nz==0 || (i/nx/ny)%nx==nz-1) {node_volume[i] = node_volume[i]/2; node_boundary[i]=true;}
  }

// create vector of nodes
int dof=0;
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
  for(int j=0; j<3; j++) idx[j]=dof++;
  DeformationNode<3>* n = new DeformationNode<3>(id,idx,temp);
  nodes.push_back( n );
  defNodes.push_back( n );
}
cout << "Number of nodes: " <<nodes.size() << endl;


//create material
double rho = 1.0;
double E = 200.0;
double nu = 0.4;
double viscosity=1.0e-2;

typedef CompNeoHookean MaterialType;
MaterialType protein( rho, E, nu );

//create body
typedef MFBody<MaterialType,MFShape> body;
body bd(protein,nodes,node_volume,support_size);
//compute the energy and forces
bd.compute(1,1,0);
  
std::cout<<"Volume of the body is "<<bd.volume()<<std::endl;
std::cout<<"Strain energy of the body is "<<bd.totalStrainEnergy()<<std::endl;
  
// create Model
Model::BodyContainer bdc;
bdc.push_back(&bd);
Model model(bdc, nodes);

//model.checkConsistency(true,false);
//model.checkRank(model.dof()-6,true);
model.normalModes_all(false,7,node_volume,rho,"try");
return 0;  
/*  
  DeformationNode<3>::PositionVector ll=defNodes[10]->position();
  std::cout<<ll(0)<<"\t"<<ll(1)<<"\t"<<ll(2)<<std::endl;
  std::cout<<"---------------------------------------------------------------"<<std::endl;
  std::cout<<x[10]<<std::endl;
  std::cout<<"---------------------------------------------------------------"<<std::endl;
//std::cout<<x<<std::endl;
MFShape::CoordinateArray p(defNodes[5000]->getPosition(0),defNodes[5000]->getPosition(1),defNodes[5000]->getPosition(2));
ss.compute(p,defNodes,80.0,10.0);
//ss.compute(x[5000],x,80,10000,10);

CoordinateArray let(0), letx(0), lety(0), letz(0);
voom::MFShape::FunctionContainer value=ss.functions(), xd=ss.xderivative(), yd=ss.yderivative(), zd=ss.zderivative();
voom::MFShape::NodeNContainer node=ss.neighb();
std::cout<<"the size of neighbour nodes is "<<value.size()<<std::endl;
double total=0,totalx=0,totaly=0,totalz=0;
for(int i=0;i<value.size();i++){
std::cout<<node[i]<<"\t"<<value[i]<<"\t"<<xd[i]<<"\t"<<yd[i]<<"\t"<<zd[i]<<"\n";
//x[node[i]]
total = total + value[i];
let += x[node[i]]*value[i];
totalx += xd[i];
letx += x[node[i]]*xd[i];
totaly += yd[i];
lety += x[node[i]]*yd[i];
totalz += zd[i];
letz += x[node[i]]*zd[i];
}
std::cout<<std::endl;
std::cout<<total<<"\t"<<totalx<<"\t"<<totaly<<"\t"<<totalz<<"\t"<<std::endl;
std::cout<<let<<"\t"<<x[5000]<<endl<<"\t"<<letx<<"\t"<<lety<<"\t"<<letz<<endl;*/
//std::cout<<ss.functions<<std::endl;
//std::cout<<ss.nodes<<std::endl;
/*tvmet::Matrix<double,4,4> M(0), invM(0), I(0);
M = 1,6,8,1,
5,2,6,4,
8,3,1,2,
7,9,2,5; 
ss.inverse(M,invM);
I = invM*M;
std::cout<<I<<std::endl;

tvmet::Matrix<double,3,1> a;
a=1,2,3;
tvmet::Matrix<double,3,3> b(0);
b=a*tvmet::trans(a);
std::cout<<b<<std::endl;*/
}
