//#include <string>
// This file reads in the configurations from MD
// either in .dat or .pdf formats
// then calculates the deformation gradient at each node
// then calculates the three invariants at each node
// and writes them to the file invariants1.dat *2.dat and *3.dat
// in the format : at each configuration a row of invariants at all the nodes
// 
//
//------>>> DIFFERENCE FROM homog.cc
//          This uses a table to create averages read from table.dat
//
//

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
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>
using namespace voom;

void readconfig_dat(char* filename, std::vector< std::vector<int> > & tableid, float* x, float* y, float* z, int npts);
void readconfig_pdb(char* filename, float* x, float* y, float* z, int npts);
void readconfig_dat(char* filename, float* x, float* y, float* z, int npts, int average);

int main(int argc, char* argv[])
{
  if( argc < 2 ) {
    cout << "Usage: homog-table inputFileName [-n start_config -N end_config -s support_size - k kernel_option]." << endl;
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
  bool datfile = true;
  bool pdbfile = false;
  bool badCommandLine=false;
  double support_size;
  int StartConfig=0, EndConfig=0;
  bool cubic_kernel=false, modified_kernel=false, singular_kernel=false;
  
  for(char option; (option=getopt(argc,argv,"n:N:s:k:")) != EOF; ) 
    switch (option) {
    default :
      badCommandLine=true;
      break;
    case 'n' :
      StartConfig = std::atoi(optarg);
      std::cout << "Starting No. of configurations: " << StartConfig << std::endl;
      break;
    case 'N' :
      EndConfig = std::atoi(optarg);
      std::cout << "Ending No. of configurations: " << EndConfig << std::endl;
      break;
    case 's' :
      support_size = std::atof(optarg);
      std::cout << "Global support size: " << support_size <<std::endl;
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
    cout << "Usage: homog-table inputFileName [-n start_configs -N end_config -s support_size - k kernel_option]." << endl;
    return(0);
  }
  
  // create input stream
  ifs.open( inputFileName.c_str(), ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << inputFileName << endl;
    exit(0);
  }

  std::ifstream ftable;
  ftable.open("table.dat",std::ios::in);
  if (!ftable) {
    cout << "Cannot open input file: table.dat " << endl;
    exit(0);
  }
  // create vector of nodes
  int dof=0;
  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;
  vector<double> node_volume,supp_size, supp_size_hat;
  double supphatfactor=0.9;
  if(modified_kernel) std::cout<<"A factor of "<<supphatfactor<<" is used for scaling the support size for primitive kernel"<<std::endl;
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
    double supphat;
    DeformationNode<3>::Point x;
    ifs >> temp >> x(0) >> x(1) >> x(2) >> weight>>supphat;
    
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
    node_volume.push_back(weight);
    supp_size.push_back(support_size);
    supp_size_hat.push_back(supphat*supphatfactor);
  }
  cout << "Number of nodes: " <<nodes.size() << endl;
  ifs.close();
  // create material, values don't matter here
  double rho = 1.0;
  double E = 200.0;
  double nu = 0.4;
  double viscosity=1.0e-2;
  typedef CompNeoHookean MaterialType;
  MaterialType protein( rho, E, nu );
  
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
  
  
  std::cout << std::endl;
  std::cout << "Body energy = " << bd.energy() << std::endl;
  std::cout << "Body volume = " << bd.volume() << std::endl;
  std::cout << std::endl;
  
  //read the table file
  std::vector< std::vector<int> > tableid;
  int N;
  ftable >> N;
  do{
    std::vector<int> id;
    for(int i=0;i<N;i++){
      int l;
      ftable >> l;
      id.push_back(l);
    }
    tableid.push_back(id);
    ftable >> N;
  }while(!ftable.eof());
  if(tableid.size()!=npts) std::cout<<"Error"<<std::endl;
  ftable.close();

  //open the output files
  ofstream file1, file2, file3;
  file1.open("inv1.dat",ios::binary);
  file2.open("inv2.dat",ios::binary);
  file3.open("inv3.dat",ios::binary);
  char tab='\t', newline='\n';
  int col=npts+1, row=EndConfig-StartConfig+1;
  file1.write((char*)&row,sizeof(int)); file1.write((char*)&col,sizeof(int));
  file2.write((char*)&row,sizeof(int)); file2.write((char*)&col,sizeof(int));
  file3.write((char*)&row,sizeof(int)); file3.write((char*)&col,sizeof(int));
  //read in the configuration files in a loop
  for(int i=StartConfig;i<=EndConfig;i++){
    char inpfile[20];
    std::cout<<"Reading "<<i<<"th configuration."<<std::endl;
    float x[npts], y[npts], z[npts];
    if(datfile){
      //std::cout<<"using .dat files"<<std::endl;
      sprintf(inpfile,"../../com-traj/%d.dat",i);
      readconfig_dat(inpfile,tableid,x,y,z,npts);}
    
    if(pdbfile){
      //std::cout<<"using .pdb files"<<std::endl;
      sprintf(inpfile,"../calpha-traj/%d.pdb",i);
      readconfig_pdb(inpfile,x,y,z,npts);} 
    std::cout<<"File has been read"<<std::endl;
    for(int inode=0;inode<npts;inode++){
      DeformationNode<3>::Point newx;
      newx(0)=x[inode]; newx(1)=y[inode]; newx(2)=z[inode];
      defNodes[inode]->setPoint(newx);
    }
    std::vector<double> I1(npts,0.), I2(npts,0.), I3(npts,0.);
    std::cout<<"Calculating the invariants"<<std::endl;
    bd.cal_invariants(I1,I2,I3);
    //bd.cal_stretch_ratios(I1,I2,I3);
    double di=double(i);
    std::cout<<"double of i is "<<di<<std::endl;
    file1.write((char*)&di,sizeof(double));
    file2.write((char*)&di,sizeof(double));
    file3.write((char*)&di,sizeof(double));
    
    //file1 << i; file2 << i; file3 << i;
    for(int inode=0;inode<npts;inode++){
      //std::cout<<defNodes[inode]->getPoint(0)<<std::endl;
      //	file1 << "\t" << I1[inode];
      //	file2 << "\t" << I2[inode];
      //	file3 << "\t" << I3[inode];
      //file1.write(&tab,sizeof(char));
      //file2.write(&tab,sizeof(char));
      //file3.write(&tab,sizeof(char));
      file1.write((char*)&I1[inode],sizeof(double)); 
      file2.write((char*)&I2[inode],sizeof(double));
      file3.write((char*)&I3[inode],sizeof(double));
    }
    //file1 << std::endl; file2 << std::endl; file3 << std::endl;
    //file1.write(&newline,sizeof(char));
    //file2.write(&newline,sizeof(char));
    //file3.write(&newline,sizeof(char));
  }
  std::cout<<"Everything done!"<<std::endl;
  return 0;
}

void readconfig_dat(char* filename, std::vector< std::vector<int> > & tableid, float* x, float* y, float* z, int npts){
  ifstream inpfile;
  inpfile.open(filename);
  if (!inpfile) {
    cout << "Cannot open input file: " << filename << endl;
    exit(0);
  }
  int id;
  double xt,yt,zt;
  inpfile >> id >> xt >> yt >> zt;
  std::vector<double> xtemp,ytemp,ztemp;
  
  do{
    xtemp.push_back(xt);
    ytemp.push_back(yt);
    ztemp.push_back(zt);
    inpfile >> id >> xt >> yt >> zt;
  }while(!inpfile.eof());

  for(int i=0; i<npts; i++){
    double xt=0.,yt=0.,zt=0.;
    for(int j=0;j<tableid[i].size();j++){
      xt += xtemp[tableid[i][j]];
      yt += ytemp[tableid[i][j]];
      zt += ztemp[tableid[i][j]];
    }
    x[i]=xt/tableid[i].size(); y[i]=yt/tableid[i].size(); z[i]=zt/tableid[i].size();
  }
  inpfile.close();
}

void readconfig_pdb(char* filename, float* x, float* y, float* z, int npts){
  ifstream inpfile;
  inpfile.open(filename);
  if (!inpfile) {
    cout << "Cannot open input file: " << filename << endl;
    exit(0);
  }
  int i=0;
  string line;
  
  while (! inpfile.eof() )
    {
      getline (inpfile,line);
      if (line.length()>0 && line.substr(0,4)=="ATOM"){
	x[i] = atof(line.substr(30,8).c_str()); 
	y[i] = atof(line.substr(38,8).c_str()); 
	z[i] = atof(line.substr(46,8).c_str()); 
	i++;
      }
    }
	inpfile.close();
}

void readconfig_dat(char* filename, float* x, float* y, float* z, int npts, int average){
  ifstream inpfile;
  inpfile.open(filename);
  if (!inpfile) {
    cout << "Cannot open input file: " << filename << endl;
    exit(0);
  }
  int id;
  
  for(int i=0; i<npts; i++){
    double xm=0., ym=0., zm=0.;
    for(int j=0; j<average; j++){
      double xj,yj,zj;
      inpfile >> id >> xj >> yj >> zj;
      xm += xj; ym += yj; zm += zj;
    }
    xm /= average; ym /= average; zm /= average;
    x[i]=xm; y[i]=ym; z[i]=zm;
  }
  inpfile.close();
}
