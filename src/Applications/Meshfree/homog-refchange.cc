//#include <string>
// This file reads in the configurations from MD
// either in .dat or .pdf formats
// then calculates the deformation gradient at each node
// then calculates the three invariants at each node
// and then calculates a functional = (sum over nodes) (I1bar_prob-3)^2 where I1bar_prob is the most probable I1bar
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
#include "Lbfgsb.h"
#include "CompNeoHookean.h"
#include "HomogMP.h"
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>
#include "ViscousRegularizer.h"
using namespace voom;

void readconfig_pdb(char* filename, float* x, float* y, float* z, int npts);
void readconfig_dat(char* filename, float* x, float* y, float* z, int npts);
void readconfig_dat(char* filename, float* x, float* y, float* z, int npts, int average);
double most_prob(blitz::Array<double,1> x);
double most_probJ(blitz::Array<double,1> x);
double most_probable(double* x); //given a distribution x finds the most probably value

int main(int argc, char* argv[])
{
  if( argc < 2 ) {
    cout << "Usage: homog inputFileName [-n start_config -N end_config -s support_size - k kernel_option -a average_number -b alpha]." << endl;
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
  double support_size, alpha;
  int StartConfig=0, EndConfig=0, avg=1;
  bool cubic_kernel=false, modified_kernel=false, singular_kernel=false;
  
  for(char option; (option=getopt(argc,argv,"n:N:s:k:a:b:")) != EOF; ) 
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
    case 'a' :
      avg=std::atoi(optarg);
      std::cout << "Average of "<<avg<<" trajectory points taken"<<std::endl;
      break;
    case 'b' :
      alpha = std::atof(optarg);
      std::cout << "alpha: " << alpha <<std::endl;
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
    cout << "Usage: homog inputFileName [-n start_configs -N end_config -s support_size - k kernel_option -a average_number -b alpha]." << endl;
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
  
  //read the displacements to be added to the reference config
  //ifstream dispfs;
  //dispfs.open("random.dat",ios::in);

 
  // read in points
  for(int i=0; i<npts; i++) {
    int id=i;
    int temp;
    double weight;
    double supphat;
    DeformationNode<3>::Point x;
    ifs >> temp >> x(0) >> x(1) >> x(2) >> weight>>supphat;
    //double dx,dy,dz;
    //dispfs >> dx >> dy >> dz;
    //x(0) += alpha*dx; x(1) += alpha*dy; x(2) += alpha*dz;
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
  //double rho = 1.0;
  //double E = 200.0;
  //double nu = 0.4;
  //double viscosity=1.0e-2;
  //typedef CompNeoHookean MaterialType;
  typedef HomogMP MaterialType;
  MaterialType protein( 1., 1., 1. );
  
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
  
  //open the output files
  //ofstream file1, file2, file3;
  //file1.open("inv1.dat",ios::binary);
  //file2.open("inv2.dat",ios::binary);
  //file3.open("inv3.dat",ios::binary);
  //char tab='\t', newline='\n';
  int col=npts, row=EndConfig-StartConfig+1;
  //file1.write((char*)&row,sizeof(int)); file1.write((char*)&col,sizeof(int));
  //file2.write((char*)&row,sizeof(int)); file2.write((char*)&col,sizeof(int));
  //file3.write((char*)&row,sizeof(int)); file3.write((char*)&col,sizeof(int));
  col/=60; row*=60; //using symmetry
  {
  blitz::Array<double,2> I1bar(row,col);
  blitz::Array<double,2> I2bar(row,col);
  blitz::Array<double,2> J(row,col);
  I1bar=0.; I2bar=0.; J=0.;


  //vary the reference configuration to see the change in J
  //random vector of length of DOF
  //std::vector<double> J(20,0.);
  //for(int irefconf=0;irefconf<20;irefconf++)
    //add disp*irefconf to reference configuration
    //recalculate the shape functions???????????


  //read in the configuration files in a loop
  for(int i=StartConfig;i<=EndConfig;i++){
    char inpfile[20];
    std::cout<<"Reading "<<i<<"th configuration."<<std::endl;
    float x[npts], y[npts], z[npts];
    if(datfile){
      //std::cout<<"using .dat files"<<std::endl;
      sprintf(inpfile,"com-traj/%d.dat",i);
      readconfig_dat(inpfile,x,y,z,npts,avg);}
    
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
    
    int idata=i-StartConfig;
    for(int inode=0;inode<npts;inode++){
      int icol=inode%col; 
      //if(icol>=col && icol<0) std::cout<<"Index out of range"<<std::endl;
      int irow=idata*60+inode/col;
      //if(irow>=row && irow <0) std::cout<<"Index out of range"<<std::endl;
      I1bar(irow,icol)=I1[inode]*pow(I3[inode]*I3[inode],-1./3.);
      I2bar(irow,icol)=I2[inode]*pow(I3[inode]*I3[inode],-2./3.);
      J(irow,icol) = I3[inode];
    }
  }
  std::cout<<"Calculated the strains for all states"<<std::endl;
  std::vector<double> I1MostProb(npts,3.), I2MostProb(npts,3.), JMostProb(npts,1.);
  std::cout<<"Calculated the most probable strain state invariants"<<std::endl;
  for(int inode=0;inode<npts/60;inode++){
    //call the function to get most probable I1bar
    std::cout<<inode<<std::endl;
    I1MostProb[inode] = 1*(most_prob(I1bar(blitz::Range::all(),inode))-3.)+3.;
    I2MostProb[inode] = 1*(most_prob(I2bar(blitz::Range::all(),inode))-3.)+3.;
    JMostProb[inode] = 1*(most_probJ(J(blitz::Range::all(),inode))-1.)+1.;
    //add to the functional (I1barMostProbable-3)*(I1barMostProbable-3)
    //J +=  (I1MostProb-3)*(I1MostProb-3) ;
  }
  for(int isym=1;isym<60;isym++){
    for(int inode=0;inode<npts/60;inode++){
      I1MostProb[isym*npts/60+inode] = I1MostProb[inode];
      I2MostProb[isym*npts/60+inode] = I2MostProb[inode];
      JMostProb[isym*npts/60+inode] = JMostProb[inode];
    }
  }
  std::cout<<"Calculated the most probable strain state invariants"<<std::endl;
  bd.setMPinv(I1MostProb,I2MostProb,JMostProb);
  }
  //reset the reference state
    for(int inode=0;inode<npts;inode++){
      defNodes[inode]->resetPoint();
    }
  bd.compute(true,true,false);
  
  
  std::cout << std::endl;
  std::cout << "Body energy = " << bd.energy() << std::endl;
  std::cout << "Body volume = " << bd.volume() << std::endl;
  std::cout << std::endl;
  return 0;
  //Creating the model to solve for the minima
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

  //add a viscous regularizer
  double viscosity = 10.;
  ViscousRegularizer vr(bd.nodes(), viscosity);
  //bd.addElement( &vr ); 

  // solve to minimize the functional
  solver.solve(&model);
  bd.printParaview("minimized2");
  //write the final result in the format of MFinput
  std::ofstream ofs("minimized.dat");
  if(!ofs)  std::cout<<"Cannot open output file for displacements."<<std::endl;
  for(int i=0;i<npts;i++){
    Vector3D x=defNodes[i]->point();
    ofs << i+1 << '\t' << x(0) << '\t' << x(1) << '\t' << x(2) << '\t' << node_volume[i] << '\t' << supp_size_hat[i]/supphatfactor << std::endl; 
  }
  ofs.close();
  std::cout<<"Everything done! "<<std::endl;//<<"J = "<<J<<std::endl;
  return 0;
}

void readconfig_dat(char* filename, float* x, float* y, float* z, int npts){
  ifstream inpfile;
  inpfile.open(filename);
  if (!inpfile) {
    cout << "Cannot open input file: " << filename << endl;
    exit(0);
  }
  int id;
  
  for(int i=0; i<npts; i++){
    inpfile >> id >> x[i] >> y[i] >> z[i];
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

double most_prob(blitz::Array<double,1> x){
  //find max, min and size of vector x (cut the limits at 3 and 6)
  double max,min;
  max=blitz::max(x(blitz::Range::all()));
  min=blitz::min(x(blitz::Range::all()));
  if(max>min+3.) {max=min+3.;}//this is only for I1bar, to prevent bad binning due to one large value of I1bar
  int size=x.rows();
  //create the bins
  int nbins=30;
  double interval=(max-min)/nbins;
  std::vector<double> bin_centers(nbins,0.);
  std::vector<int> counts(nbins,0);
 
  for(int i=0;i<nbins;i++){
    bin_centers[i] = min+(i+0.5)*interval;
  }

  //count the frequency in them
  for(int i=0;i<size;i++){
    int bin=int((x(i)-min)/interval);
    if(bin>nbins-1 || bin < 0) continue; //bin=nbins-1;
    counts[bin]++;
  }

  //find the bin with maximum counts
  int maxcount=0; int maxid=0;
  for(int i=0;i<nbins;i++){
    if(counts[i]>maxcount){
      maxcount=counts[i];
      maxid=i;
    }
  }
  //return the bin of maxid i.e. the most probable x
  return bin_centers[maxid];
}
 
double most_probJ(blitz::Array<double,1> x){
  //find max, min and size of vector x (cut the limits at 0.5 and 2)
  double max,min;
  max=blitz::max(x(blitz::Range::all()));
  min=blitz::min(x(blitz::Range::all()));
  if(max>2.) max=2.; 
  if(min<0.5) min=0.5; //this is only for J, to prevent bad binning due to one large or one small value of J
  int size=x.rows();
  //create the bins
  int nbins=30;
  double interval=(max-min)/nbins;
  std::vector<double> bin_centers(nbins,0.);
  std::vector<int> counts(nbins,0);
 
  for(int i=0;i<nbins;i++){
    bin_centers[i] = min+(i+0.5)*interval;
  }

  //count the frequency in them
  for(int i=0;i<size;i++){
    int bin=int((x(i)-min)/interval);
    if(bin>nbins-1 || bin < 0) continue; //bin=nbins-1;
    counts[bin]++;
  }

  //find the bin with maximum counts
  int maxcount=0; int maxid=0;
  for(int i=0;i<nbins;i++){
    if(counts[i]>maxcount){
      maxcount=counts[i];
      maxid=i;
    }
  }
  //return the bin of maxid i.e. the most probable x
  return bin_centers[maxid];
} 
