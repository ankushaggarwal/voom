#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <unistd.h>
#include <tvmet/Vector.h>
#include "Node.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "CGfast.h"
#include "Spring.h"
#include "GenericBody.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

using namespace voom;

// Get current date/time, format is YYYY-MM-DD-HH-mm-ss
void write_output(char* filename, const std::vector<double> &node_energy, const std::vector< DeformationNode<3>* > &defNodes, std::vector< std::vector<int> > &table);
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d-%H-%M-%S", &tstruct);

    return buf;
}

/*! This program creates an elastic network for proteins
   then solves for given BCs and, if applied, a prestrain
*/
int main(int argc, char* argv[])
{
  if(argc < 2){
    std::cout<<"Usage minimize -n NodeFileName -s SpringsFileName -b BCFileName -t TableFileName -p prestrain-value" << std::endl;
    return(0);
  }
 
  // get command-line arguments
  bool verbose=true;
  for(int i=0; i<argc; i++) {
    std::cout << std::setw(8) << i << "\t"
	      << argv[i] << std::endl;
  }
  int option;
  bool badCommandLine=false;
  std::string inputSpringsFileName, inputBCFileName, inputNodesFileName, inputTableFileName;
  std::ifstream ifs;
  float prestrain;

  while( (option=getopt(argc,argv,"s:b:n:t:p:")) != -1 ) {
    std::cout << "option = " << char(option) << std::endl;
    switch (option) {
    case 's':
      inputSpringsFileName = optarg;
      std::cout << "Springs file name: " << inputSpringsFileName << std::endl;
      break;
    case 'b':
      inputBCFileName = optarg;
      std::cout << "BC file name: " << inputBCFileName << std::endl;
      break;
    case 'n':
      inputNodesFileName=optarg;
      std::cout << "Nodes file name: " << inputNodesFileName << std::endl;
      break;
    case 't':
      inputTableFileName=optarg;
      std::cout << "Table file name: " << inputTableFileName << std::endl;
      break;
    case 'p':
      prestrain=std::atof(optarg);
      std::cout << "Prestrain: " << prestrain << " will be applied." << std::endl;
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
  // Read in the nodes
  //
  // create input stream
  ifs.open( inputNodesFileName.c_str(), std::ios::in);
  if (!ifs) {
    std::cout << "Cannot open input file: " << inputNodesFileName << std::endl;
    exit(0);
  }

  // create vector of nodes
  int dof=0,npts;
  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;
  
  // Input .vtk file containing nodes and connectivities
  ifs >> npts; 
  defNodes.reserve(npts);
  std::cout<<"No. of nodes = "<<npts<<std::endl;
  // read in points
  for(int i=0; i<npts; i++) {
    int id=i;
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  std::cout << "Number of nodes: " <<nodes.size() << std::endl;
  ifs.close();

  // default springs constant
  double k=1;

  //
  // Read in the springs
  //
  // create input stream
  ifs.open( inputSpringsFileName.c_str(), std::ios::in);
  if (!ifs) {
    std::cout << "Cannot open input file: " << inputSpringsFileName << std::endl;
    exit(0);
  }
    
  // read in triangle connectivities
  typedef Spring<3> spring_element;
  std::vector<spring_element*> spring_container;
  Body::ElementContainer elem_container;
  //std::vector< tvmet::Vector<int,2> > connectivities;
  int c1, c2;
  int nsprings=0; ifs >> nsprings;
  spring_container.reserve(nsprings);
  std::cout << "Number of springs: " << nsprings << std::endl;
  for (int i = 0; i<nsprings; i++){
    ifs >> c1 >> c2;
    Spring<3> * sp = new Spring<3>(defNodes[c1],defNodes[c2],k);
    sp->setIDA(c1);
    sp->setIDB(c2);
    spring_container.push_back( sp );
    elem_container.push_back(sp);
  }
  ifs.close();

  //
  // Read in the table for superatoms
  //
  ifs.open( inputTableFileName.c_str(), std::ios::in);
  if (!ifs) {
    std::cout << "Cannot open input file: " << inputTableFileName << std::endl;
    exit(0);
  }
  std::vector < std::vector<int> > table;
  int nsupatoms;
  ifs >> nsupatoms;
  for(int i=0;i<nsupatoms;i++){
    std::vector<int> temp;
    int np;
    ifs >> np;
    for(int j=0;j<np;j++){
      int p;
      ifs >> p;
      temp.push_back(p);
    }
    table.push_back(temp);
  }
  ifs.close();

  //
  // Read in the BCs
  //
  ifs.open( inputBCFileName.c_str(), std::ios::in);
  if (!ifs) {
    std::cout << "Cannot open input file: " << inputBCFileName << std::endl;
    exit(0);
  }

  std::vector<bool> BCindex(3*npts,false);
  std::vector<double> BCvalue(3*npts,0.);

  while(!ifs.eof()){
    int index;
    double value;
    ifs >> index >> value;
    BCindex[index]=true;
    BCvalue[index]=value;
  }

  ifs.close();
  // Crate a folder for the new output files named after current date and time
  std::string folder=currentDateTime();
  mkdir(folder.c_str(),0775);

  GenericBody bd(elem_container,nodes);

  //randomly assign pre-strain to all the springs
  srand(18756);
  float HI=prestrain,LO=-prestrain;
  for(int i=0;i<nsprings;i++){
    float random= LO + (float)rand()/((float)RAND_MAX/(HI-LO));
    spring_container[i]->resetLength(spring_container[i]->getLength()* (1+random) );
  }

  bd.compute(true,false,false);
  std::cout<<"Initial energy of the body is "<<bd.energy()<<std::endl;
  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(&bd);
  Model model(bdc,nodes);

  // consistency check
  //model.checkConsistency(true,false);
  //model.checkRank(model.dof()-6,true);
  //return 0;
  // set the solver parameters
  int m=5;
  double factr=1.0e7;
  double pgtol=1.5e-5;
  int iprint  =1;
  int maxIter=50000;
  std::ifstream lbfgsbinp("lbfgsb.inp");
  lbfgsbinp >> iprint >> factr >> pgtol >> m;
  if(verbose) 
    std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl;
  
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxIter );

  // set up bounds for solver
  blitz::Array<int,1> nbd(3*nodes.size());
  blitz::Array<double,1> lo(3*nodes.size());
  blitz::Array<double,1> hi(3*nodes.size());
  nbd = 0; //both upper and lower bound will be assigned
  lo = 0.0;
  hi = 0.0;
  double delta=0.,BCfactor;

  int MAXITER=20;
  BCfactor=5./MAXITER;
  for(int a=0; a<defNodes.size(); a++) {
    Vector3D pos;
    pos = defNodes[a]->point();
    for(int i=0;i<3;i++){
      if(BCindex[3*a+i]){
        nbd(3*a+i)=2;
        hi(3*a+i)=pos(i)+delta; //+BCfactor*BCvalue[3*a+i]+delta;
        lo(3*a+i)=pos(i)-delta; //+BCfactor*BCvalue[3*a+i]-delta;
      }
    }
  }
  ///////////////////////////////////////////////////////////////
  //                   START THE ITERATION                     //
  ///////////////////////////////////////////////////////////////
  for(int iter=0;iter<MAXITER;iter++){
    //BCfactor=iter*5./MAXITER;
    if(iter>0){
      for(int a=0; a<defNodes.size(); a++) {
        //Vector3D pos;
        //pos = defNodes[a]->point();
        for(int i=0;i<3;i++){
          if(BCindex[3*a+i]){
            //nbd(3*a+i)=2;
            hi(3*a+i) += BCfactor*BCvalue[3*a+i];
            lo(3*a+i) += BCfactor*BCvalue[3*a+i];
            //lo(3*a+i)=pos(i)-delta;
          }
        }
      }
    }

    solver.setBounds(nbd,lo,hi);

    // solve to minimize the energy
    solver.solve(&model);
    bd.compute(true,false,false);

    /////////////////////////OUTPUT////////////////////////////////////
    std::cout<<"Final energy of the body is "<<bd.energy()<<std::endl;
    std::vector<double> node_energy(npts,0.0);
    for(int i=0;i<nsprings;i++){
      double spring_energy=spring_container[i]->energy();
      node_energy[spring_container[i]->getIDA()] += spring_energy/2.;
      node_energy[spring_container[i]->getIDB()] += spring_energy/2.;
    }  
    char filename[35];
    sprintf(filename,"%s/%d.dat",folder.c_str(),iter);
    write_output(filename, node_energy, defNodes, table);
  }
  ///////////////////////////////////////////////////////////////
  //                     END THE ITERATION                     //
  ///////////////////////////////////////////////////////////////
  return 0;
}

void write_output(char* filename, const std::vector<double> &node_energy, const std::vector< DeformationNode<3>* > &defNodes, std::vector< std::vector<int> > &table){
  std::ofstream outfile;
  outfile.open(filename);
  if (!outfile) {
    std::cout << "Cannot open output file: " << filename << std::endl;
  }
  //outfile << "asdfasdf"<<std::endl;
  int npts=table.size();
  for(int i=0;i<npts;i++){
    double x=0.,y=0.,z=0.,e=0.;
    for(int j=0;j<table[i].size();j++){
      x += defNodes[table[i][j]]->getPoint(0);
      y += defNodes[table[i][j]]->getPoint(1);
      z += defNodes[table[i][j]]->getPoint(2);
      e += node_energy[table[i][j]];
    }
    x /= table[i].size();
    y /= table[i].size();
    z /= table[i].size();
    outfile << x << "\t" << y << "\t" << z << "\t" << e << std::endl;
  }
  /*int npts=node_energy.size();
  for(int i=0;i<npts;i++){
    outfile << defNodes[i]->getPoint(0) << "\t" << defNodes[i]->getPoint(1) << "\t" << defNodes[i]->getPoint(2) << std::endl;
  }*/
  outfile.close();
}
