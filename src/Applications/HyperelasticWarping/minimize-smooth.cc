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
#include "WarpingPenalty.h"
#include "LoopShellBody.h"
#include "C0MembraneBody.h"
#include "SkewHexonBody.h"
#include "WarpingBody.h"
#include "ShapeTri3.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "CGfast.h"
#include <stdlib.h>
#include <srfftw.h>
#include "DeformableContact.h"

using namespace voom;

#define idm(i,j,k) (in_data[(i)+2*(nx/2+1)*((j)+ny*(k))])
#define ndm(i,j,k) (new_data[(i)+nx*((j)+ny*(k))])
#define fm(i,j,k) (filter[(i)+2*(nx/2+1)*((j)+ny*(k))])
void print_body(Model::BodyContainer bdc);
void gaussian_smooth_fft(std::vector<double> &in_data, double sigma, int nx, int ny, int nz);

/*! This program creates a shell model with LoopShellBody for bending,
   C0MembraneBodies for stretching and one WarpingBody for penalty of
  difference between template and target image. The total energy is minimized for 
  image registration
 */
int main(int argc, char* argv[])
{
  if(argc < 2){
    std::cout<<"Usage minimize modelName -s sigma -r restartFileName -R RestartFileName -T templateFileName -S targetFileName" << std::endl;
    return(0);
  }
 
  // get command-line arguments
  bool verbose=true;
  for(int i=0; i<argc; i++) {
    std::cout << std::setw(8) << i << "\t"
	      << argv[i] << std::endl;
  }
  //double c1=0.,c2=0.,in_eta=0.;
  //const double Rcapsid = 30.;
  //int option;
  //bool badCommandLine=false;
  std::string modelName = argv[1];
  std::ifstream ifs;
  int option;
  double sigma=0.;
  bool badCommandLine=false;
  bool restart=false;
  bool restartRef=false,templateNameInput=false,targetNameInput=false;
  std::string restartFileName,targetFileName,templateFileName;

  while( (option=getopt(argc,argv,"s:r:R:S:T:")) != -1 ) {
    std::cout << "option = " << char(option) << std::endl;
    switch (option) {
      // -s sigma 
      //
      // the gaussian width for smoothing
    case 's':
      sigma = std::atof(optarg);
      std::cout << "sigma: " << sigma << std::endl;
      break;
      // -r restartfilename 
      //
      // whether restart or not with displacements 
      // added to the current position of nodes
    case 'r':
      restart=true;
      restartFileName = optarg;
      std::cout << "Restarting from file: " << restartFileName << std::endl;
      break;
      // -R restartfilename
      //
      // whether restart or not with displacements 
      // added to the reference position of nodes
    case 'R':
      restartRef=true;
      restartFileName = optarg;
      std::cout << "Restarting from file: " << restartFileName << std::endl;
      break;
      // -T templatefilename
      //
      // filename of the template image
    case 'T':
      templateNameInput=true;
      templateFileName = optarg;
      std::cout << "Template File: " << templateFileName << std::endl;
      break;
      // -S targetfilename
      //
      // filename of the target image
    case 'S':
      targetNameInput=true;
      targetFileName = optarg;
      std::cout << "Target File: " << targetFileName << std::endl;
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
  //different format connectivity for stretching body
  std::vector< std::vector<int> > s_connectivities;
  for(int i=0; i<connectivities.size(); i++) {
    std::vector<int> c(3);
    for(int j=0; j<3; j++) c[j]=connectivities[i](j);
    s_connectivities.push_back(c);
  }
/*  
  //read the two contact leaflets
  std::vector< std::vector<int> > large_leaflet_connect;
  while(token!="large_leaflet") ifs >> token;
  int datasize;
  ifs >> token >> datasize >> token;
  if(datasize!=ntri) std::cerr <<"Error while reading large leaflet component, data size not right" << std::endl;
  for(int i=0;i<ntri;i++){
    int flag;
    ifs >> flag;
    if(flag==1)large_leaflet_connect.push_back(s_connectivities[i]);
  }

  //small leaflet
  std::vector<int> small_leaflet;
  while(token!="small_leaflet") ifs >> token;
  ifs >> token >> datasize >> token;
  if(datasize!=npts) std::cerr <<"Error while reading large leaflet component, data size not right" << std::endl;
  for(int i=0;i<npts;i++){
    int flag;
    ifs >> flag;
    if(flag==1)small_leaflet.push_back(i);
  }*/
  //reading the modulus
  /*std::vector<double> modulus;
  modulus.reserve(ntri);
  while(token!="modulus") ifs >> token;
  int datasize;
  ifs >> token >> datasize >> token;
  if(datasize!=ntri) std::cerr <<"Error while reading large leaflet component, data size not right" << std::endl;
  for(int i=0;i<ntri;i++){
    int flag;
    ifs >> flag;
    if(flag==1)modulus.push_back(0.05);
    else modulus.push_back(1.);
  }*/
  std::cout << "Finished reading the input file" <<std::endl; 
  ifs.close();

  //
  // image registration penalty function
  //

  //create template and target image voxel data objects
  int nx,ny,nz;
  double xmin,ymin,zmin;
  double dx,dy,dz;
  int nVoxels;
  //first template image in vtk structured points format
  if(templateNameInput)
    inputFileName=templateFileName;
  else
    inputFileName = "template.vtk";
  // create input stream
  bool binary=false;
  if(binary)
  ifs.open( inputFileName.c_str(), std::ios::in|std::ios::binary);
  else
  ifs.open( inputFileName.c_str(), std::ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << inputFileName << endl;
    exit(0);
  }
  std::cout << "Opened file: " << inputFileName << std::endl;
  ifs >> token;
  while( token != "DIMENSIONS" ) {ifs >> token;
  //std::cout<<token<<std::endl; 
  }
  ifs >> nx >> ny >> nz
      >> token
      >> xmin >> ymin >> zmin
      >> token
      >> dx >> dy >> dz
      >> token
      >> nVoxels;
  while( token != "Density") ifs >> token;
  ifs >> token >> token >> token >> token;

  std::vector<double> voxel_data(nVoxels,0.);
  int temp1;
  for(int i=0;i<nVoxels;i++){
    //ifs.read((char*)&temp1,sizeof(int));
    //voxel_data[i]=(short)temp1;
    ifs >> voxel_data[i];
  }
  ifs.close();
  std::cout << "Finished reading file: " << inputFileName << std::endl;

  if(sigma>0.) //send for smoothing
  {
    gaussian_smooth_fft(voxel_data,sigma,nx,ny,nz);
  }
  VoxelData template_imagedata(xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,voxel_data);
  //then target image
  if(targetNameInput)
    inputFileName=targetFileName;
  else
  inputFileName = "target.vtk";
  // create input stream
  if(binary)
  ifs.open( inputFileName.c_str(), std::ios::in|std::ios::binary);
  else
  ifs.open( inputFileName.c_str(), std::ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << inputFileName << endl;
    exit(0);
  }
  std::cout << "Opened file: " << inputFileName << std::endl;

  ifs >> token; 
  while( token != "DIMENSIONS" ) ifs >> token;
  ifs >> nx >> ny >> nz
      >> token
      >> xmin >> ymin >> zmin
      >> token
      >> dx >> dy >> dz
      >> token
      >> nVoxels;
  while( token != "Density") ifs >> token;
  ifs >> token >> token >> token >> token;

  voxel_data.resize(nVoxels,0.);
  for(int i=0;i<nVoxels;i++)
    ifs >> voxel_data[i];
  ifs.close();
  std::cout << "Finished reading file: " << inputFileName << std::endl;
 
  if(sigma>0.) //send for smoothing
  {
    gaussian_smooth_fft(voxel_data,sigma,nx,ny,nz);
  }
  VoxelData target_imagedata(xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,voxel_data);
  voxel_data.clear();

  Vector3D x0(5,5,5);
  //double lambda=0.001;
  double lambda=1e-6;
  WarpingPenalty penalty(lambda,x0,&template_imagedata,&target_imagedata);
  typedef WarpingBody<TriangleQuadrature,WarpingPenalty,ShapeTri3> warping_body;
  warping_body * penalty_body = new warping_body(penalty,s_connectivities,nodes,1);
  penalty_body->compute(true,false,false);
  std::cout<<"Initial energy of the penalty body is "<<penalty_body->energy()<<std::endl;

  // Read in the starting displacement guess when restarting the simulation
  if(restart){
    std::ifstream guess(restartFileName.c_str());
    if (!guess) {
      cout << "Cannot open restart file: " << restartFileName << endl;
      return 0;
    }
    guess >> token;
    while (token != "double")
      guess >> token; // >> token >> token >> token >> token;
    for(int i=0; i<npts; i++) {
      int id=i;
      double dx=0.,dy=0.,dz=0.;
      guess >> dx >> dy >> dz;
      //defNodes[i]->addPoint(0,dx);
      //defNodes[i]->addPoint(1,dy);
      //defNodes[i]->addPoint(2,dz);
      defNodes[i]->setPoint(0,dx);
      defNodes[i]->setPoint(1,dy);
      defNodes[i]->setPoint(2,dz);
    }
    guess.close();
  }

  // Read in the starting displacement guess when restarting the simulation (reference for strain energy will be reset)
  std::vector<double> zero(3,0.0);
  std::vector< std::vector<double> > start_disp(npts,zero);
  if(restartRef){
    std::ifstream guess(restartFileName.c_str());
    if (!guess) {
      cout << "Cannot open restart file: " << restartFileName << endl;
      return 0;
    }
    guess >> token >> token >> token >> token >> token;
    for(int i=0; i<npts; i++) {
      int id=i;
      double dx=0.,dy=0.,dz=0.;
      guess >> dx >> dy >> dz;
      start_disp[i][0]=dx;
      start_disp[i][1]=dy;
      start_disp[i][2]=dz;
      defNodes[i]->addPoint(0,dx);
      defNodes[i]->addPoint(1,dy);
      defNodes[i]->addPoint(2,dz);
      defNodes[i]->resetPosition();
    }
    guess.close();
  }

  // default material constants

  // bending
  double Y=1.;
  double KC=0.1;
  double KG=-0.1;
  double nu=1./3.;
  double C0=0.;
  // stretching
  double mu=1.;
  double kS=2.;

  // now input (if file is available) user-defined material constants
  std::ifstream inp("parameters.inp");
  inp >> Y >> nu >> KC >> KG >> C0 >> mu >> kS;
  
  // print parameters out to screen
  if(verbose) 
    std::cout << " Y: " << Y << std::endl
	      << " nu: " << nu << std::endl
	      << " KC: " << KC << std::endl
	      << " KG: " << KG << std::endl
	      << " C0: " << C0 << std::endl
              << " mu: " << mu << std::endl
              << " kS: " << kS << std::endl;


  // create bending body; SCElastic only has curvature terms, no in-plane strain terms.
  int quadOrder=1;
  
  // bending material object, to be copied when body generates new elements
  typedef SCElastic BendingMaterial;
  BendingMaterial bending(KC,KG,C0);

  // Bending body will generate a mesh of Loop Subdivision shell
  // elements from the connectivity array and the container of node pointers.
  typedef LoopShellBody<BendingMaterial> LSB;
  LSB * bending_body = new LSB(bending, connectivities, nodes, quadOrder);

  // Use the body's output routines to generate vtk files with results
  bending_body->setOutput(paraview);
  bending_body->SetRefMeanCurvature();
  //bending_body->SetBendingModulus(modulus,KC);

  //create the deformable contact element
  //std::cout<<"Adding contact element"<<std::endl;
  //DeformableContact contact(defNodes,large_leaflet_connect,small_leaflet,3,0.1,1);
  //std::cout<<"Added contact element"<<std::endl;
  //std::cout<<"Initial energy of the contact = " << contact.energy() << std::endl;
  //bending_body->pushBack(&contact);
  //contact.printParaview("contact");

  // Stretching body using EvansElastic Material, which separates area strain from pure
  // shear strain (at constant area).  

  typedef EvansElastic StretchingMaterial;
  StretchingMaterial stretch(0.,0.,0.,mu,kS);

  // Use generic C0MembraneBody for stretching of pents, with the EvansElastic material.
  typedef C0MembraneBody<TriangleQuadrature,StretchingMaterial,ShapeTri3> membrane;
  membrane * stretch_body = new membrane(stretch,s_connectivities,nodes,quadOrder);
  stretch_body->setOutput(paraview);

  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(bending_body);
  bdc.push_back(stretch_body);
  bdc.push_back(penalty_body);
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
  int maxIter=15000;
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
  nbd = 2; //both upper and lower bound will be assigned
  lo = 0.0;
  hi = 0.0;
  double delta=10;

  for(int a=0; a<defNodes.size(); a++) {
    Vector3D pos;
    pos = defNodes[a]->point();
    for(int i=0;i<3;i++){
      hi(3*a+i)=pos(i)+delta;
      lo(3*a+i)=pos(i)-delta;
    }
  }

  solver.setBounds(nbd,lo,hi);


  penalty_body->compute(true,false,false);
  std::cout<<"Initial energy of the penalty body is "<<penalty_body->energy()<<std::endl;
  // solve to minimize the energy
  bdc[0]->printParaview("initial");
  solver.solve(&model);
  penalty_body->compute(true,false,false);
  //contact.printParaview("contact");
  std::cout<<"Energy of the penalty body is "<<penalty_body->energy()<<std::endl;
  for(int iter=0;iter<0;iter++){
  penalty_body->updateLagrangeMultiplier();
  //lambda *= 2;
  //penalty_body->setMaterialLambda(lambda);
  penalty_body->compute(true,false,false);
  std::cout<<"Energy of the penalty body is "<<penalty_body->energy()<<std::endl;
  solver.solve(&model);
  penalty_body->compute(true,false,false);
  std::cout<<"Final energy of the penalty body is "<<penalty_body->energy()<<std::endl;
  }
  bdc[0]->printParaview("final");
  //bdc[1]->printParaview("minimized");
  //print the postprocessing strains and stress
/*  hex_body->printParaviewEigVec("Hexbody");
  pentamer_body->printParaviewEigVec("Pentbody");

  hex_body->printParaview2("Hexbody2");
  pentamer_body->printParaview2("Pentbody2");

  //model.checkRank(model.dof()-6,true);

  //print the total energies of three bodies
  std::cout <<"Bending energy = "<<bending_body->totalStrainEnergy()<<std::endl;
  std::cout <<"Pentamer Stretching energy = "<<pentamer_body->totalStrainEnergy()<<std::endl;
  std::cout <<"Hexamer Stretching energy = "<<hex_body->totalStrainEnergy()<<std::endl;
*/

  //print the displacements
  std::ofstream ofs("disp.dat");
  if(!ofs)  std::cout<<"Cannot open output file for displacements."<<std::endl;
  ofs << "POINT_DATA " << defNodes.size() << std::endl
      << "VECTORS displacements double" << std::endl;
  //     << "LOOKUP_TABLE displacements" << endl; 
  for(int i=0; i<defNodes.size(); i++) {
    Vector3D nodalDisp;
    nodalDisp = defNodes[i]->point() - defNodes[i]->position();
    ofs << std::setprecision(16) 
        << nodalDisp(0) + start_disp[i][0]
        << '\t' <<nodalDisp(1) + start_disp[i][1]
        << '\t' <<nodalDisp(2) + start_disp[i][2]
        << std::endl;
  }
  ofs.close();
  //calculate the center
/*  tvmet::Vector<double,3> center(0,0,0);
  for(int i=0; i<defNodes.size(); i++) {
    center += defNodes[i]->point();
  }
  center /= nodes.size();
  std::cout << "Center = " << center << std::endl;

  //calculate the average radius and aspherity
  Ravg = 0;
  for(int i=0; i<defNodes.size(); i++) {
    Ravg += tvmet::norm2(defNodes[i]->point());
  }
  Ravg /= nodes.size();
  std::cout << "Ravg = " << Ravg << std::endl;

  double deltaR=0.;
  for(int i=0; i<defNodes.size(); i++) {
      double Rtemp = tvmet::norm2(defNodes[i]->point());
      deltaR += (Rtemp-Ravg)*(Rtemp-Ravg);
    }
  double aspherity = sqrt(deltaR/nodes.size())/Ravg;
  std::cout<<"Aspherity = "<<aspherity<<std::endl; */
  return 0;
}

void gaussian_smooth_fft(std::vector<double> &original_data, double sigma, int nx, int ny, int nz){
  int i,j,k,ii,jj,kk;
  fftw_real *in_data, *new_data, *filter;
  in_data = (fftw_real *)malloc(sizeof(fftw_real)*2*(nx/2+1)*ny*nz);
  new_data = (fftw_real *)malloc(sizeof(fftw_real)*nx*ny*nz);
  filter = (fftw_real *)malloc(sizeof(fftw_real)*2*(nx/2+1)*ny*nz);

  //do setup for fftw
  fftw_complex *IN, *FILTER, *NEW;
  NEW = (fftw_complex *)malloc(sizeof(fftw_complex)*(nx/2+1)*ny*nz);
  IN = (fftw_complex *)in_data;
  FILTER = (fftw_complex *)filter;
  fftw_real scale = 1.0/(nx*ny*nz);

  fftwnd_plan p, pinv;
  //pass dimensions in reverse order due to column major storage
  p = rfftw3d_create_plan(nz,ny,nx, FFTW_REAL_TO_COMPLEX,
			  FFTW_ESTIMATE|FFTW_IN_PLACE);
  pinv = rfftw3d_create_plan(nz,ny,nx, FFTW_COMPLEX_TO_REAL,
			     FFTW_ESTIMATE);

  //store values from 1D array
  int iter=0;
  for(k=0;k<nz;k++){
    for(j=0;j<ny;j++){
      for(i=0;i<2*(nx/2+1);i++){
	if(i<nx)
	  idm(i,j,k)=(fftw_real)original_data[iter++];
	else
	  idm(i,j,k) = 0.0; //pad extra space with zeros
      }
    }
  }
  double argscale = 0.5/sigma/sigma;

  //store gaussian filter in 3d real array

  // NOTE:  must translate the filter WITH WRAP-AROUND to avoid
  // imparting a translation on the output.  Put the peak at (0,0,0).
  // This amounts to using the distance from the closest corner
  // as the arguent to the gaussian.  

  double sumexp = 0.0;
  for(k=0;k<nz;k++){
    kk = k<nz/2? k:nz-k;
    for(j=0;j<ny;j++){
      jj = j<ny/2? j:ny-j;
      for(i=0;i<nx;i++){
	ii = i<nx/2? i:nx-i;
	double tmp = exp(-(ii*ii + jj*jj + kk*kk)*argscale);
	fm(i,j,k) = tmp;
	sumexp += tmp;
      }
    }
  }

  //normalize filter to add up to one
  for(k=0;k<nz;k++){
    for(j=0;j<ny;j++){
      for(i=0;i<nx;i++){
	fm(i,j,k) /= sumexp;
      }
    }
  }
  
  //transform the input data and the filter into the fourier domain
  rfftwnd_one_real_to_complex(p, in_data, NULL);
  rfftwnd_one_real_to_complex(p, filter, NULL);

  //perform convolution in fourier domain (as mult)
  for(i=0; i<(nx/2+1)*ny*nz; i++){
    NEW[i].re = (IN[i].re*FILTER[i].re - IN[i].im*FILTER[i].im)*scale;
    NEW[i].im = (IN[i].re*FILTER[i].im + IN[i].im*FILTER[i].re)*scale;
  }

  //invert for the output, new_data
  rfftwnd_one_complex_to_real(pinv, NEW, new_data);
  
  //clean up
  rfftwnd_destroy_plan(p);
  rfftwnd_destroy_plan(pinv);
  std::cout << "Finished smoothing" << std::endl << "Copying back to the original array" << std::endl;
  iter=0;
  for(int k=0;k<nz;k++){
    for(int j=0;j<ny;j++){
      for(int i=0;i<nx;i++){
        original_data[iter++] = (double)ndm(i,j,k);
      }
    }
  }
  //free the memory allocated by malloc
  free(in_data);
  free(new_data);
  free(filter);
  free(NEW);
  return;
}
