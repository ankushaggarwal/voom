
//#include <string>
// This file reads in the nodes and element (tets) data
// created by tetgen
// then calculates the sum of RK shape function at the centroid of each tet
// and output it in a vtk format to visualize in paraview
// after visualizing it, a cutoff can be decided and using cal_vol_new.py script 
// volume of the tets that pass the cutoff criterion can be calculated to form the final meshfree input file
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace std;
#include "MFShape.h"
//#include "MFShape.cc"
#include "NodeBase.h"
#include "Node.h"
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>
using namespace voom;

int main(int argc, char* argv[])
{
	if( argc < 2 ) {
		cout << "Usage: indent modelName [-c compress/no -H final_height -b boundary_conditions -n no_configs -s support_size]." << endl;
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
	ofstream ofs;
	string modelName = argv[1];
	string nodeFileName = modelName + ".1.node";
	string eleFileName = modelName + ".1.ele";
	string outFileName = modelName + "-sum-shape.vtk";
	
	bool badCommandLine=false;
	double support_size;
	int nConfig=0;
	
	for(char option; (option=getopt(argc,argv,"c:H:b:n:s:")) != EOF; ) 
		switch (option) {
			default :
			badCommandLine=true;
			break;
			case 's' :
			support_size = std::atof(optarg);
			std::cout << "Global support size: " << support_size <<std::endl;
			break;
		}
	if( badCommandLine || argc < 2 ) {
		cout << "Usage: capsid InputFileName [-c -H finalHeight -b boundaryConditionFlag -R restartFilename]." << endl;
		return(0);
	}
	
	
	// open node file
	ifs.open( nodeFileName.c_str(), ios::in);
	if (!ifs) {
		std::cout << "Cannot open node file: " << nodeFileName << std::endl;
		exit(0);
	}
	// open output file
	ofs.open(outFileName.c_str());
	if (!ofs){
		std::cout << "Cannot open output file: " <<outFileName << std::endl;
		exit(0);
	}
	// create vector of nodes
	int dof=0;

	std::vector< NodeBase* > nodes;
	std::vector< DeformationNode<3>* > defNodes;
	vector<double> node_volume,supp_size;
	
	// Input .dat file containing nodes and weights
	std::string token;
	int npts=0;
	ifs >> npts; 
	ifs >> token >> token >> token;
	defNodes.reserve(npts);
	ofs << "# vtk DataFile Version 3.0 \nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "<<npts<<" float"<<std::endl;

	// read in points
	for(int i=0; i<npts; i++) {
		int id=i;
		int temp;
		DeformationNode<3>::Point x;
		ifs >> temp >> x(0) >> x(1) >> x(2);
		x/= 10.0;
		NodeBase::DofIndexMap idx(3);
		for(int j=0; j<3; j++) idx[j]=dof++;
		DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
		nodes.push_back( n );
		defNodes.push_back( n );
		supp_size.push_back(support_size);
		ofs << x(0) << "\t" << x(1) << "\t" << x(2) << std::endl;
		}
	cout << "Number of nodes: " <<nodes.size() << endl;
	ifs.close();
	
	// open element file
	ifs.open( eleFileName.c_str(), ios::in);
	if (!ifs) {
		std::cout << "Cannot open element file: " << eleFileName << std::endl;
		exit(0);
	}
	int nElem=0;
	ifs >> nElem;
	std::cout<<"Number of elements is "<<nElem<<std::endl;
	ofs <<"\nCELLS "<<nElem<<"   "<<nElem*5<<std::endl;

	ifs >> token >>token;
	vector<double> sum;
	sum.reserve(nElem);
	std::cout<<"Calculating the sum of shape functions.\nNumber of elements completed:"<<std::endl;
	for(int i=0;i<nElem;i++){
		int node1, node2, node3, node4, temp;
		ifs >> temp >> node1 >> node2 >> node3 >> node4;
		//centeroid of the tetrahedra
		double x = (defNodes[node1-1]->getPosition(0)+defNodes[node2-1]->getPosition(0)+defNodes[node3-1]->getPosition(0)+defNodes[node4-1]->getPosition(0))/4.;
		double y = (defNodes[node1-1]->getPosition(1)+defNodes[node2-1]->getPosition(1)+defNodes[node3-1]->getPosition(1)+defNodes[node4-1]->getPosition(1))/4.;
		double z = (defNodes[node1-1]->getPosition(2)+defNodes[node2-1]->getPosition(2)+defNodes[node3-1]->getPosition(2)+defNodes[node4-1]->getPosition(2))/4.;
		//calcualte the sum of shape functions at the centroid and push back in to array "sum"
		voom::MFShape::CoordinateArray cent(x,y,z);
		voom::MFShape shape;
		sum.push_back(shape.partition_of_unity(cent,defNodes,supp_size));
		if(i%1000==0) std::cout<<i<<"\t"<<flush;
		ofs << "4\t" << node1-1 << "\t" << node2-1 << "\t" << node3-1 << "\t" << node4-1 << std::endl;
	}
	std::cout<<std::endl;
	std::cout<<"Sum of shape functions calculated. Writing to the .vtk file"<<std::endl;
	ifs.close();
	//writing the element types
	ofs << "\nCELL_TYPES   "<<nElem<<std::endl;
	for(int i=0;i<nElem;i++) ofs << "10" <<std::endl;
	
	//writing the sum of shape functions
	ofs << "\nCELL_DATA   " <<nElem << std::endl << "SCALARS density float 1" << std::endl << "LOOKUP_TABLE default\n"<<std::endl;
	for(int i=0;i<nElem;i++) ofs << sum[i] << std::endl;
	std::cout << "VTK file written."<<std::endl;
	ofs.close();
}
