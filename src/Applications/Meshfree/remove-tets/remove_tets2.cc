#include <string>
// This file reads in the nodes and element (tets) data
// created by tetgen
// then reads in the surface created by MSMS (their positions and normals)
// for each tet, it finds the closest surface point and calculates the dot product of normal and the position vector from centroid to the surface point
// if the dot product is negative, centroid is outside the surface, otherwise its inside the surface
// it outputs the tets and their "dot-product" in a .vtk file, which can be read into cal_vol_supp_threshold.py script to create MF input
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include "PointBin.h"

using namespace std;
using namespace voom;

int main(int argc, char* argv[])
{
	if( argc < 2 ) {
		cout << "Usage: indent modelName [-c compress/no -H final_height -b boundary_conditions -n no_configs -s support_size]." << endl;
		return(0);
	}
	
	bool verbose=true;
	
	/*for(int i=0; i<argc; i++) {
		std::cout << std::setw(8) << i << "\t"
		<< argv[i] << std::endl;
	}*/
	//Input file name and other parameters from run command
	string objFileName = argv[1];
	ifstream ifs;
	ofstream ofs;
	string modelName = argv[1];
	string nodeFileName = modelName + ".1.node";
	string eleFileName = modelName + ".1.ele";
	string surfaceFileName = modelName + ".vert";
	string outFileName = modelName + "-sum-shape.vtk";
	/*
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
	*/
	
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

	std::vector< std::vector<double> > nodes;
	//vector<double> node_volume,supp_size;
	
	// Input .dat file containing nodes and weights
	std::string token;
	int npts=0;
	ifs >> npts; 
	ifs >> token >> token >> token;
	nodes.reserve(npts);
	ofs << "# vtk DataFile Version 3.0 \nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS "<<npts<<" float"<<std::endl;

	// read in points
	for(int i=0; i<npts; i++) {
		int id=i;
		int temp;
		std::vector<double> x(3);
		ifs >> temp >> x[0] >> x[1] >> x[2];
	//	x/= 10.0;
	//	NodeBase::DofIndexMap idx(3);
	//	for(int j=0; j<3; j++) idx[j]=dof++;
	//	DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
		nodes.push_back( x );
		ofs << x[0] << "\t" << x[1] << "\t" << x[2] << std::endl;
		}
	cout << "Number of nodes: " <<nodes.size() << endl;
	ifs.close();

	// open surface file
	ifs.open(surfaceFileName.c_str(),ios::in);
	if(!ifs){
		std::cout<<"Error: Cannot open the file "<<surfaceFileName<<std::endl;
		exit(0);
        }
	typedef std::vector<double> pos;
        std::vector<pos> vert;
        std::vector<pos> normal;
        int nsurf;
        string line;
        getline(ifs,line);
        getline(ifs,line);
        double junk;
        ifs >> nsurf >> junk >> junk >> junk;
        double x,y,z,nx,ny,nz;
        vert.reserve(nsurf);
        normal.reserve(nsurf);
	std::cout<<"number of surface points: "<<nsurf<<std::endl;
        for(int i=0;i<nsurf;i++){
                ifs >> x >> y >> z >> nx >> ny >> nz >> junk >> junk >> junk;
		//std::cout<<"x "<<x<<" y "<<y<<" z "<<z<<std::endl;
		//std::cout<<"nx "<<nx<<" ny "<<ny<<" nz "<<nz<<std::endl;
		//getchar();
                pos temp(3); temp[0]=x; temp[1]=y; temp[2]=z;
                pos norm(3); norm[0]=nx; norm[1]=ny; norm[2]=nz;
                vert.push_back(temp);
                normal.push_back(norm);
        }
	PointBin grid(10,10,1,vert); 
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
	vector<double> dot_prod;
	dot_prod.reserve(nElem);
	std::cout<<"Calculating the dot product.\nNumber of elements completed:"<<std::endl;
	for(int i=0;i<nElem;i++){
		int node1, node2, node3, node4, temp;
		ifs >> temp >> node1 >> node2 >> node3 >> node4;
		//centeroid of the tetrahedra
		double x = (nodes[node1-1][0] + nodes[node2-1][0] + nodes[node3-1][0] + nodes[node4-1][0])/4.;
		double y = (nodes[node1-1][1] + nodes[node2-1][1] + nodes[node3-1][1] + nodes[node4-1][1])/4.;
		double z = (nodes[node1-1][2] + nodes[node2-1][2] + nodes[node3-1][2] + nodes[node4-1][2])/4.;

		double mindis=numeric_limits<double>::max();
		int minJ;
		//find the nearest surface point
		for(int j=0;j<nsurf;j++){
			//calculate distance from x,y,z
			double dis = (vert[j][0]-x)*(vert[j][0]-x) + (vert[j][1]-y)*(vert[j][1]-y) + (vert[j][2]-z)*(vert[j][2]-z);
			if(dis<mindis) {mindis=dis; minJ=j;}
		}
		//calcualte the dot product
		dot_prod.push_back( (vert[minJ][0]-x)*normal[minJ][0] + (vert[minJ][1]-y)*normal[minJ][1] + (vert[minJ][2]-z)*normal[minJ][2]);
		if(i%1000==0) std::cout<<i<<"\t"<<flush;
		ofs << "4\t" << node1-1 << "\t" << node2-1 << "\t" << node3-1 << "\t" << node4-1 << std::endl;
	}
	std::cout<<std::endl;
	std::cout<<"Dot products calculated. Writing to the .vtk file"<<std::endl;
	ifs.close();
	//writing the element types
	ofs << "\nCELL_TYPES   "<<nElem<<std::endl;
	for(int i=0;i<nElem;i++) ofs << "10" <<std::endl;
	
	//writing the sum of shape functions
	ofs << "\nCELL_DATA   " <<nElem << std::endl << "SCALARS density float 1" << std::endl << "LOOKUP_TABLE default\n"<<std::endl;
	for(int i=0;i<nElem;i++) ofs << dot_prod[i] << std::endl;
	std::cout << "VTK file written."<<std::endl;
	ofs.close();
}
