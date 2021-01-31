// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         Ankush Aggarwal
//                University of California Los Angeles
//                   (C) 2012 All Rights Reserved
//
//----------------------------------------------------------------------
// 

#if !defined(__DeformableContact_h__)
#define __DeformableContact_h__

#include <vector>
#include <cstdio>
#include <ctime>

#include "Element.h"
#include "NodeBin.h"

namespace voom {

  //! Smooth contact between two deformable surfaces
  /*! DeformableContact is an element which "simulates" contact
      by applying force on nodes close enough using a potential
      to make sure they don't come in actual contact
      Two sets of nodes on each surface of potential contact are input 
      with the "Master" nodes along with their connectivities and "Slave" nodes (without connectivity)
  */
  class DeformableContact :
    public Element {

    public:

    // typedefs
    //! 3-D deformation node
    typedef DeformationNode<3> DefNode;
    //! Container for node pointers
    typedef std::vector< DefNode* > DefNodeContainer;
    //! Iterator for node pointers
    typedef DefNodeContainer::iterator DefNodeIterator;
    //! Const iterator for node pointers
    typedef DefNodeContainer::const_iterator ConstDefNodeIterator;

    using Element::BaseNodeContainer; 
    using Element::BaseNodeIterator; 
    using Element::ConstBaseNodeIterator; 

    //! Construct from deformation nodes, connectivities of master nodes, a list of slave nodes, spacing of the grid, a cutoff distance and a "spring" constant 
    DeformableContact( const DefNodeContainer & Nodes, std::vector<std::vector<int> > MasterConnectivity, std::vector<int> SlaveNodes, double spacing, double cutoff, double k ) 
      {
	_defNodes = Nodes;
	_MasterConnectivity = MasterConnectivity;
	_SlaveNodes = SlaveNodes;
        _spacing = spacing;
	_cutoff = cutoff;
        _k = k;
	_energy = 0.0;
        std::vector<double> x(3,0.0);
        _MasterCentroids.clear();
        _MasterCentroids.resize(_MasterConnectivity.size(),x);
        _active_contacts=0;
        compute(true,true,false);
      }

    //! Compute the energy, force, and stiffness
    void compute(bool f0, bool f1, bool f2) {
      _active_contacts=0;
      if(f0) _energy = 0.0;
      
      //update the neighbor search data structure
      for(int i=0;i<_MasterConnectivity.size();i++){
        _MasterCentroids[i][0] = (_defNodes[_MasterConnectivity[i][0]]->getPoint(0) + _defNodes[_MasterConnectivity[i][1]]->getPoint(0) + _defNodes[_MasterConnectivity[i][2]]->getPoint(0) )/3.;

        _MasterCentroids[i][1] = (_defNodes[_MasterConnectivity[i][0]]->getPoint(1) + _defNodes[_MasterConnectivity[i][1]]->getPoint(1) + _defNodes[_MasterConnectivity[i][2]]->getPoint(1) )/3.;

        _MasterCentroids[i][2] = (_defNodes[_MasterConnectivity[i][0]]->getPoint(2) + _defNodes[_MasterConnectivity[i][1]]->getPoint(2) + _defNodes[_MasterConnectivity[i][2]]->getPoint(2) )/3.;
      }
      NodeBin grid(_spacing,_spacing,_spacing,_MasterCentroids);
      //Now start a check for each slave node for "possible contacts" 
      for(int i=0; i<_SlaveNodes.size();i++){
        DefNode* xstar=_defNodes[_SlaveNodes[i]];
        const Vector3D & xstar_pos = xstar->point();
        std::vector<double> temp(3,0.0); temp[0]=xstar_pos(0); temp[1]=xstar_pos(1); temp[2]=xstar_pos(2);
        //for each slave node, find the "close" neighbors from the search data structure
        std::vector<int> list;
        list=grid.GetList(temp);
        //find the closest of these neighbors
        double mindist = std::numeric_limits<double>::max();
        //double mindist = -std::numeric_limits<double>::max();
        int minid=-1;
        for(int j=0;j<list.size();j++){
          DefNode* x1=_defNodes[_MasterConnectivity[list[j]][0]]; const Vector3D & x1_pos = x1->point();
          DefNode* x2=_defNodes[_MasterConnectivity[list[j]][1]]; const Vector3D & x2_pos = x2->point();
          DefNode* x3=_defNodes[_MasterConnectivity[list[j]][2]]; const Vector3D & x3_pos = x3->point();
          //projection direction of node to the triangle is the normal of the triangle
          Vector3D dir = tvmet::cross(x1_pos-x2_pos,x1_pos-x3_pos);
          dir = dir/tvmet::norm2(dir);
          double dist=dot((xstar_pos-x1_pos),dir);  
          //closest point
          Vector3D xhat; xhat = xstar_pos - dist*dir;
          Vector3D totalarea = tvmet::cross(x2_pos-x1_pos,x3_pos-x1_pos);
          double totalarea_normsq = tvmet::dot(totalarea,totalarea);
          double N1,N2,N3;
          N1 = tvmet::dot(tvmet::cross(x2_pos-xhat,x3_pos-xhat),totalarea)/totalarea_normsq;
          N2 = tvmet::dot(tvmet::cross(x3_pos-xhat,x1_pos-xhat),totalarea)/totalarea_normsq;
          N3 = tvmet::dot(tvmet::cross(x1_pos-xhat,x2_pos-xhat),totalarea)/totalarea_normsq;
          if(N1>0 && N1<1 && N2>0 && N2<1 && N3>0 && N3<1){ //checking the barycentric coordinates
            if(mindist>std::abs(dist)) {mindist=std::abs(dist); minid=list[j];}
            //if(std::abs(mindist)>std::abs(dist)) {mindist=dist; minid=list[j];}
          }
        }
        //if the minimum distance is less than cutoff, do actual computation
        if(mindist < _cutoff){
        //if(mindist > -_cutoff){
          //three nodes of the closest triangle
          DefNode* x1=_defNodes[_MasterConnectivity[minid][0]]; const Vector3D & x1_pos = x1->point();
          DefNode* x2=_defNodes[_MasterConnectivity[minid][1]]; const Vector3D & x2_pos = x2->point();
          DefNode* x3=_defNodes[_MasterConnectivity[minid][2]]; const Vector3D & x3_pos = x3->point();
          //normal of the triangle
          Vector3D dir = tvmet::cross(x1_pos-x2_pos,x1_pos-x3_pos);
          dir = dir/tvmet::norm2(dir);
          double dist=dot((xstar_pos-x1_pos),dir);  //calculate dist again because it can be + or -, while mindist was kept absolute value
          //if(dist>0) std::cout<<"Contact the dist is positive"<<std::endl;
          //if(dist<0) std::cout<<"Contact the dist is negative"<<std::endl;
          //closest point
          Vector3D xhat; xhat = xstar_pos - dist*dir;
          if(f0 || f1) {
            Vector3D totalarea = tvmet::cross(x2_pos-x1_pos,x3_pos-x1_pos);
            double totalarea_normsq = tvmet::dot(totalarea,totalarea);
            //double totalarea=cal_area(x1_pos,x2_pos,x3_pos);
            double N1,N2,N3;
            N1 = tvmet::dot(tvmet::cross(x2_pos-xhat,x3_pos-xhat),totalarea)/totalarea_normsq;
            N2 = tvmet::dot(tvmet::cross(x3_pos-xhat,x1_pos-xhat),totalarea)/totalarea_normsq;
            N3 = tvmet::dot(tvmet::cross(x1_pos-xhat,x2_pos-xhat),totalarea)/totalarea_normsq;
            //if(N1>0 && N1<1 && N2>0 && N2<1 && N3>0 && N3<1){ //don't need to check the barycentric coordinates again
              _active_contacts++;
              //if(f0) _energy += 1./2.*_k*(_cutoff + dist)*(_cutoff + dist);
              //if(f0) _energy += 1./2.*_k*(_cutoff - std::abs(dist))*(_cutoff - std::abs(dist));
              //if(f0) _energy += _k*(1./std::abs(dist) - 1./_cutoff);
              if(f0) _energy += _k*(1./std::abs(dist)/std::abs(dist) - 1./_cutoff/_cutoff);
              if(f1){
                //double totalforce = -_k*(_cutoff + dist); //the first part of the expression is sign of dist
                //double totalforce = (dist>0?1:-1)*_k*(_cutoff - std::abs(dist)); //the first part of the expression is sign of dist
                //double totalforce = (dist>0?1:-1)*_k*(1./std::abs(dist)/std::abs(dist)); 
                double totalforce = (dist>0?1:-1)*2*_k/pow(std::abs(dist),3);
                /*if(dist<0)
                  totalforce = -_k*(_cutoff - std::abs(dist));
                else
                  totalforce = _k*(_cutoff - std::abs(dist));*/
                //total force in the normal direction on the slave node
                xstar->addForce(0,totalforce*dir(0)); xstar->addForce(1,totalforce*dir(1)); xstar->addForce(2,totalforce*dir(2));
                //three components of the total force on the master nodes in the opposite direction
                x1->addForce(0,-totalforce*N1*dir(0)); x1->addForce(1,-totalforce*N1*dir(1)); x1->addForce(2,-totalforce*N1*dir(2));
                x2->addForce(0,-totalforce*N2*dir(0)); x2->addForce(1,-totalforce*N2*dir(1)); x2->addForce(2,-totalforce*N2*dir(2));
                x3->addForce(0,-totalforce*N3*dir(0)); x3->addForce(1,-totalforce*N3*dir(1)); x3->addForce(2,-totalforce*N3*dir(2));
              }
           // }
          }
        }
      }
      //if(_active_contacts>0) {std::cerr <<"DeformableContact::The two surfaces are in contact now" << std::endl; printParaview("contact"); exit(0);}
    }

    //! print a paraview file for visualization
    void printParaview(const std::string name) 
    {
      //compute(true,true,false);
      std::vector<bool> nodes_active(_defNodes.size(),false);
      std::vector<double> nodes_mindist(_defNodes.size(),3*_spacing);
      std::vector<bool> elem_active(_MasterConnectivity.size(),false);
      //set forces on nodes equal to zero
      for(int i=0;i<_defNodes.size();i++){
        _defNodes[i]->setForce(0,0.);
        _defNodes[i]->setForce(1,0.);
        _defNodes[i]->setForce(2,0.);
      }
      //compute the active status and nodal forces
      //update the neighbor search data structure
      for(int i=0;i<_MasterConnectivity.size();i++){
        _MasterCentroids[i][0] = (_defNodes[_MasterConnectivity[i][0]]->getPoint(0) + _defNodes[_MasterConnectivity[i][1]]->getPoint(0) + _defNodes[_MasterConnectivity[i][2]]->getPoint(0) )/3.;

        _MasterCentroids[i][1] = (_defNodes[_MasterConnectivity[i][0]]->getPoint(1) + _defNodes[_MasterConnectivity[i][1]]->getPoint(1) + _defNodes[_MasterConnectivity[i][2]]->getPoint(1) )/3.;

        _MasterCentroids[i][2] = (_defNodes[_MasterConnectivity[i][0]]->getPoint(2) + _defNodes[_MasterConnectivity[i][1]]->getPoint(2) + _defNodes[_MasterConnectivity[i][2]]->getPoint(2) )/3.;
      }
      NodeBin grid(_spacing,_spacing,_spacing,_MasterCentroids);
      //Now start a check for each slave node for "possible contacts" 
      for(int i=0; i<_SlaveNodes.size();i++){
        DefNode* xstar=_defNodes[_SlaveNodes[i]];
        const Vector3D & xstar_pos = xstar->point();
        std::vector<double> temp(3,0.0); temp[0]=xstar_pos(0); temp[1]=xstar_pos(1); temp[2]=xstar_pos(2);
        //for each slave node, find the "close" neighbors from the search data structure
        std::vector<int> list;
        list=grid.GetList(temp);
        //find the closest of these neighbors
        double mindist = std::numeric_limits<double>::max();
        int minid=-1;
        for(int j=0;j<list.size();j++){
          DefNode* x1=_defNodes[_MasterConnectivity[list[j]][0]]; const Vector3D & x1_pos = x1->point();
          DefNode* x2=_defNodes[_MasterConnectivity[list[j]][1]]; const Vector3D & x2_pos = x2->point();
          DefNode* x3=_defNodes[_MasterConnectivity[list[j]][2]]; const Vector3D & x3_pos = x3->point();
          //projection direction of node to the triangle is the normal of the triangle
          Vector3D dir = tvmet::cross(x1_pos-x2_pos,x1_pos-x3_pos);
          dir = dir/tvmet::norm2(dir);
          double dist=dot((xstar_pos-x1_pos),dir);  
          //closest point
          Vector3D xhat; xhat = xstar_pos - dist*dir;
          Vector3D totalarea = tvmet::cross(x2_pos-x1_pos,x3_pos-x1_pos);
          double totalarea_normsq = tvmet::dot(totalarea,totalarea);
          double N1,N2,N3;
          N1 = tvmet::dot(tvmet::cross(x2_pos-xhat,x3_pos-xhat),totalarea)/totalarea_normsq;
          N2 = tvmet::dot(tvmet::cross(x3_pos-xhat,x1_pos-xhat),totalarea)/totalarea_normsq;
          N3 = tvmet::dot(tvmet::cross(x1_pos-xhat,x2_pos-xhat),totalarea)/totalarea_normsq;
          if(N1>0 && N1<1 && N2>0 && N2<1 && N3>0 && N3<1){ //checking the barycentric coordinates
            if(mindist>std::abs(dist)) {mindist=std::abs(dist); minid=list[j];}
          }
        }
        //if the minimum distance is less than cutoff, do actual computation
        if(mindist < _cutoff){
          nodes_mindist[_SlaveNodes[i]]=mindist;
          nodes_active[_SlaveNodes[i]]=true;
          elem_active[minid]=true;
          //three nodes of the closest triangle
          DefNode* x1=_defNodes[_MasterConnectivity[minid][0]]; const Vector3D & x1_pos = x1->point();
          DefNode* x2=_defNodes[_MasterConnectivity[minid][1]]; const Vector3D & x2_pos = x2->point();
          DefNode* x3=_defNodes[_MasterConnectivity[minid][2]]; const Vector3D & x3_pos = x3->point();
          //normal of the triangle
          Vector3D dir = tvmet::cross(x1_pos-x2_pos,x1_pos-x3_pos);
          dir = dir/tvmet::norm2(dir);
          double dist=dot((xstar_pos-x1_pos),dir);  //calculate dist again because it can be + or -, while mindist was kept absolute value
          //closest point
          Vector3D xhat; xhat = xstar_pos - dist*dir;
          Vector3D totalarea = tvmet::cross(x2_pos-x1_pos,x3_pos-x1_pos);
          double totalarea_normsq = tvmet::dot(totalarea,totalarea);
          //double totalarea=cal_area(x1_pos,x2_pos,x3_pos);
          double N1,N2,N3;
          N1 = tvmet::dot(tvmet::cross(x2_pos-xhat,x3_pos-xhat),totalarea)/totalarea_normsq;
          N2 = tvmet::dot(tvmet::cross(x3_pos-xhat,x1_pos-xhat),totalarea)/totalarea_normsq;
          N3 = tvmet::dot(tvmet::cross(x1_pos-xhat,x2_pos-xhat),totalarea)/totalarea_normsq;
          if(N1>0 && N1<1 && N2>0 && N2<1 && N3>0 && N3<1){ //checking the barycentric coordinates
          //if(f0) _energy += 1./2.*_k*(_cutoff - std::abs(dist))*(_cutoff - std::abs(dist));
          //double totalforce;
          double totalforce = (dist>0?1:-1)*_k*(1./std::abs(dist)/std::abs(dist)); 
          //if(dist<0)
          //  totalforce = -_k*(_cutoff - std::abs(dist));
          //else
          //  totalforce = _k*(_cutoff - std::abs(dist));
          //double totalforce = -_k*(_cutoff - std::abs(dist));
          //double totalforce = -2*_k/pow(std::abs(dist),3);
          //total force in the normal direction on the slave node
          xstar->addForce(0,totalforce*dir(0)); xstar->addForce(1,totalforce*dir(1)); xstar->addForce(2,totalforce*dir(2));
          //three components of the total force on the master nodes in the opposite direction
          x1->addForce(0,-totalforce*N1*dir(0)); x1->addForce(1,-totalforce*N1*dir(1)); x1->addForce(2,-totalforce*N1*dir(2));
          x2->addForce(0,-totalforce*N2*dir(0)); x2->addForce(1,-totalforce*N2*dir(1)); x2->addForce(2,-totalforce*N2*dir(2));
          x3->addForce(0,-totalforce*N3*dir(0)); x3->addForce(1,-totalforce*N3*dir(1)); x3->addForce(2,-totalforce*N3*dir(2));
          }
        }
      }

      
      std::string fileName = name + ".vtk";
      std::ofstream ofs(fileName.c_str());
      if (!ofs) {
        std::cout << "can not open output ("
	  	  << fileName
	  	  << ") file." << std::endl;
        exit(0);
      }
    
    ////////////////////////////////////////////////////////////////////
    //
    //    Node Section
    //
    ofs << "# vtk DataFile Version 2.0\n"
	<< "Test example" << std::endl
	<< "ASCII" << std::endl
	<< "DATASET POLYDATA" << std::endl
	<< "POINTS  " << _defNodes.size() << "  double" << std::endl;
    
    int nNodes=_defNodes.size();
    // output nodal postions
    for (int e=0; e<nNodes; e++){
      const Vector3D & nodalPos = _defNodes[e]->point();
      ofs << std::setprecision(16) 
	  << nodalPos(0) << "  "
	  << nodalPos(1) << "  "
	  << nodalPos(2) << std::endl;
    }

    /////////////////////////////////////////////////////////////////////
    //
    //    Element Section
    //
    int nShells = _MasterConnectivity.size();
    ofs << "POLYGONS  " << nShells << "  "
	<< 4*nShells << std::endl;
    for(int e=0; e<nShells; e++) {
      ofs << 3 << "  "
	  << std::setw(10) << _MasterConnectivity[e][0]
	  << std::setw(10) << _MasterConnectivity[e][1]
	  << std::setw(10) << _MasterConnectivity[e][2]
	  << std::endl;
    }
     	
    //////////////////////////////////////////////////////////////////////////
    //
    //  output data for each element, if it is an active contact or not
    //
    ofs << "CELL_DATA    " << nShells << std::endl;
    ofs << "SCALARS activecontact  int  1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for(int e=0; e<nShells; e++) {
      if(elem_active[e])  ofs << "1\n";
      else ofs << "0\n";
    }
    ofs << std::endl;

    ofs << endl << "POINT_DATA " << nNodes << std::endl;
    ofs << "SCALARS activecontact  int  1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for(int e=0; e<nNodes; e++) {
      if(nodes_active[e])  ofs << "1\n";
      else ofs << "0\n";
    }
    ofs << std::endl;

    ofs << "SCALARS mindist  double  1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for(int e=0; e<nNodes; e++) {
      ofs << nodes_mindist[e]<<"\n";
    }
    ofs << std::endl;
    
    ofs << endl << "VECTORS forces double" << endl;
    for (int e=0; e<nNodes; e++){
      const Vector3D & nodalForce = _defNodes[e]->force();
      ofs << std::setprecision(16) 
	  << nodalForce(0)
	  << '\t' <<nodalForce(1)
	  << '\t' <<nodalForce(2) << std::endl;
    }

    ofs.close();
    return;
    }

    //! return the viscosity proportionality factor
    double springK() const {return _k;}

    //! return the viscosity proportionality factor
    void setSpring(double k) {_k = k;}

    private:
    double cal_area(const Vector3D & x1, const Vector3D & x2, const Vector3D & x3){
      double area = 1./2.*tvmet::norm2( tvmet::cross(x2-x1,x3-x1));
      return area;
    } 
    //! spring constant, cutoff for potential and spacing for the search data structure 
    double _k, _cutoff, _spacing;
    //! active contact number
    int _active_contacts;
    //! Connectivities of master nodes
    std::vector<std::vector<int> > _MasterConnectivity;
    //! list of slave nodes
    std::vector<int> _SlaveNodes;
    //! deformation nodes
    DefNodeContainer _defNodes;
    //! Centroid positions of master nodes
    std::vector<std::vector<double> > _MasterCentroids;
  };

} // end namespace

#endif // __DeformableContact_h__
