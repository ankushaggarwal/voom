// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         Ankush Aggarwal
//                University of California Los Angeles
//                   (C) 2012 All Rights Reserved
//
//----------------------------------------------------------------------
// 

#if !defined(__NormalPressure_h__)
#define __NormalPressure_h__

#include <vector>
#include <cstdio>
#include <ctime>

#include "Element.h"
#include "NodeBin.h"

namespace voom {

  //! A normal pressure on a shell surface
  /*! NormalPressure is an element which applies a 
      normal pressure on the surface of a thin shell
  */
  class NormalPressure :
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

    //! Construct from deformation nodes and connectivity of master nodes
    NormalPressure( const DefNodeContainer & Nodes, std::vector<std::vector<int> > Connectivity) 
      {
	_defNodes = Nodes;
	_Connectivity = Connectivity;
	_energy = 0.0;
        std::vector<double> x(3,0.0);
        _normals.clear();
        _normals.resize(_Connectivity.size(),x);
        _pressure=0.;
        compute(true,true,false);
      }

    //! Compute the energy, force, and stiffness
    void compute(bool f0, bool f1, bool f2) {
      if(f0) _energy = 0.0;
      for(int i=0;i<_Connectivity.size();i++){
        DefNode* x1=_defNodes[_Connectivity[i][0]]; const Vector3D & x1_pos = x1->point();
        DefNode* x2=_defNodes[_Connectivity[i][1]]; const Vector3D & x2_pos = x2->point();
        DefNode* x3=_defNodes[_Connectivity[i][2]]; const Vector3D & x3_pos = x3->point();
        //normal and area of the triangle
        Vector3D normal=tvmet::cross(x2_pos-x1_pos,x3_pos-x1_pos);
        double area=tvmet::norm2(normal);
        normal=normal/area;
        double totalforce=area*_pressure;
        if(f0) _energy += totalforce* (tvmet::dot(normal,x1->point()-x1->position()) + tvmet::dot(normal,x2->point()-x2->position()) + tvmet::dot(normal,x3->point()-x3->position()));
        if(f1){
          x1->addForce(0,totalforce*normal(0)); x1->addForce(1,totalforce*normal(1)); x1->addForce(2,totalforce*normal(2));
          x2->addForce(0,totalforce*normal(0)); x2->addForce(1,totalforce*normal(1)); x2->addForce(2,totalforce*normal(2));
          x3->addForce(0,totalforce*normal(0)); x3->addForce(1,totalforce*normal(1)); x3->addForce(2,totalforce*normal(2));
        }
      }
    }

    //! return the current pressure
    double pressure() const {return _pressure;}

    //! return the viscosity proportionality factor
    void setPressure(double p) {_pressure = p;}

    private:
    double cal_area(const Vector3D & x1, const Vector3D & x2, const Vector3D & x3){
      double area = 1./2.*tvmet::norm2( tvmet::cross(x2-x1,x3-x1));
      return area;
    } 
    //! spring constant, cutoff for potential and spacing for the search data structure 
    double _pressure;
    //! Connectivities of master nodes
    std::vector<std::vector<int> > _Connectivity;
    //! deformation nodes
    DefNodeContainer _defNodes;
    //! Centroid positions of master nodes
    std::vector<std::vector<double> > _normals;
  };

} // end namespace

#endif // __NormalPressure_h__
