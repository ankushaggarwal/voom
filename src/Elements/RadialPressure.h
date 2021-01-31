// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         Ankush Aggarwal
//                University of California Los Angeles
//                   (C) 2011 All Rights Reserved
//
//----------------------------------------------------------------------
// 

#if !defined(__RadialPressure_h__)
#define __RadialPressure_h__

#include <vector>
#include <cstdio>
#include <ctime>

#include "Element.h"

namespace voom {
  //! Applies a pressure force at given set of nodes (in the direction of the nodes multiplied by the nodal projection area)
  // E = -p vecN.(vecx-vecX) (where vecN=vecX/|X|

  class RadialPressure :
    public Element {
   
    public:
    typedef DeformationNode<3> Node_t;

    typedef std::vector<Node_t*> NodeContainer;
    typedef NodeContainer::iterator NodeIterator;
    typedef NodeContainer::const_iterator ConstNodeIterator;

    //! Construct from a set of nodes
    RadialPressure( const NodeContainer & nodes, const std::vector<double> & node_area)
      {
        _nodes = nodes;
        _NodalArea = node_area;
        _pressure = 0;
        center();
        RefMeanRadius();
      }

    //set the pressure
    void setP(double p){_pressure=p;}

    //calculate the center
    void center(){
      NodeIterator n=_nodes.begin();
      _center.resize((*n)->dof(),0.);
      for(; n!=_nodes.end(); n++){
	for(int i=0; i<(*n)->dof(); i++) 
          _center[i] += (*n)->getPosition(i);
      }
      for(int i=0; i<_center.size(); i++) 
	_center[i] /= _nodes.size();
    }

    //calculate the mean radius in reference configuration
    double RefMeanRadius() {
      double meanR=0;
      center();
      for(NodeIterator n=_nodes.begin(); n!=_nodes.end(); n++){
        double normx=0;
	for(int i=0; i<(*n)->dof(); i++) 
          normx += ((*n)->getPosition(i)-_center[i])*((*n)->getPosition(i)-_center[i]);
        normx = sqrt(normx);
        meanR += normx;
      }
      meanR /= _nodes.size();
      _RefMeanR=meanR;
      return meanR;
      }

    //calculate the mean radius
    double MeanRadius() {
      double meanR=0;
      center();
      for(NodeIterator n=_nodes.begin(); n!=_nodes.end(); n++){
        double normx=0;
	for(int i=0; i<(*n)->dof(); i++) 
          normx += ((*n)->getPoint(i)-_center[i])*((*n)->getPoint(i)-_center[i]);
        normx = sqrt(normx);
        meanR += normx;
      }
      meanR /= _nodes.size();
      return meanR;
      }

    //! Compute energy, force and stiffness
    void compute(bool f0, bool f1, bool f2) {
      //double curR=MeanRadius();
      //double AreaRatio=curR*curR/_RefMeanR/_RefMeanR; //ratio of surface area in current to reference config
      double AreaRatio=1.;

      if(f0){
        _energy=0;
        int count=0;
        for(NodeIterator n=_nodes.begin(); n!=_nodes.end(); n++){
          //normX
          double normX=0;
	  for(int i=0; i<(*n)->dof(); i++) 
            normX += ((*n)->getPosition(i)-_center[i])*((*n)->getPosition(i)-_center[i]);
          normX = sqrt(normX);
          //X.(x-X)
          double xdot=0;
	  for(int i=0; i<(*n)->dof(); i++) 
            xdot += ((*n)->getPosition(i)-_center[i])*((*n)->getPoint(i)-(*n)->getPosition(i));
          _energy += -_pressure*AreaRatio*xdot/normX*_NodalArea[count++];
        }
      }

      if(f1){
        int count=0;
        for(NodeIterator n=_nodes.begin(); n!=_nodes.end(); n++){
          //normX
          double normX=0;
	  for(int i=0; i<(*n)->dof(); i++) 
            normX += ((*n)->getPosition(i)-_center[i])*((*n)->getPosition(i)-_center[i]);
          normX = sqrt(normX);
	  for(int i=0; i<(*n)->dof(); i++) {
	    double xi = (*n)->getPosition(i)-_center[i];
	    (*n)->addForce(i, -_pressure*AreaRatio*xi/normX*_NodalArea[count]);
          }
          count++;
        }
      }
   }
   
   private:

   //! Node container
   NodeContainer _nodes;

   //! pressure
   double _pressure;

   //! reference mean radius (to make the pressure act in current config)
   double _RefMeanR;

   //! nodal area (to be multiplied by the pressure to make it a force)
   std::vector<double> _NodalArea;

   //! center of the nodes set to apply radial force 
   std::vector<double> _center;
   };

} // end namespace

#endif // __RadialPressure_h__
