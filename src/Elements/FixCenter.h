// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         Ankush Aggarwal
//                University of California Los Angeles
//                   (C) 2011 All Rights Reserved
//
//----------------------------------------------------------------------
// 

#if !defined(__FixCenter_h__)
#define __FixCenter_h__

#include <vector>
#include <cstdio>
#include <ctime>

#include "Element.h"

namespace voom {
  //! Applies a penalty to prevent the motion of the center of given set of nodes
  // E = -p (centerx-centerX).(centerx-centerX) (where centerx=vector sum of all node positions in x(deformed) and X(reference state)

  class FixCenter :
    public Element {
   
    public:
    typedef DeformationNode<3> Node_t;

    typedef std::vector<Node_t*> NodeContainer;
    typedef NodeContainer::iterator NodeIterator;
    typedef NodeContainer::const_iterator ConstNodeIterator;

    //! Construct from a set of nodes
    FixCenter( const NodeContainer & nodes, const double k)
      {
        _nodes = nodes;
        _penalty = k;
        RefCenter();
      }

    //calculate the current/deformed center
    void center(){
      NodeIterator n=_nodes.begin();
      _center.resize((*n)->dof(),0.);
      for(; n!=_nodes.end(); n++){
	for(int i=0; i<(*n)->dof(); i++) 
          _center[i] += (*n)->getPoint(i);
      }
      for(int i=0; i<_center.size(); i++) 
	_center[i] /= _nodes.size();
    }

    //calculate the reference center
    void RefCenter(){
      NodeIterator n=_nodes.begin();
      _RefCenter.resize((*n)->dof(),0.);
      for(; n!=_nodes.end(); n++){
	for(int i=0; i<(*n)->dof(); i++) 
          _RefCenter[i] += (*n)->getPosition(i);
      }
      for(int i=0; i<_RefCenter.size(); i++) 
	_RefCenter[i] /= _nodes.size();
    }

    //! Compute energy, force and stiffness
    void compute(bool f0, bool f1, bool f2) {
      center();
      int dim=_RefCenter.size();
      if(f0){
        double normD=0.; //norm of the center displacement
        for(int i=0; i<dim;i++)
          normD += (_center[i]-_RefCenter[i])*(_center[i]-_RefCenter[i]);
        _energy=_penalty*normD;
      }

      if(f1){
        //calculate the force
        std::vector<double> force;
        force.resize(dim,0.);
        for(int i=0;i<dim;i++)
          force[i]=2*_penalty*(_center[i]-_RefCenter[i]);
        //apply the same force on all nodes
        for(NodeIterator n=_nodes.begin(); n!=_nodes.end(); n++){
	  for(int i=0; i<dim; i++)
	    (*n)->addForce(i, force[i]);
        }
      }
   }
   
   private:

   //! Node container
   NodeContainer _nodes;

   //! penalty
   double _penalty;

   //! center of the nodes in deformed and reference configuration 
   std::vector<double> _center,_RefCenter;
   };

} // end namespace

#endif // __FixCenter_h__
