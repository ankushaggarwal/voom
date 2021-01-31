// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          Ankush Aggarwal
//                University of California Los Angeles
//                   (C) 2012 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__WarpingBody_h__)
#define __WarpingBody_h__

#include <blitz/array.h>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "Body.h"
#include "WarpingElement.h"
#include "voom.h"
#include "Node.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{

  /*! 
    WarpingBody is a concrete class derived from Body, implementing
    the concept of a penalty function to align a template image to a target image using finite elements.
    Method is know as Hyperelastic Warping for image registration
    Ref: Veress, Phatak and Weiss
  */
  template < class Quadrature_t, class Material_t, class Shape_t >
  class WarpingBody : public Body
  {
  public:
    
    // typedefs
    typedef 
    WarpingElement<Quadrature_t,Material_t,Shape_t> WarpingElement_t;
    typedef typename 
    WarpingElement_t::Node_t WarpingNode_t;
    typedef typename 
    std::vector< WarpingNode_t* > WarpingNodeContainer;
    typedef typename 
    WarpingNodeContainer::iterator WarpingNodeIterator;
    typedef typename 
    WarpingNodeContainer::const_iterator ConstWarpingNodeIterator;

    typedef 
    std::vector< WarpingElement_t* > WarpingElementContainer;
    typedef typename 
    WarpingElementContainer::iterator WarpingElementIterator;
    typedef typename 
    WarpingElementContainer::const_iterator ConstWarpingElementIterator;
		
    typedef 
    std::vector< Constraint* > ConstraintContainer;
    typedef typename 
    ConstraintContainer::iterator ConstraintIterator;
    typedef typename 
    ConstraintContainer::const_iterator ConstConstraintIterator;
		
    typedef 
    std::vector<int> ElementConnectivity;
    typedef 
    std::vector<ElementConnectivity> ConnectivityContainer;
    typedef 
    ConnectivityContainer::iterator ConnectivityIterator;
    typedef 
    ConnectivityContainer::const_iterator ConstConnectivityIterator;


    //! Default Constructor
    WarpingBody() {;}

    //! Construct from stuff
    WarpingBody(
	Material_t material,
	ConnectivityContainer & connectivities,
	const NodeContainer & nodes,
	const unsigned quadOrder=1) ;
    
    //! virtual destructor
    virtual ~WarpingBody() {;}
    
     void setMaterialLambda(double lambda) {
       for(int si=0; si<_membranes.size(); si++) {
 	WarpingElement_t* s=_membranes[si];
 	s->setMaterialLambda(lambda);     	
       }
     }

     void updateLagrangeMultiplier() {
       for(int si=0; si<_membranes.size(); si++) {
 	WarpingElement_t* s=_membranes[si];
 	s->updateLagrangeMultiplier();     	
       }
     }

    //! Do mechanics on Body
    void compute( bool f0, bool f1, bool f2 );

    //! return the elements container
    const WarpingElementContainer & elements() {return _membranes;}

    //! loop over elements and add up the strain energy of each one
    double totalStrainEnergy() const {
      double totalStrainEnergy = 0.0;
      for(ConstWarpingElementIterator e=_membranes.begin(); 
	  e!=_membranes.end(); e++) {
	totalStrainEnergy += (*e)->strainEnergy();
      }	
      return totalStrainEnergy;
    }

    //! set the template voxel values
    void setTemplateVoxel(double v)
    {
      for(ConstWarpingElementIterator e=_membranes.begin(); 
	  e!=_membranes.end(); e++) {
        (*e)->setTemplateVoxel(v);
      }
    }
    
    //! loop over elements and add up the work of each one
    double work() const {
      double work = 0.0;
      for(ConstWarpingElementIterator e=_membranes.begin(); 
	  e!=_membranes.end(); e++) {
	work += -(*e)->work()+(*e)->constraintEnergy();
      }	
      return work;
    }

    virtual void pushBack( Element* e ) { 
      _elements.push_back(e); 
      _active.push_back(true);}

    void resetReference() {
      for(WarpingNodeIterator n = _membraneNodes.begin(); 
	  n != _membraneNodes.end(); n++) {
	(*n)->resetPosition();
      }

      for(WarpingElementIterator e=_membranes.begin(); 
	  e!=_membranes.end(); e++) {
	(*e)->updateRefConfiguration();
      }	
      return;
    }

    //! Query an element's activity status
    bool active(int e) { return _active[e]; }

    //! Mark an element as active so it will be computed
    void activate(int e) { _active[e] = true; }

    //! Mark an element as inactive so it will not be computed
    void deactivate(int e) { _active[e] = false; }

    //! Print paraview file 
    void printParaview(const std::string name) const; 

  private:

    //
    // data
    //
    
    //! Nodes
    WarpingNodeContainer _membraneNodes;

    //! Elements
    WarpingElementContainer 	_membranes;		

    //! vector of boolean flags to activate/deactivate elements
    std::vector<bool> _active;

#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
  };  
} // namespace voom

#include "WarpingBody.icc"

#endif // __WarpingBody_h__
