// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         Ankush Aggarwal
//                University of California Los Angeles
//                   (C) 2012 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file WarpingElement.h

  \brief Warping is a concrete class derived from element, implementing
  a C0 finite element for finite deformation membrane structures image registration using Hyperelastic Warping.

*/

#if !defined(__WarpingElement_h__)
#define __WarpingElement_h__

#include <vector>
#include <cstdio>
#include <ctime>

#include "StandardElement.h"
#include "Node.h"
#include "ShellGeometry.h"
#include "VoomMath.h"

namespace voom
{

  /*!  Concrete class for a C0 finite deformation membrane finite element.
  */
  template<class Quadrature_t, class Material_t, class Shape_t>
  class WarpingElement 
    : public StandardElement< DeformationNode<3>, Quadrature_t,
			      Material_t, Shape_t > 
  {
    
  public:

    typedef DeformationNode<3> Node_t;

    // WSK & MMG don't think this typedef should be here!!!
//     typedef TriangleQuadrature Quadrature_t;

    typedef StandardElement<Node_t, Quadrature_t, Material_t, Shape_t> Base;

    typedef typename Base::NodeContainer 		NodeContainer;
    typedef typename Base::NodeIterator 		NodeIterator;
    typedef typename Base::ConstNodeIterator 		ConstNodeIterator;
    typedef typename Base::QuadPointIterator 		QuadPointIterator;
    typedef typename Base::ConstQuadPointIterator 	ConstQuadPointIterator;

    using Base::_nodes;
    using Base::_quadPoints;
    using Base::_baseNodes;
    using Base::_energy;

    //! virtual destructor
    virtual ~WarpingElement() {;}

    //! Constructor
    WarpingElement( const Quadrature_t & quad,
		Material_t & mat,
		const NodeContainer & nodes)
    {
      //! initialize NodeContainer
      unsigned nNodes = nodes.size();
      _nodes = nodes;

      for(ConstNodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) 
	_baseNodes.push_back(*n);
      
      //! initialize materials and shape functions
      _quadPoints.clear();
      for(typename Quadrature_t::ConstPointIterator p=quad.begin(); 
	  p!=quad.end(); p++) {
	Shape_t shp( p->coords );
      
	_quadPoints.push_back( typename Base::QuadPointStruct(p->weight, mat, shp) ); 
	if( shp.functions().size() != nNodes ) {
	  std::cout << "Number of nodes: " << nNodes << std::endl
		    << "Number of functions: " << shp.functions().size()
		    << std::endl
		    << "These should be equal." << std::endl;
	  exit(0);
	}
      }
      
////////////////////////////////////////////////////////// ADDED THIS PART:::CHECK FOR CORRECTNESS
    // loop for every quadrature point and set the correct reference position in material
    for(QuadPointIterator p=_quadPoints.begin(); 
	p!=_quadPoints.end(); p++){
      Vector3D X(0.0);

      const typename Shape_t::FunctionContainer & N 
	= p->shape.functions();

      for (int b = 0; b < _nodes.size(); b++){
	const DeformationNode<3>::Point & xb = _nodes[b]->position();
	X 	       +=   N[b]     * xb;
      }
      p->material.setRefPosition(X);
    }
////////////////////////////////////////////////////////// 


      //! allocate memory for mechanics variables
      _internalForce.resize( nNodes );

      updateRefConfiguration(); 

      //! initialize values = 0 by default
      //compute(false, false, false); 
      compute(true, true, false); 
    }

  public:

    //! Set the template voxel value to a constant
    void setTemplateVoxel(double v)
    {
      for(QuadPointIterator p=_quadPoints.begin(); 
	p!=_quadPoints.end(); p++){
        p->material.setTemplateVoxel(v);
      }
    }

     void setMaterialLambda(double lambda){
       for(QuadPointIterator p=_quadPoints.begin(); p!=_quadPoints.end(); p++){
 	Material_t& material = p->material;
 	material.setLambda(lambda);
 	//material.setCytoSpring(0.0,0.0,kSpring);
       }
     }
   
     void updateLagrangeMultiplier(){
       for(QuadPointIterator p=_quadPoints.begin(); p!=_quadPoints.end(); p++){
 	Material_t& material = p->material;
 	material.updateLagrangeMultiplier();
 	//material.setCytoSpring(0.0,0.0,kSpring);
       }
     }
    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    virtual void compute(bool f0, bool f1, bool f2);
    
    double strainEnergy() const { return _strainEnergy; }

    //! compute positions
    Vector3D computePosition(const double s1, const double s2);

    //! check new positions
    void checkPositions();
		
    //! recompute the reference geometry and send it to materials
    void updateRefConfiguration();

    //! set the reference configuration in xy plane explicitly
    void SetRefConfiguration(double edgelen);

    //! set the edgelength of reference configuration triangle
//    static void setRefLength(double e){_edgelen=e;}

    //
    //   data
    //
  private:
    //! forces on the element
    blitz::Array< Vector3D, 1> _internalForce;
		
    //! strain energy
    double _strainEnergy;

   //! Reference area
   double _RefArea;

  };  // end of class
} // namespace voom

#include "WarpingElement.icc"

#endif // __WarpingElement_h__
