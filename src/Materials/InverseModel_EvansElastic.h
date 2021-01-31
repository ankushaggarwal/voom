// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         Ankush Aggarwal
//                University of California Los Angeles
//                   (C) 2012 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file InverseModel_EvansElastic.h

  \brief Interface for Inverse model with Evans material law for an elastic membrane 

*/

#if !defined(__InverseModel_EvansElastic_h__)
#define __InverseModel_EvansElastic_h__

#include<blitz/array.h>
 
 
#include<vector>
#include "voom.h"
#include "SCElastic.h"
#include "VoomMath.h"
#include "Node.h"

namespace voom
{

  /*!  Class implementing the Evans 
    model for an elastic lipid bilayer & cytoskeleton.
  */
  
  class InverseModel_EvansElastic : public SCElastic {

  protected:
    //! shear modulus
    ScalarNode* _mu;
    //double _mu; 

    //! stretching modulus
    ScalarNode* _kS;
    //double _kS;

    //! first invariant of Rt. Cauchy-Green deformation tensor
    double _trC;

    //! second invariant of Rt. Cauchy-Green deformation tensor
    double _J;

    //! stretching energy
    double _Ws;

    //! addition "stress resultants" due to material DOFs
    tvmet::Vector< Vector3D, 3 > _dnmu;  // stress resultants
    tvmet::Vector< Vector3D, 3 > _dnkS;  // stress resultants

  public:

    InverseModel_EvansElastic(const double kC, const double kG, const double C0, ScalarNode* mu, ScalarNode* kS ) 
      : SCElastic(kC, kG, C0), _mu(mu), _kS(kS), _J(0.0), _trC(0.0), _Ws(0.0)
    {};

    virtual void updateState(bool f0, bool f1, bool f2 );

    double shearModulus() const {return _mu->getPoint(0);}
    double stretchingModulus() const {return _kS->getPoint(0);}
    double J() const {return _J;}
    double trC() const {return _trC;}
    const Tensor3D DefGradient();
    const Tensor3D Cauchystress();
    virtual const tvmet::Vector< Vector3D, 3 >& mu_stressResultants() const 
    { return _dnmu; }
    virtual const tvmet::Vector< Vector3D, 3 >& kS_stressResultants() const 
    { return _dnkS; }

    void setShearModulus(double new_mu) {
      _mu->setPoint(0,new_mu);
    }

    void setStretchModulus(double new_kS) {
      _kS->setPoint(0,new_kS);
    }

    virtual double bendingEnergy() const {return _W - _Ws;}
    virtual double stretchingEnergy() const {return _Ws;}
  };
  
} //namespace voom

#endif //  !defined(__InverseModel_EvansElastic_h__)
