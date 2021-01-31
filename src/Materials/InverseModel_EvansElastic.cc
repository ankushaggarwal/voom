// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         Ankush Aggarwal
//                University of California Los Angeles
//                   (C) 2012 All Rights Reserved
//
//----------------------------------------------------------------------
//

/*! 
  \file InverseModel_EvansElastic.cc

  \brief Implementation of Evans model for an elastic cytoskeleton & bilayer.

*/
#include <iostream>
#include "InverseModel_EvansElastic.h"

namespace voom 
{
  void InverseModel_EvansElastic::updateState(bool f0, bool f1, bool f2)
  {
    // first call updateState from the SCElastic to get bending
    // contributions
    SCElastic::updateState(f0, f1, f2);

    // OK, now add the stretching response
    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;
                
    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();

    const BasisVectors & refBasis = _referenceGeometry.a();
    const BasisVectors & refDual = _referenceGeometry.aDual();
    
    Tensor3D F(0.0);
    for(int alpha=0; alpha<2; alpha++){
      for(int i=0; i<3; i++) {
	for(int J=0; J<3; J++) {
	  F(i,J) += basis(alpha)(i)*refDual(alpha)(J);	  
	}
      }
    }

    Tensor3D C(0.0);
    C = trans(F)*F;

    _trC = C(0,0)+C(1,1)+C(2,2);
    double trCSquare = 0.0;
    for(int i=0; i<3; i++) {
      for(int k=0; k<3; k++) {
	trCSquare += C(i,k)*C(k,i);
      }
    }
      
    _J = sqrt((_trC*_trC - trCSquare)/2.0); 
    
    double mu=_mu->getPoint(0);
    double kS=_kS->getPoint(0);

    Tensor3D P(0.0);
    P = (mu/_J)*F 
      + (1.0/_J)*(kS*(_J-1) - 0.5*mu*_trC/sqr(_J))*(_trC*F - F*C);

    Tensor3D dPmu(0.0);
    dPmu = (1./_J)*F 
      + (1.0/_J)*(0.-0.5*_trC/sqr(_J))*(_trC*F - F*C);
    
    Tensor3D dPkS(0.0);
    dPkS = (1.0/_J)*(_J-1)*(_trC*F - F*C);

    if( f0 ) {
      // compute strain energy
      _Ws = (kS*sqr(_J-1.0)/2.0 + 0.5*mu*(_trC/_J-2.0))/_J;
      _W += _Ws;
    }
    
    if( f1 ) {
      //
      // -------------- stress resultants --------------
      for(int alpha=0; alpha<2; alpha++){ 
	_n(alpha) += P*refDual(alpha)/_J;
	_dnmu(alpha) += dPmu*refDual(alpha)/_J;
	_dnkS(alpha) += dPkS*refDual(alpha)/_J;
      }
      _n(2) = 0.0, 0.0, 0.0;
      _dnmu(2) = 0.0, 0.0, 0.0;
      _dnkS(2) = 0.0, 0.0, 0.0;
    }
    
    return;
  }

  const Tensor3D InverseModel_EvansElastic::DefGradient(){
    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;
                
    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();

    const BasisVectors & refBasis = _referenceGeometry.a();
    const BasisVectors & refDual = _referenceGeometry.aDual();
    
    Tensor3D F(0.0);
    for(int alpha=0; alpha<2; alpha++){
      for(int i=0; i<3; i++) {
	for(int J=0; J<3; J++) {
	  F(i,J) += basis(alpha)(i)*refDual(alpha)(J);	  
	}
      }
    }
   return F;
  }
  
  const Tensor3D InverseModel_EvansElastic::Cauchystress(){
    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;
                
    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();

    const BasisVectors & refBasis = _referenceGeometry.a();
    const BasisVectors & refDual = _referenceGeometry.aDual();
    
    Tensor3D F(0.0);
    for(int alpha=0; alpha<2; alpha++){
      for(int i=0; i<3; i++) {
	for(int J=0; J<3; J++) {
	  F(i,J) += basis(alpha)(i)*refDual(alpha)(J);	  
	}
      }
    }
    Tensor3D C(0.0);
    C = trans(F)*F;

    _trC = C(0,0)+C(1,1)+C(2,2);
    double trCSquare = 0.0;
    for(int i=0; i<3; i++) {
      for(int k=0; k<3; k++) {
	trCSquare += C(i,k)*C(k,i);
      }
    }
      
    _J = sqrt((_trC*_trC - trCSquare)/2.0); 
    double mu=_mu->getPoint(0);
    double kS=_kS->getPoint(0);
    Tensor3D P(0.0);
    P = (mu/_J)*F 
      + (1.0/_J)*(kS*(_J-1) - 0.5*mu*_trC/sqr(_J))*(_trC*F - F*C);
    Tensor3D cauchy(0.0);
    cauchy=1/_J*P*trans(F);
    return cauchy;
  }
}
