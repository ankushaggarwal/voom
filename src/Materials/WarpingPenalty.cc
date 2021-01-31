// -*- C++ -*-
//----------------------------------------------------------------------
//
//                  Ankush Aggarwal, William S. Klug
//                University of California Los Angeles
//                   (C) 2012 All Rights Reserved
//
//----------------------------------------------------------------------
//

#include <iostream>
#include "VoomMath.h"
#include "WarpingPenalty.h"

//using namespace std;
namespace voom {

  // Construction/Destruction
  WarpingPenalty::WarpingPenalty(const WarpingPenalty &Input)
  {
    if (this == &Input) return;
    
    _init(Input._lambda,Input._X,Input._template_data,Input._target_data);
    //_lambda = Input._lambda;

    //_X = Input._X;
    //_TX = Input._TX;
    //_target_data=Input._target_data;
    //_template_data=Input._template_data;
  }

	
  void WarpingPenalty::_init(double lambda, Vector3D X, VoxelData * template_data, VoxelData * target_data)
  {
    _lambda = lambda;
    _X = X;
    _x = X;
    _target_data=target_data;
    _template_data=template_data;

    _TX=_template_data->GetValue(_X); //use the pointer to voxel data object here
    //std::cout<<"Debugging::MATERIAL::HERE"<<std::endl;

    // Augmented Lagrange multiplier
    _gamma=0.;
    factor=1.;

    // strain energy density
    _W = 0.0;
    //
    // first Piola stress tensor
    _force = 
      0.0, 0.0, 0.0;
  }


  // Operators


  // General methods
        
  void WarpingPenalty::updateState(bool fl0, bool fl1, bool fl2)
  {
    double Sx=_target_data->GetValue(_x);//use the pointer to voxel data object here
    Vector3D derivativeSx=_target_data->GetDerivative(_x);//use the pointer to voxel data object here
    
 
      //energy = lambda/2( T(X) - S(x))^2
    if(fl0)
      _W = _lambda/2.*(_TX-Sx)*(_TX-Sx) - factor*_gamma*(_TX-Sx);

    if(fl1)
      _force = -_lambda*(_TX-Sx)*derivativeSx + factor*_gamma*derivativeSx;
  
    if (fl2) {
    }
//    if(_W>0) std::cout<<"Sx="<<Sx
//                      <<"x="<<_x(0)<<" "<<_x(1)<<" "<<_x(2)<<std::endl;
    return;
  }

  void WarpingPenalty::updateLagrangeMultiplier()
  {
    double Sx=_target_data->GetValue(_x);
    _gamma -= _lambda*(_TX-Sx);
  }
    

  /*void WarpingPenalty::ConsistencyTest() 
  {
    std::cout << "checking consistency of force vector" << std::endl;
    updateState(false, true, false);
    Vector3D DfAna, DfNum;  DfAna = _force;                  // current 1st Piola stress
    const double eps = 1.0e-8 * max(_x);   // perturbation
    for(int i = 0; i < DfAna.rows(); i ++){
      for( int j = 0; j < DfAna.cols(); j++){
	_force(i,j) += eps;
	updateState(true, false, false);
	double W = _W;
	_x(i,j) -= 2*eps;
	updateState(true, false, false);
	W -= _W;
	_x(i,j) += eps;       // restore value
	DfNum(i,j) = W/2.0/eps;
      }
    }
    std::cout << "Analytical value of force:" << std::endl;
    std::cout << DfAna << std::endl;
    std::cout << "Numberical value of force:" << std::endl;
    std::cout << DfNum << std::endl;
    std::cout << "Difference between forces:" << std::endl;
    DfNum -= DfAna; DfNum /= max(DfAna);
    std::cout << DfNum << std::endl;
		
    std::cout << "\n\n\n\n" << std::endl;
    return;
  }*/

} // namespace voom
