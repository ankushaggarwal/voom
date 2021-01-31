#include<limits>
namespace voom{
  
  //Construction/Destruction
  
  VoxelData::VoxelData(double x0, double y0, double z0, double dx, double dy, double dz, int nx, int ny, int nz, const std::vector<double> & in_data)
  {
    //set the geometry parameters
    _dx=dx; _dy=dy; _dz=dz;
    _minx = x0-dx/2.; _miny = y0-dy/2.; _minz = z0-dz/2.;
    _nx = nx; _ny = ny; _nz = nz;
    _maxx = x0+nx*dx;
    _maxy = y0+ny*dy;
    _maxz = z0+nz*dz;

   //initialize the 3D data array
   _data.resize(nx);
   for(int i=0;i<nx;i++){
     _data[i].resize(ny);
     for(int j=0;j<ny;j++){
       _data[i][j].resize(nz);
       }
     }
     //store values from 1D array
     int iter=0;
     for(int k=0;k<nz;k++){
       for(int j=0;j<ny;j++){
	 for(int i=0;i<nx;i++)
	   _data[i][j][k]=in_data[iter++];
       }
     }
   }

  double VoxelData::GetValue(Vector3D & point){
    //calculate the box number for point
    int ix=int((point[0]-_minx)/_dx);
    int iy=int((point[1]-_miny)/_dy);
    int iz=int((point[2]-_minz)/_dz);
    //if out of range, print an error and return
    if(ix<0 || iy<0 || iz<0 || ix>_nx-1 || iy>_ny-1 || iz>_nz-1){
      std::cout<<"VoxelData::GetValue::Error: The point is out of range"<<std::endl;
      std::cout<<point[0]<<" "<<point[1]<<" "<<point[2]<<std::endl;
      exit(0);
      return 0;
    }
    //%%%%%%%%%piecewise constant%%%%%%%%%
    //return _data[ix][iy][iz];
    //for 1D array it should be _data[(ix*_ny+iy)*_nz+iz]
    
    //%%%%%%%%%linear%%%%%%%%
    double center_value=_data[ix][iy][iz];
    Vector3D deriv=GetDerivative(point);
    double deltax=point[0]-(ix+0.5)*_dx - _minx; 
    double deltay=point[1]-(iy+0.5)*_dy - _miny; 
    double deltaz=point[2]-(iz+0.5)*_dz - _minz;
    return center_value+(deltax*deriv[0]+deltay*deriv[1]+deltaz*deriv[2]); 
  }
  
  const Vector3D & VoxelData::GetDerivative(Vector3D & point){
    //_derivative.resize(3);
    //calculate the box number for point
    int ix=int((point[0]-_minx)/_dx);
    int iy=int((point[1]-_miny)/_dy);
    int iz=int((point[2]-_minz)/_dz);
    _derivative[0]=_derivative[1]=_derivative[2]=0.;
    //if out of range, print an error and return
    if(ix<1 || iy<1 || iz<1 || ix>_nx-2 || iy>_ny-2 || iz>_nz-2){
      std::cout<<"VoxelData::GetDerivative::Error: The point is out of range"<<std::endl;
      std::cout<<point[0]<<" "<<point[1]<<" "<<point[2]<<std::endl;
      return _derivative;
    }

    _derivative[0]=(_data[ix+1][iy][iz] - _data[ix-1][iy][iz])/2./_dx;
    _derivative[1]=(_data[ix][iy+1][iz] - _data[ix][iy-1][iz])/2./_dy;
    _derivative[2]=(_data[ix][iy][iz+1] - _data[ix][iy][iz-1])/2./_dz;

    return _derivative;
  }
  
}
