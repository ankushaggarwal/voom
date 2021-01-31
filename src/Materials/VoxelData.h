#ifndef __VoxelData__
#define __VoxelData__

#include<vector>
#include<iostream>
#include<limits>
namespace voom
{
  
  class VoxelData
  {
  public:
    //Constructor/destructors:
    
    //!Default constructor
    VoxelData(){}
    
    //!Construct from coordinates of "zeroth" voxel and its size three directions, number of voxels in three directions
    //!and a 1D vector containing data in the series x -> y -> z (same as vtk structured points format)
    inline VoxelData(double x0, double y0, double z0, double dx, double dy, double dz, int nx, int ny, int nz, const std::vector<double> & in_data);

    //!Get back data and its derivative at a particular point
    inline double GetValue(Vector3D & point);

    inline const Vector3D & GetDerivative(Vector3D & point);
  private:
    //!Minimum and maximum positions
    double _minx, _maxx, _miny, _maxy, _minz, _maxz;
    
    //!Size of the grid in x, y and z direction
    int _nx, _ny, _nz;
    
    //!Spacing of the grid in x, y and z direction
    double _dx, _dy, _dz;
    
    //!3D matrix of data
    std::vector< std::vector< std::vector<double> > > _data;

    //!Derivative of data returned
    Vector3D _derivative;
  };
}
#include "VoxelData.cc"
#endif //VoxelData


