#ifndef __PointBin__
#define __PointBin__

#include<vector>
#include<iostream>
#include<limits>
//namespace voom
//{
  
  class PointBin
  {
  public:
    //Constructor/destructors:
    
    //!Default constructor
    PointBin(){}
    
    //!Construct from a single grid size and a set of data
// wrong implementation    PointBin(double size, std::vector<std::vector<double> >& in_points){PointBin(size,size,size,in_points);}
    
    //!Construct from three grid sizes in three directions
    PointBin(double dx, double dy, double dz, std::vector<std::vector<double> >& in_points);

    
    //!Construct from three number of cells in three directions
    PointBin(int nx, int ny, int nz, std::vector<std::vector<double> >& in_points);

    
    //!Get back the list of int for a given point based on it's location
    std::vector<int>& GetList(std::vector<double>& point);

  private:
    //!Minimum and maximum positions
    double _minx, _maxx, _miny, _maxy, _minz, _maxz;
    
    //!Size of the grid in x, y and z direction
    int _nx, _ny, _nz;
    
    //!Spacing of the grid in x, y and z direction
    double _dx, _dy, _dz;
    
    //!std::vector of lists, where each element contains list of points
    std::vector< std::vector<int> > _grid;

    //!std::vector of list of points to be returned
    std::vector<int> _returnlist;
  };
//}
//#include "PointBin.cc"
#endif //PointBin


