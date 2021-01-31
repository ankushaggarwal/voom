#include "PointBin.h"
#include<limits>

//namespace voom{
  
  //Construction/Destruction
  
  PointBin::PointBin(double dx, double dy, double dz, std::vector<std::vector<double> >& in_points)
  {
    
    _dx=dx; _dy=dy; _dz=dz;
    //find minimum and maximum of in_points
    _minx = _miny = _minz = std::numeric_limits<double>::max();
    _maxx = _maxy = _maxz = std::numeric_limits<double>::min();
    for(int i=0;i<in_points.size();i++){
      
      if(_minx > in_points[i][0]) _minx = in_points[i][0];
      if(_miny > in_points[i][1]) _miny = in_points[i][1];
      if(_minz > in_points[i][2]) _minz = in_points[i][2];
      
      if(_maxx < in_points[i][0]) _maxx = in_points[i][0];
      if(_maxy < in_points[i][1]) _maxy = in_points[i][1];
      if(_maxz < in_points[i][2]) _maxz = in_points[i][2];
      
    }
    //increase the size by dx,dy and dz at every corner
    _minx=_minx-_dx; _miny=_miny-_dy; _minz=_minz-_dz;
    _maxx=_maxx+_dx; _maxy=_maxy+_dy; _maxz=_maxz+_dz;
    
    //calculate number of grids
    _nx = int((_maxx-_minx)/dx);
    _ny = int((_maxy-_miny)/dy);
    _nz = int((_maxz-_minz)/dz);
    
    //reserve the size
    _grid.resize(_nx*_ny*_nz);
    
    for(int i=0;i<in_points.size();i++){
      int ix = int((in_points[i][0] - _minx)/_dx);
      int iy = int((in_points[i][1] - _miny)/_dy);
      int iz = int((in_points[i][2] - _minz)/_dz);
      _grid[ix+iy*_nx+iz*_nx*_ny].push_back(i);
    }
  }

  PointBin::PointBin(int nx, int ny, int nz, std::vector<std::vector<double> >& in_points)
  {
    
    _nx=nx; _ny=ny; _nz=nz;
    //find minimum and maximum of in_points
    _minx = _miny = _minz = std::numeric_limits<double>::max();
    _maxx = _maxy = _maxz = std::numeric_limits<double>::min();
    for(int i=0;i<in_points.size();i++){
      
      if(_minx > in_points[i][0]) _minx = in_points[i][0];
      if(_miny > in_points[i][1]) _miny = in_points[i][1];
      if(_minz > in_points[i][2]) _minz = in_points[i][2];
      
      if(_maxx < in_points[i][0]) _maxx = in_points[i][0];
      if(_maxy < in_points[i][1]) _maxy = in_points[i][1];
      if(_maxz < in_points[i][2]) _maxz = in_points[i][2];
      
    }
    
    //calculate spacing of grids
    _dx = (_maxx-_minx)/_nx;
    _dy = (_maxy-_miny)/_ny;
    _dz = (_maxz-_minz)/_nz;
    
    //add 2 to nx, ny and nz
    _nx+=2; _ny+=2; _nz+=2; 
    //and increase the size by dx,dy and dz at every corner
    _minx=_minx-_dx; _miny=_miny-_dy; _minz=_minz-_dz;
    _maxx=_maxx+_dx; _maxy=_maxy+_dy; _maxz=_maxz+_dz;

    //reserve the size
    _grid.resize(_nx*_ny*_nz);
    
    for(int i=0;i<in_points.size();i++){
      int ix = int((in_points[i][0] - _minx)/_dx);
      int iy = int((in_points[i][1] - _miny)/_dy);
      int iz = int((in_points[i][2] - _minz)/_dz);
      _grid[ix+iy*_nx+iz*_nx*_ny].push_back(i);
    }
  }

  std::vector<int>& PointBin::GetList(std::vector<double>& point){
    //calculate the box number for point
    int ix=int((point[0]-_minx)/_dx);
    int iy=int((point[1]-_miny)/_dy);
    int iz=int((point[2]-_minz)/_dz);
    _returnlist.clear();
    bool flag=false;
    for(int ii=std::max(0,ix-1);ii<=std::min(ix+1,_nx-1);ii++)
      for(int jj=std::max(0,iy-1);jj<=std::min(iy+1,_ny-1);jj++)
	for(int kk=std::max(0,iz-1);kk<=std::min(iz+1,_nz-1);kk++){
          flag=true;
	  for(int mm=0;mm<_grid[ii+jj*_nx+kk*_nx*_ny].size();mm++)
	    _returnlist.push_back(_grid[ii+jj*_nx+kk*_nx*_ny][mm]);
	}
    if(flag==false) std::cout<<"Point does not lie in any of the grid cells"<<std::endl;
    return _returnlist;
  }
  
//}

