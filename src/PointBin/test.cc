#include<iostream>
#include<vector>
#include"PointBin.h"
//#include"PointBin.cc"
//using namespace voom;

int main(){

std::vector<std::vector<double> > points;

for(int i=0;i<20;i++){
double x=i/2.;
for(int j=0;j<20;j++){
double y=j/2.;
for(int k=0;k<20;k++){
double z=k/2.;
std::vector<double> p(3);
p[0]=x;p[1]=y;p[2]=z;
points.push_back(p);
}}}
PointBin grid(2.5,2.5,2.5,points);
std::vector<double> p(3);
p[0]=4.3;p[1]=6.3;p[2]=2.3;
std::vector<int> list;
std::cout<<"everything ok"<<std::endl;
list=grid.GetList(p);
}

/*class try
{
public:
try(){};
double value(){return _value;};

private:
double _value;
}*/
