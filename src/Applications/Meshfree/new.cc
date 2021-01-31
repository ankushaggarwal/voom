#include<tvmet/Vector.h>

using namespace tvmet;
using namespace std;

main(){
  tvmet::Vector<double, 3> v1(1,2,3);
  tvmet::Vector<double,3> x1,x2,x3,x4,p,xo,temp,RHS,coeff,theta;
  //centroid of the tetrahedra
  xo=(x1+x2+x3+x4)/4;
  //origin shifted to the centroid
  x1 = x1 - xo;
  x2 = x2 - xo;
  x3 = x3 - xo;
  x4 = x4 -xo;
  p = p - xo;

  //face 1 (not x4)
  temp=tvmet::cross(x1-x2,x1-x3);
temp=cross(x1-x2,x1-x3);
temp2=cross(cross(x2,x3),temp);
coeff=cross(cross(x2,p),temp)+cross(cross(p,x3),temp);
theta=temp2./coeff;
if(theta(1)==theta(2) && theta(2)==theta(3))
    if(theta(1)<1 && theta(1)>=0)
        disp('outside');
        return
    else
        disp('inside');
    end
else
    disp('error, line and plane intersection wrong');
    return
end

%face 2 (not x1)
temp=cross(x2-x3,x2-x4);
temp2=cross(cross(x3,x4),temp);
coeff=cross(cross(x3,p),temp)+cross(cross(p,x4),temp);
theta=temp2./coeff;
if(theta(1)==theta(2) && theta(2)==theta(3))
    if(theta(1)<1 && theta(1)>=0)
        disp('outside');
        return
    else
        disp('inside');
    end
else
    disp('error, line and plane intersection wrong');
    return
end

%face 3 (not x2)
temp=cross(x3-x4,x3-x1);
temp2=cross(cross(x4,x1),temp);
coeff=cross(cross(x4,p),temp)+cross(cross(p,x1),temp);
theta=temp2./coeff;
if(theta(1)==theta(2) && theta(2)==theta(3))
    if(theta(1)<1 && theta(1)>=0)
        disp('outside');
        return
    else
        disp('inside');
    end
else
    disp('error, line and plane intersection wrong');
    return
end

%face 4 (not x3)
temp=cross(x4-x1,x4-x2);
temp2=cross(cross(x1,x2),temp);
coeff=cross(cross(x1,p),temp)+cross(cross(p,x2),temp);
theta=temp2./coeff;
if(theta(1)==theta(2) && theta(2)==theta(3))
    if(theta(1)<1 && theta(1)>=0)
        disp('outside');
        return
    else
        disp('inside');
    end
else
    disp('error, line and plane intersection wrong');
    return
end
}
