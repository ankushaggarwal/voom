#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include "mkl.h"

//extern "C" void mkl_ddnscsr (const MKL_INT *job , const MKL_INT *m , const MKL_INT *n , double *adns , const MKL_INT *lda , double *acsr , MKL_INT *ja , MKL_INT *ia , MKL_INT *info );

int main(int argc, char* argv[])
{

MKL_INT job[]={0, 0, 1, 2, 9, 1};
const MKL_INT m=3, n=3;
double adns[]={1., 0., 3., 0., 0., 6., 3., 6., 9. };
MKL_INT lda[]={3};
double acsr[9];
MKL_INT ja[9], ia[4], info, eig_found;

mkl_ddnscsr (job , &m , &n , adns , lda , acsr , ja , ia , &info );
//mkl_scsrmv();
//cblas_daxpy (9, 2., adns, 1, adns, 1);

for(int i=0;i<9;i++)
std::cout << adns[i] << std::endl;
for(int i=0;i<ia[3];i++)
std::cout << ja[i] << "\t" << acsr[i] << std::endl;
ia[3]++;

double res,epsout, eig[3];
MKL_INT fpm[128], loop;
feastinit(fpm);
const char uplo[]={'F'};
const double emin=-100.,emax=100.; 
MKL_INT m0=3;
//dfeast_scsrev(uplo, n, a, ia, ja, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info)
dfeast_scsrev(uplo, &n, acsr, ia, ja, fpm,&epsout,&loop,&emin,&emax, &m0, eig, adns, &eig_found, &res, &info );

std::cout << info << std::endl;
std::cout << loop << " iterations taken to converge" << std::endl;
std::cout << eig_found << " eigenvalues found in the interval" << std::endl;
for(int i=0;i<eig_found;i++){
std::cout << i << "th-eigvenvalue=" << eig[i] << std::endl;
std::cout << "Corresponding eigenvector=" << std::endl;
std::cout << adns[n*i] << "\t" << adns[1+n*i] << "\t" << adns[2+n*i] << std::endl;
}
return 0;

}
