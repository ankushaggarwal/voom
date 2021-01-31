#include "mkl.h"

  void Model::normalModes_all_sparse(bool eigenvectors, int nmodes, std::string filename) {

    bool verbose = false;
    //for time calculation
    time_t start, end;
    double duration=0.;

    std::cout << "Entered the normal mode calculation." << std::endl;
    Storage solver;
    std::cout<<"Initializing the stiffness matrix"<<this->dof()<<std::endl;
    solver.resize(this->dof());
    std::cout<<"Initializing the stiffness matrix"<<std::endl;
    blitz::Array<double,2> &k_n = solver._DDE;
    std::cout<<"Stiffness matrix initialized"<<std::endl;
    k_n = 0.0;
    getField(solver);
    blitz::Array<double,1>& x = solver._x;

    const double h = 1.0e-8*blitz::max(blitz::abs(x));

    // numerical force and stiffness
    blitz::Range all = blitz::Range::all();
    for (int i = 0 ; i < _dof; i++) {
      if (verbose) std::cout << " i = " << i << std::endl;
      x(i) += h; 
      putField(solver);
      computeAndAssemble(solver, false, true, false);
      k_n(i,all) = solver._DE;

      x(i) -= 2.0*h; putField(solver);
      computeAndAssemble(solver, false, true, false);
      k_n(i,all) -= solver._DE;
  
      x(i) += h; putField(solver);
    }
    k_n /= 2.0*h;

    k_n = 0.5*(k_n + k_n.transpose(blitz::secondDim,blitz::firstDim));

    //Assigning zero to values below a tolerance
    double TOL=1e-12;
    nNonZeros = 0;
    double kMax = max(blitz::abs(k_n(all,all) ) );
    for (int i = 0 ; i < _dof; i++) 
      for (int j = 0 ; j < _dof; j++) 
        if( abs(k_n(i,j)/kMax) < TOL) k_n(i,j)=0.;
        else nNonZeros++;

    std::cout << "Converting to sparse format." << std::endl;

    MKL_INT job[]={0, 0, 1, 2, 9, 1};
    const MKL_INT n=_dof;
    double dens_mat[];
    dens_mat=k_n.data();
    MKL_INT lda=_dof;
    double val[nNonZeros];
    MKL_INT ja[nNonZeros], ia[_dof+1], info;

    mkl_ddnscsr (job , &n , &n , dens_mat , &lda , val , ja , ia , &info );

    std::cout << "Stiffness matrix calculated. Calculating eigenvalues." << std::endl;

    // compute Eigenvalues and Eigenvectors by calling FEAST library from MKL
    const char uplo='F';
    MKL_INT fpm[128], loop;
    feastinit(fpm);
    const double emin=0.,emax=1.e8.;
    MKL_INT m0=_dof, eig_found;
    double res,epsout, eigenvalues[_dof];

    time(&start)
    dfeast_scsrev(&uplo, &n, val, ia, ja, fpm, &epsout, &loop, &emin, &emax, &m0, eigenvalues, dens_mat, &eig_found, &res, &info );
    time(&end);

    duration=difftime(end,start);
    std::cout<<"Time taken for eigenvalue decomposition is "<<duration<<" seconds."<<std::endl;

    if (info != 0) {
      std::cout << "Something is wrong in eigensolver" << std::endl;
      exit(0);
    }

    //print out the eigenvectors
    if (eigenvectors) {
      for (int i = 0; i < nmodes; i++) { 
        double h = 1.0;//0.3/blitz::max(k(i,blitz::Range::all()));
        solver._x += h*k_n(i,blitz::Range::all());
        //solver._x += h*k_n(blitz::Range::all(),i);
        putField(solver);
        char name[20];
        std::string fname = filename;
        sprintf(name,"evec%d",i);
        fname += name;
        print(fname);
        solver._x -= h*k_n(i,blitz::Range::all());
      }
    }

    std::cout << loop << " iterations taken to converge." << std::endl;
    std::cout << eig_found << " eigenvalues found in the interval." << std::endl;
    std::ofstream evals;
    evals.open("evals.dat");
    std::cout<<"The eigenfrequencies are:"<<std::endl;
    for (int i = 0; i < eig_found; i++) 
      evals << eigenvalues[i] << std::endl;
    evals.close();
  }
