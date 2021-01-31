  void Model::normalModes_all_sparse(bool eigenvectors, int nmodes, std::string filename) {

    bool verbose = false;
    //for time calculation
    time_t start, end;
    double duration=0.;
    blitz::Range all = blitz::Range::all();
    blitz::Array<double,2> k_n(_dof,_dof);
    k_n = 0.0;

    std::cout << "Entered the normal mode calculation." << std::endl;
    { ///////////////////////////////////////////
    Storage solver;
    std::cout<<"Initializing the stiffness matrix " << _dof << std::endl;
    solver.resize(_dof);

    //blitz::Array<double,2> &k_n = solver._DDE;

    std::cout<<"Stiffness matrix initialized"<<std::endl;
    
    getField(solver);
    blitz::Array<double,1>& x = solver._x;

    const double h = 1.0e-8*blitz::max(blitz::abs(x));

    // numerical force and stiffness

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
    int nNonZeros = 0;
    double kMax = blitz::max(blitz::abs(k_n(all,all) ) );
    for (int i = 0 ; i < _dof; i++){ 
      for (int j = 0 ; j < _dof; j++) {
        double xtemp=k_n(i,j)>0?k_n(i,j):-k_n(i,j);
        if( xtemp/kMax < TOL) k_n(i,j)=0.;
        else nNonZeros++; 
      }
    }
    std::cout << "nNonZeros = " << nNonZeros << std::endl;
    std::cout << "Sparsity = " << 1.0*nNonZeros/_dof/_dof << std::endl;
    std::cout << "Converting to sparse format." << std::endl;

    const MKL_INT n=_dof;
    double val[nNonZeros];
    MKL_INT ja[nNonZeros], ia[_dof+1], info;

    int id=1;
    for (int i = 0 ; i < _dof; i++){ 
      ia[i]=id;
      for (int j = 0 ; j < _dof; j++) {
        double xtemp=k_n(i,j)>0?k_n(i,j):-k_n(i,j);
        if( xtemp/kMax < TOL) continue;
        val[id-1]=k_n(i,j);
        ja[id-1]=j+1;
        id++;
      }
    }
    solver.resize(0);
    ia[_dof]=nNonZeros+1;

    std::cout << "Stiffness matrix calculated. Calculating eigenvalues." << std::endl;

    // compute Eigenvalues and Eigenvectors by calling FEAST library from MKL
    const char uplo='F';
    MKL_INT fpm[128], loop;
    feastinit(fpm);
    fpm[0]=1;
    const double emin=-1,emax=100;
    MKL_INT m0=100, eig_found;
    double res,epsout, eigenvalues[m0];
    double eig_vec[m0*_dof];
    std::cout << "Everything ready, making call to FEAST." << std::endl;
    time(&start);
    dfeast_scsrev(&uplo, &n, val, ia, ja, fpm, &epsout, &loop, &emin, &emax, &m0, eigenvalues, eig_vec, &eig_found, &res, &info );
    time(&end);

    duration=difftime(end,start);
    std::cout<<"Time taken for eigenvalue decomposition is "<<duration<<" seconds."<<std::endl;
    if (info != 0) {
      std::cout << "Something is wrong in eigensolver. Info=" << info << std::endl;
      exit(0);
    }

    //print out the eigenvectors
    if (eigenvectors) {
    std::cout << "Printing eigenvectors." << std::endl;
    int vec_ptr=0;
    blitz::Array<double,1> eig_vec_i(_dof);
    eig_vec_i = 0.0;
    Storage solver;
    solver.resize(this->dof());
    getField(solver);
    std::cout << "getField done" << std::endl;
    blitz::Array<double,1>& x = solver._x;
      for (int i = 0; i < eig_found; i++) { 
        double h = 1.0;
        std::cout << "1" << std::endl;
        for (int j=0; j<_dof; j++) eig_vec_i(j)=eig_vec[vec_ptr++];
        std::cout << "2" << std::endl;
        for (int j=0; j<_dof; j++) {std::cout << j << std::endl; x(j) += h*eig_vec_i(j); }
        //solver._x += h*k_n(blitz::Range::all(),i);
        std::cout << "3" << std::endl;
        putField(solver);
        char name[20];
        std::string fname = filename;
        sprintf(name,"evec%d",i);
        std::cout << "4" << std::endl;
        fname += name;
        print(fname);
        x -= h*eig_vec_i;
      }
    }

    std::cout << loop << " iterations taken to converge." << std::endl;
    std::cout << eig_found << " eigenvalues found in the interval." << std::endl;
    //std::ofstream evals;
    //evals.open("evals.dat");
    std::cout<<"The eigenfrequencies are:"<<std::endl;
    for (int i = 0; i < eig_found; i++) 
      std::cout << eigenvalues[i] << std::endl;
    //evals.close();
    } ///////////////////////////////////////////
    std::cout << "Out of scope" << std::endl;
  }
