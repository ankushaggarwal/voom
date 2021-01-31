bin_PROGRAMS    = minimize normal_modes indent 
INCLUDES        =-I$(srcdir)/../../       		\
		-I$(srcdir)/../../Body/		\
		-I$(srcdir)/../../Model/		\
		-I$(srcdir)/../../Geometry/       	\
		-I$(srcdir)/../../Node/           	\
		-I$(srcdir)/../../Shape/          	\
		-I$(srcdir)/../../Materials/       	\
		-I$(srcdir)/../../Elements/       	\
		-I$(srcdir)/../../VoomMath/      	\
		-I$(srcdir)/../../Quadrature/		\
		-I$(srcdir)/../../Solvers/		\
		-I$(srcdir)/../../Mesh/		\
		-I$(blitz_includes) -I$(tvmet_includes) \
                -fopenmp
minimize_SOURCES = minimize.cc
normal_modes_SOURCES = normal-modes.cc
indent_SOURCES = indent.cc
CFLAGS = -fopenmp
AM_LDFLAGS    = -L$(blitz_libraries) \
        -Wl,--no-as-needed -fopenmp \
        -L../                          	\
        -L../../Materials/             	\
	-L../../Solvers/               	\
	-L../../Model/                 	\
	-L../../Quadrature/		\
	-L../../Geometry/		\
        -L../../Shape/			\
        -L../../Body/                   \
        -L../../Elements/		\
	-L../../VoomMath/			\
        -L../../Mesh/                   \
        -L/opt/intel/composer_xe_2015.3.187/mkl/lib/intel64/ \
        -L/opt/intel/composer_xe_2015.3.187/compiler/lib/intel64/ \
	-L/usr/lib/gcc/i486-linux-gnu/4.2/	\
	-L../../Body/	
LDADD      = -lblitz              \
        -lMaterials                    	\
	-lModel                        	\
	-lSolvers                      	\
	-lQuadrature                   	\
	-lGeometry			\
        -lShape				\
	-lElements			\
	-lBody				\
	-lVoomMath				\
        -lMesh                          \
        -lpthread                       \
        -ldl                            \
	-llapack -lblas			\
	-lgfortran                      \
        -lmkl_core        \
        -lmkl_intel_lp64 \
        -lmkl_intel_thread \
        -lmkl_sequential \
     -lpthread -liomp5 
