bin_PROGRAMS    = minimize excite_modes cal_inv homog2Dtable
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
		-I$(blitz_includes) -I$(tvmet_includes)
minimize_SOURCES = minimize.cc
excite_modes_SOURCES = excite-modes.cc
cal_inv_SOURCES = cal-inv.cc
homog2Dtable_SOURCES = homog2Dtable.cc
AM_LDFLAGS    = -L$(blitz_libraries) \
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
        -lsrfftw                         \
        -lsfftw                          \
        -lm
