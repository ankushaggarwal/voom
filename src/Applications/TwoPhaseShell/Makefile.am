bin_PROGRAMS    = minimize minimize_indent minimize_indent_restart minimize_tubule homog2D normal_modes SHenergy SHenergy_CK SHenergy_nonuniformR homog2Dtable homog2Dicossym indent calStrain calStrain2 
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
AM_CPPFLAGS = -w
minimize_SOURCES = minimize.cc
minimize_indent_SOURCES = minimize-indent.cc
minimize_indent_restart_SOURCES = minimize-indent-restart.cc
minimize_tubule_SOURCES = minimize-tubule.cc
homog2D_SOURCES = homog2D.cc
normal_modes_SOURCES = normal-modes.cc
SHenergy_SOURCES = SHenergy.cc
SHenergy_CK_SOURCES = SHenergy2.cc
SHenergy_nonuniformR_SOURCES = SHenergy-nonuniformR.cc
homog2Dtable_SOURCES = homog2Dtable.cc
homog2Dicossym_SOURCES = homog2Dicossym.cc
indent_SOURCES = indent.cc
calStrain_SOURCES = calStrain.cc
calStrain2_SOURCES = calStrain2.cc
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
	-llapack -lblas			
