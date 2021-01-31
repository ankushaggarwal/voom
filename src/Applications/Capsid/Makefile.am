bin_PROGRAMS    = relax indent refine disassemble grow growAlt
AM_CPPFLAGS =-I$(srcdir)/../../       		\
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
#
AM_LDFLAGS = -L$(blitz_libraries) \
        -L../                          	\
        -L../../Materials/             	\
	-L../../Solvers/               	\
	-L../../Model/                 	\
	-L../../Quadrature/		\
	-L../../Geometry/		\
        -L../../Shape/			\
        -L../../Body/		\
        -L../../Elements/		\
	-L../../VoomMath/		\
	-L../../Mesh/
#
LDADD	= -lblitz              		\
        -lMaterials                    	\
	-lModel                        	\
	-lSolvers                      	\
	-lQuadrature                   	\
	-lGeometry			\
        -lShape				\
	-lBody			\
	-lElements			\
	-lVoomMath				\
	-lMesh				\
	-lvtkIO 			\
	-lvtkGraphics			\
	-lvtkGenericFiltering 		\
	-lvtkFiltering 			\
	-lvtkCommon 			\
	-llapack			\
	-lblas				\
	-lgfortran
#
# capsid_lbfgsb_SOURCES    =  capsid-lbfgsb.cc 
indent_SOURCES    =  indent.cc 
#
relax_SOURCES    =  relax.cc 
#
refine_SOURCES    =  refine.cc 
#
disassemble_SOURCES    =  disassemble.cc 
#
grow_SOURCES    =  grow.cc C0MembraneGLoutput.cc
#
growAlt_SOURCES    =  grow-alt.cc C0MembraneGLoutput.cc
