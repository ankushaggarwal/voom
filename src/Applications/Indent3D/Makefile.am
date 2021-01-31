bin_PROGRAMS    = indent3D try
#CXXFLAGS= -g -ggdb -W -Wall
INCLUDES        =-I ./                 \
        -I$(blitz_includes)            \
        -I$(tvmet_includes)            \
        -I../                          \
        -I../../                       \
        -I../../../                    \
        -I../../../../                 \
        -I../../VoomMath/                  \
        -I../../Node/                  \
        -I ../../Materials/            \
        -I ../../Geometry/             \
	-I/home/software/              \
	-I../../AdsExt                 \
        -I ../../Quadrature/           \
	-I ../../Elements/             \
	-I ../../Constraints/          \
	-I ../../Potentials/           \
	-I ../../Body/                 \
	-I ../../Model/                \
	-I ../../Solvers/              \
        -I ../../Shape/
#
indent3D_SOURCES    =  indent3D.cc 
try_SOURCES = try.cc
indent3D_LDFLAGS    = -L$(blitz_libraries) \
        -L../                          	\
        -L../../Materials/             	\
	-L../../Solvers/               	\
	-L../../Model/                 	\
	-L../../Quadrature/		\
	-L../../Geometry/		\
        -L../../Shape/			\
        -L../../Elements/		\
	-L../../VoomMath/		\
	-L../../Body/
indent3D_LDADD      = -lblitz              \
        -lMaterials                    	\
	-lModel                        	\
	-lSolvers                      	\
	-lQuadrature                   	\
	-lGeometry			\
        -lShape				\
	-lBody				\
	-lElements			\
	-lVoomMath				\
	-llapack			\
	-lblas	                        \
	-lgfortran
