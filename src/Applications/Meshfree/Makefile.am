bin_PROGRAMS    = indent indent_withAL homog cons_test postProcess normal_mode pressurize calF SHenergy homog_table pressurize_load homog_icossym homog_refchange calStrain
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
		-I$(blitz_includes) -I$(tvmet_includes)
indent_SOURCES    =  indent.cc 
homog_SOURCES = homog.cc
cons_test_SOURCES = cons_test.cc
postProcess_SOURCES = postProcess.cc
indent_withAL_SOURCES = indent-withAL.cc
normal_mode_SOURCES = normal-mode.cc
pressurize_SOURCES = pressurize.cc
pressurize_load_SOURCES = pressurize-load.cc
calF_SOURCES = calF.cc
calStrain_SOURCES = calStrain.cc
SHenergy_SOURCES = SHenergy.cc
homog_table_SOURCES = homog-table.cc
homog_icossym_SOURCES = homog-icossym.cc
homog_refchange_SOURCES = homog-refchange.cc
AM_LDFLAGS    = -L$(blitz_libraries) \
        -L../                          	\
        -L../../Materials/             	\
	-L../../Solvers/               	\
	-L../../Model/                 	\
	-L../../Quadrature/		\
	-L../../Geometry/		\
        -L../../Shape/			\
        -L../../Elements/		\
	-L../../VoomMath/			\
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
	-llapack -lblas			
