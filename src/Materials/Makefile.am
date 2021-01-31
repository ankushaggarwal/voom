## Makefile.am -- Process this file with automake to produce Makefile.in
AM_CPPFLAGS=-I$(srcdir)/.. -I$(srcdir)/../Geometry -I$(srcdir)/../Node -I$(srcdir)/../VoomMath -I$(blitz_includes) -I$(tvmet_includes) 
lib_LIBRARIES=libMaterials.a
libMaterials_a_SOURCES= 	\
	SCElastic.cc 		\
	EvansElastic.cc		\
	TwoPhaseElastic.cc	\
	GLElastic.cc             \
	FVK.cc			\
	StVenant.cc 		\
	LinearizedMaterial.cc 	\
	Hookean2D.cc 		\
	Hookean.cc              \
	CompNeoHookean.cc       \
        ModEvansElastic.cc      \
        MooneyRivlin.cc         \
        HomogMP.cc         \
	EvansElastic_Skewed.cc  \
        WarpingPenalty.cc    \
        InverseModel_EvansElastic.cc
