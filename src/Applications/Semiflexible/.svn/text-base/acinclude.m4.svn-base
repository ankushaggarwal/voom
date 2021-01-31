
AC_DEFUN([VOOM_CHECK_BLITZ],
[
ac_blitz_includes=NO
ac_blitz_libraries=NO

blitz_includes="/u/home3/klug/local/blitz-gcc/include"
blitz_libraries="/u/home3/klug/local/blitz-gcc/lib"

AC_ARG_WITH(blitz-root,
    [  --with-blitz-root=DIR       Blitz root directory],
    [
	blitz_includes="$withval"/include
       	blitz_libraries="$withval"/lib
    ])

if test ! "$ac_blitz_includes" = "NO"; then
   blitz_includes="$ac_blitz_includes $blitz_includes"
fi

if test ! "$ac_blitz_libraries" = "NO"; then
   blitz_libraries="$ac_blitz_libraries $blitz_libraries"
fi

AC_MSG_CHECKING([whether the Blitz library is installed])
cat > conftest.cc <<EOF
#include <blitz/array.h>
using namespace blitz;
int main()
{
    Array<double,1> x(100);
    x = tensor::i;          // x = [ 0, 1, 2, ..., 99 ]

    Array<double,1> z(x + 150);
    Array<double,1> v(z + x * 2);

    cout << v << endl;
}
EOF

ac_link='${CXX-g++} -o conftest $CXXFLAGS -I$blitz_includes -L$blitz_libraries conftest.cc -lblitz 1>&5'

if AC_TRY_EVAL(ac_link) && test -s conftest; then
  AC_MSG_RESULT(yes)
  AC_DEFINE_UNQUOTED(HAVE_BLITZ, 1, [Define if you have a Blitz installation])
  rm -f conftest.cc
else
  AC_MSG_ERROR([
Your system fails at linking a small Blitz application!
Check if the Blitz headers and libraries can be found (you can specify include
and library directories; run ./configure --help for more info), or if your
compiler is installed correctly, or if the compiler you are using now is the
same you used to compile Blitz. For more details about this problem,
look at the end of config.log.])
fi
])

%
AC_DEFUN([VOOM_CHECK_TVMET],
[
ac_tvmet_includes=NO

tvmet_includes="/u/home3/klug/local/tvmet-gcc/include"

AC_ARG_WITH(tvmet-root,
    [  --with-tvmet-root=DIR       tvmet root directory],
    [
	tvmet_includes="$withval"/include
    ])

AC_ARG_WITH(tvmet-includes,
    [  --with-tvmet-includes=DIR  tvmet header files],
    [
       	tvmet_includes="$withval"
    ])

if test ! "$ac_tvmet_includes" = "NO"; then
   tvmet_includes="$ac_tvmet_includes $tvmet_includes"
fi

AC_MSG_CHECKING([whether the tvmet library is installed])
cat > conftest.cc <<EOF
#include <iostream>
#include <complex>
#include <tvmet/Matrix.h>

using namespace tvmet;
using std::cout;
using std::endl;

typedef Matrix<std::complex<double>,3,3>	matrix33d;

void testMM(matrix33d& res, const matrix33d& m1, const matrix33d& m2) {
  res = m1 * m2;
}

int main()
{
  matrix33d m1, m2, m3;

  m1 = 1,4,7,
       2,5,8,
       3,6,9;
  m2 = m1;

  testMM(m3, m1, m2);

  cout << m1 << "\n*\n" << m2 << "\n=";
  cout << m3 << endl;
}
EOF

ac_link='${CXX-g++} -o conftest $CXXFLAGS -I$tvmet_includes conftest.cc 1>&5'

if AC_TRY_EVAL(ac_link) && test -s conftest; then
  AC_MSG_RESULT(yes)
  AC_DEFINE_UNQUOTED(HAVE_TVMET, 1, [Define if you have a tvmet installation])
  rm -f conftest.cc
else
  AC_MSG_ERROR([
Your system fails at linking a small tvmet application!
Check if the tvmet headers and libraries can be found (you can specify include
and library directories; run ./configure --help for more info), or if your
compiler is installed correctly, or if the compiler you are using now is the
same you used to compile tvmet. For more details about this problem,
look at the end of config.log.])
fi
])

%
AC_DEFUN([VOOM_CHECK_VTK],
[
ac_vtk_includes=NO
ac_vtk_libraries=NO

vtk_includes="/u/home3/klug/local/vtk-gcc/include/vtk-5.4/"
vtk_libraries="/u/home3/klug/local/vtk-gcc/lib/vtk-5.4/"

AC_ARG_WITH(vtk-root,
    [  --with-vtk-root=DIR       vtk root directory],
    [
	vtk_includes="$withval"/include
	vtk_libraries="$withval"/lib
    ])
AC_ARG_WITH(vtk-includes,
    [  --with-vtk-includes=DIR  vtk header files],
    [
       	vtk_includes="$withval"
    ])
AC_ARG_WITH(vtk-libraries,
    [  --with-vtk-libraries=DIR  vtk lib files],
    [
       	vtk_libraries="$withval"
    ])

if test ! "$ac_vtk_includes" = "NO"; then
   vtk_includes="$ac_vtk_includes $vtk_includes"
fi

if test ! "$ac_vtk_libraries" = "NO"; then
   vtk_libraries="$ac_vtk_libraries $vtk_libraries"
fi

AC_MSG_CHECKING([whether the vtk library is installed])
cat > conftest.cc <<EOF

#include <vtkPoints.h>
int main()
{
  vtkPoints *newPts = vtkPoints::New();
  newPts->InsertNextPoint( 0.1, 0.2, 0.3 );
}
EOF

ac_link='${CXX-g++} -o conftest $CXXFLAGS -I$vtk_includes 
         -L$vtk_libraries conftest.cc 
	 -lvtkCommon -lvtksys -lpthread -ldl 1>&5'

if AC_TRY_EVAL(ac_link) && test -s conftest; then
  AC_MSG_RESULT(yes)
  AC_DEFINE_UNQUOTED(HAVE_VTK, 1, [Define if you have a vtk installation])
else
  AC_MSG_ERROR([
Your system fails at linking a small vtk application!
Check if the vtk headers and libraries can be found (you can specify include
and library directories; run ./configure --help for more info), or if your
compiler is installed correctly, or if the compiler you are using now is the
same you used to compile vtk. For more details about this problem,
look at the end of config.log.])
fi
])

