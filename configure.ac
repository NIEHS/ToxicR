# Configure Script for ToxicR
# Modified 4/27/2022
# Configure Script for ToxicR
# Modified 4/27/2022

AC_INIT([ToxicR], 1.0.1)  dnl Package version

CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_LANG(C++)
AC_PROG_CPP

#determine the OS
AC_MSG_CHECKING([Is OPENMP avaialble? ])
OS_FLAG="$(uname)"
if test x"${OS_FLAG}" == x"Darwin";then
     AC_MSG_RESULT([no])
     OPENMP="-DNO_OMP"
else
     AC_MSG_RESULT([yes])
     OPENMP="\$(SHLIB_OPENMP_CXXFLAGS) "
fi

AC_PATH_PROG([GSL_CONFIG], [gsl-config])
## If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then
    AC_MSG_RESULT([yes])
    # Use gsl-config for header and linker arguments    
    gsl_include=$(gsl-config --cflags)
    gsl_libs=$(gsl-config --libs)
else
    AC_MSG_ERROR([gsl-config not found, is GSL installed?
                 To install GSL in Ubuntu Type: 
                 sudo apt install gsl 
                 To install GSL in Fedora Type:
                 sudo yum -y install gsl
                 To iinstall GSL on macOS using homebrew type: 
                 brew install gsl
                 ])
fi

## Can we use pkg-config?
AC_PATH_PROG(have_pkg_config, pkg-config, no)
AC_MSG_CHECKING([Is pkg-config avaialble? ])
if test x"${have_pkg_config}" != x"no"; then
     AC_MSG_RESULT([yes])
     AC_MSG_CHECKING([if pkg-config knows NLopt])
     if pkg-config --exists nlopt; then
           AC_MSG_RESULT([yes])
            ## Since nlopt has been found, test for minimal version requirement
            AC_MSG_CHECKING([for pkg-config checking NLopt version])
            if pkg-config --atleast-version=2.6.0 nlopt; then
                 AC_MSG_RESULT([>= 2.6.0])
                 nlopt_include=$(pkg-config --cflags nlopt)
                 nlopt_libs=$(pkg-config --libs nlopt)
                 AC_SUBST([NLOPT_INCLUDE], "${nlopt_include}")
                 AC_SUBST([NLOPT_LIBS],    "${nlopt_libs}")   
            else
                AC_MSG_RESULT([no])
                #if it is macOS just try link and use nlopt generically
              need_to_build="yes"
            fi
    else
        AC_MSG_RESULT([no])
            #if it is macOS just try link and use nlopt generically
            if test x"${OS_FLAG}" == x"Darwin";then
                    AC_MSG_RESULT([>= 2.6.0])
                    nlopt_include=""
                    nlopt_libs="-lnlopt"
                    AC_SUBST([NLOPT_INCLUDE], "${nlopt_include}")
                    AC_SUBST([NLOPT_LIBS],    "${nlopt_libs}")   
            else
                 need_to_build="yes"
            fi
       
    fi
               
else
    AC_MSG_RESULT([no])
    need_to_build="yes"
fi

## So do we need to build
if test x"${need_to_build}" != x"no"; then
    AC_PATH_PROG(have_cmake, cmake, no)
  
    if test x"${have_cmake}" == x"no"; then
       . src/scripts/cmake_config.sh
       if test -z "${CMAKE_BIN}"; then
          ## also error to end configure here
          AC_MSG_ERROR([Could not find 'cmake'.])
       fi
    fi
    ## 'uname -m' on M1 give x86_64 which is ... not helping
    machine=`"${R_HOME}/bin/Rscript" -e 'cat(Sys.info()[["machine"]])'`
    AC_MSG_RESULT([using NLopt via local cmake build on ${machine} ])
    tools/cmake_call.sh
    ## cmake_call.sh installs into nlopt/lib, headers are copied
    nlopt_include="-I./nlopt/include"
    nlopt_libs="-L./nlopt/lib -lnlopt"
fi

SUBDIR_SOURCES="$(cd src/ && ls {code_base,polyK}/*.cpp | tr '\n' ' ')"
SRC_SOURCES="$( cd src/ && ls *.cpp | tr '\n' ' ')"

AC_SUBST(OPENMP)
AC_SUBST(SRC_SOURCES)
AC_SUBST(SUBDIR_SOURCES)
AC_SUBST([NLOPT_CPPFLAGS],["${nlopt_include}"])
AC_SUBST([NLOPT_LIBS],["${nlopt_libs}"])
AC_SUBST([GSL_CPPFLAGS],["${gsl_include}"])
AC_SUBST([GSL_LIBS],["${gsl_libs}"])
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
