#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([RNAHeliCes], [3.0.0],
        [Jiabin Huang and Bjoern Voss <jiabin.huang, bjoern.voss@biologie.uni-freiburg.de>],
        [RNAHeliCes])
AC_PREREQ([2.59])
AC_DEFINE(RELEASE_DATE, ["Apr. 23, 2013"], [Release date of current version])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR(config)    #AC_CONFIG_AUX_DIR([.])
dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AC_CONFIG_HEADERS([src/config.h])    #### combine with autoheader to generate src/config.h.in, while ./configure config.h.in ==> config.h 

AM_INIT_AUTOMAKE([1.10 -Wall no-define])


# General settings.
#AC_SUBST(VERSION)
#ISODATE=`date +%Y-%m-%d`
#AC_SUBST(ISODATE)
AC_CANONICAL_HOST
AC_GNU_SOURCE



#http://jurjenbokma.com/ApprenticesNotes/cplusplus0x_with_autotools.html
#AC_PROG_CXX
#
## Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL
#AC_C_INLINE
#AC_TYPE_SIZE_T
#
#AC_LANG([C++])
#AC_LANG_CPLUSPLUS
#AC_CHECK_HEADER(cstdint)
#AX_CXX_COMPILE_STDCXX_0X
#
#CXXFLAGS='-Wall -std=gnu++0x'

echo "o----------------------------------"
echo "| Checking for compilers"
echo "o----------------------------------"
#AC_PROG_CC([cc gcc])
AC_PROG_CC_C99
AC_ARG_VAR(FPIC,[compiler flag for generating PIC code that can be used in a shared library])
if test "$GCC" = yes; then
   AC_DEFINE(UNUSED, __attribute__ ((unused)), avoid warning about unused variables)
else
   AC_DEFINE(UNUSED,, only for gcc)
fi

#AC_PROG_F77([f77 g77 gfortran f90 xlf90 f95 g95 gfortran-dp-4.2])
AC_PROG_CXX([c++ g++ cxx CC])
AM_PROG_CC_C_O
LT_INIT([disable-static])    # AC_PROG_RANLIB
AC_PROG_LIBTOOL




#echo "o----------------------------------"
#echo "| Setting up Doxygen"
#echo "o----------------------------------"
#echo "it's standing for nothing currently."
# DX_HTML_FEATURE(ON)
# DX_CHM_FEATURE(OFF)
# DX_CHI_FEATURE(OFF)
# DX_MAN_FEATURE(OFF)
# DX_RTF_FEATURE(OFF)
# DX_XML_FEATURE(OFF)
# DX_PDF_FEATURE(OFF)
# DX_PS_FEATURE(OFF)
# DX_INIT_DOXYGEN(RNABarrier, doxygen.cfg, doc)


echo "o----------------------------------"
echo "| Checking for std headers"
echo "o----------------------------------"
#AC_CHECK_HEADERS([cstdio stdio.h], [break],)
#AC_CHECK_HEADERS([cstdarg stdarg.h], [break],)
#AC_CHECK_HEADERS([cstdlib stdlib.h], [break],)
#AC_CHECK_HEADERS([string string.h strings.h], [break],)
#AC_CHECK_HEADERS([iostream iostream.h], [break],)
#AC_CHECK_HEADERS([fstream fstream.h], [break],)
#AC_CHECK_HEADERS([stdexcept stdexcept.h], [break],)
#AC_CHECK_HEADERS([cassert assert.h], [break],)
#AC_CHECK_HEADERS([complex complex.h], [break],)
AC_CHECK_HEADERS([vector vector.h], [break],)
AC_CHECK_HEADERS([algorithm algo.h algorithm.h], [break],)
#AC_CHECK_HEADERS([map map.h], [break],)
#AC_CHECK_HEADERS([list list.h], [break],)
#AC_CHECK_HEADERS([set set.h], [break],)
#AC_CHECK_HEADERS([typeinfo], [break],)
#AC_CHECK_HEADERS([limits], [break],)
#AC_CHECK_HEADERS([iomanip iomanip.h], [break],)
#AC_CHECK_HEADERS([cmath math.h], [break],)
#AC_CHECK_HEADERS([sstream sstream.h strstream.h], [break],)

#AC_CXX_NUMERIC_LIMITS

AX_CXX_COMPILE_STDCXX_0X
#CXXFLAGS='-Wall'
#CXXFLAGS='-Wall -std=gnu++0x -Wno-deprecated'
#CXXFLAGS=-std=gnu++0x -Wno-deprecated



#### Checks for programs needed for tests
#echo "o------------------------------------------------"
#echo "| Checking for programs needed for tests"
#echo "o------------------------------------------------"
#### Check for 'diff' and get full path.
#AC_ARG_VAR([DIFF],[the 'diff' program to use for test output comparison])
#AC_PATH_PROG([DIFF],[diff],[])
#if test "x$DIFF" = "x"; then
#    AC_MSG_NOTICE([==> diff command not found!])
#    AC_MSG_NOTICE([==> Set DIFF variable if present in non-standard path!])
#    AC_MSG_ERROR([diff is mandatory to run the tests : will stop here!])
#fi

#### Check for 'cat' and get full path.
#AC_ARG_VAR([CAT],[the 'cat' program used for printing test output files])
#AC_PATH_PROG([CAT],[cat],[])
#if test "x$CAT" = "x"; then
#    AC_MSG_NOTICE([==> cat command not found!])
#    AC_MSG_NOTICE([==> Set CAT variable if present in non-standard path!])
#    AC_MSG_ERROR([cat is mandatory to run the tests : will stop here!])
#fi

#### Check for 'printf' and get full path.
#AC_ARG_VAR([PRINTF],[the 'printf' program used to print test information])
#AC_PATH_PROG([PRINTF],[printf],[])
#if test "x$PRINTF" = "x"; then
#    AC_MSG_NOTICE([==> printf command not found!])
#    AC_MSG_NOTICE([==> Set PRINTF variable if present in non-standard path!])
#    AC_MSG_ERROR([printf is mandatory to run the tests : will stop here!])
#fi




# Checks for libraries.

echo "o--------------------------------------------"
echo "| Checking for C++ Libraries BOOST"
echo "o--------------------------------------------"
# minimal reqired version of the BOOST library
#BOOST_REQUIRED_VERSION=1.36.1

##==== method 3 ====##
BOOST_REQUIRE([1.36])

BOOST_DEFUN([Bind],[BOOST_FIND_HEADER([boost/accumulators/accumulators.hpp])])  #BOOST_ACCUMULATORS
BOOST_FORMAT
BOOST_DEFUN([Integer],[BOOST_FIND_HEADER([boost/cstdint.hpp])])                 #BOOST_INTEGER
BOOST_DEFUN([Lexical_Cast],[BOOST_FIND_HEADER([boost/lexical_cast.hpp])])       #BOOST_LEXICAL_CAST
BOOST_DEFUN([Pool],[BOOST_FIND_HEADER([boost/pool/pool.hpp])])                  #BOOST_POOL
BOOST_DEFUN([Smart_Ptr],[BOOST_FIND_HEADER([boost/intrusive_ptr.hpp])])         #BOOST_SMARTPTR

# ([mt-d])
BOOST_DATE_TIME
BOOST_PROGRAM_OPTIONS

## TODO: if it fails above, it will be added with the option --with-boost 

##==== method 1 ====##
#### two methods to check libraries
#### Most modern libraries use pkg-config to provide includes and linker information
#### if pkg has not been installed, (ac_arg_with) will be used. 
#echo "o--------------------------------------------------------------"
#echo "| Checking whether pkg-config has been installed or not"
#echo "o--------------------------------------------------------------"
#AC_CHECK_PROG([HAVE_PKG_CONFIG], [pkg-config], [yes],[no])
#### check for library version using pkg-config if present
#if test "x$HAVE_PKG_CONFIG" = "xyes"; then    
#    PKG_PROG_PKG_CONFIG
#    AC_MSG_NOTICE(check if BOOST version >= $BOOST_REQUIRED_VERSION)
#    PKG_CHECK_MODULES([BOOST], [boost >= $BOOST_REQUIRED_VERSION])
#else
  #### Get Boost library and include locations
#  AC_ARG_WITH([boost-include-path],
#    [AS_HELP_STRING([--with-boost-include-path],
#      [location of the Boost headers, defaults to /usr/include])],
#    [BOOST_CFLAGS="-I$withval"],
#    [BOOST_CFLAGS='-I/usr/include'])
#  AC_SUBST([BOOST_CFLAGS])

#  AC_ARG_WITH([boost-lib-path],
#    [AS_HELP_STRING([--with-boost-lib-path], [location of the Boost libraries])],
#    [BOOST_LIBS="-L$withval -l/usr/lib64 -l/usr/lib"],
#    [BOOST_LIBS='-l/usr/lib64 -l/usr/lib'])
#  AC_SUBST([BOOST_LIBS])
#fi

##==== method 2 ====##
#### Search for specified function in specified libraries
#AC_SEARCH_LIBS( [BK_set_key],
#[boost],
#[],
#[AC_MSG_ERROR( [Cannot find the libboost library, needed to build
#this package.] )] )
#### ==> corresponding point in Makefile.am: #RNABarrier_LDFLAGS = @LDFLAGS@
#### reference page: http://www.velocityreviews.com/forums/t441986-autoconf-and-shared-library.html

# error ABORT if on of the libraries was not found
if test "$DEPENDENCYNOTFOUND" = "1"; then
    AC_MSG_NOTICE()
    AC_MSG_ERROR(library not found! Try './configure --help'.)
fi


echo "o------------------------------------------------------------"
echo "| Results are"
echo "o------------------------------------------------------------"
echo "| CC:                             $CC"
echo "| CFLAGS:                         $CFLAGS"
echo "| CPP:                            $CPP"
echo "| CXXCPP:                         $CXXCPP"
echo "| CPPFLAGS:                       $CPPFLAGS"
echo "| CXX:                            $CXX"
echo "| CXXFLAGS:                       $CXXFLAGS"
echo "| LDFLAGS:                        $LDFLAGS"
##==== method 1 ====##
#echo "| BOOST_CFLAGS:                   $BOOST_CFLAGS"
#echo "| BOOST_LIBS:                     $BOOST_LIBS"
##==== method 3 ====##
echo "| BOOST_CPPFLAGS:                 $BOOST_CPPFLAGS"
echo "| BOOST_DATE_TIME_LDFLAGS:        $BOOST_DATE_TIME_LDFLAGS"
echo "| BOOST_PROGRAM_OPTIONS_LDFLAGS:  $BOOST_PROGRAM_OPTIONS_LDFLAGS"
echo "| BOOST_DATE_TIME_LIBS:           $BOOST_DATE_TIME_LIBS"
echo "| BOOST_PROGRAM_OPTIONS_LIBS:     $BOOST_PROGRAM_OPTIONS_LIBS"
echo "| LIBS:                           $LIBS"
echo "o------------------------------------------------------------"

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(malloc.h strings.h unistd.h)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_C_INLINE([])

# Checks for functions.
AC_CHECK_FUNCS([setenv])

dnl Checks for library functions.
AC_CHECK_FUNCS(strdup strstr strchr erand48)


# And the end.  
AC_CONFIG_FILES([Makefile
                 src/librna/Makefile
                 src/rtlib/Makefile
                 src/comp/Makefile
                 src/libs/Makefile
                 src/Makefile
                 examples/Makefile
                 ])
AC_OUTPUT

AC_MSG_NOTICE([
Configuration completed. Type "make" to compile and "make install" to install.
])
