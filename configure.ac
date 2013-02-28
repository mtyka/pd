#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(pd, 0.70, mtyka@u.washington.edu )
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE( pd, 0.7 )
AM_CONFIG_HEADER([src/config.h])
AC_CONFIG_SRCDIR([])
## We want shared libraries only
AC_DISABLE_STATIC

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_RANLIB
AC_PROG_LIBTOOL

## set compiler optimisation flags - we can use conditional here later

CXXFLAGS=" -Wno-format -Wno-write-strings -O3 -funroll-loops -fomit-frame-pointer -ffast-math -felide-constructors "

## if we're using GCC then try and guess optimal architecture flags
if test "x$GCC" = "xyes"; then
  AX_GCC_ARCHFLAG(YES)
fi



############################################################
## we need python to compile and use PD. Minimum version 2.3
AM_PATH_PYTHON(2.3)
PYTHON_DEVEL_MOD 

AC_ARG_WITH(swig,
[  --without-swig      use supplied Python interface files instead of generating them using SWIG],,with_swig=yes)
if test "$with_swig" = "no"; then
  ## test for presence of pd_wrap.cxx - if exists its ok, otherwise complain and
	## explain that this option is only avalable if the PD package is from a release 
	## not from SVN. Users wanting latest SVN code MUST have SWIG >= 1.3.28 installed.

	if test -e "src/pd/pd_wrap.cxx"; then 
		
		AC_MSG_NOTICE([***********************************************************])
		AC_MSG_NOTICE([ Compiling with --without-swig])
		AC_MSG_NOTICE([ NOTE: All functions of PD will be available but you will not])
		AC_MSG_NOTICE([       able to modify the source code. ])
		AC_MSG_NOTICE([***********************************************************])
	else
		AC_MSG_NOTICE([************************* ERROR ***************************])
		AC_MSG_NOTICE([ You cannot compile with --without-swig because ])
		AC_MSG_NOTICE([ the necessary interface files (src/pd/pd.py and src/pd/pd_wrap.cxx) ])
		AC_MSG_NOTICE([ are not present. Only release packages can be built in this ])
		AC_MSG_NOTICE([ this way. To compile PD downloaded from SVN you must install ])
		AC_MSG_NOTICE([ SWIG 1.3.28 or newer. ]) 
		AC_MSG_NOTICE([***********************************************************])
		AC_MSG_ERROR([Swig 1.3.28 or newer not present.])
	fi
	
  AC_DEFINE(PD_NO_SWIG,, [use supplied Python interface files instead of generating them using SWIG])

	## We still need Python devel- set all the python paths and the like
	AC_SUBST([SWIG_PYTHON_OPT],[-python])
	AC_SUBST([SWIG_PYTHON_CPPFLAGS],[$PYTHON_CPPFLAGS])
else
	############################################################
	## we need swig to develop new code (otherwise we'll use a "precompiled" pd_wrap.cxx)
	AC_PROG_SWIG_MOD(1.3.28)
	SWIG_ENABLE_CXX_MOD

	## Python flags
	SWIG_PYTHON_MOD
fi;



############################################################
# Add some debug info: Who is building, and on what machine?
############################################################
AC_DEFINE_UNQUOTED([BUILD_TIME],["`date`"],[Date and time for build])
AC_DEFINE_UNQUOTED([BUILD_USER],["`whoami`@`hostname`"],[User doing build])
AC_DEFINE_UNQUOTED([BUILD_MACHINE],["`uname -srm`"],[Hardware and OS version for build host])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h malloc.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_STRUCT_TM
AC_C_VOLATILE

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor getcwd gethostname memset pow sqrt strchr])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/mmlib/Makefile
                 src/mmlib/manipulators/Makefile
                 src/mmlib/maths/Makefile
                 src/mmlib/maths/tntjama/Makefile
                 src/mmlib/tools/Makefile
                 src/mmlib/library/Makefile
                 src/mmlib/workspace/Makefile
                 src/mmlib/filters/Makefile
                 src/mmlib/monitors/Makefile
                 src/mmlib/fileio/Makefile
                 src/mmlib/system/Makefile
                 src/mmlib/sequence/Makefile
                 src/mmlib/forcefields/Makefile
                 src/mmlib/pickers/Makefile
                 src/mmlib/protocols/Makefile
                 src/pd/Makefile
                 src/mods/Makefile
                 src/mods/gromacs/Makefile
                 src/mods/restpermut/Makefile
                 src/mods/example/Makefile
								])

AC_OUTPUT
