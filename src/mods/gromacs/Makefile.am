MODULENAME = gromacs 
MODULE_SOURCES = ftocstr.c  gmx_system_xdr.c  gromacs.cpp  libxdrf.c 

BUILT_SOURCES = $(srcdir)/gromacs_wrap.cxx
SWIG_SOURCES = gromacs.i 

INCLUDES = -I@top_srcdir@/src/mmlib  -I@top_srcdir@/src/gromacs
pkgpython_PYTHON = gromacs.py
pkgpyexec_LTLIBRARIES = _gromacs.la
_gromacs_la_SOURCES = $(MODULE_SOURCES) $(srcdir)/gromacs_wrap.cxx $(SWIG_SOURCES)
_gromacs_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
_gromacs_la_LDFLAGS = -module
_gromacs_la_LIBADD = $(top_srcdir)/src/mmlib/libmmlib.la

$(srcdir)/gromacs_wrap.cxx : $(SWIG_SOURCES)
	$(SWIG)  $(SWIG_PYTHON_OPT) -I$(top_srcdir)/src/  -I$(top_srcdir)/src/mmlib   -o $@ $<

