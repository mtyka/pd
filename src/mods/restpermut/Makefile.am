MODULENAME = restpermut 
MODULE_SOURCES = hungarian.cpp  restpermut.cpp

BUILT_SOURCES = $(srcdir)/restpermut_wrap.cxx
SWIG_SOURCES = restpermut.i 

INCLUDES = -I@top_srcdir@/src/mmlib  -I@top_srcdir@/src/restpermut
pkgpython_PYTHON = restpermut.py
pkgpyexec_LTLIBRARIES = _restpermut.la
_restpermut_la_SOURCES = $(MODULE_SOURCES) $(srcdir)/restpermut_wrap.cxx $(SWIG_SOURCES)
_restpermut_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
_restpermut_la_LDFLAGS = -module
_restpermut_la_LIBADD = $(top_srcdir)/src/mmlib/libmmlib.la

$(srcdir)/restpermut_wrap.cxx : $(SWIG_SOURCES)
	$(SWIG)  $(SWIG_PYTHON_OPT) -I$(top_srcdir)/src/  -I$(top_srcdir)/src/mmlib   -o $@ $<

