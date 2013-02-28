BUILT_SOURCES = $(srcdir)/pd_wrap.cxx
SWIG_SOURCES = pd.i 

CXXFLAGS = -Os -fno-inline

INCLUDES = -I@top_srcdir@/src/mmlib  -I@top_srcdir@/src/pd
pkgpython_PYTHON = pd.py
pkgpyexec_LTLIBRARIES = _pd.la
_pd_la_SOURCES = accessory.cpp $(srcdir)/pd_wrap.cxx $(SWIG_SOURCES)
_pd_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
_pd_la_LDFLAGS = -module
_pd_la_LIBADD = ../mmlib/libmmlib.la

$(srcdir)/pd_wrap.cxx : $(SWIG_SOURCES)
	$(SWIG)  $(SWIG_PYTHON_OPT) -I$(top_srcdir)/src/  -I$(top_srcdir)/src/mmlib   -o $@ $<

install-exec-hook:
	echo "pd" > "$(pkgpyexecdir)/../pd.pth" 
