
lib_LTLIBRARIES = libmmlib.la
SUBDIRS = manipulators maths tools library workspace filters monitors fileio system sequence forcefields pickers protocols
libmmlib_la_SOURCES = constants.cpp constants.h exception.cpp exception.h funcgen.h global.cpp global.h interface.h mmlib.h object.cpp object.h primitives.cpp primitives.h revision.cpp revision.h traits.h typedefs.h verbosity.h
INCLUDES = -I@top_srcdir@/src/mmlib
libmmlib_la_LIBADD = manipulators/libmanipulators.la maths/tntjama/libtntjama.la maths/libmaths.la tools/libtools.la library/liblibrary.la workspace/libworkspace.la filters/libfilters.la monitors/libmonitors.la fileio/libfileio.la system/libsystem.la sequence/libsequence.la forcefields/libforcefields.la pickers/libpickers.la protocols/libprotocols.la
