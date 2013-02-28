ACLOCAL_AMFLAGS = -I m4 
AUTOMAKE_OPTIONS = foreign 
SUBDIRS = src

postbuild:
	cp src/pd/pd.py bin/
	cp src/pd/.libs/_pd.so bin/
	cp src/mmlib/.libs/libmmlib.so bin/

install-exec-hook:
	@echo ""
	@echo "PD is installed under $(prefix)."
	@echo " "
	@echo "Remember to set the PD_PARAM_PATH environment variable to point to your param directory. "
	@echo ""

mactel: all

gcc: all
	
