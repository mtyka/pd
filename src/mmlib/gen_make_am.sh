#!/bin/bash 

## remember current directory
here=$(pwd)
echo Current directory: $here

## loop over all directories (excluding .* or purgatory)
for i in $( du | grep -v "\.svn" | grep -v "\./doc" | sed 's/.*\././'  | grep -v "purgatory" | grep "./" ); do
	## print what we're working on
	echo $i/Makefile.am
	## go there (into the current directory)
	cd $i
	## blank the makefile.am
	echo > Makefile.am
	## get the current directory name (the latest branch only)
	DIR=$(pwd | sed 's/.*\///' | sed 's/ //g')
	## print the library name (mirrors the directory)
	echo noinst_LTLIBRARIES = lib$DIR.la                                          >> Makefile.am
	## do we have subdirs ? if so print them
	SUBDIRS=$(du --max-depth=1 | grep -v "\./\." | grep "\./" | sed 's/.*\.\///')
	echo SUBDIRS = $SUBDIRS                                                   >> Makefile.am 
	## find the sources in this directory
	echo lib"$DIR"_la_SOURCES = $(ls *.cpp *.h  -d 2> /dev/null | grep -v "segdistfan\|segbonded\|segrejoin" 2> /dev/null)            >> Makefile.am
	## add other shit
	echo 'INCLUDES = -I@top_srcdir@/src/mmlib'                                 >> Makefile.am


	## no go back
	cd $here
done


## write top level makefile for libmmlib.so
## as above for toplevel make file
echo > Makefile.am
DIR=$(pwd | sed 's/.*\///' | sed 's/ //g')

## except this time its a proper library (shared)
echo lib_LTLIBRARIES = libmmlib.la     >> Makefile.am

SUBDIRS=$(du --max-depth=1 | grep -v "\./doc" | grep -v "\./\." | grep "\./" | sed 's/.*\.\///' | grep -v _purg  )
echo SUBDIRS = $SUBDIRS                                                   >> Makefile.am 
echo lib"$DIR"_la_SOURCES = $(ls *.cpp *.h  -d 2> /dev/null | grep -v "segdistfan\|segbonded\|segrejoin" 2> /dev/null)            >> Makefile.am
echo 'INCLUDES = -I@top_srcdir@/src/mmlib'                                 >> Makefile.am

## and we need to find and add the names of all the convenience libraries
echo lib"$DIR"_la_LIBADD = $(du  | grep -v "\./doc" | grep -v "/\." | grep "\./" | grep -v _purg | sed 's/.*\.\///' | sed 's/\w*$/&\/lib&.la/')    >> Makefile.am

