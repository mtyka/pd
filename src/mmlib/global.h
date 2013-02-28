// PD is a free, modular C++ library for biomolecular simulation with a 
// flexible and scriptable Python interface. 
// Copyright (C) 2003-2013 Mike Tyka and Jon Rea
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __GLOBAL_H
#define __GLOBAL_H

// #define NO_API_FLAGS // is set using the Debug_NOAPI config
// This define for a '~2x' increase in compile time performance
// when just compiling MMLib for the 'Scratch' project.
// '__declspec(dllexport)' is required for correct DLL symbol export for pd
// and similar DLL-dependent applications, has an adverse effect on compilation 
// times, and is unnececary for scratch-based MMLib library debugging...

#ifdef WIN32
	// Exporting functions and classes to DLLs requires 
	// this 'silly' define under windows. 
	#ifdef NO_API_FLAGS
		#define PD_API
	#else
		#ifndef PD_EXPORTS
			#define PD_API __declspec(dllexport) 
		#else
			#define PD_API __declspec(dllimport) 
		#endif
	#endif

	// Heavy use of the STL e.g. std::vector is really SLOW is VisStudio '05 debug mode due to
	// 'Checked Iterators' ... we can disable this to speed up debug code at the expense of safety.
	#ifndef ULTRA_DEBUG
		#define _HAS_ITERATOR_DEBUGGING 0
	#endif

#else
	// Shared libraries under linux auto-export everything 
	// so no define is needed
	#define PD_API
#endif
// define classes you'd like to be available from modules etc. like this
// class PD_API MyClass {...}
//
// define function you'd like to be available from modules etc. like this
// PD_API int MyFunc(...){...}

// Other Defines
#define DLINE printf("%s %d\n",__FILE__,__LINE__);

// Cut some crap out of VS headers (we do use winsock.h somewhere)
#define VC_EXTRALEAN
#define WIN32_LEAN_AND_MEAN

// Damn you microsoft ;-)
#define _CRT_SECURE_NO_WARNINGS

// global.h: include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently. This file should contain just system headers
// and the bare-bones of the mmLib minimal core

// ---------------------------
//  Code Structure Philosophy
// ---------------------------
// The code tree:
//  1) Each code level should depend soley on subordinate level(s)
//  2) Each code level should have a conceptual purpose
//  3) Each code level should be include complete namespace(s)

// --------------------------
//  Level-0: System Includes
// --------------------------

// General System
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <string.h>

// C-Style
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

// CPP-Style
#include <ctime>
#include <string>
#include <vector>
#include <valarray>
#include <fstream>
#include <iostream>

// ------------------------------
//  Level-1: MMLib-Base Includes
// ------------------------------

// MMLib Root Members
#include "exception.h"
#include "object.h"
#include "typedefs.h"
#include "interface.h"
#include "constants.h"
#include "verbosity.h"

// Maths Namespace
#include "maths/maths.h"
#include "maths/rms.h"
#include "maths/fastrandom.h"

// MMLib Root members requiring dvector
#include "primitives.h"

// Tools Namespace
#include "tools/io.h"
#include "tools/streamtool.h"
#include "tools/quicksort.h"
#include "tools/stringtool.h"
#include "tools/safeformat.h"
#include "tools/stringbuilder.h"

#endif


