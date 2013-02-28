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

#ifndef __TYPE_DEFS_H
#define __TYPE_DEFS_H

// -------------------------------------------------------------
// TypeDefs
// In this file we define types that are used throughout
// various sections of PD e.g. the loopbuilder. These types can
// then be easily changed as the needs of the code change.
// e.g. a larger storage counter is required
// -------------------------------------------------------------

// -----------------------------------
// System Types
// -----------------------------------
// SIZE_MAX is the maximum value of a 'size_t'
// Unix doesnt seem to define SIZE_MAX - I tried to include <stdint.h>,
// which is where I got the define below, but it didn't work for some reason...
#ifndef SIZE_MAX
	#if __WORDSIZE == 64
		#define SIZE_MAX (18446744073709551615UL)
	#else
		#define SIZE_MAX (4294967295U)
	#endif
#endif

// -----------------------------------
// Commonly Used Tokens
// -----------------------------------
const std::string TOKEN_WHITESPACE = "\t ";  ///< tokens representing whitespace
#ifdef WIN32
const std::string TOKEN_LINEFEED   = "\r\n"; ///< tokens representing linefeeds to remove newline characters
#else
const std::string TOKEN_LINEFEED   = "\n";   ///< tokens representing linefeeds to remove newline characters
#endif

// -----------------------------------
// Typical Conformer Types
// i.e. Counting very large numbers
// -----------------------------------
typedef size_t conformer_count_type; ///< When counting conformers, you need a big storage container size_t should be fine at the moment.
typedef unsigned char conformer_type; ///< Conformers are stored as this typedef in all classes. unsigned char gives plenty of capacity for even the largest angleset.

// -----------------------------------
// Fail Bits
// -----------------------------------
// const size_t SIZE_T_FAIL = SIZE_MAX is oddly an overflow warning on 64bit windows as the MAX definition is _UI64_MAX (64bit), different to size_t which oddly seems to remain 32bit on x86_64
const size_t SIZE_T_FAIL = UINT_MAX; ///< Use the maximum value of a 'size_t' as a 'fail flag' as -1 is not available in an unsigned value Type
const conformer_type conformer_type_fail = UCHAR_MAX;

#endif

