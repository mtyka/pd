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

/*
*
* Template Numerical Toolkit (TNT)
*
* Mathematical and Computational Sciences Division
* National Institute of Technology,
* Gaithersburg, MD USA
*
*
* This software was developed at the National Institute of Standards and
* Technology (NIST) by employees of the Federal Government in the course
* of their official duties. Pursuant to title 17 Section 105 of the
* United States Code, this software is not subject to copyright protection
* and is in the public domain. NIST assumes no responsibility whatsoever for
* its use by other parties, and makes no guarantees, expressed or implied,
* about its quality, reliability, or any other characteristic.
*
*/


#ifndef TNT_SUBSCRPT_H
#define TNT_SUBSCRPT_H


//---------------------------------------------------------------------
// This definition describes the default TNT data type used for
// indexing into TNT matrices and vectors. The data type should
// be wide enough to index into large arrays. It defaults to an
// "int", but can be overriden at compile time redefining TNT_SUBSCRIPT_TYPE,
// e.g.
//
// c++ -DTNT_SUBSCRIPT_TYPE='unsigned int' ...
//
//---------------------------------------------------------------------
//

#ifndef TNT_SUBSCRIPT_TYPE
#define TNT_SUBSCRIPT_TYPE int
#endif

namespace TNT
{
	typedef TNT_SUBSCRIPT_TYPE Subscript;
} /* namespace TNT */


// () indexing in TNT means 1-offset, i.e. x(1) and A(1,1) are the
// first elements. This offset is left as a macro for future
// purposes, but should not be changed in the current release.
//
//
#define TNT_BASE_OFFSET (1)

#endif
