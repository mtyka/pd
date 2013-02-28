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

#if _DEBUG

#ifndef DEBUG_TOOLS_H
#define DEBUG_TOOLS_H

#include "maths/maths.fwd.h"

class PD_API Particle;

//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API DebugTools
{
public:
	static void DebugTools::printDebugAxis( int &atomNumberStart, int &molNumber );
	static void DebugTools::printDebugOrigin( int &atomNumberStart, int &molNumber );
	static void DebugTools::printDebugAtomLine( Particle *atom, int index );

	template <class T>
	static void printDebugAtomLine( int AtomNumber, char element, int MolNumber, Maths::Tvector<T> *v )
	{
		printDebugAtomLine( AtomNumber, element, MolNumber, v->x, v->y, v->z );
	}

	template <class T>
	static void printDebugAtomLine( int AtomNumber, char element[4], int MolNumber, Maths::Tvector<T> *v )
	{
		printDebugAtomLine( AtomNumber, element, MolNumber, v->x, v->y, v->z );
	}

	template <class T>
	static void printDebugAtomLine(int AtomNumber, char element[4], int MolNumber, T x, T y, T z )
	{
		printf("ATOM %5d %c%c%c%c ALA %3d %8.3lf%8.3lf%8.3lf%\n",
			AtomNumber,
			element[0],
			element[1],
			element[2],
			element[3],
			MolNumber,
			x, y, z);
	}

	template <class T>
	static void printDebugAtomLine(int AtomNumber, char element, int MolNumber, T x, T y, T z )
	{
		printf("ATOM  %5d  %c   MOL   %3d    %8.3f%8.3f%8.3f%\n",
			AtomNumber,
			element,
			MolNumber,
			x, y, z);
	}


private:
	DebugTools()
	{
	}
};

#endif // DEBUG_TOOLS_H

#endif // _DEBUG
