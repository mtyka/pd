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

// Example Module
//
// Shows how functionality can be added to mmlib through a module. Use this as a template
// for creating real modules.
//
#ifndef __EXAMPLE_H
#define __EXAMPLE_H

#include "forcefields/forcefield.h"

namespace Physics
{
	void exampleStandAloneFunction();






	//-------------------------------------------------
	//
	/// \brief  This is an example class which derives from a base class defined in PD 
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	///
	/// \author You! 
	///
	/// \todo STATE OF DEVELOPMENT?
	///
	/// \bug BUGS?
	///
	class ExampleModuleForcefield: public Physics::ForcefieldBase
	{
	public:
		ExampleModuleForcefield( WorkSpace &newwspace ): 
			ForcefieldBase( newwspace ) 
		{

		}

		virtual ~ExampleModuleForcefield()
		{
		}

		virtual ExampleModuleForcefield *clone() const 
		{
			return new ExampleModuleForcefield(*this);
		}

		virtual int setup(){ return 0; }
		virtual void calcForces(){}

	protected:
	};

}

#endif

