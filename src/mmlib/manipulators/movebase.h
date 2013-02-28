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

#ifndef __MOVEBASE_H
#define __MOVEBASE_H

#include "object.h"
#include "workspace/workspace.fwd.h"


namespace Manipulator
{
	//-------------------------------------------------
	//
	/// \brief Interface Executive Base class for Monte Carlo Move-Type 
	/// changes to current workspace structure
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	///
	/// \author Mike Tyka & Jon Rea 
	class PD_API MoveBase : public Object
	{
	public:
		MoveBase(WorkSpace& _wspace);
		virtual ~MoveBase();
		virtual MoveBase* clone() const = 0;

		/// \brief This method must be overloaded by implementations of Monte Carlo
		/// moves and should execute the move it implements.
		virtual int apply() = 0;

		/// \brief This is a convenient little function that simply calls apply() rounds
		/// times and saves it in a trajectory (you must load the trajectory into 
		/// the workspace previously)
		void test(int rounds);

		WorkSpace& getWspace() { return *wspace; }
		const WorkSpace& getWspace() const { return *wspace; }

	protected:
		WorkSpace* wspace;
	};
}

#ifdef SWIG
%template(ObjectContainer_Move) ObjectContainer<Manipulator::MoveBase>;
#endif

namespace Manipulator
{
	//-------------------------------------------------
	//
	/// \brief This is a container of moves which can be used to execute multiple moves
	/// at once in a Monte Carlo simulation
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
	class PD_API MoveSet: public MoveBase, public ObjectContainer<MoveBase> 
	{
	public:
		MoveSet(WorkSpace& newwspace);
		virtual ~MoveSet();
		virtual MoveSet* clone() const;
		virtual int apply();
	};
}

#endif


