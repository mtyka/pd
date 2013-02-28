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

// -------------------------------------------------------------------------------------
// Description: A minimisation protocol that is very tollerant of initial steric clash.
// --------------------------------------------------------------------------------------

#ifndef __DUAL_FF_MINIMISER_H
#define __DUAL_FF_MINIMISER_H

// Essential Headers
#include "protocols/minimise.h"
#include "workspace/workspace.fwd.h"
#include "forcefields/forcefield.fwd.h"


namespace Protocol
{





//-------------------------------------------------
//
/// \brief A combine protocol which does two minimisations on different forcefields 
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
///

	class PD_API DualFFMinimiser: public Minimisation
	{
	public:
		
		DualFFMinimiser( Physics::Forcefield & _ff );
		DualFFMinimiser( Physics::Forcefield & _ff, Physics::Forcefield & _stericff);
		DualFFMinimiser( Physics::Forcefield & _ff, const PickBase& _Picker );
		DualFFMinimiser( Physics::Forcefield & _ff, Physics::Forcefield & _stericff, const PickBase& _Picker );

		virtual DualFFMinimiser* clone() const 
		{ 
			return new DualFFMinimiser(*this); 
		}


		virtual int runcore();


		/// Perform a quick Steepest Descent Minimisation first?
		int SDPreMinSteps; 

		/// Perform a steric conjugate gradients Minimisation?
		int StericMinSteps; 

		/// If the steric minimisation is enabled and doesnt get below this energy, the full minim is not performed.
		double StericKillFull; 

		/// The size of the step size in the steric forcefield minimisation
		double StericStepSize; 

		/// The energy gradient cutoff uses in the steric forcefield minimisation
		double StericSlopeCutoff; 

	protected:
		void setDefaults();
	private:
		Physics::Forcefield* stericff; // steric forcefield if one was given
	};
}

#endif

