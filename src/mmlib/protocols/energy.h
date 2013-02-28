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

#ifndef __ENERGY_H
#define __ENERGY_H

// Essential Headers
#include "protocols/protocolbase.h" // Provides a base class
#include "workspace/snapshot.h" // Provides a class member
#include "workspace/workspace.fwd.h"
#include "forcefields/forcefield.fwd.h"

namespace Protocol
{

//-------------------------------------------------
//
/// \brief  A passive protocol that calculates the energy 
///
/// \details  runcore does nothing except for calling calcEnergy of the supplied Forcefield
///
/// \author Mike Tyka 
///
	class PD_API Energy: public ProtocolBase {
	public:
		Energy(	Physics::Forcefield & _ff):
		ProtocolBase(_ff)
		{ };

		virtual Energy* clone() const 
		{ 
			return new Energy(*this); 
		}

		virtual ~Energy(){};

		/// runs quietly, calculates energy and leaves
		virtual int runcore(); 

	private:

		/// prints a little block of parameter information
		virtual void info() const{}; 

		/// prints a line of current energies
		virtual void infoLine() const{}; 

		/// prints the headers for the above function
		virtual void infoLineHeader() const{}; 
	};


} // namespace 'Protocol'

#endif

