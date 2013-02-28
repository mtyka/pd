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

#ifndef __FFRESTRAINT_ATOMDIST_H
#define __FFRESTRAINT_ATOMDIST_H

#include "forcefields/forcefield.h" // provides a base class
#include "forcefields/restraintbase.h" // provides a base class
#include "forcefields/ffcustom.h"   // provides a member variable

namespace Physics
{

	//-------------------------------------------------
	//
	/// \brief Simple Harmonic Restrant between two atoms 
	///
	/// \details 
	///
	/// \author Mike Tyka  
	///
	/// \todo 
	///
	/// \bug 
	///
	class PD_API FF_Restraint_AtomDistance : public RestraintForcefieldBase
	{
	public:
		FF_Restraint_AtomDistance(WorkSpace &newwspace);

		virtual FF_Restraint_AtomDistance* clone() const 
		{ 
			return new FF_Restraint_AtomDistance(*this); 
		}

		/// prints a little block of parameter information
		virtual void info(); 
		
		/// prints a little block of parameter information
		virtual void detail() { info(); }; 

		/// restraint Power (default = 2, i.e. harmonic)
		int Power;  

		/// atom index 1
		int Atom_i;  

		/// atom indices;
		int Atom_j; 

		/// equilibrium distance
		double Dist_ij; 

	private:
		virtual void setup();

		virtual void calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level);
		virtual void calcEnergies();
		virtual void calcForces();
		virtual double calcEnergyAtQ(double newQ);

		virtual void infoLine() const;
		virtual void infoLineHeader() const; // prints the headers for the above function

		FF_Custom cff;
	};


} // namespace Physics


#endif

