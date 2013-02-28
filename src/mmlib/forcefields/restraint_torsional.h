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

#ifndef __FFRESTRAINT_TORSIONAL_H
#define __FFRESTRAINT_TORSIONAL_H

#include "forcefields/forcefield.h" // provides a base class
#include "forcefields/restraintbase.h" // provides a base class

#include "forcefields/ffbonded.h"   // provides a member variable

namespace Physics
{
	class PD_API Torsion;

	//-------------------------------------------------
	//
	/// \brief Restrains all or some of the internal torsions within a molecule 
	///
	/// \details 
	///    
	/// \author Mike Tyka  
	///
	/// \todo HALFBAKED 
	///
	/// \bug 
	///
	class PD_API FF_Restraint_Torsional : public AtomRestraintForcefieldBase
	{
	public:
		FF_Restraint_Torsional(  WorkSpace &newwspace );
		virtual FF_Restraint_Torsional* clone() const { return new FF_Restraint_Torsional(*this); }

		virtual double calcEnergyAtQ(double newQ);
		
		/// prints a little block of parameter information
		virtual void info(); 
		
		/// prints a little block of parameter information
		virtual void detail(); 

		bool OneRestraintPerBond; 
	protected:
		virtual void setup();
		virtual void calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level);
		virtual void calcEnergies();
		virtual void calcForces();

	private:
	 std::vector< Torsion > torsion;
	};


} // namespace Physics


#endif

