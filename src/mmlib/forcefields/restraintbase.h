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

#ifndef __RESTRAINTBASE_H
#define __RESTRAINTBASE_H

#include <vector>

#include "verbosity.h"
#include "tools/cloneholder.h"
#include "pickers/pickbase.h"          // provides a member variable
#include "workspace/pospointer.h"   // provides a member variable
#include "workspace/componentbase.h"   // provides a member variable
#include "workspace/workspace.fwd.h"
#include "workspace/snapshot.fwd.h"
#include "forcefields/forcefield.h"

namespace Protocol
{
	class ProtocolBase;
}

class Object;
class FFParamSet;
class PickBase;

namespace Physics
{


	//-------------------------------------------------
	//
	/// \brief Baseclass for all restraint-type forcefields 
	///
	/// \details 
	/// 
	/// A generic class of forcefield components that has a force constant - usually
	/// these are restraint of some sort, where the general function f(..) 
	/// takes a premultiplier
	/// k
	///
	/// epot = k*deviat = k*f(Q); where f() is the function of the restraint
	/// e.g. a harmonic restraint: f(x) = x^2;
	///
	///
	///
	/// \author Mike Tyka & Jon Rea 
	///
	/// \todo 
	///
	/// \bug 
	///
	class PD_API RestraintForcefieldBase: public ForcefieldBase 
	{
	public:
		
		RestraintForcefieldBase(WorkSpace &newwspace, double _k = 0.0);
		virtual RestraintForcefieldBase* clone() const = 0;
		virtual ~RestraintForcefieldBase(){}

		/// setup routines. You MUST call this from derived functions like
		/// this: DerivedClass::setup() { RestraintForcefieldBase::setup(); }
		virtual void setup();

		/// This returns the restraint energy at a given order parameter Q
		/// this functionality is needed by umbrella sampling approaches
		virtual double calcEnergyAtQ(double newQ) = 0;

		/// Every restraint forcefield is of the form E = kf(X)
		/// where k is the force constant.
		double k;

		// Current data (to supplement epot, inherited from ForcefieldBase)

		/// order parameter
		double Q; 

		/// deviation	
		double deviat; 

		/// set the structure to which the restraint will act. If this is not called
		/// before the restraint is used in action, the structure will be set to 
		/// what ever workspace is at the time setup() gets called.
		virtual void setRestraintStructure(const SnapShot &psp);     

		/// Save the current positions in the WorkSpace to the restraint position.
		/// After this call the restraint energy should be zero.
		virtual void saveCurrentAtomPositions();

	protected:

		const SnapShot & getRestraintStructure() const {
			if( m_RestraintStruc_set ) return m_RestraintStruc;
			else                       return getWSpace().cur;
		}
		/// Internal store of restraint positions
		SnapShot m_RestraintStruc;
		bool m_RestraintStruc_set;
	};


	
	class PD_API AtomRestraintForcefieldBase: public RestraintForcefieldBase 
	{
	public:
		
		AtomRestraintForcefieldBase(WorkSpace &newwspace, double _k = 0.0);
		virtual AtomRestraintForcefieldBase* clone() const = 0;

		/// setup routines. You MUST call this from derived functions like
		/// this: DerivedClass::setup() { RestraintForcefieldBase::setup(); }
		virtual void setup();

		/// Set the selection of atoms to be restrained by providing a Picker
		/// or combination of pickers.
		/// Example:
		/// myrest.setSelection( PickResidue( 3 ) )
		/// will restrain only residue no 3
		virtual void setSelection(const PickBase &_Picker);


		/// Save the current positions in the WorkSpace to the restraint position.
		/// After this call the restraint energy should be zero.
		virtual void saveCurrentAtomPositions();
	
	
		virtual void  detail();

	protected:

		/// Internal store for the Picker used 
		CloneHolder<PickBase> m_Picker;

		/// Internal store of restraint positions (only the picked ones)
		PosStore m_Selection;

	};



}


#endif


