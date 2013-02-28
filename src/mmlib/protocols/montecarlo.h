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

#ifndef __MC_H
#define __MC_H

#include "filters/filterbase.h" // Provides 'FilterContainer'
#include "protocols/protocolbase.h" // ProtocolBase base class
#include "protocols/temperature.h"  // Providec a class member
#include "workspace/workspace.fwd.h"
#include "forcefields/forcefield.fwd.h"

namespace Manipulator
{
	class MoveSet;
}

namespace Protocol
{
	class Temp; // Temperature

	//-------------------------------------------------
	//
	/// \brief  General Metropolis Monte Carlo Protocol 
	///
	/// \details 
	///
	/// \author Mike Tyka & Jon Rea 
	///
	/// \todo 
	///
	/// \bug 
	///
	class PD_API MonteCarlo: public ProtocolBase
	{
	public:

		MonteCarlo(
			ProtocolBase&  _evaluator,
			Manipulator::MoveSet& _moveset
			):
		ProtocolBase(_evaluator.getFF()),    // Initiate our own wspace/ff from the passed protocol base
			evaluator(&_evaluator),         // save the pointer to the evaluator protocol
			moveset(&_moveset)
		{
			PreFilters.name = "MC Pre-Filters";
			PostFilters.name = "MC Post-Filters";
			UpdateScr=1;                 // These setting will print every step that is accepted
			UpdateScrAcc = true;
			UpdateScrRej = false;
			UpdateTraAcc = false;  
			UpdateTraRej = false;
			DoZeroGeometry = false;
			ReportFilterFailReasons = false;
			Temperature = new ConstantProfile( 300 ); 
			FinalState = LastAcc;
			reset();
		};

		virtual MonteCarlo* clone() const 
		{ 
			return new MonteCarlo(*this); 
		}

		/// The main run function
		virtual int runcore();

		/// prints a little block of parameter information
		virtual void info() const; 

		/// After each move has been made, but before evaluation of the accept() criterion, 
		/// the pre-filter collection is asked if the strcuture passes. If false, the structure is rejected outright.
		/// This allows some rapid and computationally cheap checks to be made on generated structrues.
		/// The evaluator is then run, and accept() called. If this passes, we use the post-filters, which should, i guess, ususally pass the structure,
		/// bit its is possible to have low energy structures which are known to be poor via some non-energetic criterion.

		/// Filters run before the evaluator
		FilterContainer PreFilters; 

		/// Filters run after the evaluator
		FilterContainer PostFilters; 

		/// Add pre monitors to the simulation, just pass a reference to any monitor
		void addPreFilter(FilterBase &filter);
		
		/// Add post monitors to the simulation, just pass a reference to any monitor
		void addPostFilter(FilterBase &filter);

		/// The temperature
		void setTemperature( const ProfileBase & _Temperature ){ Temperature = _Temperature.clone(); };
		void setTemperature( unsigned int ){ Temperature = new ConstantProfile( 300 ); };
private:
		CloneHolder< ProfileBase > Temperature;							 
public:

		/// print acceptances to screen ?
		bool UpdateScrAcc;  
		/// print rejections to screen ?
		bool UpdateScrRej;
		/// save acceptances to trajectory ?
		bool UpdateTraAcc;  
		/// save rejections to trajectory ?
		bool UpdateTraRej;


		/// Should the geometry be zeroed after every move?
		bool DoZeroGeometry; 

		/// Report why filters have failed
		bool ReportFilterFailReasons; 

		enum FinalStateType
		{
			Last, 
			LastAcc, 
			LowestEpot
		} FinalState;

		void reset()
		{
			acceptances = 0;
			acceptances_block = 0;
			rejections = 0;
			nonacceptances = 0;
			evaluations = 0;
			lowestEne = DBL_MAX;
		}
	protected:
		
		/// prints a line of current state
		virtual void infoLine() const;        
		
		/// prints the headers for the above function
		virtual void infoLineHeader() const;   

		/// acceptance function ( Metropolis function by default ) 
		virtual bool accept(
			double enenew, 
			double eneold
			) const;

		/// lowest energy state so far
		SnapShot lowstate;  
	
		/// last accepted state 
		SnapShot oldstate;  

	private:
		int Step;

		enum ValidationClass
		{
			Accept, 
			Reject, 
			PreFilterFail,
			PostFilterFail
		} validity;

		Protocol::ProtocolBase * evaluator;
		Manipulator::MoveSet * moveset;

		double lowestEne;

		int acceptances;
		mutable int	acceptances_block;
		int rejections;
		int nonacceptances;
		int evaluations;
	};

} // namespace 'Protocol'

#endif

