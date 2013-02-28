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

#ifndef __PROTOCOL_BASE_H
#define __PROTOCOL_BASE_H

// Essential Headers
#include "tools/cloneholder.h" // Provides a class member
#include "monitors/monitorbase.h" // Provides a class member
#include "pickers/basicpickers.h" // Provides a class member
#include "workspace/pospointer.h" // Provides a class member
#include "workspace/componentbase.h" // Provides a class member
#include "workspace/workspace.fwd.h"
#include "forcefields/forcefield.fwd.h"
#include "protocols/temperature.h"

namespace Protocol
{
	/// \class Protocol
	/// \brief Base class to all protocols
	/// \author Mike Tyka
	///
	/// \details
	/// A Protocol is a functor class that performs some sort simulation or
	/// calculation on a WorkSpace and is thus the central active class PD_API in pd.
	/// The base class is ProtocolBase which defines a common interface - most
	/// importantly the function run core() which is to be overloaded by derived
	/// classes and implements the simulation. Protocol also contains a container
	/// for Monitors which plug in to any Protocol and provide functionality to
	/// perform measurements through the simulation. Finally Protocol has timing
	/// functions that allow benchmarking. Due to the standardised interface
	/// Protocols can be combined and cascaded to create new types of Protocols.
	/// For example, the protocol MonteCarlo takes another Protocol as an argument,
	/// which can be a protocol performing a single energy calculation (Energy)
	/// leading to simple Metropolis Monte Carlo, a protocol performing a
	/// (Minimisation) leading to Monte Carlo with Minimisation (REF) or even
	/// a MolecularDynamics protocol resulting in Heterogeneous Monte Carlo
	/// (REF). In this way complex new protocols can be created at the level
	/// of the user interface without writing a single line of C++ . Because 99%
	/// of the compute time is spent in inner loops of the energy calculation
	/// itself this high-level recombination of protocols does not cause any
	/// performance loss.

	class PD_API ProtocolBase: 
		public Object,
		public WorkSpaceOperatorBase
	{
	public:
		ProtocolBase( Physics::Forcefield & _ff );
		virtual ProtocolBase* clone() const = 0;

		/// NOTE!!: Run should return either the number of Steps performed, or -1 on failure.
		virtual int run(); 
		virtual int runcore() = 0;

		/// prints a little block of parameter information	
		virtual void info() const = 0; 

		/// List the runtime types of the forcefields that this protocol uses.
		void listForcefields(); 

		virtual void printFinalStatistics(size_t steps=0);

		// -----------------------------------------------------------------------
		//  Wrapper functions for MonitorContainer 'mon' (mainly for SWIG export)
		// -----------------------------------------------------------------------

		/// Add monitors to the simulation, just pass a reference to any monitor
		virtual void addMonitor(Monitors::MonitorBase &monitor);

		/// Remove all monitors from this simulation
		virtual void clearMonitors();

		/// remove the last monitor
		virtual void popMonitor();

		/// Monitor container analysing the run - this is public for convenient access from C++.
		Monitors::MonitorContainer mon;

		inline Physics::Forcefield& getFF(){ return *ff; }

		// --------------------
		//  Public Member Data
		// --------------------
		
		/// screen/stdout update frequency (0=off)
		int UpdateScr;   
		
		/// trajectory update frequency (0=off)
		int UpdateTra;   
		
		/// monitor update frequency  (0=off)
		int UpdateMon;   
		
		/// number of Steps between full neighborlist updates;
		int UpdateNList; 
		
		/// Controls screen output
		Verbosity::Type OutputLevel;     

		/// Total Number of Steps to be done
		size_t Steps; 

		/// Temperature (only if Thermostats are active) - note this does not take a number -
		/// it takes a Temperature type, i.e. Temp, TempLinear or TempExp. These allow the
		/// user to define different temperature protocols.

		virtual void setTargetTemp( const ProfileBase & _TargetTemp ){ TargetTemp = _TargetTemp.clone(); };
		virtual void setTargetTemp( double _Temperature ){ TargetTemp = new ConstantProfile( _Temperature ); };
		double getTargetTemp( double ratio ){ return TargetTemp->get(ratio); };
		
		/// prints a line of info
		virtual void infoLine() const = 0;       

		/// prints the headers for the above function
		virtual void infoLineHeader() const = 0; 

	protected:
		CloneHolder< ProfileBase > TargetTemp;							 
	
	protected:

		/// sets default parameter values
		void settodefault();

		/// Forcefield this protocol works on
		Physics::Forcefield *ff;                 

        /// Ensures that the primary forcefield 'ff' is setup correctly
		virtual int ensureFFSetup();

        /// Allows derived classes to ensure that additional forcefields are setup correctly.
        int ensureFFSetup(Physics::Forcefield& ff);

		/// check energies are still real numbers or throw an exception
		void assertStability();                  
		
		/// run all monitors that have been aded to this protocol
		void runmonitors();
		
		/// do not rerun but re-supply the last set of values to the monitors
		void repeatmonitors();

		/// refresh the Neighbourlist in the Workspace.
		void refreshNeighborList();

		int mon_counter;
		int nlist_counter;

		// common statistics

		/// computer time at start of simulation (in secs since 1980)
		long starttime;                          
		/// computer time at end of simulation (in secs since 1980)
		long endtime;
	};

	class RangedProtocolBase : public ProtocolBase
	{
	public:
		RangedProtocolBase( Physics::Forcefield& _ff );

		/// A second base related base constructor for those protocols that can function on a sub-range of atoms
		RangedProtocolBase( Physics::Forcefield& _ff, const PickAtomRange& _newRange ); 

		inline const PickAtomRange& getRange() const { return m_Range; }
		inline void setRange(const PickAtomRange& _newRange) { m_Range = _newRange; }
		void setFullRange();

		inline size_t getStartAtom() const { return m_Range.getStartAtomIndex(); }
		inline size_t getEndAtom() const { return m_Range.getEndAtomIndex(); }
		inline size_t getNAtoms() const { return m_Range.getNAtoms(); }

	private:
		PickAtomRange m_Range;
	};

	class RangesProtocolBase : public ProtocolBase
	{
	public:
		RangesProtocolBase( Physics::Forcefield& _ff );

		/// A second base related base constructor for those protocols that can function on a sub-range of atoms
		RangesProtocolBase( Physics::Forcefield& _ff, const PickAtomRange& _newRange ); 
		RangesProtocolBase( Physics::Forcefield& _ff, const PickAtomRanges& _newRange );

		inline const PickAtomRanges& getRanges() const { return m_Range; }
		inline void setRange(const PickAtomRange& _newRange) { m_PickChangeSerial++; m_Range.clear(); m_Range.addRange(_newRange); }
		void setRange(const PickAtomRanges& _newRange);
		void setFullRange();

		inline size_t getStartAtom(size_t i) const { return m_Range.getStartAtomIndex(i); }
		inline size_t getEndAtom(size_t i) const { return m_Range.getEndAtomIndex(i); }
		inline size_t getNAtoms(size_t i) const { return m_Range.getNAtoms(i); }
		inline bool getReversed(size_t i) const { return m_Range.getReversed(i); }

		size_t getPickerSerial() const { return m_PickChangeSerial; }

	private:
		size_t m_PickChangeSerial;
		PickAtomRanges m_Range;
	};

	class PickedProtocolBase : public ProtocolBase
	{
	public:
		PickedProtocolBase( Physics::Forcefield& _ff );

		/// A second base related base constructor for those protocols that can function on an arbritrary sub-set of atoms
		PickedProtocolBase( Physics::Forcefield& _ff, const PickBase& _Picker ); 

		void setPicking(const PickBase& _Picker);
		void setPickEverything();

		const PickBase& getPicker() const { return m_Picker.data(); }
		const PosPointer getPosPointer() const { return m_Atoms; }
		size_t getPickerSerial() const { return m_PickerSerial; }

	protected:
		PosPointer& getPosPointer() { return m_Atoms; }

	private:
		CloneHolder<PickBase> m_Picker;
		PosPointer m_Atoms;	

		/// A variable that is incremented each time the picker is changed
		size_t m_PickerSerial; 
	};


	/// Helper function which determines wether an action should happen, given the current step and an action interval.
	bool every(int _Step, int _Interval); 
}



#ifdef SWIG
%template(ObjectContainer_ProtocolBase) ObjectContainer<Protocol::ProtocolBase>;
#endif

namespace Protocol 
{
	//-------------------------------------------------
	//
	/// \brief   
	///
	/// \details 
	///    
	///
	/// \author Mike Tyka & Jon Rea 
	///
	
	
	class PD_API ProtocolSet: 
		public ProtocolBase, 
		public ObjectContainer<ProtocolBase> 
	{
	public:

		ProtocolSet( Physics::Forcefield & _ff );
		
		virtual ProtocolSet* clone() const { return new ProtocolSet(*this); }


		/// NOTE!!: Run should return either the number of Steps performed, or -1 on failure.
		virtual int run(); 
		virtual int runcore();

		/// prints a little block of parameter information	
		virtual void info() const;

		virtual void printFinalStatistics(size_t steps=0);

		/// Add monitors to the simulation, just pass a reference to any monitor
		virtual void addMonitor(Monitors::MonitorBase &monitor);

		/// Remove all monitors from this simulation
		virtual void clearMonitors();

		/// remove the last monitor
		virtual void popMonitor();

		/// Temperature (only if Thermostats are active) - note this does not take a number -
		/// it takes a Temperature type, i.e. Temp, TempLinear or TempExp. These allow the
		/// user to define different temperature protocols.

		virtual void setTargetTemp( const ProfileBase & _TargetTemp ); 
		virtual void setTargetTemp( double _Temperature );
		
		/// prints a line of info
		virtual void infoLine() const {};

		/// prints the headers for the above function
		virtual void infoLineHeader() const {}; 

	protected:
		virtual void preAdditionHook( ProtocolBase & newelement )
		{
		};

		virtual void postAdditionHook( ProtocolBase & newelement )
		{
			generateName();
		};
		
		virtual void generateName();
	private:
	};




} // namespace Protocol






#endif


