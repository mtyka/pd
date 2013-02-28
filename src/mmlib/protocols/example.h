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

#ifndef __EXAMPLE_PROTOCOL_H
#define __EXAMPLE_PROTOCOL_H

// Essential Headers
#include "protocols/protocolbase.h" // Provides the base class
#include "workspace/workspace.fwd.h"


namespace Protocol
{
	/// \brief Exemplifies how to create a new protocol and acts as a template
	///
	/// \details
	/// Detailed description
	///
	/// \author You!

	class PD_API ExampleProtocol : public ProtocolBase
	{
	public: 

		/// Protocols always take at least a pointer to a forcefield
		/// at creation time - that paramter is passed on to the 
		/// base class Protocol base, you dont need to worry about it
		ExampleProtocol( Physics::Forcefield& _ff );
		
		/// Destructor - clean up your mess in here
		virtual ~ExampleProtocol();

		/// This a so called clone function. It create a copy of this object.
		/// You don't need to worry about it, all you need to do to replace the class name.
		virtual ExampleProtocol* clone() const 
		{ 
			return new ExampleProtocol(*this); 
		}

		/// prints a little block of parameter information, the class
		/// should print its entire "state" through this call
		virtual void info() const; 
		
		
		/// this function is the main function which carries out the
		/// main loop of the simulation. It's outline implementation 
		/// is in example.cpp
		virtual int runcore(); 

	public:

		/// public data i.e. parameters. 
		
		/// Note that your Base class ProtcolBase already
		/// provides a number of public parameters which this class will inherit
		/// and that you can use directly (see protocol.h for more details)
		/// 
		/// int UpdateScr;   // screen/stdout update frequency (0=off)
		/// int UpdateTra;   // trajectory update frequency (0=off)
		/// int UpdateMon;   // monitor update frequency ( (0=off)
		/// int UpdateNList; // number of Steps between full neighborlist updates;
		/// Verbosity::Type OutputLevel;     // if true then completely silent mode (no screen output at all)
		/// size_t Steps;    // number of steps to be done

		/// Examples: (Names should be CapitalisedCamelCase !)
		/// Note that you should provide sensible default paramters !

		enum ExampleMode { Mode1, Mode2, Mode3 } RunMode;
		int ExampleParam1;
		double ExampleParam2;

	protected:

		/// Put private data here such as internal stores of current step number, state,
		/// Library of structures, etc...
		/// 
		/// e.g.:

		size_t Step;

		/// Put private member function such as helpers and bits of algorithms
		/// etc.. here. You can use virtual functions to allow other people
		/// to modify your algorithm by overloading specific functional aspects 
		/// of your algorithm. Note that calling virtual functions does carry 
		/// a small overhead.

		/// prints a line of current energies/information/stepnumber etc..
		virtual void infoLine() const;       
		
		/// prints the headers for the above function
		virtual void infoLineHeader() const; 
	};
}

#endif

