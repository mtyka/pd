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

// you must always include the precompiled header global.h
#include "global.h"

// This is needed in virtually all protocols - workspaces are after all what
// protocols work on
#include "workspace/workspace.h"
#include "forcefields/forcefield.h"

// further additional headers here
// ...

// include your "own" .h file last
#include "protocols/example.h"

// Declare any namespaces you might need to use (optional)
//using namespace Physics;
//using namespace Maths;
//using namespace Monitors;
// ...

// Protocols are all in one common namespace "Protocol"
namespace Protocol
{
	ExampleProtocol::ExampleProtocol( Physics::Forcefield & _ff):
		ProtocolBase(_ff)
	{
		// Make sure you assign sensible default parameters in the constructor!
		// e.g.
		RunMode = Mode1;
		ExampleParam1 = 100;
		ExampleParam2 = 1E15;

		// You should also initilise your private/protected member variables
		Step = 0;
	}

	ExampleProtocol::~ExampleProtocol()
	{
		// If you need to create dynamically allocated memory, deallocate it here.
		// otherwise leave destructor blank.
	}

	void ExampleProtocol::info() const
	{
		printf("ExampleProtocol: \n");
		printf("  RunMode:    " );
		switch (RunMode)
		{
			case Mode1: printf("Mode1\n"); break;
			case Mode2: printf("Mode2\n"); break;
			case Mode3: printf("Mode3\n"); break;
			default: throw(CodeException("Unknown Mode in ExampleProtocol"));
		}
		printf("  ExampleParam1: %d Units \n", ExampleParam1);
		printf("  ExampleParam2: %f Units \n", ExampleParam2);
	}

	int ExampleProtocol::runcore()
	{
		// print a header line or some other information
		info();
		infoLineHeader();

		// When the execution gets here the base class ProtocolBase will
		// have already set up various things for you:
		//
		// wspace will be a valid pointer to the workspace you'll be working on
		//        its elements are thus accssible as getWSpace().data etc..
		// getWSpace().cur  is a class of type SnapShot which contains the current
		//        atom positions of the simulation. During the course of the protocol
		//        you will be changing these according to energies/forces claculated by:
		// ff     A pointer to a forcefield that works on wspace. You will issue calls
		//        to ff->calcForces() or ff->calcEnergies() to calculate potential energies
		//        and their derivatives

		// Main simulation loop, e.g. loops of the number of steps
		// to be simulated. The protocol here is a dummy protocol and serves
		// exemplary function only.
		for(Step = 0; Step < Steps; Step++)
		{
			// refresh the neighbor list - this function is from the base class
			// ProtocolBase and will make sure the neighbor list in wspace
			// is uptodate according to the update frequence (UpdateNList)
			// that the user has set
			refreshNeighborList();

			// Calculate the forces or energies of the system
			ff->calcForces();
			// If you dont need the forces but just the energies then just
			// call ff->calcEnergies(); instead, that is sometimes faster

			// The total potential energy will we available in
			// getWSpace().ene.epot

			// Call any monitors that might have been added to you as a protocol
			runmonitors();

			// Now apply the forces/change the system/etc..
			// The current atom positions  are in getWSpace().cur.atom[AtomNumber].p.{x/y/z}
			// The current atom forces     are in getWSpace().cur.atom[AtomNumber].f.{x/y/z}
			// The current atom velocities are in getWSpace().cur.atom[AtomNumber].v.{x/y/z}
			// getWSpace().old.atom[AtomNumber].{p/f/v} can be used to store these quantities
			// from the last time step if required. Alternatively a local isntance of
			// SnapShot can be used for the purpose of storing intermediate information

			// If your protocol is to support periodic boundary conditions call this:
			getWSpace().cleanSpace(); // move any stray molecules back into simulation

			// NOTE: It is important that getWSpace().outtra.append(); is called prior to
			// infoLine(); as appending to a tra triggers the cRMS calculation, which is
			// reported by inforLine().

			// Save the coordinates in the trajectory as required
			if((OutputLevel)&&(UpdateTra > 0)&&
				 (((Step) % UpdateTra) == 0))
			{
				getWSpace().outtra.append();
			}

			// Display your infoline every so often (UpdateScr)
			if((OutputLevel)&&(UpdateScr > 0)&&
				 (((Step) % UpdateScr) == 0))
			{
				infoLine();
			}
		}

		// print some statistics if you want
		if(OutputLevel)
		{
			printFinalStatistics();
		}

		// runcore should return the number of force/energy evaluations
		// done
		return Step;
	}


	// This function prints a single line of information that is
	// printed every so often (every UpdateScr) steps to inform the user
	// about the state of the simulation.
	void ExampleProtocol::infoLine() const
	{
		// prints a line of current state of system, monitors and energies

		// Information on the state of the protocol, i.e. step number, temperature,
		// etc...
		printf("%5d",Step );

		// Information on the current data in the monitors
		mon.printCurData();

		// Information on the current energies
		ff->infoLine();

		// Newline (the previous calls don't print a newline character)
		printf("\n");
	}

	// this function should "match up" and label the columns produced by the
	// above function infoLine() and is usually called at the beginning of a run
	void ExampleProtocol::infoLineHeader() const
	{
		// prints the headers for the above function
		printf("%5s", "NStep");

		// Headers of the monitors
		mon.printHeader();

		// Headers of the forcefield components
		ff->infoLineHeader();

		// New Line
		printf("\n");
	}
}

