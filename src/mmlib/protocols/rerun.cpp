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

#include "global.h"

#include "workspace/snapshot.h"
#include "fileio/intra.h"

#include "protocols/rerun.h"

namespace Protocol
{
	Rerun::Rerun( 
		Physics::Forcefield &_ff,
		InputTrajectory &_inputtra
		)
		: ProtocolBase( _ff ),
		m_DummyFF ( _ff.getWSpace() ),  
		m_InputTra ( &_inputtra ),
		m_Evaluator(NULL)
	{
		m_RecalculateEnergies = true;
		setDefaults();
	}

	Rerun::Rerun( 
		WorkSpace &_wspace,	
		InputTrajectory &_inputtra
		)
		: m_DummyFF ( _wspace ),  // initialise this **first**
		ProtocolBase( m_DummyFF ),    // now pass m_DummyFF to the BaseClass 
		m_InputTra ( &_inputtra ),
		m_Evaluator(NULL)
	{
		m_RecalculateEnergies = false; // we dont have a forcefield so what's the point ?
		setDefaults();
	}

	Rerun::Rerun( 
		Physics::Forcefield &_ff,
		InputTrajectory &_inputtra, 
		ProtocolBase& _evaluator
		)
		: ProtocolBase( _ff ),
		m_DummyFF ( _ff.getWSpace() ),  
		m_InputTra ( &_inputtra ),
		m_Evaluator(&_evaluator)
	{
		m_RecalculateEnergies = true;
		setDefaults();
	}

	Rerun::Rerun( 
		WorkSpace &_wspace,	
		InputTrajectory &_inputtra, 
		ProtocolBase& _evaluator
		)
		: m_DummyFF ( _wspace ),  // initialise this **first**
		ProtocolBase( m_DummyFF ),    // now pass m_DummyFF to the BaseClass 
		m_InputTra ( &_inputtra ),
		m_Evaluator(&_evaluator)
	{
		m_RecalculateEnergies = false; // we dont have a forcefield so what's the point ?
		setDefaults();
	}

	void Rerun:: setDefaults()
	{
		Steps = 1;       // this is just to fool the base class which requires that Steps>0.
		Step = 0;
		
		UpdateScr = 1;   // screen/stdout update frequency (0=off)
		UpdateTra = 1;   // trajectory update frequency (0=off)
		UpdateMon = 1;   // monitor update frequency ( (0=off)
		UpdateNList = 1; // number of Steps between full neighborlist updates;
		OutputLevel = Verbosity::Silent;  // if true then completely silent mode (no screen output at all)
	}

	Rerun* Rerun::clone() const 
	{ 
		return new Rerun(*this); 
	}

	void Rerun::info() const 
	{
	} 

	int Rerun::runcore()
	{
		Step=0;
		while( !m_InputTra->isEndOfFile() )
		{
			SnapShot ss;

			// Attempt to read an entry and check if the end of file was reached.
			if( m_InputTra->readNext(ss) ) 
				break; // eof if not false

			// attempt to load the coordinates into workspace
			try
			{
				getWSpace().load( ss );
			}
			catch(ArgumentException) 
			{
				throw(IOException("The trajectory provided has a different number of atoms (" + int2str(ss.nAtoms()) +
					") then the workspace (" + int2str( getWSpace().nAtoms() ) + ")." ) );
			}
			catch(...)
			{ 
				throw; 
			}

			// now do what you have to do:
			getWSpace().Step = Step; 
			refreshNeighborList();
			ff->calcForces();

			// This allows us to cache states in a trajectory and do something to them in a second protocol
			if( m_Evaluator )
			{
				m_Evaluator->run();
			}

			runmonitors();
			// Save the coordinates in the trajectory as required
			if(every(Step,UpdateTra)) getWSpace().outtra.append();

			if((OutputLevel)&&every(Step,UpdateScr)) infoLine();

			Step++;
		}

		return Step;
	}

	void Rerun::infoLine() const 
	{ 
		printf("%7d\t", Step);	

		if(getWSpace().getVolume()>=0)
		{
			printf("%5.0lf\t",getWSpace().getVolume());
		}
		mon.printCurData();
		ff->infoLine();
		printf("\n");
	}

	void Rerun::infoLineHeader() const 
	{
		printf("%7s", "Step");
		if(getWSpace().getVolume()>=0)
		{
			printf("%9s","Volume");
		}		
		mon.printHeader();
		ff->infoLineHeader();
		printf("\n");
	}

}  //namespace Protocol


