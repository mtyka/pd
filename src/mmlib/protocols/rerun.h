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

#ifndef __RERUN_H
#define __RERUN_H

// Essential Headers
#include "protocols/protocolbase.h" // Provides the base class
#include "fileio/intra.fwd.h" // Provides the forwards to InputTrajectory class; 
#include "forcefields/forcefield.h" // Provides the Physics::Forcefield class
#include "workspace/workspace.fwd.h"


namespace Protocol
{
	/// \brief Replays a coordinate trajectory through a workspace 
	///
	/// \details 
	///  Detailed description
	///
	/// \author Mike Tyka 

	class PD_API Rerun : public ProtocolBase
	{
	public:
		/// This constructor takes a forcefield and a trajectory and
		/// and can thus recalculate energies
		Rerun( Physics::Forcefield &_ff, InputTrajectory &_inputtra );

		/// This constructor takes a workspace and a trajectory and
		/// and can thus not calculate energies
		Rerun( WorkSpace &_wspace, InputTrajectory &_inputtra );

		/// This constructor takes a forcefield and a trajectory and
		/// and can thus recalculate energies
		Rerun( Physics::Forcefield &_ff, InputTrajectory &_inputtra, ProtocolBase& _evaluator );

		/// This constructor takes a workspace and a trajectory and
		/// and can thus not calculate energies
		Rerun( WorkSpace &_wspace, InputTrajectory &_inputtra, ProtocolBase& _evaluator );

		virtual ~Rerun(){};

		virtual Rerun* clone() const;

		virtual void info() const; 

		virtual int runcore();

	protected:
		InputTrajectory *m_InputTra;

		bool m_RecalculateEnergies;

		/// prints a line of current energies/information/stepnumber etc..
		virtual void infoLine() const;       

		/// prints the headers for the above function
		virtual void infoLineHeader() const; 

	private:
		Physics::Forcefield m_DummyFF;
		ProtocolBase* m_Evaluator;

		void setDefaults();

		int Step;
	};
}

#endif

