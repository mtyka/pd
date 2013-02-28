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
#include "forcefields/forcefield.h"
#include "minimise.h"
#include "dualffminimiser.h"

using namespace Physics;

namespace Protocol
{
	DualFFMinimiser::DualFFMinimiser( Physics::Forcefield & _ff ) 
		: Minimisation( _ff ), stericff(NULL)
	{
		setDefaults();
	}

	DualFFMinimiser::DualFFMinimiser( Physics::Forcefield & _ff, Physics::Forcefield & _stericff)  
		: Minimisation( _ff ), stericff(&_stericff)
	{
		setDefaults();
	}

	DualFFMinimiser::DualFFMinimiser( Physics::Forcefield & _ff, const PickBase& _Picker )  
		: Minimisation( _ff, _Picker ), stericff(NULL)
	{
		setDefaults();
	}

	DualFFMinimiser::DualFFMinimiser( Physics::Forcefield & _ff, Physics::Forcefield & _stericff, const PickBase& _Picker )
		: Minimisation( _ff, _Picker ), stericff(&_stericff)
	{
		setDefaults();
	}

	void DualFFMinimiser::setDefaults()
	{
		SDPreMinSteps = 0;
		StericMinSteps = 0;
		StericStepSize = 2E5;
		StericSlopeCutoff = -1.0;
		StericKillFull = DBL_MAX;
	}

	int DualFFMinimiser::runcore()
	{
		double unitConvStericKillFull = StericKillFull * (Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na);

		int Steps = 0;

		if( SDPreMinSteps > 0 )
		{
			// Perform a steepest descent preminimisation using the steric ff if 
			// available and the normal ff if not.
			Physics::Forcefield* ffSend = (stericff!=NULL) ? stericff : ff;
			ASSERT( ffSend != NULL, CodeException, "DualFFMinimiser internal NULL reference!");
			if( OutputLevel ) printf("Steep Pre-Minimisation:\n");
			Minimisation prepreminimise( *ffSend, getPicker() );
			prepreminimise.Algorithm = Minimisation::SteepestDescent;
			prepreminimise.Steps = SDPreMinSteps;

			// Mirror parent update params
			prepreminimise.UpdateScr = UpdateScr;
			prepreminimise.UpdateTra = UpdateTra;
			prepreminimise.UpdateMon = UpdateMon;
			prepreminimise.UpdateNList = UpdateNList;
			prepreminimise.OutputLevel = OutputLevel;

			Steps += prepreminimise.run();
		}

		if( StericMinSteps > 0 )
		{
			ASSERT( stericff != NULL, ArgumentException, "StericMinSteps can only be used if a steric forcefield has been provided");
			if( OutputLevel ) printf("Steric Minimisation:\n");
			Minimisation stericMinimise( *stericff, getPicker() );
			stericMinimise.Algorithm = Minimisation::ConjugateGradients;
			stericMinimise.Steps = StericMinSteps;
			stericMinimise.StepSize = StericStepSize;
			stericMinimise.SlopeCutoff = StericSlopeCutoff;

			// Mirror parent update params
			stericMinimise.UpdateScr = UpdateScr;
			stericMinimise.UpdateTra = UpdateTra;
			stericMinimise.UpdateMon = UpdateMon;
			stericMinimise.UpdateNList = UpdateNList;
			stericMinimise.OutputLevel = OutputLevel;

			Steps += stericMinimise.run();
		}
		else if( stericff != NULL )
		{
			printf("CODE WARNING: DualFFMinimiser stericff is defined, but StericMinSteps is '0'. No steric minimisation has been performed!\n");
		}

		ff->calcEnergies(); // trigger a single energy calc on the full forcefield
		// Then test to see we have done what we have resolved steric issues
		if( getWSpace().ene.epot < unitConvStericKillFull )
		{
			if( OutputLevel ) 
			{
				printf("Full Minimisation:\n");
			}
			Steps += Minimisation::runcore();
		}
		else
		{
			if( OutputLevel ) 				
			{
				printf("Full Minimisation: Killed as Steric did not succeed. Ene:'%8.3f', Cutoff:'%8.3f'\n",
					getWSpace().ene.epot * Physics::PhysicsConst::Na / Physics::PhysicsConst::kcal2J, StericKillFull);
			}
		}

		return Steps;
	}
} // namespace




