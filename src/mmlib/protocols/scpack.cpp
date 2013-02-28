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
#include "workspace/workspace.h"
#include "forcefields/forcefield.h"
#include "library/rotamerlib.h"
#include "manipulators/movebase.h"
#include "manipulators/basicmoves.h"
#include "protocols/dualffminimiser.h"
#include "protocols/montecarlo.h"
#include "scpack.h"

void PD_API MCPackSideChains( WorkSpace& wspace, Physics::Forcefield& ffs, Physics::Forcefield& ff, const Library::RotamerLibrary& rotLib )
{
	// Our moveset
	Manipulator::MoveSet moves( wspace );
	Manipulator::SidechainRotamerLibMove* rotamer = new Manipulator::SidechainRotamerLibMove( wspace, rotLib, 1.0, 1.0 );
	moves.addWithOwnership( rotamer );

	// Minimise only the sidechains
	Protocol::DualFFMinimiser eval(ff, ffs, PickSidechains() );
	eval.SDPreMinSteps = 31;
	eval.StericMinSteps = 501;
	eval.Steps = 501;
	eval.StepSize = 2E1;
	//eval.!OutputLevel = true;
	eval.UpdateScr = 50;
	//eval.UpdateTra = 1;
	//eval.UpdateMon = 10;
	eval.run();

	// Perform a montecarlo procedure
	Protocol::MonteCarlo mc( eval, moves );
	mc.Steps = 10;
	mc.UpdateTra = 1;
	mc.UpdateScr = 1;
	mc.UpdateTraAcc = true;
	//mc.UpdateTraRej = true;
	mc.FinalState = Protocol::MonteCarlo::LowestEpot;
	mc.run(); // Pack sidechains using the rotamers	
}

void PD_API MCPackSideChains( WorkSpace& wspace, Physics::Forcefield& ffs, Physics::Forcefield& ff, const std::string& rotLibPath )
{
	Library::RotamerLibrary rotLib(wspace.ffps());
	rotLib.readLib( rotLibPath );
	MCPackSideChains( wspace, ffs, ff, rotLib );
}

