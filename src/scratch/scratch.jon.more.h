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

#ifndef __SCRACTH_JON_MORE_H
#define __SCRACTH_JON_MORE_H

#include "global.h"

#include "forcefields/forcefield.h"
#include "protocols/montecarlo.h"
#include "workspace/workspace.fwd.h"

#include "scratch.jon.basics.h"

namespace Protocol
{
	class MonteCarloLowKeep : public MonteCarlo
	{
	public:
		MonteCarloLowKeep(ProtocolBase&  _evaluator,
		Manipulator::MoveSet& _moveset);
		virtual int runcore();

	protected:
		// Redefine acceptance function
		virtual bool accept(
			double enenew,
			double eneold
		) const;

	private:
		double m_StartEne;
		mutable double m_BestEne;
	};
}

void Test_Conformer();
//void TestSoftSteric( WorkSpace& wspace );
//void TestSoftSteric();
//void Test_CoiledCoil();
//void TestPops();
//void Test_Mike();
//void main_CraigSim();
//int Tweaky(int argc, char** argv);
void main_GentleHarmonicRestraintMinimisation();
//void TempFunc();
void TestProximityGrid();
void randomRotamerApplication( bool _testSterics, bool importForeignLib );
void TheWandersOfRotamers();
void testGraphTheory();
void doMD();
void coiledCoilMake();
void zincFingerPrimaryModelBuilder();

#endif

