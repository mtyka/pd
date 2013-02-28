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

#include "forcefields/ffbonded.h"
#include "workspace/workspace.h"
#include "workspace/snapshot.h"
#include "workspace/neighbourlist.h"
#include "workspace/space.h"

#include "forcefields/restraint_atomdist.h"

using namespace IO;
using namespace Maths;

namespace Physics
{

	FF_Restraint_AtomDistance::FF_Restraint_AtomDistance(WorkSpace &newwspace)
		: RestraintForcefieldBase(newwspace), 
		cff(newwspace, Verbosity::Silent) 
	{
		cff.name_verbose = "Distance Rest:";
		cff.name_infoheader = "distrest";
		Dist_ij = -1;
		Atom_i = -1;
		Atom_j = -1;
	}

	void FF_Restraint_AtomDistance::info() 
	{ 
		printf(" -- Atom Distance Restraint ------------ \n");
		printf(" k:       %8.4lf kcal/mol/A^2\n", k);
		printf(" Dist_ij:  %.2f\n", Dist_ij);
		printf(" atom i:   %d\n", Atom_i);
		printf(" atom j:   %d\n", Atom_j);
	}

	void FF_Restraint_AtomDistance::calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level)
	{
		cff.calcEnergiesVerbose(level);
		epot = cff.epot;
		Q = getWSpace().calcAtomDistance(Atom_i,Atom_j);
		deviat = cff.epot/(k/PhysicsConst::Na*PhysicsConst::kcal2J);
	}

	void FF_Restraint_AtomDistance::calcEnergies()
	{
		cff.calcEnergies();
		epot = cff.epot;
		Q = getWSpace().calcAtomDistance(Atom_i,Atom_j);
		deviat = cff.epot/(k/PhysicsConst::Na*PhysicsConst::kcal2J);
	}

	void FF_Restraint_AtomDistance::calcForces()
	{
		cff.calcForces();
		epot = cff.epot;
		Q = getWSpace().calcAtomDistance(Atom_i, Atom_j);
		deviat = cff.epot/(k/PhysicsConst::Na*PhysicsConst::kcal2J);
	}

	double FF_Restraint_AtomDistance::calcEnergyAtQ( double newQ )
	{
		Q = newQ;
		epot = 0.5*(k/PhysicsConst::Na*PhysicsConst::kcal2J)*Maths::sqr(Q - Dist_ij);
		deviat = Maths::sqr(Q);
		epot = (k/PhysicsConst::Na*PhysicsConst::kcal2J)*deviat;
		return epot;
	}

	void FF_Restraint_AtomDistance::infoLine() const 
	{ 
		// prints a line of current energies
		cff.infoLine();
	}
	void FF_Restraint_AtomDistance::infoLineHeader() const 
	{ 
		// prints the headers for the above function
		cff.infoLineHeader();
	}

	void FF_Restraint_AtomDistance::setup()
	{
		// Proxies
		printf("setting up FF_Restraint_AtomDistance \n");
		if(Dist_ij < 0)
		{
			throw(ProcedureException("FF_Restraint_AtomDistance: Dist_ij must be set and larger than 0"));
		}
		if(Atom_j < 0)
		{
			throw(ProcedureException("FF_Restraint_AtomDistance: Atom_i must set and larger than 0"));
		}
		if(Atom_j < 0)
		{
			throw(ProcedureException("FF_Restraint_AtomDistance: Atom_j must set and larger than 0"));
		}
		cff.clear();
		CustomForce cf(Atom_i, Atom_j);
		cf.setToHarmonic(Dist_ij, k);

		cff.addForce(cf);
	}


} // namespace Physics

