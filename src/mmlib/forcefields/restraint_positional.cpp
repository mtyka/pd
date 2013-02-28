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

#include "forcefields/restraint_positional.h"

using namespace IO;
using namespace Maths;

namespace Physics
{
	void FF_Restraint_Positional::info() const
	{ 
		// prints a little block of parameter information
		ForcefieldBase::info(); 
		printf(" k:                  %8.4lf kcal/mol/A^2\n", k);
		printf(" Power:              %.0f\n", Power);
		printf(" Atoms restrained:   %d\n", m_Selection.size() );
	}



	void FF_Restraint_Positional::setup()
	{
		AtomRestraintForcefieldBase::setup();
	}

	void FF_Restraint_Positional::calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level)
	{
		calcForces();
		printf(" Cartesian Rest:%10.3lf kcal/mol\n", epot * PhysicsConst::J2kcal * PhysicsConst::Na);
	}

	void FF_Restraint_Positional::calcEnergies()
	{
		calcForces();
	}

	void FF_Restraint_Positional::calcForces()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		const double k_SI = k / (PhysicsConst::J2kcal * PhysicsConst::Na); // convert to SI units for calculation
		const double half_k_SI = k_SI * 0.5;

		// Assign epot of this class
		epot = 0.0;

		for(size_t i = 0; i < m_Selection.size(); i++)
		{	
			dvector force;
			force.diff(m_Selection.p(i), m_Selection.sp(i)); // The difference between the current and stores positions
			wspace.boundary().getClosestImage(force);
			double dist = force.mag();

			double forcemag = half_k_SI;
			for(int p = 0; p < (Power - 2); p++) 
			{
				forcemag *= dist;
			}

			double ene = forcemag;
			ene *= dist * dist;
			forcemag *= double (Power);

			force.mul(-forcemag / PhysicsConst::Angstrom);

			epot += ene;
			if(!Passive) 
			{
				wspace.ene.epot += ene;
				wspace.atom[m_Selection.iat(i)].epot += ene;
				m_Selection.f(i).add(force);
			}
		}

		deviat = 2.0 * epot / k_SI;
	}

	double FF_Restraint_Positional::calcEnergyAtQ(double newQ)
	{
		return 0.0;
	}



} // namespace Physics

