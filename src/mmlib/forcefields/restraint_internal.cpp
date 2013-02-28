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

#include "forcefields/restraint_internal.h"

using namespace IO;
using namespace Maths;

namespace Physics
{

	FF_Restraint_Internal::FF_Restraint_Internal( WorkSpace &newwspace )
		: AtomRestraintForcefieldBase(newwspace) 
	{
		name = "InternalRest";
		ShortName = "IntRest";

		RestCutoff = 1000.0; // Angstrom
		DivByNumber = false; // divide k by number of pairs ??
	}

	void FF_Restraint_Internal::setup()
	{
		AtomRestraintForcefieldBase::setup();

		WorkSpace& wspace = getWSpace();
		SnapShot savewspace;
		savewspace = wspace.save();
		SnapShot reststruc = getRestraintStructure();
		getWSpace().load(reststruc);   // if a specific restraint structure was specified then load this first
		
		createRestraints();            // create restraints as defined by parameters
		
		wspace.load(savewspace);       // restore workspace
	}

	void FF_Restraint_Internal::calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level)
	{
		calcForces();
		printf(" IntRestr.      %10.3lf kcal/mol\n", epot * PhysicsConst::J2kcal * PhysicsConst::Na);
	}

	void FF_Restraint_Internal::calcEnergies()
	{
		calcForces();
	}

	void FF_Restraint_Internal::calcForces()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		SnapShotAtom *atom = wspace.cur.atom;

		int i, j; // i & j are the atom indices
		int ib; // ib is the bond counting variable
		double force_magnitude, fk, d; // force magnitude, force constant
		dvector vb; // bond dvector  j-->i
		double potential; // potential energy of the individual bond
		epot = 0;
		// calculate bonds
		for(ib = 0; ib < rest.size(); ib++) 
		{
			i = rest[ib].i; // get the atom indices of the bond atoms
			j = rest[ib].j;
			fk = rest[ib].k; // get force constant in J A-2

			vb.setTo(atom[i].p); // create the bond dvector  j-->i
			vb.sub(atom[j].p);
			d = vb.mag();        // get distance
			force_magnitude =  fk * (d - rest[ib].l) / PhysicsConst::Angstrom; //calculate force
			potential = 0.5 * fk * sqr(d - rest[ib].l);
			epot += potential;

			vb.mul(force_magnitude / d); // divide by length (unify) and times by force_magnitude(forcemagnitude)

			atom[i].f.sub(vb); // add force and reaction force to the two bond atoms
			atom[j].f.add(vb);
		}

		wspace.ene.epot += epot;

		Q = sqrt(2.0*epot/(k/PhysicsConst::Na*PhysicsConst::kcal2J));
		deviat = 2.0*epot/(k/PhysicsConst::Na*PhysicsConst::kcal2J);
	}

	double FF_Restraint_Internal::calcEnergyAtQ(double newQ)
	{
		Q = newQ;
		epot = k*Maths::sqr(Q);
		epot = epot;
		deviat = epot/k;
		return epot;
	}

	void FF_Restraint_Internal::info() const
	{
		printf(" k:                                          %8.4lf kcal/mol/A^2\n", k);
		printf(" Atoms restrained:                           %d\n", m_Selection.size() );
		printf(" Number of restraints:                       %4d \n", nrest);
		printf(" RestCutoff:                                 %8.4lf \n", RestCutoff);
		printf(" divide effective k by number of restraints?: %s \n", DivByNumber ? "yes" : "no");
	}

	int FF_Restraint_Internal::createRestraints()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		int i, j, n = 0;
		double dist;
		// delete all previous forces
		rest.clear();
		for(i = 0; i < wspace.atom.size(); i++) {
			if(!m_Picker->matches(wspace.atom[i])) continue;
			for(j = i + 1; j < wspace.atom.size(); j++) {
				if(!m_Picker->matches(wspace.atom[j])) continue;
				dist = wspace.calcAtomDistance(i, j);
				if(dist > RestCutoff) continue;

				rest.push_back( Bond(i,j,dist,k * PhysicsConst::kcal2J / PhysicsConst::Na) );
				n++;
			}
		}
		nrest = n;

		if(DivByNumber){
			for(i = 0; i < rest.size(); i++) {
				rest[i].k /= (double)nrest;
			}
		}
		return n;
	}


} // namespace Physics

