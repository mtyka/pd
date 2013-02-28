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

#include "forcefields/restraint_torsional.h"

using namespace IO;
using namespace Maths;

namespace Physics
{

	FF_Restraint_Torsional::FF_Restraint_Torsional( WorkSpace &newwspace )
		: AtomRestraintForcefieldBase(newwspace) 
	{
		k = 0;
		OneRestraintPerBond = false; 
	}

	double FF_Restraint_Torsional::calcEnergyAtQ(double newQ)
	{
		return 0;
	}

	void FF_Restraint_Torsional::setup()
	{
		// Call your parent's setup function first !
		AtomRestraintForcefieldBase::setup();
		// Proxies
		WorkSpace& wspace = getWSpace();

		FF_Bonded ffbonded(wspace);
		ffbonded.DoBonds = false;
		ffbonded.DoAngles = false;
		ffbonded.setup();

		torsion.clear();
		int i;

		// look up the torsions relevant to the backbone
		for(i = 0; i < ffbonded.torsion.size(); i++) 
		{

			// see if all the atoms of each possible torsion match the 
			// selection
			if( !
			( 
				m_Selection.getPicker().matches( wspace.atom[ffbonded.torsion[i].j] )  &&
				m_Selection.getPicker().matches( wspace.atom[ffbonded.torsion[i].b] )  &&
				m_Selection.getPicker().matches( wspace.atom[ffbonded.torsion[i].a] )  &&
				m_Selection.getPicker().matches( wspace.atom[ffbonded.torsion[i].i] )  
			) ) continue;

			// avoid duplications
			if( OneRestraintPerBond )
			{				
				printf("OneRestraintPerBond: Removing duplicate torsions ... \n");
				size_t ti;
				for(ti=0;ti<torsion.size();ti++)
				{
					if((torsion[ti].a == ffbonded.torsion[i].a)&&
						(torsion[ti].b == ffbonded.torsion[i].b)) break;
					if((torsion[ti].a == ffbonded.torsion[i].b)&&
						(torsion[ti].b == ffbonded.torsion[i].a)) break;
				}
				if(ti<torsion.size()) continue;
			}


			// Create Torsion
			Torsion Temptors;

			Temptors = ffbonded.torsion[i];
			Temptors.gamma[0] = 
				calcTorsionAngle(
					wspace.cur.atom[ffbonded.torsion[i].i].p,
					wspace.cur.atom[ffbonded.torsion[i].a].p,
					wspace.cur.atom[ffbonded.torsion[i].b].p,
					wspace.cur.atom[ffbonded.torsion[i].j].p);
			Temptors.Vn[0] = k/(PhysicsConst::J2kcal * PhysicsConst::Na);
			Temptors.n[0] = 0;
			Temptors.terms = 1;

			torsion.push_back(Temptors);
		}

		printf("Found %d torsional restraints",torsion.size());

		needsetup = false;
	}

	void FF_Restraint_Torsional::calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level)
	{
		calcForces();
		printf(" Tors. Restraint:    %10.3lf kcal/mol\n", epot * PhysicsConst::J2kcal * PhysicsConst::Na);
		
		if( level > Summary ) { 	
			printf("\n Breakdown by torsional restraint: \n"); 
			WorkSpace& wspace = getWSpace();

			epot = 0;
			if(!Passive)
			{
				for( int i = 0; i < torsion.size(); i++) 
				{ 
					calcDihedralForcesVerbose(wspace,torsion[i],epot);
				}
				wspace.ene.epot += epot;
			}
			else
			{
				for( int i = 0; i < torsion.size(); i++) 
				{ 
					calcDihedralForcesVerbosePassive(wspace,torsion[i],epot);
				}
			}
		}

	}

	void FF_Restraint_Torsional::calcEnergies()
	{
		calcForces();
	}

	void FF_Restraint_Torsional::calcForces()
	{
		WorkSpace& wspace = getWSpace();

		epot = 0;
		if(!Passive)
		{
			for( int i = 0; i < torsion.size(); i++) 
			{ 
				calcDihedralForcesNonVerbose(wspace,torsion[i],epot);
			}
			wspace.ene.epot += epot;
		}
		else
		{
			for( int i = 0; i < torsion.size(); i++) 
			{ 
				calcDihedralForcesNonVerbosePassive(wspace,torsion[i],epot);
			}
		}
	}

	// prints a little block of parameter information
	void FF_Restraint_Torsional::info() 
	{
		printf(" k:                                %8.4lf kcal/mol/rad\n", k);
		printf(" OneRestraintPerBond:              % 4s\n", OneRestraintPerBond ? "yes" : "no" );
		printf(" Number of Torsional restraints:   %4d\n", torsion.size());
	}


	// prints a little block of parameter information
	void FF_Restraint_Torsional::detail() 
	{
		setup();
		info();
		printf("\n");
		WorkSpace& wspace = getWSpace();
		printf("\n  No  ResName  #(Atom)[Type] #(Atom)[Type] #(Atom)[Type] #(Atom)[Type] \n");
		for( int i = 0; i < torsion.size(); i++) 
		{ 
			printf("%5d %4s %5d(%4s)[%3s] %5d(%4s)[%3s] %5d(%4s)[%3s] \n ",
				i,
				wspace.atom[torsion[i].a].parentl3name.c_str(),
				torsion[i].i, wspace.atom[torsion[i].i].pdbname.c_str(),
				wspace.ffps().AtomType[wspace.atom[torsion[i].i].FFType].name.c_str(),
				torsion[i].a, wspace.atom[torsion[i].a].pdbname.c_str(),
				wspace.ffps().AtomType[wspace.atom[torsion[i].a].FFType].name.c_str(),
				torsion[i].b, wspace.atom[torsion[i].b].pdbname.c_str(),
				wspace.ffps().AtomType[wspace.atom[torsion[i].b].FFType].name.c_str());
		}

	}

} // namespace Physics

