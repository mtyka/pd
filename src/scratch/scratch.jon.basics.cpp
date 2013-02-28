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
#include "forcefields/breakablebonded.h"
#include "forcefields/nonbonded.h"
#include "forcefields/gbff.h"
#include "forcefields/lcpo.h"
#include "forcefields/ffsoftvdw.h"
#include "workspace/workspace.h"

#include "scratch.jon.basics.h"

using namespace Physics;

Forcefield createffs( WorkSpace& wspace, bool useBreakableFF, bool summary )
{
	Forcefield ff = Forcefield(wspace);

	if( useBreakableFF )
	{
		FF_BreakableBonded* bonds = new FF_BreakableBonded(wspace);
		ff.addWithOwnership( bonds ) ;
	}
	else
	{
		FF_Bonded* bonds = new FF_Bonded(wspace);
		ff.addWithOwnership( bonds ) ;
	}

	FF_SoftVDW *sff = new FF_SoftVDW(wspace);
	ff.addWithOwnership(sff);

	if( summary ) ff.printEnergySummary();

	return ff;
}

Forcefield createffts( WorkSpace& wspace, bool useBreakableFF, bool summary)
{
	Forcefield ff = Forcefield(wspace);

	if( useBreakableFF )
	{
		FF_BreakableBonded* bonds = new FF_BreakableBonded(wspace);
		bonds->DoBonds = false;
		bonds->DoAngles = false;
		ff.addWithOwnership( bonds ) ;
	}
	else
	{
		FF_Bonded* bonds = new FF_Bonded(wspace);
		bonds->DoBonds = false;
		bonds->DoAngles = false;
		ff.addWithOwnership( bonds ) ;
	}

	FF_SoftVDW *sff = new FF_SoftVDW(wspace);
	ff.addWithOwnership(sff);

	if( summary ) ff.printEnergySummary();

	return ff;
}

Forcefield createffVac(WorkSpace& wspace, bool useBreakableFF, bool summary)
{
	Forcefield ff = Forcefield(wspace);

	if( useBreakableFF )
	{
		FF_BreakableBonded* bonds = new FF_BreakableBonded(wspace);
		ff.addWithOwnership( bonds ) ;
	}
	else
	{
		FF_Bonded* bonds = new FF_Bonded(wspace);
		ff.addWithOwnership( bonds ) ;
	}

	FF_NonBonded* nb = new FF_NonBonded(wspace);
	nb->Cutoff = 12.0;
	nb->InnerCutoff = 6.0;
	ff.addWithOwnership( nb );

	if( summary ) ff.printEnergySummary();

	return ff;
}

Forcefield createff(WorkSpace& wspace, bool useBreakableFF, double dielec, bool summary )
{
	Forcefield ff = Forcefield(wspace);

	if( useBreakableFF )
	{
		FF_BreakableBonded* bonds = new FF_BreakableBonded(wspace);
		ff.addWithOwnership( bonds ) ;
	}
	else
	{
		FF_Bonded* bonds = new FF_Bonded(wspace);
		ff.addWithOwnership( bonds ) ;
	}

	FF_GeneralizedBorn_Still* gbsa = new FF_GeneralizedBorn_Still( wspace ); // used to take nb
	gbsa->FastMode = true;
	gbsa->DielectricSolute = dielec;
	ff.addWithOwnership( gbsa );

	FF_SASA_LCPO* sasa = new FF_SASA_LCPO(wspace);
	sasa->GlobalASP = 0.009;
	ff.addWithOwnership( sasa );

	if( summary ) ff.printEnergySummary();

	return ff;
}

double getMeRMS( const std::vector<Maths::dvector>& native, const std::vector<Maths::dvector>& conformer )
{
	double sumSqr = 0.0;
	for( size_t i = 0; i < native.size(); i++ )
	{
		sumSqr += native[i].sqrdist( conformer[i] );
	}
	return std::sqrt( sumSqr / (double)native.size() );
}

