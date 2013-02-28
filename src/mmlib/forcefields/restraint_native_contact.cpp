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

#include "forcefields/restraint_native_contact.h"

using namespace IO;
using namespace Maths;

namespace Physics
{


	FF_Restraint_NativeContact::FF_Restraint_NativeContact( WorkSpace &newwspace )
		: AtomRestraintForcefieldBase( newwspace ) 
	{
		name = "NatContactRestraint";
		ShortName = "NCRest";
		k = 0;
		NativeDist = 7.0;
		Steepness = 5.0;
		p0 = 1.0;
	}

	void FF_Restraint_NativeContact::setup()
	{
		AtomRestraintForcefieldBase::setup();

		WorkSpace& wspace = getWSpace();
		SnapShot savewspace;
		savewspace = wspace.save();
		SnapShot reststruc = getRestraintStructure();
		getWSpace().load(reststruc);   // if a specific restraint structure was specified then load this first
		
		createContacts(); // create restraints as defined by parameters
		
		wspace.load(savewspace);       // restore workspace
	}

	double FF_Restraint_NativeContact::calcEnergyAtQ(double newQ)
	{
		return 0.0;
	}

	void FF_Restraint_NativeContact::info() const
	{
		AtomRestraintForcefieldBase::info();
		printf(" k:            %8.4lf kcal/mol/A^2\n", k);
		printf(" NativeDist:   %4lf \n", NativeDist);
		printf(" Steepness:    %4lf \n", Steepness);
		printf(" p0:           %4lf \n", p0);
		printf(" contacts:     %4d \n",  contact.size());
	}

	void FF_Restraint_NativeContact::detail() const
	{
		info();
		const WorkSpace& wspace = RestraintForcefieldBase::getWSpace();
		for(size_t icontact=0;icontact<contact.size();icontact++)
		{
			printf("NC: %d  %4d %4s %4d %4s \n", icontact,
				contact[icontact].i,
				wspace.atom[contact[icontact].i].rawname.c_str(),
				contact[icontact].j,
				wspace.atom[contact[icontact].j].rawname.c_str());
		}
	}

	int FF_Restraint_NativeContact::createContacts()
	{
		printf("Searching for native contacts\n");

		// Proxies
		WorkSpace& wspace = RestraintForcefieldBase::getWSpace();

		int i, j = 0;
		double dist;
		// delete all previous forces
		contact.clear();
		for(i = 0; i < wspace.atom.size(); i++) {
			
			// Does this atom need to be part of the restraint set ?
			if(!m_Picker->matches(wspace.atom[i])) continue;
			
			for(j = i + 1; j < wspace.atom.size(); j++) 
			{
				// Does this atom need to be part of the restraint set ?
				if(!m_Picker->matches(wspace.atom[j])) continue;

				// get the distance between atoms i and j - this will be assumed to be their "native" distance.
				dist = wspace.calcAtomDistance(i, j);
				if(dist > (NativeDist - (log(1.0/0.1 - 1)/Steepness)) ) continue;
				// only accept contacts that are will at least score 0.9 out of 1.0 on
				// the sigmoidal function

				//printf("Native Contact: %d %d %d \n",contact.size(), i ,j );
				contact.push_back(IndexPair(i, j));
			}
		}
		if(contact.size() == 0){
			throw(ProcedureException("No native contacts found! Are you sure there isnt an error somewhere ?"));
		}
		return 0;
	}

	void FF_Restraint_NativeContact::calcEnergies()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		size_t icontact;
		double p;
		p = 0;
		for(icontact=0;icontact<contact.size();icontact++){
			double coni;
			double Dist_ij;
			Dist_ij = wspace.cur.atom[contact[icontact].i].p.dist(
				wspace.cur.atom[contact[icontact].j].p);
			coni = 1.0/(1.0 + exp( Steepness * ( Dist_ij - NativeDist ) ) );
			p += coni;
		}
		Q = p/double(contact.size());// normalise such that the native state is at Q=1 and
		// and with 0 native contacts Q=0
		// the restraint acts as a harmonic on Q
		epot = 0.5 * k * sqr(p0 - Q) / PhysicsConst::Na * PhysicsConst::kcal2J; // k is in units kcal/mol
		wspace.ene.epot += epot;
	}

	void FF_Restraint_NativeContact::calcForces()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		size_t icontact;
		double p= 0;
		double invN = 1.0/double(contact.size());
		for(icontact=0;icontact<contact.size();icontact++)
		{
			double Dist_ij;
			Dist_ij = wspace.cur.atom[contact[icontact].i].p.dist(
				wspace.cur.atom[contact[icontact].j].p);
			double t = exp( Steepness * ( Dist_ij - NativeDist ) );
			p += 1.0/(1.0 + t);
		}
		Q = p*invN;// normalise such that the native state is at Q=1 and
		// and with 0 native contacts Q=0
		// the restraint acts as a harmonic on Q
		epot = 0.5 * k * sqr(p0 - Q) / PhysicsConst::Na * PhysicsConst::kcal2J;	
		//printf(" %f  %f  %e \n",Q,k,epot * PhysicsConst::J2kcal * PhysicsConst::Na);
		wspace.ene.epot += epot;
		double depotdQ = k * (p0 - Q) / PhysicsConst::Na * PhysicsConst::kcal2J / PhysicsConst::Angstrom;

		dvector fv;
		for(icontact=0;icontact<contact.size();icontact++)
		{
			double Dist_ij;
			unsigned i = contact[icontact].i;
			unsigned j = contact[icontact].j;
			Dist_ij = wspace.cur.atom[i].p.dist(
				wspace.cur.atom[j].p);
			double t = exp( Steepness * ( Dist_ij - NativeDist ) );
			double dQdd = invN * Steepness * t / sqr(1.0 + t);
			double force_magnitude = depotdQ * dQdd;
			fv.setTo(wspace.cur.atom[i].p);
			fv.sub(wspace.cur.atom[j].p);
			fv.mul(force_magnitude / Dist_ij);
			wspace.cur.atom[i].f.sub(fv);
			wspace.cur.atom[j].f.add(fv);
		}
	}



} // namespace Physics



