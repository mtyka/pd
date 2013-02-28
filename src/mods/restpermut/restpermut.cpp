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

// Pre-compiled Header
#include "global.h"

#include "tools/stringtool.h"
#include "maths/maths.h"
#include "maths/maths_vector.h"
#include "hungarian.h"
#include "workspace/space.h"
#include <valarray>

// Self Header Include Should Be Last
#include "restpermut.h"

void restpermut_cflush()
{
	fflush(stdout);
}


namespace Physics
{

	using namespace Maths;

	void Restraint_PermuteSolvent::setup()
	{
		saveCurrentAtomPositions();
	}

	void Restraint_PermuteSolvent::saveCurrentAtomPositions()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		size_t i,imol;
		atom_index.clear();
		nmolsize = -1;
		for(imol=0;imol<wspace.mol.size();imol++){
			if(wspace.mol[imol].name != MolName) continue;

			unsigned foundatoms = 0;

			// find reference atom
			for(i=wspace.mol[imol].ifirst; i<=wspace.mol[imol].ilast; i++ ){
				if(wspace.atom[i].rawname == RefAtomName) break;
			}
			if(i>wspace.mol[imol].ilast){
				throw(ArgumentException("Found a molecule with mathcing name but found no ref atom in it"));
			}
			foundatoms++;
			atom_index.push_back(i);
			// collect the rest of the atoms
			for(i=wspace.mol[imol].ifirst; i<=wspace.mol[imol].ilast; i++ ){
				if(wspace.atom[i].rawname == RefAtomName) continue; 
				atom_index.push_back(i);
				foundatoms++;
			}
			if(nmolsize < 0){
				nmolsize = foundatoms;
			}else{
				if(foundatoms != nmolsize){
					throw(ArgumentException("Found a molecule with mathching name but differing atom number"));
				}
			}
		}

		printf("Found %d molecules called %s. Molsize = %d  \n",
			atom_index.size()/nmolsize,
			MolName.c_str(),
			nmolsize );

		for(i=0;i<nmolsize;i++){
    	kmul_store.push_back(1.0);
		}

		// add the constraints positions
		atom_constraint_pos.clear();
		for(i=0;i<atom_index.size();i++){
			atom_constraint_pos.push_back(wspace.cur.atom[atom_index[i]].p);
		}

		// set 1 to 1 assignement
		resetAssignment();

		costm.resize(sqr(ref_atom_2_real_atom.size()));
		printf("Ref atoms: %d\n",ref_atom_2_real_atom.size());

		calcNewCostMatrix();
		calcForces();
		updateAssignment();
		calcForces();
	};

	void Restraint_PermuteSolvent::calcEnergies()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		int imol;
		int i;
		dvector force;
		double k_SI = k / (PhysicsConst::J2kcal * PhysicsConst::Na); // convert to SI units for calculation
		double forcemag;
		double ene,dist;

		if( (wspace.Step%max(UpdateCostmatrix,(unsigned)1))==0){
			calcNewCostMatrix();
			updateAssignment();
		}
		else if( (wspace.Step%max(UpdateAssignment,(unsigned)1))==0){
			updateCostMatrix();
			updateAssignment();
		}

		stepCounter++;
		// This is a hack that only works on molecules with two symetry swap atoms
		// like water.
		bool doSymetrySwap=true;
		//bool doSwap;
		if(Sym_2_i<0) doSymetrySwap=false;
		if(Sym_2_j<0) doSymetrySwap=false;
		if(Sym_2_i>=nmolsize) doSymetrySwap=false;
		if(Sym_2_j>=nmolsize) doSymetrySwap=false;

		deviat = 0;

		epot = 0;
		for(imol=0;imol<real_atom_2_ref_atom.size();imol++){
			if(doSymetrySwap){
				int real_index1   = imol*nmolsize + Sym_2_i;
				int assign_index1 = real_atom_2_ref_atom[imol]*nmolsize + Sym_2_i;
				int real_index2   = imol*nmolsize + Sym_2_j;
				int assign_index2 = real_atom_2_ref_atom[imol]*nmolsize + Sym_2_j;

				double cost1 = sqrdist(wspace.cur.atom[atom_index[real_index1]].p,
					atom_constraint_pos[assign_index1]) +
					sqrdist(wspace.cur.atom[atom_index[real_index2]].p,
					atom_constraint_pos[assign_index2]);
				double cost2 = sqrdist(wspace.cur.atom[atom_index[real_index1]].p,
					atom_constraint_pos[assign_index2]) +
					sqrdist(wspace.cur.atom[atom_index[real_index2]].p,
					atom_constraint_pos[assign_index1]);
				// if the swapped cost is less swap the atom positions
				// this can be done on every step since it's an operation of order N
				if(cost2 < cost1){
					dvector temp(atom_constraint_pos[assign_index2]);
					atom_constraint_pos[assign_index2] = atom_constraint_pos[assign_index1];
					atom_constraint_pos[assign_index1] = temp;	
				}

			}
			for(i=0;i<Maths::min(EneRestrictToFirst,nmolsize);i++){
				int real_index = imol*nmolsize + i;
				int assign_index = real_atom_2_ref_atom[imol]*nmolsize + i;

		//		if( wspace.old.atom[atom_index[real_index]].p == wspace.cur.atom[atom_index[real_index]].p ) continue;

				force.diff(wspace.cur.atom[atom_index[real_index]].p,
					atom_constraint_pos[assign_index]);

				wspace.boundary().getClosestImage(force);
				dist = force.mag();

				forcemag = 0.5 * k_SI * kmul_store[i];  // set the appropriate k multiplier (all 1.0 by default)
				//printf("%d %s %d %f\n",atom_index[real_index],wspace.atom[atom_index[real_index]].pdbname.c_str(),i,kmul_store[i]);
				for(int p = 0; p < (Power - 2); p++) {
					forcemag *= dist;
				}
				ene = forcemag * dist * dist;
				epot += ene;
				if(!Passive) {
					wspace.ene.epot += ene;
				}
			}

		}
		deviat = 2.0*epot/k_SI;
	}

	void Restraint_PermuteSolvent::calcForces()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		int imol;
		int i;
		dvector force;
		double k_SI = k / (PhysicsConst::J2kcal * PhysicsConst::Na); // convert to SI units for calculation
		double forcemag;
		double ene,dist;

		if( (wspace.Step%max(UpdateCostmatrix,(unsigned)1))==0){
			calcNewCostMatrix();
			updateAssignment();
		}
		else if( (wspace.Step%max(UpdateAssignment,(unsigned)1))==0){
			updateCostMatrix();
			updateAssignment();
		}

		stepCounter++;

		// This is a hack that only works on molecules with two symetry swap atoms
		// like water.
		bool doSymetrySwap=true;
		//bool doSwap;
		if(Sym_2_i<0) doSymetrySwap=false;
		if(Sym_2_j<0) doSymetrySwap=false;
		if(Sym_2_i>=nmolsize) doSymetrySwap=false;
		if(Sym_2_j>=nmolsize) doSymetrySwap=false;

		deviat = 0;

		epot = 0;
		for(imol=0;imol<real_atom_2_ref_atom.size();imol++){
			if(doSymetrySwap){
				int real_index1   = imol*nmolsize + Sym_2_i;
				int assign_index1 = real_atom_2_ref_atom[imol]*nmolsize + Sym_2_i;
				int real_index2   = imol*nmolsize + Sym_2_j;
				int assign_index2 = real_atom_2_ref_atom[imol]*nmolsize + Sym_2_j;

				double cost1 = sqrdist(wspace.cur.atom[atom_index[real_index1]].p,
					atom_constraint_pos[assign_index1]) +
					sqrdist(wspace.cur.atom[atom_index[real_index2]].p,
					atom_constraint_pos[assign_index2]);
				double cost2 = sqrdist(wspace.cur.atom[atom_index[real_index1]].p,
					atom_constraint_pos[assign_index2]) +
					sqrdist(wspace.cur.atom[atom_index[real_index2]].p,
					atom_constraint_pos[assign_index1]);
				// if the swapped cost is less swap the atom positions
				// this can be done on every step since it's an operation of order N
				if(cost2 < cost1){
					dvector temp(atom_constraint_pos[assign_index2]);
					atom_constraint_pos[assign_index2] = atom_constraint_pos[assign_index1];
					atom_constraint_pos[assign_index1] = temp;	
				}

			}
			for(i=0;i<Maths::min(EneRestrictToFirst,nmolsize);i++){
				int real_index = imol*nmolsize + i;
				int assign_index = real_atom_2_ref_atom[imol]*nmolsize + i;

				force.diff(wspace.cur.atom[atom_index[real_index]].p,
					atom_constraint_pos[assign_index]);

				wspace.boundary().getClosestImage(force);
				dist = force.mag();

				forcemag = 0.5 * k_SI * kmul_store[i];  // set the appropriate k multiplier (all 1.0 by default)
				//printf("%d %s %d %f\n",atom_index[real_index],wspace.atom[atom_index[real_index]].pdbname.c_str(),i,kmul_store[i]);
				for(int p = 0; p < (Power - 2); p++) {
					forcemag *= dist;
				}
				ene = forcemag;
				ene *= dist * dist;
				forcemag *= double (Power);

				force.mul(-forcemag / PhysicsConst::Angstrom);

				epot += ene;
				if(!Passive) {
					wspace.ene.epot += ene;
					wspace.atom[atom_index[real_index]].epot += ene;
					wspace.cur.atom[atom_index[real_index]].f.add(force);
				}
			}

		}
		deviat = 2.0*epot/k_SI;
	}

	// calculate a new costmatrix from scratch
	void Restraint_PermuteSolvent::calcNewCostMatrix()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		dvector dc;
		size_t i,j,iat;
		int count = 0;
		for(i=0;i<ref_atom_2_real_atom.size();i++){
			for(j=0;j<ref_atom_2_real_atom.size();j++){
				int costij=0;
				for(iat=0;iat<Maths::min(VoroRestrictToFirst,nmolsize);iat++){

					dc.diff(wspace.cur.atom[atom_index[i*nmolsize]+iat].p,
						atom_constraint_pos[j*nmolsize+iat]);
					wspace.boundary().getClosestImage(dc);
					costij += ( int( double(dc.innerdot() * kmul_store[iat]*100000.0) ) );  // costs are square distances
				}	
				costm[count] = costij;
				count++;

			}
		}
	}

	// calculate a new costmatrix from scratch
	void Restraint_PermuteSolvent::updateCostMatrix()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		dvector dc;
		size_t i,j,iat;
		int count = 0;
		for(i=0;i<ref_atom_2_real_atom.size();i++){
			for(j=0;j<ref_atom_2_real_atom.size();j++){
				if(costm[count] > 40.0*100000.0){ count++; continue; }
				int costij=0;
				for(iat=0;iat<Maths::min(VoroRestrictToFirst,nmolsize);iat++){

					dc.diff(wspace.cur.atom[atom_index[i*nmolsize]+iat].p,
						atom_constraint_pos[j*nmolsize+iat]);
					wspace.boundary().getClosestImage(dc);
					costij += ( int( double(dc.innerdot()* kmul_store[iat] *100000.0) ) );  // costs are square distances
				}	
				costm[count] = costij;
				count++;

			}
		}
	}

	// reassign the atoms to their anchor points such as to end up with the smallest amount 
	// of displacement
	void Restraint_PermuteSolvent::updateAssignment()
	{
		// do assignement
		LinearAssignmentJVC(
			costm, 
			ref_atom_2_real_atom,	
			real_atom_2_ref_atom
			);
	}
	
	// reassign the atoms to their anchor points such as to end up with the smallest amount 
	// of displacement
	void Restraint_PermuteSolvent::updateAssignment_Munkres()
	{
		// do assignement
		LinearAssignmentHungarian(
			costm, 
			ref_atom_2_real_atom,	
			real_atom_2_ref_atom
			);
	}

	// reassign the atoms to their anchor points such as to end up with the smallest amount 
	// of displacement
	void Restraint_PermuteSolvent::resetAssignment()
	{
		// create refatom array by recording the indices fo the refernce atoms in the 
		// serial atom array - these atom indices will be reassigned using the hungarian algorithm
		// to assing the closest assignemnt
		size_t imol = 0;
		ref_atom_2_real_atom.clear();
		real_atom_2_ref_atom.clear();

		for(size_t i=0;i<atom_index.size();i+=nmolsize){
			ref_atom_2_real_atom.push_back(imol);
			real_atom_2_ref_atom.push_back(imol);
			imol ++;
		}	
	}




} 



