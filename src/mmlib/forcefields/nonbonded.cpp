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

#include "forcefields/nonbonded.h" // provides base class
#include "workspace/workspace.h"
#include "workspace/neighbourlist.h"
#include "workspace/space.h"
#include "traits.h"
#include "exception.h"
#include "funcgen.h"

namespace Physics
{



	/// \details calculates both vacuum electrostatics as well as VdW forces & energies
	/// in one function for an arbitrary boundary (no boundary, periodic boundary, etc..)
	/// Basic structure:
	/// For every atom
	///   For every atom in neighbor list
	///     calculate potential energy
	///     calculate Forces along i<->j
	///
	/// VdW: A = 2.0*Maths::Maths::sqr(radius^6)/d12;
	/// B = radius^6/d6;
	/// potential energy: 4.0*epsilon*(0.5*A-B);
	/// force: -(24.0*epsilon/Dist_ij) * (A-B);
	///
	/// Electrostatics:
	/// potential energy: (1/(PhysicsConst::4pi_e0*Dielectric)) * (qi*qj)/(Dist_ij);
	/// force: -elec_potential/Dist_ij;
	///
	/// This function is a "generator" function - it provides a blue print or template
	/// vor various specialized and generic versions of it
	/// This is implemented using c++ templates and traits.
	/// The idea is when an instantiation of a particular template, certain parts
	/// of the generator function are excluded/included in the compilation process.
	/// this is achieved using traits:
	/// In the function body there will be if-block  like this one:
	///
	///	if(is_same_type<T_Space,PeriodicBox>::value){
	///   ...
	/// }
	///
	/// is_same_type<A,B> returns true when the two types are identical or false
	/// otherwise. However this result is known at **compile time** which means
	/// that when the compiler gets to the optimisation stage it sees a
	/// if(true){ ... } or if(false){ ... } 
	/// Any reaasonable compiler will thus optimise the if statement away leaving only
	/// the content of the if statement bare (if the test evaluated to true) or 
	/// nothing if the test evaluated to false.
	///
	/// This mechanism enables one to write hand-inlined code to do certain operations
	/// much faster than calling a function (as does the GenericSpace version)




	template <typename T_Space, bool verbosemode>
	void FF_NonBonded_CalcEnergies_Update_T(
		FF_NonBonded &ff, 
		WorkSpace& wspace,
		std::vector <size_t> &changed_atom,
		std::vector <size_t> &old_changed_atom,
		double &epot_reference,
		ForcefieldBase::AtomicVerbosity verbose_level = ForcefieldBase::Summary)
	{
		using namespace Maths;	

		// set up proxies to workspace to make code more readable.
		size_t natom = wspace.atom.size();             // number of atoms in workspace
		const ParticleStore& atomparam = wspace.atom;  // atom parameter array
		SnapShotAtom *atom = wspace.cur.atom;						// atom coordinate array
		const NeighbourData *fnbor = wspace.nlist().getData(); // neighborlist

		// Basic stuff
		int i, j, nj;              // i&j are atom indices, nj count through neighbor list
		dvector fv;     // force dvector
		double force_magnitude= 0;             // force magnitude
		double intvir = 0;         // part of internal virial
		double vdw_force = 0;      // individual force magnitudes
		double elec_force = 0;
		double vdw_potential = 0;  // individual potential energy contributions
		double elec_potential = 0;
		double vdw14scale = 1.0;
		double elec14scale = 1.0;

		double sqrdistij,Dist_ij;             // distance in Angstrom
		double invdistij; //,  invd6;
		const double sqrcutoff = sqr(ff.Cutoff);
		const double sqrinnercutoff = sqr(ff.InnerCutoff);

		double radiusij, epsilon, A, B;
		double qi, qj;
		const double invdielectric = 1.0 / ff.Dielectric;

		// periodic spaces
		unsigned image;
		const unsigned maximages = wspace.boundary().ncells();
		unsigned nimages;
		PeriodicBox *periodic_box;

		if(is_same_type<T_Space,PeriodicBox>::value)
		{
			periodic_box = (PeriodicBox *) &wspace.boundary();
		}

		// local temporary variables
		dvector dc;
		double atomi_radius;
		double atomi_epsilon;
		dvector fv_iatom;
		dvector fv_jatom;

		// elec switching
		const double Swidth = ff.Cutoff - ff.InnerCutoff;
		const double invSwidth = 1.0 / Swidth;
		double S;
		//double dSdd;

		// vdw switching
		const double vdwSwidth = ff.VdwCutoff - ff.VdwInnerCutoff;
		const double vdwinvSwidth = 1.0 / vdwSwidth;
		double vdwS;
		//double vdwdSdd;

		// precalculated stuff for force switching
		const double sA = 1.0/cube( sqr(ff.Cutoff) - sqr(ff.InnerCutoff) );
		const double sB = -(  cube(sqr(ff.Cutoff)) - 3.0* sqr(ff.Cutoff) * sqr(ff.Cutoff) * sqr(ff.InnerCutoff));
		const double sC = 6.0* sqr(ff.Cutoff) * sqr(ff.InnerCutoff);
		const double sD = -(sqr(ff.Cutoff) + sqr(ff.InnerCutoff));
		const double sE = 2.0/5.0;
		const double fswitch_innerV = sA * (sB * (1.0/ff.InnerCutoff) + sC * ff.InnerCutoff + sD * cube(ff.InnerCutoff) + sE * cube(ff.InnerCutoff) * sqr(ff.InnerCutoff));
		const double fswitch_cutoffV = sA * (sB * (1.0/ff.Cutoff)     + sC * ff.Cutoff      + sD * cube(ff.Cutoff)      + sE * cube(ff.Cutoff)      * sqr(ff.Cutoff));
		double eshift;

		if(ff.InnerCutoff < ff.Cutoff) eshift = 1/ff.InnerCutoff + (fswitch_innerV - fswitch_cutoffV);
		else                     eshift = 0;

		// statistics
		int totalpairs = 0;
		int pairs14 = 0;
		int vdwpairs = 0;
		int elecpairs = 0;

		// DONT set epot to zero - we need its old value
		// epot = 0; 

		//if(has changed 

		double refmul=-1;
		int    twoloopstart = 0;

		wspace.ene.epot_vdw = 0;
		wspace.ene.epot_elec = 0;


		old_changed_atom = changed_atom;
		changed_atom.clear();

		// loop over all particles
		for(i = 0; i < natom; i++){
			if( !(wspace.old.atom[i].p == wspace.cur.atom[i].p) ) changed_atom.push_back(i);
		}

		// check if THE SAME particles have move since last time
		if(old_changed_atom == changed_atom){
			// if so skip the first par tof the inner loop
			twoloopstart=1;
			ff.epot = epot_reference;
		}else{
			// else do both parts but record a new reference energy for next time
			epot_reference=ff.epot;
		}


		// loop over all changed particles
		size_t ichanged;
		for(ichanged = 0; ichanged < changed_atom.size(); ichanged++) 
		{ 
			i = changed_atom[ichanged];

			//for(i = 0; i < natom; i++){
			//	if( wspace.old.atom[i].p == wspace.cur.atom[i].p ) continue;

			atomi_radius = atomparam[i].radius;
			atomi_epsilon = atomparam[i].epsilon;
			qi = atomparam[i].charge;

			fv_iatom.zero();

			// loop over all it's neighbors
			for(nj = 0; nj < fnbor[i].n; nj++) 
			{ 
				j = fnbor[i].i[nj];
				fv_jatom.zero();
				if( (fnbor[i].Type[nj] & 127) > 1)	continue;// only get non-bonded and 1-4 neighbors but include shadows! 

				if( ff.IgnoreIntraResidue )
					if( atomparam[i].ir == atomparam[j].ir ) continue; 
				// if j has also moved, only calculate interaction once
				if( !(wspace.old.atom[j].p == wspace.cur.atom[j].p) ){
					if( (fnbor[i].Type[nj]) > 1)	continue;
				}

				vdw14scale = 1.0;
				elec14scale = 1.0;

				if(fnbor[i].Type[nj] == 1) {
					pairs14++;
					vdw14scale = ff.Vdw14Scaling;
					elec14scale = ff.Elec14Scaling;
				}

				//force_magnitude= 0; 
				//intvir = 0;
				double mul;

				mul=-1.0;
				refmul=-1.0;
				int twoloop;
				dc.diff(wspace.old.atom[j].p,wspace.old.atom[i].p);
				for(twoloop=twoloopstart;twoloop<2;twoloop++){

					if(twoloop == 1){
						dc.diff(wspace.cur.atom[j].p,wspace.cur.atom[i].p); 
						mul *= -1.0;
						refmul = 0.0;
					}

					nimages=maximages;
					if(is_same_type<T_Space,PeriodicBox>::value){
						if(maximages > 1)
							if( dc.innerdot() < sqr( periodic_box->savedist*2.0 - ff.Cutoff )){
								nimages = 1;
							}
					}

					image = 0;
imageloop:

					fv.setTo(dc);

					if(is_same_type<T_Space,PeriodicBox>::value){
						// Specialised code for PeriodicBox space boundaries
						fv.add(periodic_box->celloffset[image]);
						while(fv.x >  periodic_box->halfCellSize.x) fv.x -= periodic_box->cellSize.x;
						while(fv.x < -periodic_box->halfCellSize.x) fv.x += periodic_box->cellSize.x;
						while(fv.y >  periodic_box->halfCellSize.y) fv.y -= periodic_box->cellSize.y;
						while(fv.y < -periodic_box->halfCellSize.y) fv.y += periodic_box->cellSize.y;
						while(fv.z >  periodic_box->halfCellSize.z) fv.z -= periodic_box->cellSize.z;
						while(fv.z < -periodic_box->halfCellSize.z) fv.z += periodic_box->cellSize.z;

						sqrdistij = fv.innerdot();
						if(sqrdistij > sqrcutoff) goto imageloop_continue;			
						Dist_ij = sqrt(sqrdistij);
					}
					else if(is_same_type<T_Space,InfiniteSpace>::value){
						// Specialised code for InfiniteSpace space boundaries
						// I.e.: do nothing in infinite space - the real and imaginary 
						// images are identical. 
						sqrdistij = fv.innerdot();
						Dist_ij = sqrt(sqrdistij);
					}
					else					
					{
						// Specialised code for PeriodicBox space boundaries
						wspace.boundary().getImage(fv,image);
						sqrdistij = fv.innerdot();
						if(sqrdistij > sqrcutoff) goto imageloop_continue;			
						Dist_ij = sqrt(sqrdistij);
					}



					totalpairs++;

					invdistij = 1.0 / Dist_ij;

					// Van der Waals Force
					vdw_potential = 0;
					vdw_force = 0;
					elec_potential = 0;
					elec_force = 0;
					if(ff.DoVdw){
						if(Dist_ij < ff.VdwCutoff) {
							vdwpairs++;
							// first calculate well depth (epsilon) and zero-force-separation (radiusij)
							// Use Lorentz-Berthelot mixing rules
							// arthemic mean for radiusij (=sum of radii) and geometric mean for epsilon
							radiusij = atomi_radius + atomparam[j].radius;
							// radiusij *= 0.890898718; // div 2^(1/6)

							// two different mixing rules - at the moment only geometric mean is implemented
							epsilon = vdw14scale * atomi_epsilon * atomparam[j].epsilon;

							B = sqr(radiusij*invdistij); // * sqr(radiusij) * sqr(radiusij) * invd6;
							B *=B*B;
							A = sqr(B); 
							vdw_potential = (2.0 * epsilon) * (0.5 * A - B);

							if(Dist_ij > ff.VdwInnerCutoff) {
								vdwS = (1.0 - sqr(vdwinvSwidth * (Dist_ij - ff.VdwInnerCutoff))); // Dimension less
								vdwS = sqr(vdwS);
								vdw_potential *= vdwS;
							}

						}
					}
					// electrostatics --------------------------------

					if(ff.DoElec) {
						if(Dist_ij < ff.Cutoff) {
							qj = atomparam[j].charge;
							elecpairs++;
							if(ff.EnergySwitch){
								if(ff.DDDielectric) {   // distance dependent ff.Dielectric 
									double incdistij = ((ff.Dielectric + ff.DDDielectricAlpha * Dist_ij) * Dist_ij);
									elec_potential = PhysicsConst::econv_joule * elec14scale * (qi * qj) / incdistij;
								} else {
									elec_potential = PhysicsConst::econv_joule * invdielectric * elec14scale * (qi * qj) * invdistij;
								}

								if(Dist_ij > ff.InnerCutoff) {
									S = (1.0 - sqr(invSwidth * (Dist_ij - ff.InnerCutoff))); // Dimension less
									S = sqr(S);
									elec_potential *= S;
								}
							} else
								if(ff.ForceSwitch){
									elec_potential = PhysicsConst::econv_joule * invdielectric * elec14scale * (qi * qj);
									if(Dist_ij > ff.InnerCutoff) {
										elec_potential *=  -  ( sA * Dist_ij *
											(sB * sqr(invdistij) + 
											sC + 
											sD * sqrdistij + 
											sE * sqr(sqrdistij))
											-fswitch_cutoffV); 
									}else{
										elec_potential   *=  (invdistij - eshift);		
									}
								} else {
									elec_potential = PhysicsConst::econv_joule * invdielectric * elec14scale * (qi * qj) * invdistij;
								}
						}
					}
					// add up potentials -----------------------------


					// update energy
					ff.epot += mul*(vdw_potential + elec_potential);
					epot_reference += refmul*(vdw_potential + elec_potential);

					// an "explicit" for loop 
imageloop_continue:

					if(!is_same_type<T_Space,InfiniteSpace>::value){
						image++;
						if(image<nimages) goto imageloop;
					}

				} // twoloop (in the first sweep 'remove' the old interactions, in the second add the new ones)

			}

		}

		wspace.ene.epot += ff.epot; 
	}





















	void FF_NonBonded::calc_LJ_Coulomb_Force(
		// general
		double Dist_ij,
		double invdistij,

		// vdw
		double radiusij,
		double epsilon,

		// elec
		double qi, 
		double qj,

		// result	
		double &vdw_potential,
		double &vdw_force,
		double &elec_potential,
		double &elec_force
		){

			using namespace Maths;

			double A, B;
			const double invdielectric = 1.0 / Dielectric;

			const double sqrcutoff = sqr(Cutoff);
			const double sqrinnercutoff = sqr(InnerCutoff);

			// elec switching
			const double Swidth = Cutoff - InnerCutoff;
			const double invSwidth = 1.0 / Swidth;
			double S;
			double dSdd;

			// vdw switching
			const double vdwSwidth = VdwCutoff - VdwInnerCutoff;
			const double vdwinvSwidth = 1.0 / vdwSwidth;
			double vdwS;
			double vdwdSdd;

			// precalculated stuff for force switching
			const double sA = 1.0/cube( sqr(Cutoff) - sqr(InnerCutoff) );
			const double sB = -(  cube(sqr(Cutoff)) - 3.0* sqr(Cutoff) * sqr(Cutoff) * sqr(InnerCutoff));
			const double sC = 6.0* sqr(Cutoff) * sqr(InnerCutoff);
			const double sD = -(sqr(Cutoff) + sqr(InnerCutoff));
			const double sE = 2.0/5.0;
			const double fswitch_innerV = sA * (sB * (1.0/InnerCutoff) + sC * InnerCutoff + sD * cube(InnerCutoff) + sE * cube(InnerCutoff) * sqr(InnerCutoff));
			const double fswitch_cutoffV = sA * (sB * (1.0/Cutoff)     + sC * Cutoff      + sD * cube(Cutoff)      + sE * cube(Cutoff)      * sqr(Cutoff));
			double eshift;


			double sqrdistij = sqr(Dist_ij);

			if(InnerCutoff < Cutoff) eshift = 1/InnerCutoff + (fswitch_innerV - fswitch_cutoffV);
			else                     eshift = 0;

			// Van der Waals Force
			vdw_potential = 0;
			vdw_force = 0;
			elec_potential = 0;
			elec_force = 0;
			if(Dist_ij < VdwCutoff) {
				B = sqr(radiusij*invdistij); // * sqr(radiusij) * sqr(radiusij) * invd6;
				B *=B*B;
				A = sqr(B); 
				vdw_potential = (2.0 * epsilon) * (0.5 * A - B);
				vdw_force = -(12.0E10 * epsilon * invdistij) * (A - B);

				if(Dist_ij > VdwInnerCutoff) {
					vdwS = (1.0 - sqr(vdwinvSwidth * (Dist_ij - VdwInnerCutoff))); // Dimension less
					vdwdSdd = -4.0E10 * sqr(vdwinvSwidth) * (Dist_ij - VdwInnerCutoff) * vdwS; // Units of per meter (not per angstrom!)
					vdwS = sqr(vdwS);
					vdw_force = vdwS * vdw_force + vdw_potential * vdwdSdd;

					vdw_potential *= vdwS;
				}

			}
			// electrostatics --------------------------------

			if(Dist_ij < Cutoff) {
				if(EnergySwitch){
					if(DDDielectric) {   // distance dependent Dielectric 
						double incdistij = ((Dielectric + DDDielectricAlpha * Dist_ij) * Dist_ij);
						elec_potential = PhysicsConst::econv_joule * (qi * qj) / incdistij;
						elec_force = -elec_potential * (Dielectric + 2.0 * DDDielectricAlpha * Dist_ij) / (incdistij * PhysicsConst::Angstrom);

					} else {
						elec_potential = PhysicsConst::econv_joule * invdielectric * (qi * qj) * invdistij;
						elec_force = -elec_potential * invdistij * 1E10;
					}

					if(Dist_ij > InnerCutoff) {
						S = (1.0 - sqr(invSwidth * (Dist_ij - InnerCutoff))); // Dimension less
						dSdd = -4.0E10 * sqr(invSwidth) * (Dist_ij - InnerCutoff) * S ; // Units of per meter (not per angstrom!)
						S = sqr(S);
						elec_force = S * elec_force + elec_potential * dSdd;
						elec_potential *= S;
					}
				} else
					if(ForceSwitch){

						elec_potential = PhysicsConst::econv_joule * invdielectric * (qi * qj);
						if(Dist_ij > InnerCutoff) {
							elec_force = -1E10 * elec_potential* 
								sA * sqr(invdistij) * 
								sqr(sqrcutoff - sqrdistij)*
								(sqrcutoff  - 3.0*sqrinnercutoff + 2.0*sqrdistij);

							elec_potential *=  -  ( sA * Dist_ij *
								(sB * sqr(invdistij) + 
								sC + 
								sD * sqrdistij + 
								sE * sqr(sqrdistij))
								-fswitch_cutoffV); 
						}else{
							elec_force =  -1E10 * elec_potential* sqr(invdistij);
							elec_potential   *=  (invdistij - eshift);		
						}
					} else {
						elec_potential = PhysicsConst::econv_joule * invdielectric * (qi * qj) * invdistij;
						elec_force = -elec_potential * invdistij * 1E10;
					}
			}


	}






	/// Latest Fast & Tempalated NonBonded Calculation


	const int T_SqrtFPU = 0;                       // use built-in sqrt function
	const int T_SqrtTable = 1;                     // use sqrttable + 1 Newton-Raphson iteration (single precision)

	const int T_VdwMode_None = 0;               	 // No VDW
	const int T_VdwMode_Normal = 1;             	 // Normal VDW, no switching
	const int T_VdwMode_EnergySwitch = 2;       	 // VDW + Energy Switch

	const int T_ElecMode_None = 0;                 // No electrostatics
	const int T_ElecMode_Normal = 1;               // Normal electrostatics, no switching
	const int T_ElecMode_DDDielectric_Levitt = 2;  // As above but with DDDielectric (a la Levitt)
	const int	T_ElecMode_EnergySwitch = 3;         // Energy Switching
	const int	T_ElecMode_ForceSwitch = 4;          // Force switching

	const int T_VerboseMode_False = 0;
	const int T_VerboseMode_True = 1;
	template <
		int T_SqrtMode,
		int T_VdwMode, 
		int T_ElecMode, 
		int T_VerboseMode
	>
	void FF_NonBonded_CalcForces_T_fast2(
	FF_NonBonded      &ff, 
	NonBonded_Pack    *local_atomparam,
	Maths::dvector    *basisvector,
	WorkSpace         &wspace, 
	ForcefieldBase::AtomicVerbosity verbose_level = ForcefieldBase::Summary)
	{ 
 
		using namespace Maths;	
		// set up proxies to workspace to make code more readable.
		size_t natom = wspace.atom.size();             // number of atoms in workspace
		const ParticleStore& atomparam = wspace.atom;  // atom parameter array
		SnapShotAtom *atom = wspace.cur.atom;						// atom coordinate array
		const NeighbourData *fnbor = wspace.nlist().getData(); // neighborlist

		// Basic stuff
		int j, nbor_type, i , nj;              // i&j are atom indices, nj count through neighbor list
		dvector fv;     // force dvector
		double force_magnitude= 0.0;             // force magnitude
		double intvir = 0.0;         // part of internal virial
		double vdw_force = 0.0;      // individual force magnitudes
		double elec_force = 0.0;
		double vdw_potential = 0.0;  // individual potential energy contributions
		double elec_potential = 0.0;
		double vdw14scale = 1.0;
		double elec14scale = 1.0;

		double sqrdistij,Dist_ij;             // distance in Angstrom
		double invdistij; 
		const double sqrcutoff = sqr(ff.Cutoff);
		const double sqrinnercutoff = sqr(ff.InnerCutoff);

		double radiusij=0.0, epsilon=0.0, A, B;
		double qi, qj = DBL_MAX;
		const double invdielectric = 1.0 / ff.Dielectric;


		// local temporary variables
		double atomi_radius=0;
		double atomi_epsilon=0;
		dvector fv_iatom;
		dvector fv_jatom;

		// elec switching
		const double Swidth = ff.Cutoff - ff.InnerCutoff;
		const double invSwidth = 1.0 / Swidth;
		double S;
		double dSdd;

		// vdw switching
		const double vdwSwidth = ff.VdwCutoff - ff.VdwInnerCutoff;
		const double vdwinvSwidth = 1.0 / vdwSwidth;
		double vdwS;
		double vdwdSdd;

		// precalculated stuff for force switching
		const double sA = 1.0/cube( sqr(ff.Cutoff) - sqr(ff.InnerCutoff) );
		const double sB = -(  cube(sqr(ff.Cutoff)) - 3.0* sqr(ff.Cutoff) * sqr(ff.Cutoff) * sqr(ff.InnerCutoff));
		const double sC = 6.0* sqr(ff.Cutoff) * sqr(ff.InnerCutoff);
		const double sD = -(sqr(ff.Cutoff) + sqr(ff.InnerCutoff));
		const double sE = 2.0/5.0;
		const double fswitch_innerV = sA * (sB * (1.0/ff.InnerCutoff) + 
			sC * ff.InnerCutoff + 
			sD * cube(ff.InnerCutoff) + 
			sE * cube(ff.InnerCutoff) * sqr(ff.InnerCutoff));
		const double fswitch_cutoffV = sA * (sB * (1.0/ff.Cutoff)     + 
			sC * ff.Cutoff      + 
			sD * cube(ff.Cutoff)      + 
			sE * cube(ff.Cutoff)      * sqr(ff.Cutoff));
		double eshift;
		if(ff.InnerCutoff < ff.Cutoff) eshift = 1/ff.InnerCutoff + (fswitch_innerV - fswitch_cutoffV);
		else                     eshift = 0;

		// statistics
		int totalpairs = 0;
		int pairs14 = 0;
		int vdwpairs = 0;
		int elecpairs = 0;

		// stuff for the verbose modes
		double lastenergy_vdw = 0;
		double lastenergy_elec = 0;


		// low bondorder scaling
		double tabVdw14Scaling[8] = {0.0, 0.0, 0.0, ff.Vdw14Scaling, 1.0, 1.0, 1.0, 1.0};
		double tabElec14Scaling[8] = {0.0, 0.0, 0.0, ff.Elec14Scaling, 1.0, 1.0, 1.0, 1.0};


		// inv sqrt tables precalcs
		int     invsqrttable_bins = 256;
		double   invsqrttable_binsize;
		double   invinvsqrttable_binsize; 
		double *invsqrttable; 
		if( T_SqrtMode != T_SqrtFPU ){
			invsqrttable_binsize = sqr(ff.Cutoff+4.0)/double(invsqrttable_bins);
			invinvsqrttable_binsize = 1.0/invsqrttable_binsize; 
			invsqrttable = new double [invsqrttable_bins];
			for(i=0;i<invsqrttable_bins;i++){
				invsqrttable[i] = 1.0/sqrt(invsqrttable_binsize* (i+1) );
			}
		}



		if(T_VerboseMode==T_VerboseMode_True)
		{ 
			// only in the verbose instantiation
			printf("\nPairwise simple forces \n\n");
			if(verbose_level == ForcefieldBase::Detailed)
				printf("Pair: i nr(Name) nj nr(Name) Vdw: dist radiusij Vvdw Elec: chgi chgj Dist_ij Velec \n\n");
		}


		// initialise energies
		ff.epot = 0;
		wspace.ene.epot_vdw = 0;
		wspace.ene.epot_elec = 0;

		// make a little proxy
		int *nlistptr=&fnbor[0].i[0];
		// loop over all particles
		for(i = 0; i < natom; i++) 
		{
			// get properties of atom i
			atomi_radius  = local_atomparam[i].radius;
			atomi_epsilon = local_atomparam[i].epsilon;
			qi            = local_atomparam[i].charge;

			fv_iatom.zero();
			int const fnborn = fnbor[i].n;
			nlistptr=&fnbor[i].i[0];
			// loop over all its neighbors
			for(nj = 0; nj < fnborn; nj++) 
			{ 
				// work out the atom number and extra info such as bondorder and space vector.
				// the atom number is in the lower 24 bits, the bondorder and space vector index in the remainder
				//j         = *(nlistptr); nlistptr++; 
				//nbor_type = *(nlistptr); nlistptr++;
				j         = *(nlistptr); nlistptr++; 
				nbor_type = j>>24;
				j        &= 0x00FFFFFF; 
				if( j >= i ) break;  // ignore shadow neighbors !!

				// determine scaling factors - grab them out of a table using the bornorder (which is
				// now in the lower 4 bits of the variable nbor_type 
				vdw14scale  = tabVdw14Scaling[(nbor_type&7)]; 
				elec14scale = tabElec14Scaling[(nbor_type&7)];

				intvir = 0;
				// obtain the vector from atom i to atom j and add the periodic space shift vector
				// also read out of an array 
				fv.diff(wspace.cur.atom[j].p,wspace.cur.atom[i].p);
				fv.add( basisvector[ (nbor_type>>3)&31] );

				// get the squared distance and determine if we need to consider this atom pair
				sqrdistij = fv.innerdot();
				if( sqrdistij > sqrcutoff ){
					continue;
				}

				totalpairs++;

				// get the distance and inverse distance
				if( T_SqrtMode == T_SqrtFPU ){
					Dist_ij = sqrt(sqrdistij);
					invdistij = 1.0 / Dist_ij;
				}else{
					invdistij = invsqrttable[ int(sqrdistij *  invinvsqrttable_binsize) ];
					invdistij = 0.5*(3.0 - sqr(invdistij)*sqrdistij)*invdistij;
					Dist_ij   = sqrdistij*invdistij; 
				}

				// Van der Waals Force
				vdw_potential = 0;
				vdw_force = 0;
				if( T_VdwMode > T_VdwMode_None){
					if(Dist_ij < ff.VdwCutoff){
						radiusij = atomi_radius + local_atomparam[j].radius;
						epsilon = vdw14scale * atomi_epsilon * local_atomparam[j].epsilon;

						B = sqr(radiusij*invdistij); 
						B *=B*B;
						A = sqr(B); 
						vdw_potential = (2.0 * epsilon) * (0.5 * A - B);
						vdw_force = -(12.0E10 * epsilon * invdistij) * (A - B);

						if( T_VdwMode == T_VdwMode_EnergySwitch ){	
							if(Dist_ij > ff.VdwInnerCutoff) {
								vdwS = (1.0 - sqr(vdwinvSwidth * (Dist_ij - ff.VdwInnerCutoff))); 
								vdwdSdd = -4.0E10 * sqr(vdwinvSwidth) * (Dist_ij - ff.VdwInnerCutoff) * vdwS; 
								vdwS = sqr(vdwS);
								vdw_force = vdwS * vdw_force + vdw_potential * vdwdSdd;

								vdw_potential *= vdwS;
							}
						}
					}
				}



				if( T_ElecMode > T_ElecMode_None){
					qj = local_atomparam[j].charge;
					if( T_ElecMode == T_ElecMode_Normal )
					{
						elec_potential = PhysicsConst::econv_joule * invdielectric * elec14scale * (qi * qj) * invdistij;
						elec_force = -elec_potential * invdistij * 1E10;
					}

					if( T_ElecMode == T_ElecMode_DDDielectric_Levitt )
					{
						double incdistij = ((ff.Dielectric + ff.DDDielectricAlpha * Dist_ij) * Dist_ij);
						elec_potential = PhysicsConst::econv_joule * elec14scale * (qi * qj) / incdistij;
						elec_force = -elec_potential * (ff.Dielectric + 2.0 * ff.DDDielectricAlpha * Dist_ij) / (incdistij * PhysicsConst::Angstrom);
					}

					if( T_ElecMode == T_ElecMode_EnergySwitch )
					{

						elec_potential = PhysicsConst::econv_joule * invdielectric * elec14scale * (qi * qj) * invdistij;
						elec_force = -elec_potential * invdistij * 1E10;

						if(Dist_ij > ff.InnerCutoff) {
							S = (1.0 - sqr(invSwidth * (Dist_ij - ff.InnerCutoff))); // Dimension less
							dSdd = -4.0E10 * sqr(invSwidth) * (Dist_ij - ff.InnerCutoff) * S ; // Units of per meter (not per angstrom!)
							S = sqr(S);
							elec_force = S * elec_force + elec_potential * dSdd;
							elec_potential *= S;
						}
						//printf("%e   %e   %e \n",Dist_ij, elec_potential / qi / qj, elec_force /qi /qj );
					}

					if( T_ElecMode == T_ElecMode_ForceSwitch )
					{
						elec_potential = PhysicsConst::econv_joule * invdielectric * elec14scale * (qi * qj);
						if(Dist_ij > ff.InnerCutoff) {
							elec_force = -1E10 * elec_potential* 
								sA * sqr(invdistij) * 

								sqr(sqrcutoff - sqrdistij)*

								(sqrcutoff  - 3.0*sqrinnercutoff + 2.0*sqrdistij);

							elec_potential *=  -  ( sA * Dist_ij *
								(sB * sqr(invdistij) + 
								sC + 
								sD * sqrdistij + 
								sE * sqr(sqrdistij))
								-fswitch_cutoffV); 
						}else{
							elec_force =  -1E10 * elec_potential* sqr(invdistij);
							elec_potential   *=  (invdistij - eshift);		
						}
					}

				}else{
					elec_potential = 0;
					elec_force = 0;
				}

				wspace.ene.epot += vdw_potential + elec_potential;;
				wspace.ene.epot_vdw += vdw_potential;
				wspace.ene.epot_elec += elec_potential;

				force_magnitude= vdw_force + elec_force;

				// for intermolecular forces add up virial components
				wspace.ene.InternalVirial += Dist_ij * force_magnitude* PhysicsConst::Angstrom;

				fv.mul(invdistij * force_magnitude);
				fv_iatom.add(fv);

				if( T_VerboseMode == T_VerboseMode_True){
					if(verbose_level == ForcefieldBase::Detailed) {
						printf("Pair:%5d(%4s)%5d(%4s) Vdw: %5.2lf %5.2lf %5.2lf %8.3lf Elec: % 4.2lf % 4.2lf %5.2lf % 8.3lf\n",
							i, atomparam[i].pdbname.c_str(),
							j, atomparam[j].pdbname.c_str(),
							Dist_ij,
							radiusij,
							epsilon * PhysicsConst::J2kcal * PhysicsConst::Na,
							vdw_potential * PhysicsConst::J2kcal * PhysicsConst::Na,
							qi, 
							qj, 
							Dist_ij, 
							elec_potential * PhysicsConst::J2kcal * PhysicsConst::Na);
					}
				}


				atom[j].f.sub(fv);

			}
			atom[i].f.add(fv_iatom);


			// add a vdw correction if required
			if(ff.VdwCor){
				// this correction is the analytical integral of the attractive portion of the LJ potential
				// from ff.VdwCorCutoff to infinity. usually you;d want to set ff.VdwCorCutoff == ff.VdwCutoff;

				vdw_potential = -0.5*ff.VdwCorDensity * ff.VdwCorEpsilon / (PhysicsConst::J2kcal * PhysicsConst::Na)
					* (8.0/3.0) * MathConst::PI * sqr(ff.VdwCorRadius) * 
					sqr(ff.VdwCorRadius) * sqr(ff.VdwCorRadius) /
					cube(ff.VdwCorCutoff);
				// for the purpose of vdw path assume it's non hydrogen/non hydorgen
				wspace.ene.epot += vdw_potential;
				wspace.ene.epot_vdw += vdw_potential;
				wspace.ene.epot_vdw_att += vdw_potential;
			}


			if( T_VerboseMode == T_VerboseMode_True){
				printf("Atom: %12.7lf %12.7lf %12.7lf %5d(%4s) %5d Rad: %6.4lf Eps: %6.4lf Vdw: % 8.3lf Elec: % 4.4lf % 8.3lf\n",
					atom[i].p.x, atom[i].p.y, atom[i].p.z,
					i, atomparam[i].pdbname.c_str(),
					fnbor[i].n,
					atomparam[i].radius,
					sqr(atomparam[i].epsilon) * PhysicsConst::J2kcal * PhysicsConst::Na,
					(wspace.ene.epot_vdw - lastenergy_vdw) * PhysicsConst::J2kcal * PhysicsConst::Na, 
					atomparam[i].charge, 
					(wspace.ene.epot_elec - lastenergy_elec) * PhysicsConst::J2kcal * PhysicsConst::Na);
				lastenergy_vdw = wspace.ene.epot_vdw;
				lastenergy_elec = wspace.ene.epot_elec;
			}
		}

		if(T_VerboseMode == T_VerboseMode_True){
			printf("\n");
			printf(" ----------------- -----------------\n");
			printf(" Total: %8.3lf Total: %8.3lf \n\n\n",
				wspace.ene.epot_vdw * PhysicsConst::J2kcal * PhysicsConst::Na, 
				wspace.ene.epot_elec * PhysicsConst::J2kcal * PhysicsConst::Na);
			printf(" Total number of interactions:\n");
			printf("    Total:         %8d\n",totalpairs);
			printf("    1-4:           %8d\n",pairs14);
			printf("    Van d. Waals:  %8d\n",vdwpairs);
			printf("    Electrostatic: %8d\n",elecpairs);
		}

		ff.epot = wspace.ene.epot_vdw + wspace.ene.epot_elec;
		if( T_SqrtMode != T_SqrtFPU ){
			delete [] invsqrttable;
		}

	}


	template <int a, int b, int c, int d>
	struct FF_NonBonded_CalcForces_T_fast2_wrap{
		FF_NonBonded_CalcForces_T_fast2_wrap():ptr(&FF_NonBonded_CalcForces_T_fast2<a, b, c, d>){}
		typedef void (*Tptr)(
			FF_NonBonded      &ff, 
			NonBonded_Pack    *local_atomparam,
			Maths::dvector    *basisvector,
			WorkSpace         &wspace, 
			ForcefieldBase::AtomicVerbosity verbose_level);
		Tptr ptr;
		const static int limita=T_SqrtTable; 
		const static int limitb=T_VdwMode_EnergySwitch;
		const static int limitc=T_ElecMode_ForceSwitch;
		const static int limitd=T_VerboseMode_True;
		static void overflow(){
			throw(CodeException("CODE ERROR: FF_NonBonded_CalcForces_T_fast2_wrap template requested that is out of bounds.") ); 
		}
	};




	FF_NonBonded::FF_NonBonded( WorkSpace &newwspace ): 
	ForcefieldBase( newwspace )
	{
		name = "Non-Bonded Forcefield";
		ShortName = "NonBonded";
		local_atomparam = NULL;
		settodefault();
	}

	FF_NonBonded* FF_NonBonded::clone() const 
	{ 
		return new FF_NonBonded(*this); 
	}

	FF_NonBonded::~FF_NonBonded()
	{
		delete [] local_atomparam;
	}


	/// prints a little block of parameter information
	void FF_NonBonded::info() const
	{ 
		ForcefieldBase::info(); // standard info 
		printf("Solute Dielectric:        %4.1lf\n", Dielectric);
		printf("VDW 1-4 scaling factor:   %6.3lf\n", Vdw14Scaling);
		printf("Elec 1-4 scaling factor:  %12.8lf\n", Elec14Scaling);
		printf("Distance dep. Dielectric: %s\n", DDDielectric == true ? "yes" : "no");
		if(DDDielectric)
		{
			printf("DDD form:              e(r) = %3.2lf + %3.2lf*r\n",Dielectric, DDDielectricAlpha);
		}
		printf("VDW Cutoff:               %4.1lf A\n", VdwCutoff);
		printf("VDW inner Cutoff:         %4.1lf A\n", VdwInnerCutoff);
		printf("Elec.static Cutoff:       %4.1lf A\n", Cutoff);
		printf("Elec.static inner Cutoff: %4.1lf A\n", InnerCutoff);
		printf("Switching:                %s\n", ForceSwitch ? "Force switching" : (EnergySwitch ? "Potential switch" : "none"));
		if(VdwCor)
		{
			printf("VDW correction          yes\n");
			printf(" +- density             %9.6lf N/A^3\n",VdwCorDensity); 
			printf(" +- epsilon             %7.3lf kcal/mol\n",VdwCorEpsilon);
			printf(" +- radius              %4.2lf A\n",VdwCorRadius);
			printf(" `- Cutoff              %4.2lf A\n",VdwCorCutoff);  
		}
		printf("UsePartialRecalc:         %s\n", UsePartialRecalc ? "Yes" : "No");
		printf("IgnoreIntraResidue:       %s\n", IgnoreIntraResidue ? "Yes" : "No");
	}


	void FF_NonBonded::settodefault()
	{
		DoVdw = true;
		DoElec = true;
		Cutoff = 15;
		InnerCutoff = 12;

		VdwCutoff = Cutoff;
		VdwInnerCutoff = InnerCutoff;
		VdwAttenuator = 1.0;
		ElecAttenuator = 1.0;

		AutoScaling14 = true;
		Vdw14Scaling = 1.0;
		Elec14Scaling = 1.0;

		Dielectric = 1.0;
		DDDielectric = false;
		DDDielectricAlpha = 1.0;

		ForceSwitch=true;    
		EnergySwitch=false;

		// longrange vdw correction
		VdwCor = false;
		VdwCorDensity = 0.0;
		VdwCorEpsilon = 0.0;
		VdwCorRadius = 0.0;
		VdwCorCutoff = 0.0;

		fullrecalc = true;
		UsePartialRecalc = false;

		IgnoreIntraResidue = false;

		fast = 0;
	}


	// The setup function sets up some general stuff before runtime
	void FF_NonBonded::setup()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		Active = true;
		wspace.nlist().requestCutoff(Maths::max(Cutoff,VdwCutoff));
		if(AutoScaling14)
		{
			Vdw14Scaling = wspace.ffps().Vdw14Scaling;
			Elec14Scaling = wspace.ffps().Elec14Scaling;
		}


		// setting up atom energy array
		fullrecalc = true;  // next calculation must be a full one


		// load local parameters (for faster access)
		local_atomparam = new NonBonded_Pack[wspace.nAtoms()];

		for(int i = 0; i < wspace.nAtoms(); i ++)
		{
			local_atomparam[i].epsilon = wspace.atom[i].epsilon;
			local_atomparam[i].radius  = wspace.atom[i].radius;
			local_atomparam[i].charge  = wspace.atom[i].charge;
		}

		setupBasisVectors();
	}

	void  FF_NonBonded::setupBasisVectors(){
		// periodic space make local basisvector table
		WorkSpace& wspace = getWSpace();
		Space *space = &wspace.boundary();
		int nbasisvecs = space->nBasisVectors();
		if( nbasisvecs > 32 )
		{
			throw(ProcedureException(
				"FF_NonBonded: Cannot use more than 32 periodic boundary image vectors (" 
				+ int2str(nbasisvecs) + " basis vectors). \nThis may be because you have chosen " +
				"a non-bonded cutoff which is greater than half the periodic box width.\n"
				"Your settings are: Cutoff=" + double2str(Cutoff) + "\n"+
				"                   Half Smallest Box width=" + double2str(space->getSmallestChord()/2.0) + "\n"
				"Such a setup is not supported by this forcefield. \n"));
		}
		for(int b=0;b<nbasisvecs;b++)
		{
			basisvector[b] = space->getBasisVector(b);
		}
	}


	// Include all variables who's change should trigger a resetup
	unsigned long FF_NonBonded::calcCheckSum()
	{
		unsigned long sum = ForcefieldBase::calcCheckSum();
		sum += (unsigned long)( Cutoff*VdwCutoff*100000.0 ) + 
			1354231*(int)AutoScaling14;;
		return sum;
	}

	void FF_NonBonded::infoLine() const 
	{ 
		// prints a line of current energies
		const WorkSpace& wspace = getWSpace();
		const Hamiltonian *printene = &wspace.ene;
		printf("% 8.1lf% 8.1lf",
			double(printene->epot_vdw) * PhysicsConst::J2kcal * PhysicsConst::Na,
			double(printene->epot_elec) * PhysicsConst::J2kcal * PhysicsConst::Na);
	}

	void FF_NonBonded::infoLineHeader() const
	{ 
		// prints the headers for the above function
		printf("%8s%8s", "EVdw", "Ecoul");
	}


	void FF_NonBonded::calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level)
	{ 
		WorkSpace& wspace = getWSpace();
		if(level == Summary) 
		{
			calcForces();
			printf(" VdW:           %10.3lf kcal/mol\n", double(wspace.ene.epot_vdw) * PhysicsConst::J2kcal * PhysicsConst::Na);
			printf(" Electrostatic: %10.3lf kcal/mol\n", double(wspace.ene.epot_elec) * PhysicsConst::J2kcal * PhysicsConst::Na);
// // DEBUG STUFF TO TEST VIRIAL
//		wspace.ene.InternalVirial = 0;
//		calcForces(); 
//		double ene0 = wspace.ene.epot_elec + wspace.ene.epot_vdw;
//		printf("Ene:    %e Virial:  %e\n", ene0, wspace.ene.InternalVirial );
//		ClosedSpace *cs = dynamic_cast<ClosedSpace *> (&wspace.boundary());
//		double v1 = cs->volume();
//		double dV = 1.0000001;
//		wspace.scaleSystem(dV);
//		wspace.ene.InternalVirial = 0;
//		calcForces(); 
//		double ene1 = wspace.ene.epot_elec + wspace.ene.epot_vdw;
//
//		printf("Ene:    %e Virial:  %e   V1  %e\n", ene1, wspace.ene.InternalVirial, v1 );
//		printf("dEnedV = %e    %e  \n ", (ene1 - ene0)/(v1*dV- v1), wspace.ene.InternalVirial / (3.0 * v1));
//		printf("dEne = %e \n", ene1 - ene0 );
//		printf("V1 V2 dV : %e  %e   %e  \n",v1,cs->volume(), (v1*dV- v1) );
//		//wspace.
			return;
		}
		verbose_level = level;	
 
		setupBasisVectors();
 
		int T_SqrtMode    = T_SqrtFPU;
		int T_VdwMode     = T_VdwMode_EnergySwitch; 
		int T_ElecMode    = T_ElecMode_Normal;
		int T_VerboseMode = T_VerboseMode_True;
 
		if(EnergySwitch)T_ElecMode    =  T_ElecMode_EnergySwitch;
		if(ForceSwitch) T_ElecMode    =  T_ElecMode_ForceSwitch;

		if(!DoElec)T_ElecMode         =  T_ElecMode_None;
		if(!DoVdw)T_VdwMode           =  T_VdwMode_None;
 
		tcall4<FF_NonBonded_CalcForces_T_fast2_wrap>
			(T_SqrtMode, T_VdwMode, T_ElecMode, T_VerboseMode)
			( *this,local_atomparam, &basisvector[0], wspace, level );
 
		epot_elec = wspace.ene.epot_elec;
		epot_vdw  = wspace.ene.epot_vdw; 	
	};

	void FF_NonBonded::calcEnergies(){  
		setupBasisVectors();

		// force a full recalculation of the energy if the neighbor list has been updated
		unsigned currentCount = getWSpace().nlist().getFullUpdateCount();
		if ( currentCount != m_Nlist_FullUpdateCount){
			fullrecalc = true;
			m_Nlist_FullUpdateCount = currentCount;  // remember latest count for next time
		}

		if( (!fullrecalc) && (UsePartialRecalc) ){ 
			//getWSpace().old = getWSpace().cur;
			//calcEnergies_T<false>(); 

			calcEnergies_Update(); 
			//printf("A-Partial: %f \n", epot *  PhysicsConst::J2kcal * PhysicsConst::Na );
			//double epot_s = epot;
			//calcEnergies_T<false>(); 
			//printf("Echeck:    %f \n", (epot - epot_s) *  PhysicsConst::J2kcal * PhysicsConst::Na );

			// save the coordinates for the fast energy evaluation to work next time 
			getWSpace().old = getWSpace().cur;
		}else{
			calcEnergies_Full(); 
			// save the coordinates for the fast energy evaluation to work next time 
			getWSpace().old = getWSpace().cur;
			fullrecalc = false; // full recalc not needed anymore
		}
	};

	void FF_NonBonded::calcEnergies_Full(){ 
		WorkSpace& wspace = getWSpace();

		int T_SqrtMode    = T_SqrtFPU;
		int T_VdwMode     = T_VdwMode_EnergySwitch; 
		int T_ElecMode    = T_ElecMode_Normal;
		int T_VerboseMode = T_VerboseMode_False;

		if(EnergySwitch)T_ElecMode    =  T_ElecMode_EnergySwitch;
		if(ForceSwitch) T_ElecMode    =  T_ElecMode_ForceSwitch;
		if(!DoElec)T_ElecMode         =  T_ElecMode_None;
		if(!DoVdw)T_VdwMode           =  T_VdwMode_None;

		tcall4<FF_NonBonded_CalcForces_T_fast2_wrap>
			(T_SqrtMode, T_VdwMode, T_ElecMode, 0)
			( *this,local_atomparam, &basisvector[0], wspace, Forcefield::Summary );

		epot_elec = wspace.ene.epot_elec;
		epot_vdw  = wspace.ene.epot_vdw; 	

	};

	void FF_NonBonded::calcEnergies_Update(){ 
		WorkSpace& wspace = getWSpace();	
		if( dynamic_cast<const PeriodicBox*>(&wspace.boundary() ) != NULL ){
			FF_NonBonded_CalcEnergies_Update_T<PeriodicBox, false> ( *this, wspace, changed_atom, old_changed_atom,  epot_reference );  
		}	else 
			if( dynamic_cast<const InfiniteSpace*>(&wspace.boundary() ) != NULL ){
				FF_NonBonded_CalcEnergies_Update_T<InfiniteSpace, false> ( *this, wspace, changed_atom, old_changed_atom,  epot_reference ); 
			}	else {
				FF_NonBonded_CalcEnergies_Update_T<int , false>  ( *this, wspace, changed_atom, old_changed_atom,  epot_reference ); 
			}

			epot_elec = wspace.ene.epot_elec;
			epot_vdw  = wspace.ene.epot_vdw; 	

	};

	void FF_NonBonded::calcForces()
	{ 
		setupBasisVectors();

		WorkSpace& wspace = getWSpace();

		int T_SqrtMode    = T_SqrtFPU;
		int T_VdwMode     = T_VdwMode_EnergySwitch; 
		int T_ElecMode    = T_ElecMode_Normal;
		int T_VerboseMode = T_VerboseMode_False;

		if(EnergySwitch)T_ElecMode    =  T_ElecMode_EnergySwitch;
		if(ForceSwitch) T_ElecMode    =  T_ElecMode_ForceSwitch;

		if(!DoElec)T_ElecMode         =  T_ElecMode_None;
		if(!DoVdw)T_VdwMode           =  T_VdwMode_None;

		tcall4<FF_NonBonded_CalcForces_T_fast2_wrap>
			(T_SqrtMode, T_VdwMode, T_ElecMode, 0)
			( *this,local_atomparam, &basisvector[0], wspace, Forcefield::Summary );

		epot_elec = wspace.ene.epot_elec;
		epot_vdw  = wspace.ene.epot_vdw; 	

	};


}



