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

// OpenMP headers for multi-core parallelisation
#ifdef HAVE_OPENMP
	#include <omp.h>
#endif

#include "traits.h"
#include "funcgen.h"

#include "forcefields/nonbonded.h" // provides base class
#include "forcefields/nonbonded_ti_linear_openmp.h" // provides base class

#include "workspace/workspace.h"
#include "workspace/neighbourlist.h"
#include "workspace/space.h"

namespace Physics
{

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

	const int T_VerboseMode_False = 1;
	const int T_VerboseMode_True = 1;

	template <
		int T_SqrtMode,
		int T_VdwMode, 
		int T_ElecMode, 
		int T_VerboseMode
	>
	void FF_NonBonded_TI_Linear_OPENMP_CalcForces_T_fast2(
	FF_NonBonded_TI_Linear_OPENMP      &ff, 
	NonBonded_Pack    *local_atomparam,
	Maths::dvector    *basisvector,
	WorkSpace         &wspace, 
	double 						 current_lambda,
	double 						&dEdlambda,
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

#ifdef HAVE_OPENMP
			// with OPENMP enabled, still use only one thread when 
			// calculating stuff verbosely
			omp_set_dynamic(0);
			omp_set_num_threads(1);
#endif

			// only in the verbose instantiation
			printf("\nPairwise simple forces \n\n");
			if(verbose_level == ForcefieldBase::Detailed)
				printf("Pair: i nr(Name) nj nr(Name) Vdw: dist radiusij Vvdw Elec: chgi chgj Dist_ij Velec \n\n");
		}else{
#ifdef HAVE_OPENMP
			omp_set_dynamic(1);  // let system decide on number of threads
#endif
		}


		// initialise energies
		ff.epot = 0;
		wspace.ene.epot_vdw = 0;
		wspace.ene.epot_elec = 0;
		dEdlambda = 0.0;

		// make a little proxy
		int *nlistptr=&fnbor[0].i[0];
		// loop over all particles

#ifdef HAVE_OPENMP
		// blocks of 25 seem reasonable. could be made larger.
		#pragma omp parallel for schedule(guided, 25)
#endif

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

				// E = lambda * Eff
				// dEdlambda = Eff
				//
				if(  ff.DecoupleVdw ){
					dEdlambda += vdw_potential;
					vdw_potential *= current_lambda;
					vdw_force     *= current_lambda;
				}

				if(  ff.DecoupleElec ){
					dEdlambda += elec_potential;
					elec_potential *= current_lambda;
					elec_force     *= current_lambda;
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

		dEdlambda *=  PhysicsConst::J2kcal * PhysicsConst::Na;

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
	struct FF_NonBonded_TI_Linear_OPENMP_CalcForces_T_fast2_wrap{
		FF_NonBonded_TI_Linear_OPENMP_CalcForces_T_fast2_wrap():ptr(&FF_NonBonded_TI_Linear_OPENMP_CalcForces_T_fast2<a, b, c, d>){}
		typedef void (*Tptr)(
			FF_NonBonded_TI_Linear_OPENMP      &ff, 
			NonBonded_Pack    *local_atomparam,
			Maths::dvector    *basisvector,
			WorkSpace         &wspace, 
			double 						 current_lambda,
			double 						&dEdlambda,
			ForcefieldBase::AtomicVerbosity verbose_level);
		Tptr ptr;
		const static int limita=T_SqrtTable; 
		const static int limitb=T_VdwMode_EnergySwitch;
		const static int limitc=T_ElecMode_ForceSwitch;
		const static int limitd=T_VerboseMode_True;
		static void overflow()
		{
			throw(CodeException("CODE ERROR: FF_NonBonded_TI_Linear_OPENMP_CalcForces_T_fast2_wrap template requested that is out of bounds.") ); 
		}
	};

	void FF_NonBonded_TI_Linear_OPENMP::calcEnergies()
	{ 
		setupBasisVectors();
		calcForces();
	}

	void FF_NonBonded_TI_Linear_OPENMP::calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level)
	{ 
		WorkSpace& wspace = getWSpace();
		if(level == Summary) {
			calcEnergies();
			printf(" VdW:           %10.3lf kcal/mol\n", wspace.ene.epot_vdw * PhysicsConst::J2kcal * PhysicsConst::Na);
			printf(" Electrostatic: %10.3lf kcal/mol\n", wspace.ene.epot_elec * PhysicsConst::J2kcal * PhysicsConst::Na);
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

		tcall4<FF_NonBonded_TI_Linear_OPENMP_CalcForces_T_fast2_wrap>
			(T_SqrtMode, T_VdwMode, T_ElecMode, 0)
			( *this,local_atomparam, &basisvector[0], wspace, 
				lambda, dEdlambda, level );
	}



	void FF_NonBonded_TI_Linear_OPENMP::calcForces(){ 
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

		tcall4<FF_NonBonded_TI_Linear_OPENMP_CalcForces_T_fast2_wrap>
			(T_SqrtMode, T_VdwMode, T_ElecMode, 0)
			( *this,local_atomparam, &basisvector[0], wspace, 
			lambda, dEdlambda, Forcefield::Summary );
	}
}
