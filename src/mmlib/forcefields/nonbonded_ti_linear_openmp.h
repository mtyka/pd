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

#include "workspace/space.h"
#include "forcefields/nonbonded.h"

namespace Physics
{




	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	///
	/// \author Mike Tyka  
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class PD_API FF_NonBonded_TI_Linear_OPENMP : public FF_NonBonded,
																				public FF_Extension_TI
	{
	public:
		FF_NonBonded_TI_Linear_OPENMP( WorkSpace &newwspace ):
				FF_NonBonded(newwspace),
				FF_Extension_TI(newwspace)
		{
			settodefault();
		}

		~FF_NonBonded_TI_Linear_OPENMP(){
		};


		bool DecoupleElec;

		bool DecoupleVdw;

		/// prints a little block of parameter information
		void info() const{ 
			FF_NonBonded::info();    // call base class info function (the nonboned part)
			FF_Extension_TI::info(); // call base class info function (the TI extension )
			printf("INFO: DecoupleElec:             %s\n", DecoupleElec ? "Yes\n" : "No\n");
			printf("INFO: DecoupleVdw:              %s\n", DecoupleVdw  ? "Yes\n" : "No\n");
		}

	protected:
		void settodefault()
		{
			FF_NonBonded::settodefault(); // call base class info function
			DecoupleElec = true;
			DecoupleVdw  = true;
		}

		virtual void calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level);
		virtual void calcEnergies();
		virtual void calcForces();
	};


















/*






	/// A specialised forcefield for Thermodynamic Integration. Basically
	/// as well as calculating the nonbonded interaction it allows scaling 
	/// of the energy components (seperately) and it calculates the derivative of
	/// the energy wrt to the lambda parameter
	/// it can be used to determine the free energy of liquids for example (by turning of
	/// nonbonded interactions turns the the system into a perfect gas whose's precise
	/// free energy is known analytically)

	template <typename T_Space>





//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka  
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API NonBonded_TI: public NonBondedForcefield_TemplateBoundary< T_Space >, public FF_Extension_TI
	{
	public:
		NonBonded_TI( WorkSpace &newwspace ) : 
			NonBondedForcefield_TemplateBoundary< PeriodicBox >( newwspace ),
			FF_Extension_TI()
		{
			settodefault();
		}
		~NonBonded_TI(){};

		/// lambda power for hydrogen-hydrogen                             
		int lpower_elec_HH;   

		/// lambda power for hydrogen-non-hydrogen 
		int lpower_elec_HnH;  

		/// lambda power for non-hydrogen-non-hydrogen interactions
		int lpower_elec_nHnH; 

		/// lambda power for hydrogen-hydrogen 
		int lpower_vdw_HH;    

		/// lambda power for hydrogen-non-hydrogen 
		int lpower_vdw_HnH;     

		/// lambda power for non-hydrogen-non-hydrogen interactions
		int lpower_vdw_nHnH;    

		void info() const
		{ 
			// prints a little block of parameter information
			NonBondedForcefield_TemplateBoundary< T_Space >::info(); // call base class info function
			// add some more info
			printf("INFO: Lambda:                   %6.3lf\n", lambda);
			printf("INFO: lpower_elec_HH            %d\n",lpower_elec_HH    );
			printf("INFO: lpower_elec_HnH           %d\n",lpower_elec_HnH    );
			printf("INFO: lpower_elec_nHnH          %d\n",lpower_elec_nHnH   );
			printf("INFO: lpower_vdw_HH             %d\n",lpower_vdw_HH      );
			printf("INFO: lpower_vdw_HnH            %d\n",lpower_vdw_HnH     );
			printf("INFO: lpower_vdw_nHnH           %d\n",lpower_vdw_nHnH    );
			int i;
			double lambda_elec_HH  =1.0;for(i=0;i<lpower_elec_HH  ;i++) lambda_elec_HH   *= lambda;
			double lambda_elec_HnH =1.0;for(i=0;i<lpower_elec_HnH ;i++) lambda_elec_HnH  *= lambda;
			double lambda_elec_nHnH=1.0;for(i=0;i<lpower_elec_nHnH;i++) lambda_elec_nHnH *= lambda;
			double lambda_vdw_HH   =1.0;for(i=0;i<lpower_vdw_HH   ;i++) lambda_vdw_HH   *= lambda;
			double lambda_vdw_HnH  =1.0;for(i=0;i<lpower_vdw_HnH  ;i++) lambda_vdw_HnH   *= lambda;
			double lambda_vdw_nHnH =1.0;for(i=0;i<lpower_vdw_nHnH ;i++) lambda_vdw_nHnH  *= lambda;
			printf("INFO: lambda_elec_HH            %.8lf\n", lambda_elec_HH   );
			printf("INFO: lambda_elec_HnH           %.8lf\n", lambda_elec_HnH  );
			printf("INFO: lambda_elec_nHnH          %.8lf\n", lambda_elec_nHnH );
			printf("INFO: lambda_vdw_HH             %.8lf\n", lambda_vdw_HH    );
			printf("INFO: lambda_vdw_HnH            %.8lf\n", lambda_vdw_HnH   );
			printf("INFO: lambda_vdw_nHnH           %.8lf\n", lambda_vdw_nHnH  );
		}

	protected:

		void settodefault()
		{
			NonBondedForcefield_TemplateBoundary< T_Space >::settodefault(); // call base class info function
			// default the values to those used in Chem Phys Lett 189, 273 (1992)
			lpower_elec_HH = 3; 
			lpower_elec_HnH = 3;
			lpower_elec_nHnH = 3;

			lpower_vdw_HH = 5;  
			lpower_vdw_HnH = 11; 
			lpower_vdw_nHnH = 5;
		}

		virtual void calcEnergiesVerbose(ForcefieldBase::Verbosity level)
		{ 
			typedef  NonBondedForcefield_TemplateBoundary< T_Space > Base;
			WorkSpace& wspace = Base::getWSpace();
			if(level == Base::Summary) 
			{
				calcEnergies();
				printf(" VdW:           %10.3lf kcal/mol\n", wspace.ene.epot_vdw * PhysicsConst::J2kcal * PhysicsConst::Na);
				printf(" Electrostatic: %10.3lf kcal/mol\n", wspace.ene.epot_elec * PhysicsConst::J2kcal * PhysicsConst::Na);
				return;
			}
			Base::verbose_level = level;		
			calcForces_T<true>(); 
		};

		virtual void calcEnergies(){ calcForces_T<false>(); };
		virtual void calcForces(){	 calcForces_T<false>(); };

		/// Generating function:
		template <bool verbosemode>
		void calcForces_T();
	};



	template <typename T_Space>
	template <bool verbosemode>
	void NonBonded_TI<T_Space>::calcForces_T()
	{
		using namespace Maths;	
		typedef  NonBondedForcefield_TemplateBoundary< T_Space > Base;

		// set up proxies to workspace to make code more readable.
		WorkSpace& wspace = Base::getWSpace();
		size_t natom = wspace.atom.size();             // number of atoms in workspace
		const ParticleStore& atomparam = wspace.atom;  // atom parameter array
		SnapShotAtom *atom = wspace.cur.atom;						// atom coordinate array
		const NeighbourData *fnbor = wspace.nlist().getData(); // neighborlist

		// Basic stuff
		int i, j, nj;              // i&j are atom indices, nj count through neighbor list
		dvector fv;     // force dvector
		double force_magnitude= 0;             // force magnitude
		double vdw_force = 0;      // individual force magnitudes
		double elec_force = 0;
		double vdw_potential = 0;  // individual potential energy contributions
		double elec_potential = 0;
		double vdw14scale = 1.0;
		double elec14scale = 1.0;

		double sqrdistij,Dist_ij; // distance in Angstrom
		double invdistij, invd6;
		double sqrcutoff = sqr(Base::Cutoff);

		double radiusij, epsilon, A, B;
		double qi, qj;
		double invdielectric = 1.0 / Base::Dielectric;

		// periodic spaces
		unsigned image;
		unsigned maximages = wspace.boundary().ncells();
		unsigned nimages;
		PeriodicBox *periodic_box;
		if(is_same_type<T_Space,PeriodicBox>::value){
			periodic_box = (PeriodicBox *) &wspace.boundary();
		}

		// local temporary variables
		dvector dc;
		double atomi_radius;
		double atomi_epsilon;
		dvector fv_iatom;
		dvector fv_jatom;

		// elec switching
		double Swidth = Base::Cutoff - Base::InnerCutoff;
		double invSwidth = 1.0 / Swidth;
		double S;
		double dSdd;

		// vdw switching
		double vdwSwidth = Base::VdwCutoff - Base::VdwInnerCutoff;
		double vdwinvSwidth = 1.0 / vdwSwidth;
		double vdwS;
		double vdwdSdd;

		// precalculated stuff for force switching
		double sA = 1.0/cube( sqr(Base::Cutoff) - sqr(Base::InnerCutoff) );
		double sB = -(  cube(sqr(Base::Cutoff)) - 3.0* sqr(Base::Cutoff) * sqr(Base::Cutoff) * sqr(Base::InnerCutoff));
		double sC = 6.0* sqr(Base::Cutoff) * sqr(Base::InnerCutoff);
		double sD = -(sqr(Base::Cutoff) + sqr(Base::InnerCutoff));
		double sE = 2.0/5.0;
		double fswitch_innerV = sA * (sB * (1.0/Base::InnerCutoff) + sC * Base::InnerCutoff + sD * cube(Base::InnerCutoff) + sE * cube(Base::InnerCutoff) * sqr(Base::InnerCutoff));
		double fswitch_cutoffV = sA * (sB * (1.0/Base::Cutoff) + sC * Base::Cutoff + sD * cube(Base::Cutoff) + sE * cube(Base::Cutoff) * sqr(Base::Cutoff));

		// statistics
		int totalpairs = 0;
		int pairs14 = 0;
		int vdwpairs = 0;
		int elecpairs = 0;

		// stuff for the verbose modes
		double lastenergy_vdw = 0;
		double lastenergy_elec = 0;

		// stuff for the lambda values

		double invlambda = 1.0/lambda;
		double lambda_elec_HH  =1.0;for(i=0;i<lpower_elec_HH  ;i++) lambda_elec_HH   *= lambda;
		double lambda_elec_HnH =1.0;for(i=0;i<lpower_elec_HnH ;i++) lambda_elec_HnH  *= lambda;
		double lambda_elec_nHnH=1.0;for(i=0;i<lpower_elec_nHnH;i++) lambda_elec_nHnH *= lambda;

		double lambda_vdw_HH   =1.0;for(i=0;i<lpower_vdw_HH   ;i++) lambda_vdw_HH   *= lambda;
		double lambda_vdw_HnH  =1.0;for(i=0;i<lpower_vdw_HnH  ;i++) lambda_vdw_HnH   *= lambda;
		double lambda_vdw_nHnH =1.0;for(i=0;i<lpower_vdw_nHnH ;i++) lambda_vdw_nHnH  *= lambda;

		if(verbosemode){ // only in the verbose instantiation
			printf("\nPairwise simple forces \n\n");
			if(Base::verbose_level == Base::Detailed)
				printf("Pair: i nr(Name) nj nr(Name) Vdw: dist radiusij Vvdw Elec: chgi chgj Dist_ij Velec \n\n");
		}

		// initialise stuff
		wspace.ene.epot_vdw = 0;
		wspace.ene.epot_elec = 0;

		dEdlambda = 0.0;

		for(i = 0; i < natom; i++) 
		{ 
			// loop over all particles
			atomi_radius = atomparam[i].radius;
			atomi_epsilon = atomparam[i].epsilon;
			qi = atomparam[i].charge;

			fv_iatom.zero();

			for(nj = 0; nj < fnbor[i].n; nj++) { // and all their neighbours
				j = fnbor[i].i[nj];
				fv_jatom.zero();
				if(i > j) break;
				if(fnbor[i].Type[nj] > 1)	continue;// only get non-bonded and 1-4 neighbors (half matrix, i.e. ignore shadows)
				vdw14scale = 1.0;
				elec14scale = 1.0;

				if(fnbor[i].Type[nj] == 1) {
					pairs14++;
					vdw14scale = Base::Vdw14Scaling;
					elec14scale = Base::Elec14Scaling;
				}

				force_magnitude= 0; // total force/total pot. energy = 0
				dc.diff(wspace.cur.atom[j].p,wspace.cur.atom[i].p);

				nimages=maximages;
				if(is_same_type<T_Space,PeriodicBox>::value){
					if(maximages > 1)
						if( dc.innerdot() < sqr( periodic_box->savedist*2.0 - Base::Cutoff )){
							nimages = 1;
						}
				}

				image = 0;
imageloop:
				fv.setTo(dc);

				if(is_same_type<T_Space,PeriodicBox>::value){
					// Specialised code for PeriodicBox space boundaries

					fv.add(periodic_box->celloffset[image]);

					while(fv.x >  periodic_box->halfCellSize.x){ fv.x -= periodic_box->cellSize.x; }
					while(fv.x < -periodic_box->halfCellSize.x){ fv.x += periodic_box->cellSize.x; }
					while(fv.y >  periodic_box->halfCellSize.y){ fv.y -= periodic_box->cellSize.y; }
					while(fv.y < -periodic_box->halfCellSize.y){ fv.y += periodic_box->cellSize.y; }
					while(fv.z >  periodic_box->halfCellSize.z){ fv.z -= periodic_box->cellSize.z; }
					while(fv.z < -periodic_box->halfCellSize.z){ fv.z += periodic_box->cellSize.z; }

					sqrdistij = fv.innerdot();
					if(sqrdistij > sqrcutoff) goto imageloop_continue;			
					Dist_ij = sqrt(sqrdistij);
				}
				else if(is_same_type<T_Space,InfiniteSpace>::value){
					// do nothing in infinite space - the real and imaginary 
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
				invd6 = cube(invdistij) * cube(invdistij);

				// Van der Waals Force
				vdw_potential = 0;
				vdw_force = 0;
				elec_potential = 0;
				elec_force = 0;
				if(Base::DoVdw){
					if(Dist_ij < Base::VdwCutoff) {
						vdwpairs++;
						// first calculate well depth (epsilon) and zero-force-separation (radiusij)
						// Use Lorentz-Berthelot mixing rules
						// arthemic mean for radiusij (=sum of radii) and geometric mean for epsilon
						radiusij = atomi_radius + atomparam[j].radius;
						// radiusij *= 0.890898718; // div 2^(1/6)

						// two different mixing rules - at the moment only geometric mean is implemented
						epsilon = vdw14scale * atomi_epsilon * atomparam[j].epsilon;

						B = sqr(radiusij) * sqr(radiusij) * sqr(radiusij) * invd6;
						A = sqr(B); 

						vdw_potential = (2.0 * epsilon) * (0.5 * A - B);
						vdw_force = -(12.0 * epsilon * invdistij) * (A - B) / PhysicsConst::Angstrom;


						if(Dist_ij > Base::VdwInnerCutoff) {
							vdwS = (1.0 - sqr(vdwinvSwidth * (Dist_ij - Base::VdwInnerCutoff))); // Dimension less
							vdwdSdd = -4.0 * sqr(vdwinvSwidth) * (Dist_ij - Base::VdwInnerCutoff) * vdwS / PhysicsConst::Angstrom; // Units of per meter (not per angstrom!)
							vdwS = sqr(vdwS);
							vdw_force = vdwS * vdw_force + vdw_potential * vdwdSdd;

							vdw_potential *= vdwS;
						}

					}
				}
				// electrostatics --------------------------------

				if(Base::DoElec) {
					if(Dist_ij < Base::Cutoff) {
						qj = atomparam[j].charge;
						if(Base::EnergySwitch){
							if(Base::DDDielectric) {   // distance dependent Base::Dielectric 
								double incdistij = ((Base::Dielectric + Base::DDDielectricAlpha * Dist_ij) * Dist_ij);
								elec_potential = PhysicsConst::econv_joule * elec14scale * (qi * qj) / incdistij;
								elec_force = -elec_potential * (Base::Dielectric + 2.0 * Base::DDDielectricAlpha * Dist_ij) / (incdistij * PhysicsConst::Angstrom);

							} else {
								elecpairs++;
								elec_potential = PhysicsConst::econv_joule * invdielectric * elec14scale * (qi * qj) * invdistij;
								elec_force = -elec_potential * invdistij / PhysicsConst::Angstrom;
							}

							if(Dist_ij > Base::InnerCutoff) {
								S = (1.0 - sqr(invSwidth * (Dist_ij - Base::InnerCutoff))); // Dimension less
								dSdd = -4.0 * sqr(invSwidth) * (Dist_ij - Base::InnerCutoff) * S / PhysicsConst::Angstrom; // Units of per meter (not per angstrom!)
								S = sqr(S);
								elec_force = S * elec_force + elec_potential * dSdd;
								elec_potential *= S;
							}
						}
						if(Base::ForceSwitch){
							if(Dist_ij > Base::InnerCutoff) {
								elec_force =   -sA * sqr(invdistij) * 
									( sqr(sqr(Base::Cutoff) - sqrdistij)*
									(sqr(Base::Cutoff) + 2.0*sqrdistij - 3.0*sqr(Base::InnerCutoff)) );

								elec_potential =  -  ( sA * Dist_ij *
									(sB * sqr(invdistij) + 
									sC + 
									sD * sqrdistij + 
									sE * sqr(sqrdistij))
									-fswitch_cutoffV);
							}else{
								elec_force =  -sqr(invdistij);
								elec_potential   =  (invdistij - 1/Base::InnerCutoff - (fswitch_innerV - fswitch_cutoffV) );		//
							}

							elec_force *=  PhysicsConst::econv_joule / PhysicsConst::Angstrom * (qi * qj);
							elec_potential *=    PhysicsConst::econv_joule * elec14scale * (qi * qj);
						}
					}
				}

				// do TI bits and bobs - i.e. scale energy and forces by seperate lambda values
				// and determine the derivate of the energy wrt to lambda;
				if((wspace.atom[i].isHydrogen())&&
					(wspace.atom[j].isHydrogen())){
						vdw_potential *=lambda_vdw_HH;
						elec_potential*=lambda_elec_HH;
						vdw_force     *=lambda_vdw_HH;
						elec_force    *=lambda_elec_HH;
						dEdlambda += invlambda*vdw_potential * double(lpower_vdw_HH);
						dEdlambda += invlambda*elec_potential* double(lpower_elec_HH);
				}else
					if((!wspace.atom[i].isHydrogen())&&
						(!wspace.atom[j].isHydrogen())){
							vdw_potential *=lambda_vdw_nHnH;
							elec_potential*=lambda_elec_nHnH;
							vdw_force     *=lambda_vdw_nHnH;
							elec_force    *=lambda_elec_nHnH;
							dEdlambda += invlambda*vdw_potential * double(lpower_vdw_nHnH);
							dEdlambda += invlambda*elec_potential* double(lpower_elec_nHnH);
					}else{
						vdw_potential *= lambda_vdw_HnH;
						elec_potential*= lambda_elec_HnH;
						vdw_force     *=lambda_vdw_HnH;
						elec_force    *=lambda_elec_HnH;
						dEdlambda += invlambda*vdw_potential * double(lpower_vdw_HnH);
						dEdlambda += invlambda*elec_potential* double(lpower_elec_HnH);
					}

					// add up potentials -----------------------------

					wspace.ene.epot += vdw_potential;
					wspace.ene.epot_vdw += vdw_potential;
					if(vdw_potential > 0.0) wspace.ene.epot_vdw_rep += vdw_potential;
					else wspace.ene.epot_vdw_att += vdw_potential;

					wspace.ene.epot += elec_potential;
					wspace.ene.epot_elec += elec_potential;

					// now apply the force
					force_magnitude= 0;
					force_magnitude+= vdw_force;
					force_magnitude+= elec_force;

					// for intermolecular forces add up virial components
					if( atomparam[i].imol != atomparam[j].imol ){
						wspace.ene.InternalVirial += Dist_ij * force_magnitude* PhysicsConst::Angstrom;
					}
					fv.mul(invdistij * force_magnitude);
					fv_iatom.add(fv);
					fv_jatom.add(fv);

					if(verbosemode){
						if(Base::verbose_level == Base::Detailed) {
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

					// an "explicit" for loop 
imageloop_continue:

					if(!is_same_type<T_Space,InfiniteSpace>::value){
						image++;
						if(image<nimages) goto imageloop;
					}
					atom[j].f.sub(fv_jatom);

			}
			atom[i].f.add(fv_iatom);

			// add a vdw correction if required
			if(Base::VdwCor){
				// this correction is the analytical integral of the attractive portion of the LJ potential
				// from VdwCorCutoff to infinity. usually you;d want to set VdwCorCutoff == VdwCutoff;

				vdw_potential = -Base::VdwCorDensity * Base::VdwCorEpsilon / (PhysicsConst::J2kcal * PhysicsConst::Na)
					* (4.0/3.0) * MathConst::PI * sqr(Base::VdwCorRadius) * 
					sqr(Base::VdwCorRadius) * sqr(Base::VdwCorRadius) /
					cube(Base::VdwCorCutoff);
				// for the purpose of vdw path assume it's non hydrogen/non hydorgen
				vdw_potential *=lambda_vdw_nHnH;
				dEdlambda += invlambda*vdw_potential * double(lpower_vdw_nHnH);
				wspace.ene.epot += vdw_potential;
				wspace.ene.epot_vdw += vdw_potential;
				wspace.ene.epot_vdw_att += vdw_potential;
			}

			if( verbosemode){
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

		if(verbosemode){
			printf("\n");
			printf(" ----------------- -----------------\n");
			printf(" Total: %8.3lf Total: %8.3lf \n\n\n",
				wspace.ene.epot_vdw * PhysicsConst::J2kcal * PhysicsConst::Na, 
				wspace.ene.epot_elec * PhysicsConst::J2kcal * PhysicsConst::Na);
			printf(" dEdLambda         %10.5f\n", dEdlambda);
			printf(" Total number of interactions:\n");
			printf("    Total:         %8d\n",totalpairs);
			printf("    1-4:           %8d\n",pairs14);
			printf("    Van d. Waals:  %8d\n",vdwpairs);
			printf("    Electrostatic: %8d\n",elecpairs);
		}

	}

#ifdef SWIG
	%template(NonBondedFF_TI_PeriodicBox)  NonBonded_TI<PeriodicBox>;
#endif


















	/// VDW ti after zacharias et al. doesnt do electrostatics !

	template <typename T_Space>





//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka  
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API NonBonded_TI_Zacharias: public NonBondedForcefield_TemplateBoundary< T_Space >,
		public FF_Extension_TI
	{
	public:
		NonBonded_TI_Zacharias( WorkSpace &newwspace ):
				NonBondedForcefield_TemplateBoundary< PeriodicBox >(newwspace),
				FF_Extension_TI()
		{
			settodefault();
			m_LambdaDependentPick				 = new PickNothing();   // applies to all atoms by default
			m_LambdaInverseDependentPick = new PickNothing();   // applies to all atoms by default

		}

		void setLambdaDependent(const PickBase &_Picker){
			delete m_LambdaDependentPick;
			m_LambdaDependentPick = _Picker.clone();	
		}

		void setLambdaInverseDependent(const PickBase &_Picker){
			delete m_LambdaInverseDependentPick;
			m_LambdaInverseDependentPick = _Picker.clone();
		}

		~NonBonded_TI_Zacharias(){
			delete m_LambdaDependentPick;
			delete m_LambdaInverseDependentPick;
		};

		double TI_Delta;

		void info() const{ // prints a little block of parameter information
			NonBondedForcefield_TemplateBoundary< T_Space >::info(); // call base class info function
			// add some more info
			printf("INFO: Lambda:                   %6.3lf\n", lambda);
			printf("INFO: Delta:                    %6.3lf\n", TI_Delta  );
		}

	protected:
		void settodefault()
		{
			NonBondedForcefield_TemplateBoundary< T_Space >::settodefault(); // call base class info function
			TI_Delta = 9;
		}

		virtual void setup(){
			typedef  NonBondedForcefield_TemplateBoundary< T_Space > Base;
			NonBondedForcefield_TemplateBoundary< T_Space >::setup();

			const WorkSpace &wspace = Base::getWSpace();

			setupAtomGroups(wspace.natom());

			for(size_t i = 0; i < wspace.natom(); i++){
				if( m_LambdaDependentPick->matches( wspace.atom[i] ) ){
					mAtomGroup[i] = LambdaDependent;
				}
				if( m_LambdaInverseDependentPick->matches( wspace.atom[i] ) ){
					mAtomGroup[i] = LambdaInverseDependent;
				}
			}
			return 0;
		};

		void printLambdaDependence(){
			typedef  NonBondedForcefield_TemplateBoundary< T_Space > Base;
			const WorkSpace &wspace = Base::getWSpace();
			for(size_t i = 0; i < wspace.natom(); i++){
				printf("LambdaState:  %d \n", int(mAtomGroup[i]) );
			}
		}

		PickBase     *m_LambdaDependentPick;
		PickBase     *m_LambdaInverseDependentPick;


		virtual void calcEnergiesVerbose(ForcefieldBase::Verbosity level)
		{ 
			typedef  NonBondedForcefield_TemplateBoundary< T_Space > Base;
			WorkSpace& wspace = Base::getWSpace();
			if(level == Base::Summary) {
				calcEnergies();
				printf(" VdW:           %10.3lf kcal/mol\n", wspace.ene.epot_vdw * PhysicsConst::J2kcal * PhysicsConst::Na);
				return;
			}
			Base::verbose_level = level;		
			calcForces_T<true>(); 
		}

		virtual void calcEnergies(){ calcForces_T<false>();   };
		virtual void calcForces(){	 calcForces_T<false>(); };

		/// Generating function:
		template <bool verbosemode>
		void calcForces_T();

	};



	template <typename T_Space>
	template <bool verbosemode>
	void NonBonded_TI_Zacharias<T_Space>::calcForces_T()
	{
		using namespace Maths;	
		typedef  NonBondedForcefield_TemplateBoundary< T_Space > Base;

		// set up proxies to workspace to make code more readable.
		WorkSpace& wspace = ForcefieldBase::getWSpace();
		size_t   natom = wspace.atom.size();             // number of atoms in workspace
		const ParticleStore& atomparam = wspace.atom;  // atom parameter array
		SnapShotAtom *atom = wspace.cur.atom;						// atom coordinate array
		const NeighbourData *fnbor = wspace.nlist().getData(); // neighborlist

		// Basic stuff
		int i, j, nj;              // i&j are atom indices, nj count through neighbor list
		dvector fv;     // force dvector
		double force_magnitude= 0;             // force magnitude
		double vdw_force = 0;      // individual force magnitudes
		double elec_force = 0;
		double vdw_potential = 0;  // individual potential energy contributions
		double elec_potential = 0;
		double vdw14scale = 1.0;
		double elec14scale = 1.0;

		double sqrdistij,Dist_ij;             // distance in Angstrom
		double invdistij, invd6;
		double sqrcutoff = sqr(Base::Cutoff);
		const double sqrinnercutoff = sqr(Base::InnerCutoff);

		double radiusij, epsilon, A, B;
		double qi, qj;
		double invdielectric = 1.0 / Base::Dielectric;

		// periodic spaces
		unsigned image;
		unsigned maximages = wspace.boundary().ncells();
		unsigned nimages;
		PeriodicBox *periodic_box;
		if(is_same_type<T_Space,PeriodicBox>::value){
			periodic_box = (PeriodicBox *) &wspace.boundary();
		}

		// local temporary variables
		dvector dc;
		double atomi_radius;
		double atomi_epsilon;
		dvector fv_iatom;
		dvector fv_jatom;

		// elec switching
		double Swidth = Base::Cutoff - Base::InnerCutoff;
		double invSwidth = 1.0 / Swidth;
		double S;
		double dSdd;

		// vdw switching
		double vdwSwidth = Base::VdwCutoff - Base::VdwInnerCutoff;
		double vdwinvSwidth = 1.0 / vdwSwidth;
		double vdwS;
		double vdwdSdd;

		// precalculated stuff for force switching
		const double sA = 1.0/cube( sqr(Base::Cutoff) - sqr(Base::InnerCutoff) );
		const double sB = -(  cube(sqr(Base::Cutoff)) - 3.0* sqr(Base::Cutoff) * sqr(Base::Cutoff) * sqr(Base::InnerCutoff));
		const double sC = 6.0* sqr(Base::Cutoff) * sqr(Base::InnerCutoff);
		const double sD = -(sqr(Base::Cutoff) + sqr(Base::InnerCutoff));
		const double sE = 2.0/5.0;
		const double fswitch_innerV = sA * (sB * (1.0/Base::InnerCutoff) + sC * Base::InnerCutoff + sD * cube(Base::InnerCutoff) + sE * cube(Base::InnerCutoff) * sqr(Base::InnerCutoff));
		const double fswitch_cutoffV = sA * (sB * (1.0/Base::Cutoff)     + sC * Base::Cutoff      + sD * cube(Base::Cutoff)      + sE * cube(Base::Cutoff)      * sqr(Base::Cutoff));
		double eshift;

		// statistics
		int totalpairs = 0;
		int pairs14 = 0;
		int vdwpairs = 0;
		int elecpairs = 0;

		// stuff for the verbose modes
		double lastenergy_vdw = 0;
		double lastenergy_elec = 0;

		double rsh,rsh_2,rsh_3,rsh_4,rsh_6,rsh_7;
		double dpotdlambda;





		// stuff for the lambda values

		if(verbosemode){ // only in the verbose instantiation
			printf("\nPairwise simple forces \n\n");
			if(Base::verbose_level == Base::Detailed)
				printf("Pair: i nr(Name) nj nr(Name) Vdw: dist radiusij Vvdw Elec: chgi chgj Dist_ij Velec \n\n");
		}

		// initialise stuff
		wspace.ene.epot_vdw = 0;
		wspace.ene.epot_elec = 0;

		dEdlambda = 0.0;


		for(i = 0; i < natom; i++) { // loop over all particles
			atomi_radius = atomparam[i].radius;
			atomi_epsilon = atomparam[i].epsilon;
			qi = atomparam[i].charge;

			fv_iatom.zero();

			for(nj = 0; nj < fnbor[i].n; nj++) { // and all their neighbours
				j = fnbor[i].i[nj];
				fv_jatom.zero();
				if(i > j) break;
				if(fnbor[i].Type[nj] > 1)	continue;// only get non-bonded and 1-4 neighbors (half matrix, i.e. ignore shadows)
				vdw14scale = 1.0;
				elec14scale = 1.0;

				if(fnbor[i].Type[nj] == 1) {
					pairs14++;
					vdw14scale = Base::Vdw14Scaling;
					elec14scale = Base::Elec14Scaling;
				}

				force_magnitude= 0; // total force/total pot. energy = 0
				dc.diff(wspace.cur.atom[j].p,wspace.cur.atom[i].p);

				nimages=maximages;
				if(is_same_type<T_Space,PeriodicBox>::value){
					if(maximages > 1)
						if( dc.innerdot() < sqr( periodic_box->savedist*2.0 - Base::Cutoff )){
							nimages = 1;
						}
				}

				image = 0;
imageloop:
				fv.setTo(dc);

				if(is_same_type<T_Space,PeriodicBox>::value){
					// Specialised code for PeriodicBox space boundaries

					fv.add(periodic_box->celloffset[image]);

					while(fv.x >  periodic_box->halfCellSize.x){ fv.x -= periodic_box->cellSize.x; }
					while(fv.x < -periodic_box->halfCellSize.x){ fv.x += periodic_box->cellSize.x; }
					while(fv.y >  periodic_box->halfCellSize.y){ fv.y -= periodic_box->cellSize.y; }
					while(fv.y < -periodic_box->halfCellSize.y){ fv.y += periodic_box->cellSize.y; }
					while(fv.z >  periodic_box->halfCellSize.z){ fv.z -= periodic_box->cellSize.z; }
					while(fv.z < -periodic_box->halfCellSize.z){ fv.z += periodic_box->cellSize.z; }

					sqrdistij = fv.innerdot();
					if(sqrdistij > sqrcutoff) goto imageloop_continue;			
					Dist_ij = sqrt(sqrdistij);
				}
				else if(is_same_type<T_Space,InfiniteSpace>::value){
					// do nothing in infinite space - the real and imaginary 
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
				invd6 = cube(invdistij) * cube(invdistij);

				// Van der Waals Force
				vdw_potential = 0;
				vdw_force = 0;
				elec_potential = 0;
				elec_force = 0;
				if(Base::DoVdw){
					if(Dist_ij < Base::VdwCutoff) {
						vdwpairs++;
						// first calculate well depth (epsilon) and zero-force-separation (radiusij)
						// Use Lorentz-Berthelot mixing rules
						// arthemic mean for radiusij (=sum of radii) and geometric mean for epsilon
						radiusij = atomi_radius + atomparam[j].radius;
						// radiusij *= 0.890898718; // div 2^(1/6)

						// two different mixing rules - at the moment only geometric mean is implemented
						epsilon = vdw14scale * atomi_epsilon * atomparam[j].epsilon;

						B = sqr(radiusij) * sqr(radiusij) * sqr(radiusij) ;
						A = sqr(B); 


						rsh = 1.0/(sqrdistij + TI_Delta * (1.0 - lambda) );
						rsh_2 = sqr(rsh);
						rsh_3 = cube(rsh);
						rsh_4 = sqr(rsh_2);
						rsh_6 = rsh_4*rsh_2;
						rsh_7 = rsh_4*rsh_3;


						vdw_potential = (2.0 * epsilon) * (0.5 * A*rsh_6 - B*rsh_3);
						vdw_force = -(2.0 * epsilon) * (3.0*A*rsh_7 - 3.0*B*rsh_4);
						dpotdlambda = vdw_potential - TI_Delta*lambda*vdw_force;
						//printf(" --> %e  %e   %e    %e %e  %e   :   %e   %e \n",
						//	Dist_ij, epsilon, A, B, lambda, TI_Delta,  
						//	vdw_potential, - TI_Delta*lambda*vdw_force );

						vdw_force *= lambda*2.0*Dist_ij / PhysicsConst::Angstrom;
						vdw_potential *= lambda;

						if(Dist_ij > Base::VdwInnerCutoff) {
							vdwS = (1.0 - sqr(vdwinvSwidth * (Dist_ij - Base::VdwInnerCutoff))); // Dimension less
							vdwdSdd = -4.0 * sqr(vdwinvSwidth) * (Dist_ij - Base::VdwInnerCutoff) * vdwS / PhysicsConst::Angstrom; // Units of per meter (not per angstrom!)
							vdwS = sqr(vdwS);
							vdw_force = vdwS * vdw_force + vdw_potential * vdwdSdd;

							vdw_potential *= vdwS;
							dpotdlambda *= vdwS;
						}

						dEdlambda += dpotdlambda;

					}
				}
		// electrostatics --------------------------------

				if(Base::DoElec) {
				if(Dist_ij < Base::Cutoff) {
					qj = atomparam[j].charge;
					elecpairs++;
					if(Base::EnergySwitch){
						if(Base::DDDielectric) {   // distance dependent Dielectric 
							double incdistij = ((Base::Dielectric + Base::DDDielectricAlpha * Dist_ij) * Dist_ij);
							elec_potential = PhysicsConst::econv_joule * elec14scale * (qi * qj) / incdistij;
							elec_force = -elec_potential * (Base::Dielectric + 2.0 * Base::DDDielectricAlpha * Dist_ij) / (incdistij * PhysicsConst::Angstrom);
							
						} else {
							elec_potential = PhysicsConst::econv_joule * invdielectric * elec14scale * (qi * qj) * invdistij;
							elec_force = -elec_potential * invdistij * 1E10;
						}

						if(Dist_ij >Base:: InnerCutoff) {
							S = (1.0 - sqr(invSwidth * (Dist_ij -Base:: InnerCutoff))); // Dimension less
							dSdd = -4.0E10 * sqr(invSwidth) * (Dist_ij -Base:: InnerCutoff) * S ; // Units of per meter (not per angstrom!)
							S = sqr(S);
							elec_force = S * elec_force + elec_potential * dSdd;
							elec_potential *= S;
						}
					} else
					if(Base::ForceSwitch){
						//elec_force *=     PhysicsConst::econv_joule * invdielectric * elec14scale * (qi * qj) * 1E10 ;

						elec_potential = PhysicsConst::econv_joule * invdielectric * elec14scale * (qi * qj);
						if(Dist_ij >Base:: InnerCutoff) {
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
						elec_potential = PhysicsConst::econv_joule * invdielectric * elec14scale * (qi * qj) * invdistij;
						elec_force = -elec_potential * invdistij * 1E10;
					}
				}
				}


				// add up potentials -----------------------------

				wspace.ene.epot += vdw_potential;
				wspace.ene.epot_vdw += vdw_potential;
				if(vdw_potential > 0.0) wspace.ene.epot_vdw_rep += vdw_potential;
				else wspace.ene.epot_vdw_att += vdw_potential;

				wspace.ene.epot += elec_potential;
				wspace.ene.epot_elec += elec_potential;

				// now apply the force
				force_magnitude= 0;
				force_magnitude+= vdw_force;

				// now apply the force
				force_magnitude= 0;
				force_magnitude+= vdw_force;
				force_magnitude+= elec_force;

				// for intermolecular forces add up virial components
				if( atomparam[i].imol != atomparam[j].imol ){
					wspace.ene.InternalVirial += Dist_ij * force_magnitude* PhysicsConst::Angstrom;
				}
				fv.mul(invdistij * force_magnitude);
				fv_iatom.add(fv);
				fv_jatom.add(fv);

				if(verbosemode){
					if(Base::verbose_level == Base::Detailed) {
						printf("Pair:%5d(%4s)%5d(%4s) Vdw: %5.2lf %5.2lf %5.2lf %8.3lf\n",
							i, atomparam[i].pdbname.c_str(),
							j, atomparam[j].pdbname.c_str(),
							Dist_ij,
							radiusij,
							epsilon * PhysicsConst::J2kcal * PhysicsConst::Na,
							vdw_potential * PhysicsConst::J2kcal * PhysicsConst::Na);
					}
				}

				// an "explicit" for loop 
imageloop_continue:

				if(!is_same_type<T_Space,InfiniteSpace>::value){
					image++;
					if(image<nimages) goto imageloop;
				}
				atom[j].f.sub(fv_jatom);

			}
			atom[i].f.add(fv_iatom);

			// add a vdw correction if required
			if(Base::VdwCor){
				// this correction is the analytical integral of the attractive portion of the LJ potential
				// from VdwCorCutoff to infinity. usually you;d want to set VdwCorCutoff == VdwCutoff;

				vdw_potential = -Base::VdwCorDensity * Base::VdwCorEpsilon / (PhysicsConst::J2kcal * PhysicsConst::Na)
					* (4.0/3.0) * MathConst::PI * sqr(Base::VdwCorRadius) * 
					sqr(Base::VdwCorRadius) * sqr(Base::VdwCorRadius) /
					cube(Base::VdwCorCutoff);
				// for the purpose of vdw path assume it's non hydrogen/non hydorgen
				vdw_potential *=lambda;
				dEdlambda += vdw_potential;
				wspace.ene.epot += vdw_potential;
				wspace.ene.epot_vdw += vdw_potential;
				wspace.ene.epot_vdw_att += vdw_potential;
			}

			if( verbosemode){
				printf("Atom: %12.7lf %12.7lf %12.7lf %5d(%4s) %5d Rad: %6.4lf Eps: %6.4lf Vdw: % 8.3lf \n",
					atom[i].p.x, atom[i].p.y, atom[i].p.z,
					i, atomparam[i].pdbname.c_str(),
					fnbor[i].n,
					atomparam[i].radius,
					sqr(atomparam[i].epsilon) * PhysicsConst::J2kcal * PhysicsConst::Na,
					(wspace.ene.epot_vdw - lastenergy_vdw) * PhysicsConst::J2kcal * PhysicsConst::Na);
				lastenergy_vdw = wspace.ene.epot_vdw;
				lastenergy_elec = wspace.ene.epot_elec;
			}
		}

		dEdlambda *=  PhysicsConst::J2kcal * PhysicsConst::Na;

		if(verbosemode){
			printf("\n");
			printf(" ----------------- -----------------\n");
			printf(" Total: %8.3lf Total: %8.3lf \n\n\n",
				wspace.ene.epot_vdw * PhysicsConst::J2kcal * PhysicsConst::Na, 
				wspace.ene.epot_elec * PhysicsConst::J2kcal * PhysicsConst::Na);
			printf(" dEdLambda         %10.5f\n", dEdlambda);
			printf(" Total number of interactions:\n");
			printf("    Total:         %8d\n",totalpairs);
			printf("    1-4:           %8d\n",pairs14);
			printf("    Van d. Waals:  %8d\n",vdwpairs);
		}

	}

#ifdef SWIG
	%template(NonBondedFF_TI_Zacharias_PeriodicBox)  NonBonded_TI_Zacharias<PeriodicBox>;
#endif

*/

}


