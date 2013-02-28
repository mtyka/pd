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

#ifndef __FF_NONBONDED_H
#define __FF_NONBONDED_H

// OpenMP headers for multi-core parallelisation
#ifdef _OPENMP
#include <omp.h>
#endif

#include "forcefields/forcefield.h" // provides base class
#include "workspace/workspace.h"
#include "workspace/neighbourlist.h"
#include "monitors/monitorbase.h"

class PD_API Particle;

namespace Physics
{
	
	// This holds a local copy of the parameters for (much) faster access.
	struct NonBonded_Pack{
		 double radius;
 		 double epsilon;
     double charge;
		 double fill;
	};






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
	class FF_NonBonded: public ForcefieldBase{
	public:
		FF_NonBonded(WorkSpace &newwspace);

		virtual FF_NonBonded* clone() const;
		virtual ~FF_NonBonded();

		// -----------------
		// public parameters
		// -----------------
		int fast; // 0 normal 1 fast 2 faster etc..

		/// Obtain 1-4 scaling parameters from forcefield parameter set (ffps) or use defaults 
		bool AutoScaling14;
		
		/// Interaction cutoff distance [Angstrom]
		double Cutoff;

		/// Inner cutoff. If switching functions are used, the interaction energies will
		/// start to taper of at the InnerCutoff and reach 0 at 'Cutoff'. [Angstrom]
		double InnerCutoff;

		/// calculate forces/energies for vdw forces ? (12-6 Lennard-Jones Potential)
		bool DoVdw;

		/// Interaction cutoff distance for VdW energies [Angstrom]
		double VdwCutoff;

		/// Inner cutoff distance for VdW energies [Angstrom]
		double VdwInnerCutoff;

		/// Factor for scaling 1-4 VdW interactions (many forcefields use factor < 1.0)
		double Vdw14Scaling;

		/// Factor for scaling all VdW interactions. Default = 1.0
		double VdwAttenuator;

		/// Use Longrange VdW correction by integrating forces analytically over 
		/// average particle density space outside the cutoff distance
		bool   VdwCor;

		/// VdW correction density of particles [Particles/Angstrom^3]
		double VdwCorDensity;

		/// VdW correction Lennard Jones Well depth of particles [kcal/mol]
		double VdwCorEpsilon;
		/// VdW correction Lennard Jones Well radius of interaction [Angstrom]
		double VdwCorRadius;

		/// VdW correction: From which distance to start intergrating (set this equal to VdwCutoff) [Angstrom]
		double VdwCorCutoff; 

		/// Calculate electrostatics
		bool DoElec;
		/// Factor for scaling 1-4 electrostatic interactions (many forcefields use factor < 1.0)
		double Elec14Scaling;

		/// Factor for scaling all electrostatic interactions. Default = 1.0
		double ElecAttenuator;

		/// Dielectric constant 
		double Dielectric;

		/// Use a Distance dependent dielectric constant ? Form is e(r) = Dielectric + alpha * r
		bool   DDDielectric;

		/// Slope of distance dependent dielectric
		double DDDielectricAlpha;

		/// do force shifting ?
		///
		/// The method used here (only for electrostatic forces) is from
		/// P.J. Steinbach & B.R. Brooks, New Spherical-Cutoff Methods for Long-Range Forces
		/// in Macromolecular Simulation, (1994) J. Comp. Chem. 15, 667-683.
		bool   ForceSwitch;

		/// do potential shifting ?
		bool   EnergySwitch;

		bool   UsePartialRecalc;

		bool   IgnoreIntraResidue;

		/// prints a little block of parameter information
		void info() const;

		void printStandard_LJ_Coulomb_Force(){
			double vdw_potential;
			double vdw_force;
			double elec_potential;
			double elec_force;
			for( double Dist_ij = 2.5; Dist_ij < Cutoff ; Dist_ij += 0.025){
				calc_LJ_Coulomb_Force( Dist_ij, 1/Dist_ij , 3.0, 1.0 / (PhysicsConst::J2kcal * PhysicsConst::Na), 1, 1, 
						vdw_potential, vdw_force,
						elec_potential, elec_force );

				printf(" %8.4f  %8.4f  %8.4f  %8.4f  %8.4f \n",
					Dist_ij,
					vdw_potential  * PhysicsConst::J2kcal * PhysicsConst::Na ,
					vdw_force      * PhysicsConst::J2kcal * PhysicsConst::Na * 1E-10, 
					elec_potential * PhysicsConst::J2kcal * PhysicsConst::Na,
					elec_force     * PhysicsConst::J2kcal * PhysicsConst::Na * 1E-10
				);

			}

		}

		void calc_LJ_Coulomb_Force(
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
		);

    double getEElec() const { return epot_elec; };
    double getEVdw()  const { return epot_vdw; };

	protected:

		double epot_elec;
		double epot_vdw; 	

		// variables to do with partial updating of the energies
		bool fullrecalc;
		unsigned m_Nlist_FullUpdateCount;
		double epot_reference;
		std::vector <size_t> changed_atom;
		std::vector <size_t> old_changed_atom;

		NonBonded_Pack *local_atomparam;

		// a local store for the basis vectors
		Maths::dvector basisvector[32];

		ForcefieldBase::AtomicVerbosity verbose_level;

		void settodefault();

		// Include all variables who's change should trigger a resetup
		virtual unsigned long calcCheckSum();

		void infoLine() const; 

		void infoLineHeader() const;
	
		virtual void setup();
		virtual void setupBasisVectors();

		virtual void calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level);

		virtual void calcEnergies();
		virtual void calcEnergies_Full(); 
		virtual void calcEnergies_Update();

		virtual void calcForces();
	};







//-------------------------------------------------
//
/// \brief Monitors Vdw Energy in an FF_NonBonded Forcefield
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
	class PD_API Monitor_FF_NonBonded_EVdw: public Monitors::MonitorBase
	{
	public:

		/// Obtain a pointer to the FF_Bonded in the constructor
		Monitor_FF_NonBonded_EVdw(const FF_NonBonded &_ffnonbonded)
		{
			name = "Monitor_FF_NonBonded_EVdw";  // set a recognisable name
			ffnonbonded = &_ffnonbonded;         // remember the pointer to the bonded forcefield
		}

		virtual Monitor_FF_NonBonded_EVdw* clone() const { return new Monitor_FF_NonBonded_EVdw(*this); }

	protected:

		/// curdata just returns the current bonded energy to the monitor machinery
		void setcurdata()
		{
			// convert to natural units instead of SI units
			addData(  ffnonbonded->getEVdw() * Physics::PhysicsConst::J2kcal * Physics::PhysicsConst::Na );
		}

	private:
		const FF_NonBonded *ffnonbonded;
	};








//-------------------------------------------------
//
/// \brief  Monitors Electrostatic Energy in an FF_NonBonded Forcefield
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
	class PD_API Monitor_FF_NonBonded_EElec: public Monitors::MonitorBase
	{
	public:

		/// Obtain a pointer to the FF_Bonded in the constructor
		Monitor_FF_NonBonded_EElec(const FF_NonBonded &_ffnonbonded)
		{
			name = "Monitor_FF_NonBonded_EElec";  // set a recognisable name
			ffnonbonded = &_ffnonbonded;         // remember the pointer to the bonded forcefield
		}

		virtual Monitor_FF_NonBonded_EElec* clone() const { return new Monitor_FF_NonBonded_EElec(*this); }

	protected:

		/// curdata just returns the current bonded energy to the monitor machinery
		void setcurdata()
		{
			// convert to natural units instead of SI units
			addData(  ffnonbonded->getEElec() * Physics::PhysicsConst::J2kcal * Physics::PhysicsConst::Na );
		}

	private:
		const FF_NonBonded *ffnonbonded;
	};
}

#endif


