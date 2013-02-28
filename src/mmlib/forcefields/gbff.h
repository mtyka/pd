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

#ifndef __GBFF_H
#define __GBFF_H

#include <vector>

#include "forcefields/forcefield.h" // provides base class
#include "forcefields/nonbonded.h"  // provides constructor reference argument 
#include "monitors/monitorbase.h"   // provides constructor reference argument 

#include "maths/maths.fwd.h"        // forward declarations

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
	class PD_API ContinuumElectrostatic
	{
	public:
		ContinuumElectrostatic()
		{
			DielectricSolvent = 80.0;
			DielectricSolute = 1.0;
		}
		virtual ~ContinuumElectrostatic(){};
	public:
		double DielectricSolvent;
		double DielectricSolute;
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
	class PD_API FF_GeneralizedBorn :
		public FF_NonBonded,
		public ContinuumElectrostatic
	{
	public:
		FF_GeneralizedBorn( WorkSpace &newwspace ): 
		  FF_NonBonded( newwspace ) {}
		  FF_GeneralizedBorn(  const FF_NonBonded &_clone ):
		  FF_NonBonded( _clone ) {}
		  virtual ~FF_GeneralizedBorn() {}

		  virtual FF_GeneralizedBorn* clone() const { return new FF_GeneralizedBorn(*this); }

		  double getEPol() const { return epot_pol; }

	protected:
		double epot_pol; // born energy
	};



	//-------------------------------------------------
	//
	/// \brief Generalised Born Forcefield (Implicit Solvation)
	/// \details
	/// Description: Generalised Born / (Surface Area is separate) Implicit Solvent
	/// Implements a "fastgbsa" function which calculates vdw,
	/// vacuum electrostatics and gb forces at once for highest performance
	/// Born Radii are calculated using the Clarke Still Approximation [2]
	/// \version 0.1
	/// \author Michael Tyka
	///
	/// References:
	///
	/// These classes/functions are based on the following papers:
	/// [1] W. Clark Still, Anna Tempczyk, Ronald C. Hawley and Thomas
	/// Hendrickson, Semianalytical Treatment of Solvation for Molecular
	/// Mechanics and Dynamics, J Am. Chem. Soc. 1990, 112, 6127-6129
	/// [2] Di Qui, Peter S. Shenkin, Frank P. Hollinger and
	/// W. Clark Still The GB/SA Continuum Model for Solvation. A fast
	/// Analytical Method for the calculation of Approximate Born Radii
	/// J. Phys. Chem A 1997, 101, 3005-3014
	/// [3] Vickie Tsui & David A. Case Biopolymers (Nucleic Acid Sciences),
	/// Theory and Applications of the generalized Born Solvation Model
	/// in Macromolecular Simulations Vol 56, 275-291 (2001)
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class PD_API FF_GeneralizedBorn_Still : public FF_GeneralizedBorn
	{
	public:
		FF_GeneralizedBorn_Still(	WorkSpace &newwspace );
		virtual FF_GeneralizedBorn_Still* clone() const { return new FF_GeneralizedBorn_Still(*this); }
		virtual ~FF_GeneralizedBorn_Still(){}

		bool FastMode; // ultrafast mode doing vdw and estat as well in one function
		double DielectricSolvent;
		double DielectricSolute;
		double GbsaStillCutoff;
		double ExpApproxThreshold;
		double DielectricOffset;
		double BornRadiusOffset;

	protected:
		virtual void setup();

		virtual void calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level);
		virtual void calcEnergies();
		virtual void calcForces();

		virtual void info() const; // prints a little block of parameter information
		virtual void infoLine() const; // prints a line of current energies
		virtual void infoLineHeader() const; // prints the headers for the above function

		struct GB_AtomType
		{
			double radius;
		};

		// General Atom Parameters and a little scratch space.
		struct GB_Atom_Param
		{
			double bornradius;
			double charge;
			double radiusij;
			double epsilon;
			double deda;
			Maths::dvector position;
			Maths::dvector dedi;
		};

		// this is stuff for the Di Qiu, Peter S. Shenkin, Frank P. Hollinger, and W. Clark Still
		// analytical born radii pairwise approximation
		struct GB_Atom_Param_Still
		{
			double terms123;
			double Vj;
		};

		std::vector<GB_AtomType> GBtype;
		std::vector<GB_Atom_Param_Still> GB_atom_param_still;
		std::vector<GB_Atom_Param> GB_atom_param; 

		double bornpremul;
		double P1, P2, P3, P4, P5;

		double epot_covalent;
		double epot_pol_self; // components of the above
		double epot_pol_cross;

		int readGeneralisedBornSolvationSection();

		// Born Radii calculations
		int calcFixedBornRadiiTerms();
		int calcBornRadii_constant(double constBornRadius);
		int calcBornRadii_PairwiseApprox();

		void calcBornEnergy_covalentTerms();
		void calcBornEnergy();
		void calcBornEnergy_verbose(ForcefieldBase::AtomicVerbosity verboselevel);
		void calcBornForces();
		void calcBornForces_numerical(double dc, int dtype);

		// Specialized versions of the above that assume the bornradii are constant
		void calcBornForces_constantBornRadii();
		void calcForces_constantBornRadii_IncludingVacuo();

		void calcForcesIncludingVacuo();
	};


	//-------------------------------------------------
	//
	/// \brief  Monitors Electrostatic Energy in an FF_GeneralizedBorn Forcefield
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
	class PD_API Monitor_FF_GeneralizedBorn_EPol: public Monitors::MonitorBase
	{
	public:

		/// Obtain a pointer to the FF_Bonded in the constructor
		Monitor_FF_GeneralizedBorn_EPol(const FF_GeneralizedBorn &_ffgb)
		{
			name = "Monitor_FF_GeneralizedBorn_EPol"; // set a recognisable name
			ffgb = &_ffgb; // remember the pointer to the bonded forcefield
		}

		virtual Monitor_FF_GeneralizedBorn_EPol* clone() const { return new Monitor_FF_GeneralizedBorn_EPol(*this); }

	protected:

		/// curdata just returns the current bonded energy to the monitor machinery
		void setcurdata()
		{
			// convert to natural units instead of SI units
			addData(  ffgb->getEPol() * Physics::PhysicsConst::J2kcal * Physics::PhysicsConst::Na );
		}

	private:
		const FF_GeneralizedBorn *ffgb;
	};

} // namespace Physics

#endif

