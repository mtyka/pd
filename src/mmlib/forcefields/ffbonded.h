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

#ifndef __FFBONDED_H
#define __FFBONDED_H

#include "system/fundamentals.h"
#include "forcefields/forcefield.h" // Provides the base class
#include "monitors/monitorbase.h" // Provides the base class
#include "workspace/workspace.fwd.h"

class PD_API Particle;

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
	class PD_API Bond: public IndexPairDistance{ /// bond i-j
	public:
		Bond(){};
		Bond(int _i, int _j, double _l, double _k):
		IndexPairDistance(_i,_j,_l),k(_k) {}
		double k; ///< force constant
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
	class PD_API Angle: public IndexTriplet{
	public:
		double k; /// force constant
		double theta0; /// ideal value
		double theta; /// double value
		double l; /// distance I<->K
		double d; /// distance I<->K
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
	class PD_API Torsion: public IndexQuartet
	{
	public:
		char Type; ///< 0 for torsion, 1 for improper
		int terms; ///< how many fourier terms
		double phi; ///< current angle
		double Vn[4]; ///< the individual parameters for each term
		double n[4];
		double gamma[4];
	};

	enum ForcefieldScope
	{
		EntireSystem,
		OnlyBackbone
	};

	class PD_API TorsionalRestraint;

	//-------------------------------------------------
	//
	/// \brief calculates Bonded Forcefield Components: Bonds Stretching, Angle Bending, Torsion & Impropers
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
	class PD_API FF_Bonded : public ForcefieldBase
	{
	public:

		// This is baaad - can we do this some other way ?
		friend class PD_API FF_Restraint_Torsional;

		FF_Bonded( WorkSpace &newwspace ): 
		ForcefieldBase( newwspace ) 
		{
			name = "Bonded Forcefield";
			ShortName = "Bonds";
			Scope = EntireSystem;
			DoBonds = true;
			DoAngles = true;
			DoTorsions = true;
			DoImpropers = true;
			BondAttenuator = 1.0;
			AngleAttenuator = 1.0;
			TorsionAttenuator = 1.0;
			ImproperAttenuator = 1.0;
			IgnoreForcefieldParams = false;
		}

		virtual FF_Bonded* clone() const { return new FF_Bonded(*this); }

		inline int getNumBonds() { return (int)bond.size(); }
		inline int getNumAngles() { return (int)angle.size(); }
		inline int getNumTorsions() { return (int)torsion.size(); }
		inline int getNumImpropers(){ return (int)improper.size(); }

		Physics::Bond getBond( int i, int j ) const;

		bool isBonded( int i, int j ) const; ///< Returns true if atoms i and j are bonded in the workspace this forcefield is setup on

		int findBond(int i, int j) const; ///< Finds the index of the bond between atoms i and j in the workspace this forcefield is setup on, or returns -1 if there is no such bond.
		int findAngle(int i, int a, int j) const;
		int findTorsion(int i, int a, int b, int j) const;
		int findImproper(int i, int a, int b, int j) const;

		// Output/Debug Functions
		void printBondedParameters();
		int printBondedParameters(char *filename);
		void printBondedParameters(FILE * target);

		int printPSFfile_bondedparams(FILE *file);

		// Parameters
		ForcefieldScope Scope;

		bool DoBonds;
		bool DoAngles;
		bool DoTorsions;
		bool DoImpropers;

		double BondAttenuator;
		double AngleAttenuator;
		double TorsionAttenuator;
		double ImproperAttenuator;

		void setIgnoreForcefieldParams(bool _val = true){ IgnoreForcefieldParams = _val; }

		// takes the current bondlengths/angles as equilibrium values;
		void resetEquilibriumState() { resetEquilibriumBonds(); resetEquilibriumAngles(); }

		void resetEquilibriumBonds();
		void resetEquilibriumAngles();
		void resetEquilibriumHarmonicDihedrals();

		virtual void info() const; // prints a little block of parameter information

		double getEBond() const    { return epot_bond; }
		double getEAngle() const   { return epot_angle; }
		double getETorsion() const { return epot_torsion; }
		double getEImproper() const{ return epot_improper; }

	protected:
		bool IgnoreForcefieldParams;
		virtual void setup();

		virtual void calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level);
		virtual void calcEnergies();
		virtual void calcForces();

		virtual void infoLine() const; // prints a line of current energies
		virtual void infoLineHeader() const; // prints the headers for the above function

		int assembleBondList();
		int assembleAngleList();
		int assembleTorsionList();
		int assembleImproperList();

		std::vector<Bond> bond;
		std::vector<Angle> angle;
		std::vector<Torsion> torsion;
		std::vector<Torsion> improper;

		double epot_bond;
		double epot_angle;
		double epot_torsion;
		double epot_improper;

		void calcBondForces();
		void calcBondEnergies_Verbose();

		void calcAngleForces();
		void calcAngleEnergies_Verbose();

		void calcTorsionForces();
		void calcTorsionEnergies_Verbose();

		void calcImproperForces();
		void calcImproperEnergies_Verbose();

	};

#ifndef SWIG
	void calcDihedralForcesNonVerbose(WorkSpace &wspace, Torsion &dihedral, double &epot_dihedral);
	void calcDihedralForcesVerbose(WorkSpace &wspace, Torsion &dihedral, double &epot_dihedral);
	void calcDihedralForcesNonVerbosePassive(WorkSpace &wspace, Torsion &dihedral, double &epot_dihedral);
	void calcDihedralForcesVerbosePassive(WorkSpace &wspace, Torsion &dihedral, double &epot_dihedral);
#endif

	// --------------------------------------------------------------------------------------------------
	// Monitors and other auxillary classes -------------------------------------------------------------

	//-------------------------------------------------
	//
	/// \brief  Monitors Bond Energy in an FF_Bonded Forcefield
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
	class PD_API Monitor_FF_Bonded_EBond: public Monitors::MonitorBase
	{
	public:

		/// Obtain a pointer to the FF_Bonded in the constructor
		Monitor_FF_Bonded_EBond(FF_Bonded &_ffbonded)
		{
			name = "mon_ff_bonded_bonds";  // set a recognisable name
			ffbonded = &_ffbonded;         // remember the pointer to the bonded forcefield
		}

		virtual Monitor_FF_Bonded_EBond* clone() const { return new Monitor_FF_Bonded_EBond(*this); }

	protected:

		/// curdata just returns the current bonded energy to the monitor machinery
		void setcurdata()
		{
			// convert to natural units instead of SI units
			addData(  ffbonded->getEBond() * Physics::PhysicsConst::J2kcal * Physics::PhysicsConst::Na );
		}

	private:
		FF_Bonded *ffbonded;
	};


	//-------------------------------------------------
	//
	/// \brief Monitors Angle Energy in an FF_Bonded Forcefield
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
	class PD_API Monitor_FF_Bonded_EAngle: public Monitors::MonitorBase
	{
	public:

		/// Obtain a pointer to the FF_Bonded in the constructor
		Monitor_FF_Bonded_EAngle(FF_Bonded &_ffbonded)
		{
			name = "mon_ff_bonded_bonds";  // set a recognisable name
			ffbonded = &_ffbonded;         // remember the pointer to the bonded forcefield
		}

		virtual Monitor_FF_Bonded_EAngle* clone() const { return new Monitor_FF_Bonded_EAngle(*this); }

	protected:

		/// curdata just returns the current bonded energy to the monitor machinery
		void setcurdata()
		{
			// convert to natural units instead of SI units
			addData(  ffbonded->getEAngle() * Physics::PhysicsConst::J2kcal * Physics::PhysicsConst::Na );
		}

	private:
		FF_Bonded *ffbonded;
	};


	//-------------------------------------------------
	//
	/// \brief  Monitors Torsion Energy in an FF_Bonded Forcefield
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
	class PD_API Monitor_FF_Bonded_ETorsion: public Monitors::MonitorBase
	{
	public:

		/// Obtain a pointer to the FF_Bonded in the constructor
		Monitor_FF_Bonded_ETorsion(FF_Bonded &_ffbonded)
		{
			name = "mon_ff_bonded_bonds";  // set a recognisable name
			ffbonded = &_ffbonded;         // remember the pointer to the bonded forcefield
		}

		virtual Monitor_FF_Bonded_ETorsion* clone() const { return new Monitor_FF_Bonded_ETorsion(*this); }

	protected:

		/// curdata just returns the current bonded energy to the monitor machinery
		void setcurdata()
		{
			// convert to natural units instead of SI units
			addData(  ffbonded->getETorsion() * Physics::PhysicsConst::J2kcal * Physics::PhysicsConst::Na );
		}

	private:
		FF_Bonded *ffbonded;
	};


	//-------------------------------------------------
	//
	/// \brief  Monitors Improper Energy in an FF_Bonded Forcefield
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
	class PD_API Monitor_FF_Bonded_EImproper: public Monitors::MonitorBase
	{
	public:

		/// Obtain a pointer to the FF_Bonded in the constructor
		Monitor_FF_Bonded_EImproper(FF_Bonded &_ffbonded)
		{
			name = "mon_ff_bonded_bonds";  // set a recognisable name
			ffbonded = &_ffbonded;         // remember the pointer to the bonded forcefield
		}

		virtual Monitor_FF_Bonded_EImproper* clone() const { return new Monitor_FF_Bonded_EImproper(*this); }

	protected:

		/// curdata just returns the current bonded energy to the monitor machinery
		void setcurdata()
		{
			// convert to natural units instead of SI units
			addData(  ffbonded->getEImproper() * Physics::PhysicsConst::J2kcal * Physics::PhysicsConst::Na );
		}

	private:
		FF_Bonded *ffbonded;
	};
} // namespace Physics

#endif

