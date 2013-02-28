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

#ifndef __LCPO_H
#define __LCPO_H

#include "sasabase.h" // provides base class

// Maximal number of neighbors in the SASA Neighbor_list
#define Nlistmax 80

namespace Physics
{
	//-------------------------------------------------
	//
	/// \brief Description: Surface Area calculations Using Numerical Counting, the linear
	/// approximate LCPO algorithm [1] 
	/// Energies/Forces are calculated using Atomic Solvation Paraemters
	///
	/// \details 
	///    
	/// References:
	/// [1] Jorg Weiser, Peter S. Shenkin, W. Clark Still
	/// Approximate Atomic Surfaces from Linear Combinations of Pairwise Overlaps (LCPO)
	/// J. Comp. Chem., Vol. 20, No. 2, 217-230 (1999)
	///
	/// \author Mike Tyka  
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///

	// **NOTE**: Should derive from public 'FF_SurfaceArea' when re-fitted
	class PD_API FF_SASA_LCPO : public ForcefieldBase 		
	{
	public:
		FF_SASA_LCPO(WorkSpace &newwspace) 
			: ForcefieldBase(newwspace)
		{
			name = "LCPO SASA";
			settodefault();
		}

		virtual FF_SASA_LCPO* clone() const { return new FF_SASA_LCPO(*this); }

		std::string ASPsection_name;
		double GlobalASP;
		int UpdateSasa;

		virtual void settodefault()
		{
			ASPsection_name = "";
			GlobalASP = 0.0;
			UpdateSasa = 5;
		};

	protected:
		virtual void setup();
		virtual void calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level);
		virtual void calcEnergies();
		virtual void calcForces();

		virtual void info() const;           ///< prints a little block of parameter information
		virtual void infoLine() const;       ///< prints a line of current energies
		virtual void infoLineHeader() const; ///< prints the headers for the above function

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
		class PD_API SASA_Atom
		{
		public:		
			SASA_Atom() : sigma(0)
			{
			}
			int attype;               ///< link to the atom type
			double P1, P2, P3, P4;    ///< LCPO Parameters
			double radius;            ///< Vdw Radius + 1.4
			double SASA;              ///< SASA of individual atom
			Maths::dvector SASAderiv; ///< derivative of SASA (SASAderiv * sigma = force)
			double sigma;             ///< atomic solvation parameter
			signed char neighbors;    ///< Nr of neighbors |???
			char use;                 ///< include in calculation
		};

		class PD_API member
		{
		public:
			int size;
			int i[Nlistmax];          ///< index of neighbor
			double s[Nlistmax];       ///< this atoms overlap with its neighbor i
		};

	private:
		std::vector<SASA_Atom> SASAtype;
		std::vector<SASA_Atom> SASAatom;
		std::vector<member> Nlist;

		double epot_cav;
		double totalSASA;

		/// Can be used to implement forcefields like Eisenberg et al., Ooi et al etc..
		int readASPSection(const std::string &sectionname);

		/// This reads in the algorithm parameters (LCPO), called by setup();
		int readLCPOSection();

		// Solvent Accesable Surface Area calculations
		int calcNumericalSASA();
		int benchmarkNumericalSASA();

		int calcLCPOSASA();
		int calcLCPOSasaEnergies(bool recalc_nlist=true);
		int calcLCPOSasaForces_num();
		int calcLCPOSasaForces(bool dofullcalc);

		void testDerivatives();
	};

} // namespace Physics

#endif

