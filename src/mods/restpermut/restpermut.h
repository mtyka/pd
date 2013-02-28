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

#ifndef __RESTPERMUT_H
#define __RESTPERMUT_H

#include "forcefields/forcefield.h"
#include "forcefields/restraintbase.h"

#include <valarray>

void restpermut_cflush();

namespace Physics
{
	// Highly special class 





//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API Restraint_PermuteSolvent:public RestraintForcefieldBase{
	public:

		Restraint_PermuteSolvent(WorkSpace &newwspace): 
			RestraintForcefieldBase(newwspace) 
		{
			name = "VoroSolventRest";
			ShortName = "Voro";
			stepCounter = 0;
			UpdateAssignment = 1;
			UpdateCostmatrix = 16;
			k = 0;
			Power = 2;
			MolName = "TIP3";
			RefAtomName = "O";
			
			Sym_2_i=-1;
			Sym_2_j=-1;
			EneRestrictToFirst = 1000000;
			VoroRestrictToFirst = 1000000;
		}

		virtual Restraint_PermuteSolvent* clone() const { return new Restraint_PermuteSolvent(*this); }

		double Power;
		unsigned UpdateAssignment;
		unsigned UpdateCostmatrix;
		std::string MolName;
		std::string RefAtomName;

		// Hacks for water
		int Sym_2_i,Sym_2_j;
		int EneRestrictToFirst;
		int VoroRestrictToFirst;

		virtual void info() const{};					///< prints a little block of parameter information

		void saveCurrentAtomPositions();

		void setKMul(size_t i, double kmul){ kmul_store[i] = kmul; }
		/// resets the current assignement to the original 1 to 1
		void resetAssignment();   

		void benchAssignment(size_t nsteps){
			resetAssignment();
			printf(" ---> %d ",ref_atom_2_real_atom.size() );
			for(size_t i=0;i<nsteps;i++){
				updateAssignment();
			}
		}

		void benchAssignment_Munkres(size_t nsteps){
			resetAssignment();
			printf(" ---> %d ",ref_atom_2_real_atom.size() );
			for(size_t i=0;i<nsteps;i++){
				updateAssignment_Munkres();
			}
		}
	protected:
		virtual void setup();
		virtual void calcForces();
		virtual void calcEnergies();
		virtual double calcEnergyAtQ(double newQ){ return 0.0;};

		void calcNewCostMatrix();
		void updateCostMatrix();
		void updateAssignment();
		void updateAssignment_Munkres();

		int nmolsize;
		std::valarray<int>	 costm;			// cost matrix
		std::vector<int>		 atom_index;		
		std::vector<int>		 atom_index_reassigned;
		std::vector<Maths::dvector> atom_constraint_pos;
		std::vector<int>		 ref_atom_2_real_atom;
		std::vector<int>		 real_atom_2_ref_atom;
	
		std::vector<double>	 kmul_store;	// stores the k	multipliers
		unsigned stepCounter;
	};



}

#endif



