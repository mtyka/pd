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

#ifndef __TORSIONAL_MINIMISATION_H
#define __TORSIONAL_MINIMISATION_H

// Essential Headers
#include "protocols/protocolbase.h" // Provides a base class
#include "system/atomrangestore.h" // Provides a class member
#include "workspace/rotbond.h" // Provides a base class
#include "maths/maths.fwd.h"
#include "workspace/workspace.fwd.h"

class PD_API RotatableBond;

namespace Protocol
{
	/// \brief A helper class for TorsionalMinimisation
	/// \details Not intended for use by the user
	/// \author Jon Rea 
	class PD_API RotationDefinition_TM: public RotationDefinition
	{
	public:
		RotationDefinition_TM(): RotationDefinition()
		{
		}

		// extra required parameters, set during initialisation in the parent class
		double Value;
		double PrevValue;
		int A; ///< index of atom2
		int B; ///< index of atom3
		double ReciprocalLen;
		size_t rangeLink;
		inline void rotateByValue() { perturbRotation( Value ); }
	};

	//-------------------------------------------------------------------
	/// \brief Performs a torsional steepest-descent energy minimisation
	/// \details 
	//// Implementation of a Torsional Minimisation that works within the confines of small section 
	//// (allowing for breaks).
	/// \author Jon Rea 
	class PD_API TorsionalMinimisation: public RangesProtocolBase
	{
	public: // public constructors and functions

		enum RotationExcludeMode // flagged enumeration
		{
			NothingAtAll = 0,
			ImproperSingles = 1,
			Omega = 2,
			AllValid = ImproperSingles | Omega,
			Default = AllValid
		};

		TorsionalMinimisation( Physics::Forcefield &_ff );

		/// Minimise only rotatable bonds which involve atoms contained by this atomic range
		TorsionalMinimisation( Physics::Forcefield &_ff, const PickAtomRange& _newRange ); 
		TorsionalMinimisation( Physics::Forcefield &_ff, const PickAtomRanges& _newRange );
		virtual TorsionalMinimisation* clone() const 
		{ 
			return new TorsionalMinimisation(*this); 
		}
		virtual void settodefault();


		/// the main core call
		virtual int runcore(); 


		/// prints a little block of parameter information
		virtual void info() const; 

		/// prints a line of current energies
		virtual void infoLine() const; 

		/// prints the headers for the above function
		virtual void infoLineHeader() const; 

		inline int getStepNumber()const { return Step; }


		/// Change in energy below which to exit the minimisation. <= 0.0 to disable.
		double SlopeCutoff; 

		/// the initial AngleCap and SCAngleCap are set to this fraction of the IntendedMaxAngleCap and IntendedMaxSCAngleCap respectively
		double InitialCapFactor; 

	private:

		// Private static constants.
		static const double RESET_PREV_SQUARE;

		// internal rotation definitions and initialisation function

		/// initialise the arrays below
		void initRotations(); 
		void setupRotDefPhiPsi( RotatableBond &rotBond, size_t forRange );
		void setupRotDefChi( RotatableBond &rotBond );


		/// The atomic range that was last assigned to this instance. runCore checks this to see if reinitialisation needs to occur
		size_t m_InitPickerSerial; 

		/// Used to store the best conformation during the minimisation.
		PosStore m_BestStore; 


		/// Array of available chi rotation definitions
		std::vector<RotationDefinition_TM> ChiRotDefs; 

		/// Array of available phi/psi rotation definitions
		std::vector<RotationDefinition_TM> PhiPsiRotDefs; 

		// Internal main loop functions

		/// Cconverts the current forces to angle change factors
		void ApplyXYZAbsoluteDotVector(); 

		/// Modifies the previous factors with respect to the previous gradients, 
		/// then normalises and scales the anglechange to obtain the valuephis, valuepsis and valuechis
		void ActionConjugateGradient(); 

		// Step info

		/// the current Step number
		int Step; 

		/// count of successive uphill Steps
		int Uphill; 

		// Energy info

		/// the potential energy of the previous Step
		double previous_epot; 

		/// values for conjugate gradient calc
		double PrevSquare; 

		// simulation parameters		

		/// Torsional mode: which torsions can the minimisation move
		RotationExcludeMode RotationExcldMode; 

		/// starting and max value of angle Step sizes
		double IntendedMaxAngleCap; 

		/// ditto above
		double IntendedMaxSCAngleCap; 

		/// the max angle in radians that a backbone torsion can move per-Step
		double SCAngleCap; 

		/// the normalisation factor for sidechain rotations
		double SCAngleFactor; 

		/// the max angle in radians that a sidechain torsion can move per-Step
		double AngleCap; 

		/// the factor to increase caps by when going downhill
		double CapRise; 

		/// the factor to decrease caps by when going uphill		
		double CapDrop; 
	};
}

#endif

