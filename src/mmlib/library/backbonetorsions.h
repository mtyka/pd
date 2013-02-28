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

#ifndef BackboneTorsionSet_H
#define BackboneTorsionSet_H

#include <string.h>
#include "typedefs.h"

namespace Library
{
	//-------------------------------------------------
	//
	/// \brief  
	/// This class is used to house an angle set for a given residue. It is used as the base class PD_API for both the
	/// AngleSet and ConformerSet classes as a holder
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
	class PD_API BackboneTorsionSet
	{
	public:
		// constructor destructor and initialisation functions
		BackboneTorsionSet();
		virtual ~BackboneTorsionSet(){};

		virtual bool finalise() = 0; // pure virtual function, makes the class abstract

		// Inline public accessors
		inline bool isInitialised() const { return m_ResSingleID != '\0'; }
		inline bool isFinalised()const { return m_finalised; }

		// public per class PD_API properties
		inline size_t size() const { return m_AngleCount; }
		inline std::string getResName() const { return m_ResName; }
		inline char getResSingleID ()const { return m_ResSingleID; }

		// public per angle properties

		// individual
		// index should be of Type 'conformer_type' as they are indexed by a conformer made of 'conformer_type's
		inline char getAngleClass ( conformer_type index ) const { return m_AngleClass [ index ]; }
		inline int getAngleSetID ( conformer_type index ) const { return m_AngleSetID [ index ]; }
		inline double getPhi ( conformer_type index ) const { return m_Phis [ index ]; }
		inline double getPsi ( conformer_type index ) const { return m_Psis [ index ]; }
		inline double getOmega ( conformer_type index ) const { return m_Omegas [ index ]; }
		inline double getPropensity ( conformer_type index )const { return m_Propensities [ index ]; }

		// whole arrays
		inline const std::vector<double> getPhis () const { return m_Phis; }
		inline const std::vector<double> getPsis () const { return m_Psis; }
		inline const std::vector<double> getOmegas()const { return m_Omegas; }

		// query functions
		int getClosestAngleGroup(
			double currentPhi, double &closestPhi,
			double currentPsi, double &closestPsi,
			double currentOmega, double &closestOmega
			) const;

		int getClosestAngleGroup(
			double currentPhi,
			double currentPsi,
			double currentOmega
			) const;

		// User output functions
		void info() const;

	protected:
		bool m_finalised; ///< Flagged once finalisation has occured
		size_t m_IntendedMaxSize; ///< The eventual number of definitions - finalise() will check this

		// residue data
		char m_ResSingleID; ///< single letter ResidueID
		std::string m_ResName; ///< 3 letter string for that residueID
		size_t m_AngleCount; ///< number of angles held in this set

		std::vector<int>    m_AngleSetID; ///< number representing the index in the original angleset library (should be 0,1,2,3,etc. in the library itself, but doesnt have to be in subsets)
		std::vector<char>   m_AngleClass; ///< alpha / beta / Left-handed-turn / ( '_' == null )
		std::vector<double> m_Phis; ///< phi torsion angle data
		std::vector<double> m_Psis; ///< phi torsion angle data
		std::vector<double> m_Omegas; ///< omega torsion angle data
		std::vector<double> m_Propensities; ///< angle propensity information

		void ScalePropensities(); ///< If the propensities dont add up to 1.00, scale them so that they do ... this is mainly for use in angleset sub-sets which will not include all the original angles, but also so that we can ensure that logic requiring a total probability of 1.00 will actually work.
		void SortAngleSet(); ///< The internal angle sets should be sorted by order of propensity. Slow sort algorithm, but easy, and not a limiting procedure

	private:
		// swapping functions for SortAngleSet()
		void SwapAll( int indexI, int indexJ );
	};


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
	class PD_API BackboneTorsionLibrary: public BackboneTorsionSet
	{
	public:
		BackboneTorsionLibrary();

		bool finalise();
		bool init( char resSingleID, const std::string& resName, int maxAngleContent );
		bool addAngleGroup( char angleClass, double phi, double psi, double omega, double propensity );

	private:
		bool ValidateLibraryPropensities();
	};


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
	class PD_API BackboneTorsionSubSet: public BackboneTorsionSet
	{
	public:
		BackboneTorsionSubSet();

		bool finalise();
		bool init( const BackboneTorsionLibrary &sourceLibrary ); // maxAngleContent is set as the count in the donor angle set
		bool addAngleGroup( int libraryIndex );

		// public accessors
		inline size_t libSize() const { return m_SourceLibrary->size(); } ///< Return the number of definitions in the parent library
	
	private:
		const BackboneTorsionLibrary * m_SourceLibrary;
	};
} // namespace Library

#endif

