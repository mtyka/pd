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

#ifndef ANGLESET_H
#define ANGLESET_H

#include <string>
#include <vector>
#include <fstream>

#include "backbonetorsions.h"

namespace Library
{
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
	class PD_API AngleSet
	{
	public:
		// constructor destructor logic
		AngleSet( const std::string& fileName );
		AngleSet();

		// FileIO
		void loadFromFile( const std::string& _FileName ); ///< initialise the angleset

		// Class cleanup
		void clear(); // return the angleset to an uninitialised state

		//public functions
		void info() const;
		bool getIsDefined( char resID ) const;

		// public accessors
		inline bool isInitialised() const { return m_Initialised; } ///< return if the angleset is usable
		inline size_t size() const { return m_BackboneTorsionSets.size(); }
		inline const std::string& getResidueIDs() const { return m_StoredIDs; }
		const BackboneTorsionLibrary& getBackboneTorsionSet( size_t index ) const;
		const BackboneTorsionLibrary& getBackboneTorsionSet( char resID ) const;

		/// getClosestAngleGroup returns the index of the closest angle definition for the given residue ID.
		/// -1 is returned if there was an error, for example the residue is not present in the angleset.
		/// An overloaded version will also return the angles of the closest member
		int getClosestAngleGroup( char molID,
			double currentPhi, double &closestPhi,
			double currentPsi, double &closestPsi,
			double currentOmega, double &closestOmega
			) const;

		/// getClosestAngleGroup returns the index of the closest angle definition for the given residue ID.
		/// -1 is returned if there was an error, for example the residue is not present in the angleset.
		int getClosestAngleGroup( char molID,
			double currentPhi,
			double currentPsi,
			double currentOmega
			) const;

		// helper functions
		static void EnsureRadianRange( double &angle );

	private:
		// Internal initialisation flag
		bool m_Initialised; ///< internal flag that the set is ready for use

		// Internal member data
		std::vector<BackboneTorsionLibrary> m_BackboneTorsionSets;

		// File IO data
		std::string m_Filename;
		std::string m_Descriptor; ///< The third line of the input file, info on what the angleset is ...
		std::string m_StoredIDs; ///< Store the IDs of the residues in the angleset in the order we have them

		// File IO functions
		bool AssertHeader( std::ifstream& inputFile ); ///< Asserts that the file header is as expected and is the correct version
		bool AssertNotEoF( std::ifstream& inputFile ); ///< Asserts we are not at the end of the current file passed in ...
		int getResidueCountFromFile( std::ifstream& inputFile ); ///< Returns the number of residues in the angle set. -1 on error
		bool readAngleSetBlock( std::ifstream& inputFile, int index );
		bool readAngleSetBlocks( std::ifstream& inputFile, int blockCount ); ///< Calls readAngleSetBlock() per block, the number of blocks is defined in the header
	}; // class PD_API AngleSet
} // namespace Library

#endif

