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
#include "angleset.h"

// namespace includes
using namespace Maths;

namespace Library
{
	AngleSet::AngleSet( const std::string& fileName )
		: m_Initialised(false)
	{
		loadFromFile( fileName );
	}

	AngleSet::AngleSet()
		: m_Initialised(false)
	{
	}

	void AngleSet::clear()
	{
		m_BackboneTorsionSets.clear();
		m_Filename.clear();
		m_Descriptor.clear();
		m_StoredIDs.clear();
		m_Initialised = false;
	}

	void AngleSet::loadFromFile( const std::string& _FileName )
	{
		m_Filename = _FileName;

		if( m_Initialised )
		{
			THROW(ProcedureException,"AngleSet is already initialised, call clear() to reset.");
		}
		else
		{
			printf("Loading angleset from file: '%s'...", m_Filename.c_str());
		}

		// File initialisation
		std::ifstream inputFile(m_Filename.c_str(), std::ifstream::in);
		if( !inputFile.is_open() ) throw IOException(".AngleSet file not found: '" + m_Filename + "'!" );

		// load the 'inputFile' contents
		if( AssertHeader(inputFile) )
		{
			// passed the header OK ...
			int resCount = getResidueCountFromFile(inputFile);
			if( resCount != -1 )
			{
				m_StoredIDs.resize(resCount,'\0'); // Allocate memory to hold the desired residues
				// the obtained residue count was valid and >= 0
				if( readAngleSetBlocks( inputFile, resCount ) )
				{
					printf(" Done! %s\n",m_Filename.c_str());
					m_Initialised = true; // set the initialised flag ...

					inputFile.close(); // we are done with the 'inputFile', close the file handle

					return; // We have sucessfully read the whole file
				}
			}
		}

		// make sure we clean up ...
		inputFile.close(); // we are done with the 'inputFile', close the file handle

		// Run past the if statements, there was a failure during parsing...
		THROW(ParseException,"Critical error in parsing the '.angleset' input file");
	}

	void AngleSet::info() const
	{
		printf("Verbose Angleset Dataprint:\n\n");
		printf(" File Descriptor: %s\n", m_Descriptor.c_str());
		printf(" AngleSet import produced '%d' residue angle set definitions: %s\n\n", m_BackboneTorsionSets.size(), m_StoredIDs.c_str());
		for( int i = 0; i < m_BackboneTorsionSets.size(); i++ )
		{
			const BackboneTorsionLibrary& res = m_BackboneTorsionSets[i];
			res.info();
		}
	}

	bool AngleSet::AssertNotEoF( std::ifstream& inputFile )
	{		
		if( inputFile.eof() )
		{
			printf( "ERROR: Unexpected end of input file duing parsing: <unknown filename>\n");
			return false;
		}
		else
		{
			return true;
		}
	}

	// This function is used to look at the 1st two lines of an angle set 'inputFile' termed the
	// "header". These lines should be %ANGLE_SET, followed by %VERSION 1.0. If otherwise then
	// do some moaning ...

	bool AngleSet::AssertHeader(std::ifstream& inputFile)
	{
		std::string lineBuffer;

		//%ANGLE_SET		
		std::getline(inputFile,lineBuffer);
		if( lineBuffer.compare(0,10,"%ANGLE_SET") != 0 )
		{
			printf("ERROR: The input file specified is not and angleset!\n");
			return false;
		}
		if( !AssertNotEoF(inputFile) ) return false;


		//%VERSION 1.0
		std::getline(inputFile,lineBuffer);
		double version;
		if( lineBuffer.size() < 10 )
		{
			printf("ERROR: Parsed line is too short!");
			return false;
		}
		int n = sscanf(&lineBuffer.c_str()[8],"%lf",&version);
		if( 0 != lineBuffer.compare(0,8,"%VERSION") || n != 1 )
		{
			printf("ERROR: The second input line of the angle set is not a valid version descriptor!\n");
			return false;
		}
		if( version != 1.0 )
		{
			printf("ERROR: The angle set version given is not 1.0!\n");
			return false;
		}
		if( !AssertNotEoF(inputFile) ) return false;


		//%DESCRIPTOR <- some description ->
		std::getline(inputFile,lineBuffer);
		if( 0 != lineBuffer.compare(0,12,"%DESCRIPTOR ") || n != 1 )
		{
			printf("ERROR: The third input line of the angle set is not a valid file \"%%DESCRIPTOR \"!\n");
			return false;
		}
		m_Descriptor = lineBuffer.substr(12,lineBuffer.size()-12);
		if( !AssertNotEoF(inputFile) ) return false;


		// and if all that was OK ...
		return true; // parsing the header was successful!
	}

	int AngleSet::getResidueCountFromFile( std::ifstream& inputFile)
	{
		std::string lineBuffer;

		//%RESIDUE_COUNT 20
		std::getline(inputFile,lineBuffer);
		if( lineBuffer.size() < 16 )
		{
			printf("ERROR: Parsed line is too short!");
			return -1;
		}
		int resCount;
		int n = sscanf(&lineBuffer.c_str()[14],"%d",&resCount);
		if( 0 != lineBuffer.compare(0,14,"%RESIDUE_COUNT") || n != 1 )
		{
			printf("ERROR: The fourth input line of the angle set is not a valid residuecount descriptor!\n");
			return -1;
		}
		if( resCount < 0 )
		{
			printf("ERROR: The third input line of the angle set contained an residuecount descriptor with a count less than 0: %d \n", resCount);
			return -1;
		}
		if( !AssertNotEoF(inputFile) ) return -1;

		// and if all that was OK ...
		return resCount; // parsing the RESIDUE_COUNT was successful!
	}

	bool AngleSet::readAngleSetBlocks( std::ifstream&inputFile, int blockCount )
	{
		for( int i = 0; i < blockCount; i++ )
		{
			if( !readAngleSetBlock(inputFile,i) )
			{
				return false;
			}
		}
		return true;
	}

	bool AngleSet::readAngleSetBlock( std::ifstream&inputFile, int index )
	{
		std::string lineBuffer;
		BackboneTorsionLibrary resAng;

		// A Ala 3
		// -140 153 180 0.130
		// -72 145 180 0.800
		// -122 117 180 0.070

		// get an residue angleset header line...
		std::getline(inputFile,lineBuffer);
		if( !AssertNotEoF(inputFile) ) return false;

		int n = 0;
		char resSingleID;
		char resID[4]; // 3 plus '\0' string termination
		int angCount = -1;

		// extract the variables
		if( lineBuffer.size() < 5 )
		{
			printf("\nERROR: Parsed line is too short!");
			return false;
		}
		n = sscanf(lineBuffer.c_str(),"%c %3s %d", &resSingleID, &resID[0], &angCount);
		if( n != 3 )
		{
			printf("\nERROR: The residue angleset header line did not parse correctly, should show format e.g. A Ala 3\n");
			return false;
		}

		resAng.init( toupper(resSingleID), resID, angCount );
		m_StoredIDs[ index ] = resSingleID; // assign this too for the external accessor and getIsDefined()

		double phi;
		double psi;
		double omega;
		double propensity;
		char angClass;

		for( int i = 0; i < angCount; i++ )
		{
			if( !AssertNotEoF(inputFile) ) return false;

			std::getline(inputFile,lineBuffer);

			if( lineBuffer.size() < 9 )
			{
				printf("\nERROR: Parsed line is too short!");
				return false;
			}
			n = sscanf( lineBuffer.c_str(), "%lf %lf %lf %lf %c", &phi, &psi, &omega, &propensity, &angClass );
			if( n != 5 )
			{
				printf("\nERROR: The residue angleset angle line did not parse correctly!\n");
				return false;
			}

			// ANGLE MESSING:
			// our conformer builder internally uses different torsion definitions to the angleset file.
			//
			// Omega is nommally the torsion of the peptide group, however the 'H' atom isnt always present.
			// (especially as they are not moved by our builder even if they are in the forcefield!)
			// Therefore the CA on res+1 is used instead of the H and the desired angle is flipped by 180 degrees.
			//
			// In addition to the omega definition change, the Psi definition is also altered to use the O of the
			// current residue, rather than the N of the next residue. This is because in the special case of
			// anchor-2 the N of the next residue is available, but its position is broken, and therefore not usable.
			// The O will however still be valid as it is rotated with the loop, and therefore should be used instead
			// to give a valid Psi angle.
			//
			// These two modifications are allowed as idealised geometry is used and the peptide group is planar:
			// i.e. these two torsion pairs are always related by a 180 degree difference when ideal...
			// also ... keep the angles in the range -Maths::MathConst::PI < angle < +Maths::MathConst::PI

			phi = DegToRad( phi );
			EnsureRadianRange( phi );

			psi = DegToRad( psi );
			psi -= MathConst::PI;
			EnsureRadianRange( psi );

			omega = DegToRad( omega );
			omega -= MathConst::PI;
			EnsureRadianRange( omega );

			if( !SigFigEquality( omega, 0.0, 4 ) && // using SigFigEquality as they only need to be roughly idealised to be good enough
				!SigFigEquality( omega, MathConst::PI, 4 ) &&
				!SigFigEquality( omega, -MathConst::PI, 4 ) // both +and- Maths::MathConst::PI values are fine ...
				)
			{
				printf("\nERROR: Due to some simplifying code assumptions used in the conformer builder code, omega angles are only allowed to be 0.0 or 180.0 degrees (within tollrance). i.e. only simple idealised Cis and Trans values are allowed\n");
				return false;
			}

			// toupper() : force upper case in internal data
			if( !resAng.addAngleGroup( toupper(angClass), phi, psi, omega, propensity ) )
			{
				return false;
			}
		}

		ASSERT(resAng.isFinalised(), CodeException, "\nERROR: BackboneTorsionLibrary is incomplete!");
		m_BackboneTorsionSets.push_back(resAng);

		// all done :-D
		return true;
	}

	void AngleSet::EnsureRadianRange( double &angle )
	{
		// ensure our angle is within the range +-Maths::MathConst::PI
		while( angle < -MathConst::PI )
		{
			angle += MathConst::TwoPI;
		}
		while( angle > MathConst::PI )
		{
			angle -= MathConst::TwoPI;
		}
	}

	bool AngleSet::getIsDefined( char resID )const
	{
		for( int i = 0; i < m_BackboneTorsionSets.size(); i++ )
		{
			if( m_StoredIDs[i] == resID )
			{
				return true;
			}
		}
		return false;
	}

	const BackboneTorsionLibrary& AngleSet::getBackboneTorsionSet( size_t index )const
	{
		if( index >= size() )
		{
			THROW(OutOfRangeException,"AngleSet::getBackboneTorsionSet()");
		}
		else
		{
			// all ok ... return it ...
			return m_BackboneTorsionSets[ index ];
		}
	}

	const BackboneTorsionLibrary& AngleSet::getBackboneTorsionSet( char resID )const
	{
		// search for that resID in the stored ID list ...
		for( int i = 0; i < m_BackboneTorsionSets.size(); i++ )
		{
			if( m_StoredIDs[i] == resID )
			{
				return m_BackboneTorsionSets[i];
			}
		}
		THROW(ArgumentException,"AngleSet::getBackboneTorsionSet() could not find the required resID within the angleset.");
	}

	int AngleSet::getClosestAngleGroup(
		char molID,
		double currentPhi,
		double currentPsi,
		double currentOmega
		)const
	{
		const BackboneTorsionLibrary& angleLib = getBackboneTorsionSet( molID );
		return angleLib.getClosestAngleGroup(
			currentPhi,
			currentPsi,
			currentOmega );
	}

	int AngleSet::getClosestAngleGroup(
		char molID,
		double currentPhi, double &closestPhi,
		double currentPsi, double &closestPsi,
		double currentOmega, double &closestOmega
		)const
	{
		const BackboneTorsionLibrary& angleLib = getBackboneTorsionSet( molID );
		return angleLib.getClosestAngleGroup(
			currentPhi, closestPhi,
			currentPsi, closestPsi,
			currentOmega, closestOmega );
	}
} // namespace Library


