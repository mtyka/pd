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

#include "backbonetorsions.h"

// Namespace includes
using namespace Maths;

namespace Library
{
	// Poxy bubble-sort, but its not rate limiting, so we dont care...
	void BackboneTorsionSet::SortAngleSet()
	{
		int i,j;
		bool changed;
		while(true) // swap higher to lower
		{
			changed = false;
			for( i = 0; i < m_AngleCount; i++ )
			{
				for( j = i + 1; j < m_AngleCount; j++ )
				{
					if( m_Propensities[j] > m_Propensities[i] )
					{
						SwapAll( j, i );
						changed = true;
					}
				}
			}
			if( !changed ) break;
		}
	}

	void BackboneTorsionSet::SwapAll( int i, int j )
	{
		// std::swap(m_AngleSetID[j], m_AngleSetID[i]) : NOOOOO!
		std::swap(m_AngleClass[j], m_AngleClass[i]);
		std::swap(m_Phis[j], m_Phis[i]);
		std::swap(m_Psis[j], m_Psis[i]);
		std::swap(m_Omegas[j], m_Omegas[i]);
		std::swap(m_Propensities[j], m_Propensities[i]);
	}

	BackboneTorsionLibrary::BackboneTorsionLibrary() 
		: BackboneTorsionSet()
	{
	}

	bool BackboneTorsionLibrary::finalise()// if it is a Library then this set should be initialised with residueangle propensities that have a sum of 1.00
	{
		ASSERT( m_IntendedMaxSize == m_AngleCount, ProcedureException, "Not all angles have been added to the BackboneTorsionLibrary. Finalise must only be called internally after all angles have been added.");
		if( !ValidateLibraryPropensities() ) return false; // fail, the sum is not very near 1.00 (we tollerate small rounding deviations to 0.99/1.01, and scale to 1.00 below)
		ScalePropensities();
		SortAngleSet();
		m_finalised = true; // flag this
		return true;
	}

	bool BackboneTorsionLibrary::init( char resSingleID, const std::string& resName, int numberOfAnglesToBeAdded )
	{
		if( isFinalised() )
		{
			printf("ResiduesAngles cannot be initialised, it has already been finalised!\n");
			return false;
		}
		if( isInitialised() )
		{
			printf("ResiduesAngles cannot be initialised more than once!\n");
			return false;
		}
		// check the anglecount
		if( numberOfAnglesToBeAdded < 0 )
		{
			printf("ERROR: BackboneTorsionSet cannot be initialised with a count less than 0! Residue: %s\n",resName.c_str());
			return false;
		}

		m_IntendedMaxSize = numberOfAnglesToBeAdded;

		// set the internal members
		m_AngleCount = 0;
		m_ResName = resName;
		m_ResSingleID = resSingleID;

		// allocate memory for the residue angleset on the heap
		m_AngleSetID.resize(m_IntendedMaxSize,-1);
		m_AngleClass.resize(m_IntendedMaxSize,'\0');

		m_Phis.resize(m_IntendedMaxSize,DBL_MAX);
		m_Psis.resize(m_IntendedMaxSize,DBL_MAX);
		m_Omegas.resize(m_IntendedMaxSize,DBL_MAX);
		m_Propensities.resize(m_IntendedMaxSize,DBL_MAX);

		return true;
	}

	bool BackboneTorsionLibrary::addAngleGroup( char angleClass, double phi, double psi, double omega, double propensity )
	{
		// error check
		if( !isInitialised() )
		{
			printf("CODE SEQUENCE ERROR: ResiduesAnglesLibrary has not been initialised, addAngleGroup() cannot be called!\n");
			return false;
		}
		if( isFinalised() )
		{
			printf("CODE SEQUENCE ERROR: ResiduesAnglesLibrary has been finalised, no further calls to addAngleGroup() are allowed!\n");
			return false;
		}
		if( m_IntendedMaxSize <= m_AngleCount )
		{
			printf("ERROR: in ResiduesAnglesLibrary::addAngleGroup(): Maximum angle contents overflow!\n");
			return false;
		}

		// make the assignments
		m_AngleSetID[ m_AngleCount ] = m_AngleCount;
		m_AngleClass[ m_AngleCount ] = angleClass;
		m_Phis[ m_AngleCount ] = phi;
		m_Psis[ m_AngleCount ] = psi;
		m_Omegas[ m_AngleCount ] = omega;
		m_Propensities[ m_AngleCount ] = propensity;

		m_AngleCount++;
		if( m_IntendedMaxSize == m_AngleCount )
		{
			finalise();
		}

		return true;
	}

	bool BackboneTorsionLibrary::ValidateLibraryPropensities()
	{
		double sum = 0.0;

		for( int i = 0; i < m_AngleCount; i++ )
		{
			sum += m_Propensities[i];
		}

		//if( sum != 1.00 ) double equality testing is evil >;-)
		if( !SigFigEquality( sum, 1.00,3 ) )
		{
			double sumDiff = 1.00 - sum;
			if( sumDiff < -0.02 || sumDiff > 0.02 )
			{
				printf("ERROR: Propensity sum for '%c' (%3s) is '%8.6lf' not 1.00!\n", m_ResSingleID,m_ResName.c_str(),sum);
				return false; // the difference is too big, we will tollerate 0.99 ish due to rounding, no more
			}
			else
			{
				printf("WARNING: Propensity sum for '%c' (%3s) is '%8.6lf' not 1.00!\n", m_ResSingleID,m_ResName.c_str(),sum);
			}
		}
		return true;
	}








	BackboneTorsionSubSet::BackboneTorsionSubSet()
		: BackboneTorsionSet()
		, m_SourceLibrary( 0 ) // null pointer
	{
	}

	bool BackboneTorsionSubSet::finalise()
		// this subset should be finalised to hold propensities that have a sum of 1.00, scale this to make it true.
	{
		ScalePropensities();
		m_finalised = true; // flag this
		return true;
	}

	bool BackboneTorsionSubSet::init( const BackboneTorsionLibrary &sourceLibrary )
	{
		if( isFinalised() )
		{
			printf("ResiduesAngles cannot be initialised, it has already been finalised!\n");
			return false;
		}
		if( isInitialised() )
		{
			printf("ResiduesAngles cannot be initialised more than once!\n");
			return false;
		}

		// set the internal members
		m_SourceLibrary = &sourceLibrary;
		m_IntendedMaxSize = m_SourceLibrary->size();
		m_AngleCount = 0;

		m_ResName = m_SourceLibrary->getResName();
		m_ResSingleID = m_SourceLibrary->getResSingleID();

		m_AngleSetID.clear();
		m_AngleClass.clear();

		m_Phis.clear();
		m_Psis.clear();
		m_Omegas.clear();
		m_Propensities.clear();

		return true;
	}

	bool BackboneTorsionSubSet::addAngleGroup( int sourceLibraryIndexToAdd )
	{
		// error check
		if( !isInitialised() )
		{
			THROW(CodeException,"ResiduesAngles has not been initialised, addAngleGroup() cannot be called!");
		}
		if( isFinalised() )
		{
			THROW(CodeException,"ResiduesAngles has been finalised, no further calls to addAngleGroup() are allowed!");
		}
		if( m_AngleCount >= m_IntendedMaxSize )
		{
			THROW(CodeException,"Overflow in BackboneTorsionSet::addAngleGroup(): Maximum theoretical angle contents overflow!");
		}
		if( sourceLibraryIndexToAdd < 0 || sourceLibraryIndexToAdd >= m_SourceLibrary->size() )
		{
			THROW(CodeException,"BackboneTorsionSet::addAngleGroup(): AngleLib index is invalid!");
		}
		for( int i = 0; i < m_AngleCount; i++ )
		{
			if( m_AngleSetID[i] == sourceLibraryIndexToAdd )
			{
				THROW(CodeException,"BackboneTorsionSet::addAngleGroup(): The library already contains that index, it cannot be added twice!");
			}
		}

		// make the assignments
		m_AngleSetID.push_back( m_SourceLibrary->getAngleSetID(sourceLibraryIndexToAdd) );
		m_AngleClass.push_back( m_SourceLibrary->getAngleClass(sourceLibraryIndexToAdd) );
		m_Phis.push_back( m_SourceLibrary->getPhi(sourceLibraryIndexToAdd) );
		m_Psis.push_back( m_SourceLibrary->getPsi(sourceLibraryIndexToAdd) );
		m_Omegas.push_back( m_SourceLibrary->getOmega(sourceLibraryIndexToAdd) );
		m_Propensities.push_back( m_SourceLibrary->getPropensity(sourceLibraryIndexToAdd) );

		m_AngleCount++;

		return true;
	}












	BackboneTorsionSet::BackboneTorsionSet()
		: m_ResSingleID('\0') // IMPORNANT : initialise to '\0' to flag the class as uninitialised
		, m_finalised( false ),
		m_AngleCount(0)
	{
	}

	void BackboneTorsionSet::ScalePropensities()
	{
		double scalingFactor = 0.0;
		for( int i = 0; i < m_AngleCount; i++ )
		{
			scalingFactor += m_Propensities[i];
		}
		scalingFactor = 1.0000 / scalingFactor;
		for( int i = 0; i < m_AngleCount; i++ )
		{
			m_Propensities[i] *= scalingFactor;
		}

		scalingFactor = 0.0;
		for( int i = 0; i < m_AngleCount; i++ )
		{
			scalingFactor += m_Propensities[i];
		}
		if( !SigFigEquality( scalingFactor, 1.00, 10 ) )
		{
			THROW(CodeException,"ScalePropensities() failed!");
		}
	}

	void BackboneTorsionSet::info() const
	{
		if( !isInitialised() )
		{
			THROW(CodeException,"Cannot ResiduesAngles::info(), the class has not been initialised!");
		}
		if( !isFinalised() )
		{
			THROW(CodeException,"Cannot ResiduesAngles::info(), the class has not been finalised!");
		}

		printf(" %c %3s AngleCount:%d\n",m_ResSingleID,m_ResName.c_str(),m_AngleCount);
		printf(" #     Phi     Psi   Omega  Propens Class\n");

		for( int i = 0; i < m_AngleCount; i++ )
		{
			printf( " %d %7.3lf %7.3lf %7.3lf %8.5lf %c\n",m_AngleSetID[i],m_Phis[i],m_Psis[i],m_Omegas[i],m_Propensities[i],m_AngleClass[i]) ;
		}

		printf("\n");
	}


	int BackboneTorsionSet::getClosestAngleGroup(
		double currentPhi,
		double currentPsi,
		double currentOmega
		)const
	{
		int bestID = -1;
		double bestDistance = DBL_MAX;
		double phiDiff, psiDiff, omegaDiff;
		for( int i = 0; i < m_AngleCount; i++ )
		{
			// calculate the distances
			// keep in mind that angle space is circular, and therefore we have to check the
			// distances accross the 180 degree boundary ...
			phiDiff = currentPhi - m_Phis[i];
			if( phiDiff < 0.0 ) phiDiff = -phiDiff;
			if( phiDiff > Maths::MathConst::PI ) phiDiff = Maths::MathConst::TwoPI - phiDiff;

			psiDiff = currentPsi - m_Psis[i];
			if( psiDiff < 0.0 ) psiDiff = -psiDiff;
			if( psiDiff > Maths::MathConst::PI ) psiDiff = Maths::MathConst::TwoPI - psiDiff;

			omegaDiff = currentOmega - m_Omegas[i];
			if( omegaDiff < 0.0 ) omegaDiff = -omegaDiff;
			if( omegaDiff > Maths::MathConst::PI ) omegaDiff = Maths::MathConst::TwoPI - omegaDiff;

			// use square distance to avoid an expensive square root...
			double sqrDistance = sqr( phiDiff ) + sqr( psiDiff ) + sqr( omegaDiff );
			if( sqrDistance < bestDistance )
			{
				bestDistance = sqrDistance;
				bestID = i;
			}
		}

		return bestID; // all done :-D
	}

	int BackboneTorsionSet::getClosestAngleGroup(
		double currentPhi, double &closestPhi,
		double currentPsi, double &closestPsi,
		double currentOmega, double &closestOmega
		)const
	{
		int bestID = -1;
		double bestDistance = DBL_MAX;
		double phiDiff, psiDiff, omegaDiff;
		for( int i = 0; i < m_AngleCount; i++ )
		{
			// calculate the distances
			// keep in mind that angle space is circular, and therefore we have to check the
			// distances accross the 180 degree boundary ...
			phiDiff = currentPhi - m_Phis[i];
			if( phiDiff < 0.0 ) phiDiff = -phiDiff;
			if( phiDiff > Maths::MathConst::PI ) phiDiff = Maths::MathConst::TwoPI - phiDiff;

			psiDiff = currentPsi - m_Psis[i];
			if( psiDiff < 0.0 ) psiDiff = -psiDiff;
			if( psiDiff > Maths::MathConst::PI ) psiDiff = Maths::MathConst::TwoPI - psiDiff;

			omegaDiff = currentOmega - m_Omegas[i];
			if( omegaDiff < 0.0 ) omegaDiff = -omegaDiff;
			if( omegaDiff > Maths::MathConst::PI ) omegaDiff = Maths::MathConst::TwoPI - omegaDiff;

			// use square distance to avoid an expensive square root...
			double sqrDistance = sqr( phiDiff ) + sqr( psiDiff ) + sqr( omegaDiff );
			if( sqrDistance < bestDistance )
			{
				bestDistance = sqrDistance;
				bestID = i;
			}
		}

		closestPhi = m_Phis[bestID];
		closestPsi = m_Psis[bestID];
		closestOmega = m_Omegas[bestID];

		return bestID; // all done :-D
	}
} // namespace Library


