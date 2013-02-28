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

#include "alignment.h"

#include "tools/stringtool.h"
#include "tools/stringbuilder.h"
#include "library/residues.h"

using namespace Library;
using namespace std;

namespace Sequence
{
	BestPathBase::BestPathBase( int gapPenalty, size_t xDimension, size_t yDimension ) : m_GapPenalty(gapPenalty),
		m_XDimension(xDimension),
		m_YDimension(yDimension),
		m_ScoreMatrix(xDimension,yDimension),
		m_PathStoreMatrix((int)xDimension-1,(int)yDimension-1,2) // has to be 1 smaller than the scoreMartix, the 3rd dimension holds the 2 coordinates i and j of the desired next cell in the scoreMatrix
	{
	}

	int* BestPathBase::getEquiv( AlignmentDef &_AlignDef )
	{
		return _AlignDef.m_Equiv;
	}

	void BestPathBase::getBestPath()
	{
		int IPlusOne, JPlusOne, bestI, bestJ;
		double bestValue;

		FillScoreMatrix();

		for( int i = m_XDimension - 2; i >= 0; i-- ) // counting backwards from the far corner of the table and up one diagonally
		{
			for( int j = m_YDimension - 2; j >= 0; j-- )
			{
				IPlusOne = i + 1; // the plus ones are the next cell down diagonally, i.e. the first evaluation for best path ..
				JPlusOne = j + 1;
				bestI = IPlusOne; // we will initially asume that the diagonal is the best initial path, to be edited in a moment if not
				bestJ = JPlusOne;
				bestValue = m_ScoreMatrix[bestI][bestJ]; // so what is the assumptions' value

				for( int countI = i + 2; countI < m_XDimension; countI++ )
				{
					double additionScore = m_ScoreMatrix[countI][JPlusOne] - m_GapPenalty;
					// gap penalty is invoked on anything other than the direct diagonal;
					if( additionScore > bestValue )
					{
						bestI = countI;
						bestJ = JPlusOne;
						bestValue = additionScore;
					}
				}

				for( int countJ = j + 2; countJ < m_YDimension; countJ++ )
				{
					double additionScore = m_ScoreMatrix[IPlusOne][countJ] - m_GapPenalty;
					if( additionScore > bestValue )
					{
						bestI = IPlusOne;
						bestJ = countJ;
						bestValue = additionScore;
					}
				}

				m_PathStoreMatrix[i][j][0] = bestI; // this records an ID for the cell we found for following the path backwards
				m_PathStoreMatrix[i][j][1] = bestJ;

				m_ScoreMatrix[i][j] += bestValue; // the best possible scoring total from the path from the current i,j pairing
			}
		}

		// we now need to find the best starting cell and set m_BestPathStartCellIDX, and m_BestPathStartCellIDY
		double bestScoreMatrixValue = m_ScoreMatrix[0][0]; // begin with this cell
		m_BestPathStartCellIDX = 0;
		m_BestPathStartCellIDY = 0;
		// is it the best, probably not...

		for( int i = 1; i < m_XDimension; i++ )
		{
			if( m_ScoreMatrix[i][0] > bestScoreMatrixValue )
			{
				m_BestPathStartCellIDX = i;
				bestScoreMatrixValue = m_ScoreMatrix[i][0];
			}
		}
		for( int j = 1; j < m_YDimension; j++ )
		{
			if( m_ScoreMatrix[0][j] > bestScoreMatrixValue )
			{
				m_BestPathStartCellIDX = 0; // we neeed to be in the 1st row or the 1st column
				m_BestPathStartCellIDY = j;
				bestScoreMatrixValue = m_ScoreMatrix[0][j];
			}
		}
		// The m_PathStoreMatrix has been found and the starting cell for the path has been defined.
		// Our work here is done!
	}



	void AlignmentDef::ResetAlignment()
	{	
		m_Aligned = false;
		if( m_Equiv == NULL || m_AllocatedMemory < m_Seq1->size() )
		{
			delete[] m_Equiv;
			m_AllocatedMemory = (size_t) m_Seq1->size();
			m_Equiv = new int[m_AllocatedMemory];
			memset(m_Equiv,-1,sizeof(int)*m_AllocatedMemory);			
		}
		else
		{
			memset(m_Equiv,-1, (sizeof(int) * m_AllocatedMemory));
		}
	}

	void AlignmentDef::CloneAlignment( const AlignmentDef& _Clone )
	{
		if( _Clone.m_Equiv != NULL )
		{
			m_Aligned = _Clone.m_Aligned;
			m_AllocatedMemory = _Clone.m_AllocatedMemory;
			m_Equiv = new int[m_AllocatedMemory];
			memcpy(m_Equiv,_Clone.m_Equiv,sizeof(int)*m_AllocatedMemory);
		}
		else
		{
			m_Aligned = false; // if the equiv array is NULL, how can this have ever been aligned?
			m_AllocatedMemory = 0;
			m_Equiv = NULL;
		}
	}

	int AlignmentDef::nextEquivFromPoint( int lookFrom )const
	{
		size_t length = m_Seq1->size();
		for( size_t i = lookFrom; i < length; i++ )
		{
			if( m_Equiv[i] != -1 ) return i;
		}
		return -1;
	}

	void AlignmentDef::printToScreen()const
	{
		const int lineHeaderOverhead = 6; // Being the: "S1: \'" and "'"
		const int totalWidth = 75;
		const int chopSize = totalWidth - lineHeaderOverhead; // the resultant width of the output!

		const Sequence::BioSequence& seq1 = getSeq1();
		const Sequence::BioSequence& seq2 = getSeq2();

		// Initilaise the builder arrays
		StringBuilder seq1Builder;
		StringBuilder seq2Builder;

		// setup the starts to be in-line
		int firstIDFoundAt = nextEquivFromPoint(0);
		if( firstIDFoundAt == -1 )
		{
			cout << endl << "Alignment Identity: ???%" << endl << endl;
			cout << "S1: \'" << seq1Builder << '\'' << endl;
			cout << "FAIL : The equivalency list is empty, no alignment can be shown!" << endl;
			cout << "S2: \'" << seq2Builder << '\'' << endl;
			return;
		}

		// the first pair deemed to be structurally equivelent
		int mol1Index = firstIDFoundAt;
		int mol2Index = m_Equiv[mol1Index];

		seq1Builder.append( seq1[mol1Index] );
		seq2Builder.append( seq2[mol2Index] );

		mol1Index--; // set the cursor to one backwards
		mol2Index--;
		// now backtrack to fill the start
		while( mol1Index >= 0 || mol2Index >= 0 )
		{
			if( mol1Index >= 0 )
			{
				seq1Builder.insert(0, seq1[mol1Index] );
			}
			else
			{
				seq1Builder.insert(0,'-');
			}
			if( mol2Index >= 0 )
			{
				seq2Builder.insert(0, seq2[mol2Index] );
			}
			else
			{
				seq2Builder.insert(0,'-');
			}

			mol1Index--; // work our way back to the start of the sequences
			mol2Index--;
		}

		// now fill to the end
		mol1Index = firstIDFoundAt;
		mol2Index = m_Equiv[ mol1Index ];
		int mol1NextIndex;
		int mol2NextIndex;

		while( (mol1NextIndex = nextEquivFromPoint( mol1Index + 1 ) ) != -1 ) // -1 is returned once the equivs have ended
		{
			// get mol2's index
			mol2NextIndex = m_Equiv[ mol1NextIndex ];

			// somehow fill between the equivs if there are residues in between
			mol2Index = m_Equiv[ mol1Index ]; // the next residues

			// POTENTIAL ERROR: Found this in later debug, whats is for ???
			// mol1Index = mol1Index; // the next residues

			while( true )
			{
				mol2Index++;
				mol1Index++;
				if( mol1Index < mol1NextIndex && mol2Index < mol2NextIndex )
				{
					seq1Builder.append( seq1[mol1Index] );
					seq2Builder.append( seq2[mol2Index] );
				}
				else if( mol1Index < mol1NextIndex )
				{
					seq1Builder.append( seq1[mol1Index] );
					seq2Builder.append('-');
				}
				else if( mol2Index < mol2NextIndex )
				{
					seq1Builder.append('-');
					seq2Builder.append( seq2[mol2Index] );
				}
				else
				{
					break; // we are done here, exit condition
				}
			}

			// now add the equiv add an equiv
			seq1Builder.append( seq1[mol1NextIndex] );
			seq2Builder.append( seq2[mol2NextIndex] );

			mol1Index = mol1NextIndex;
			mol2Index = m_Equiv[ mol1Index ];
		}

		// now fill to the end
		while( mol1Index < seq1.size() || mol2Index < seq2.size() )
		{
			mol1Index++; // work our way back to the start of the sequences
			mol2Index++;
			if( mol1Index < seq1.size() && mol2Index < seq2.size() )
			{
				seq1Builder.append( seq1[mol1Index] );
				seq2Builder.append( seq2[mol2Index] );
			}
			else if( mol1Index < seq1.size() )
			{
				seq1Builder.append( seq1[mol1Index] );
				seq2Builder.append('-');
			}
			else if( mol2Index < seq2.size() )
			{
				seq1Builder.append('-');
				seq2Builder.append( seq2[mol2Index] );
			}
		}

		// Now calculate the sequence identity!
		StringBuilder identityBuilder;
		double identity = 0;
		double corespondance = 0;
		char sb1, sb2;
		for( size_t i = 0; i < seq1Builder.size(); i++ )
		{
			sb1 = seq1Builder[i];
			if( sb1 != '-' )
			{
				sb2 = seq2Builder[i];
				if( sb2 != '-' )
				{
					if( sb1 == sb2 )
					{
						identityBuilder.append('*');
						corespondance += 1.0;
						identity += 1.0;
					}
					else
					{
						if( HasSameResidueProperty(sb1,sb2) )
						{
							identityBuilder.append('.');
							identity += 1.0;
						}
						else
						{
							identityBuilder.append(' ');
						}
					}
				}
				else
				{
					identityBuilder.append(' ');
				}
			}
			else
			{
				identityBuilder.append(' ');
			}
		}

		identity = ((identity / (double)seq1Builder.size()) * 100.0);
		corespondance = ((corespondance / (double)seq1Builder.size()) * 100.0);
		cout << "Alignment Identity     : " << identity << "%" << endl;
		cout << "Alignment Corespondance: " << corespondance << "%" << endl;
		cout << "S1 Length              : " << seq1.size() << endl;
		cout << "S2 Length              : " << seq2.size() << endl;
		cout << endl;

		size_t length = seq1Builder.size();
		int itterations = length / chopSize;

		for( int i = 0; i < itterations; i++ )
		{
			cout << "S1: \'" << seq1Builder.toString(i*chopSize,chopSize) << '\'' << endl;
			cout << "    \'" << identityBuilder.toString(i*chopSize,chopSize) << '\'' << endl;
			cout << "S2: \'" << seq2Builder.toString(i*chopSize,chopSize) << '\'' << endl << endl;
		}

		int remainder = length % chopSize;
		cout << "S1: \'" << seq1Builder.toString(itterations*chopSize,remainder) << '\'' << endl;
		cout << "    \'" << identityBuilder.toString(itterations*chopSize,remainder) << '\'' << endl;
		cout << "S2: \'" << seq2Builder.toString(itterations*chopSize,remainder) << '\'' << endl << endl;
	}




	bool ExpSeqPair::isValid() const
	{
		// This should check that the biological sequence (if you take gaps into account)
		// has a 1:1 residue correspondance with the structural sequence - otherwise something fishy
		// is going on ...

		//THROW(NotImplementedException); - it's not, but an exception is undesirable, cos we want the code to work...

		return true;
	}

	void AlignerBase::Align( AlignmentDef &_AlignDef ) const
	{
		const Sequence::BioSequence &seq1 = _AlignDef.getSeq1();
		const Sequence::BioSequence &seq2 = _AlignDef.getSeq2();
		int len1 = seq1.size();
		int len2 = seq2.size();

		if( _AlignDef.m_Aligned == true )
		{
			printf("Warning 'AlignerBase::Align()' called on a definition with a previous alignment. This will be overridden!\n");
		}
		_AlignDef.ResetAlignment(); // This also ensures that we have enough memory allocation in the equiv array.
		if( len1 == 0 || len2 == 0 )
		{
			// Do sweet feck all.
		}
		else if( len1 == 1 || len2 == 1 )
		{
			// length of 2 is the minimum limit for best_path
			PoxyAlign(_AlignDef); 
		}
		else
		{
			AlignCore(_AlignDef);
		}
		_AlignDef.m_Aligned = true;
	}

	int* AlignerBase::getEquiv( AlignmentDef &_AlignDef  ) const
	{
		return _AlignDef.m_Equiv;
	}
	
	void AlignerBase::PoxyAlign( AlignmentDef &_AlignDef ) const
	{
		// Used if one of the sequences is so short that the 'BestPathBase' internal matrix would crash.
		// Just do a direct correspondnce alignment - nothing even remotely fancy.

		int* equiv = getEquiv(_AlignDef);
		const Sequence::BioSequence &seq1 = _AlignDef.getSeq1();
		const Sequence::BioSequence &seq2 = _AlignDef.getSeq2();
		int len1 = seq1.size();
		int len2 = seq2.size();
		
		// Sequence lengths are asserted in the align() call
		ASSERT( len1 != 0 && len2 != 0, CodeException, "Sequence length must be greater than 0 for alignment.");
		ASSERT( len1 == 1 || len2 == 1, CodeException, "One sequence has to be of length 1 to use PoxyAlign()");
		
		if( len1 == 1 )
		{
			for( int i = 0; i < len2; i++ )
			{
				if( seq1[0] == seq2[i] )
				{
					equiv[0] = i;
					break;
				}
			}
		}
		else
		{
			for( int i = 0; i < len1; i++ )
			{
				if( seq1[i] == seq2[0] )
				{
					equiv[i] = 0;
					break;
				}
			}
		}
	}

	SimpleBestPath::SimpleBestPath( int gapPenalty, AlignmentDef &_AlignDef ) : 
		m_AlignDef(&_AlignDef),
		BestPathBase( gapPenalty, _AlignDef.getSeq1().size(), _AlignDef.getSeq2().size() )
	{
	}

	void SimpleBestPath::FillEquivList()
	{
		// now we need to fill the _equiv from the 'best path' through the newly summed m_ScoreMatrix
		int* _equiv = getEquiv( *m_AlignDef );

		// init locals
		int iPointer = m_BestPathStartCellIDX; // these were recorded in the getBestPath() in the baseclass call...
		int jPointer = m_BestPathStartCellIDY;

		m_AlignDef->ResetAlignment();

		// lets move through the path array to set the Equivelencies in '_equiv'
		// while() i.e. until we fall off the end of the table
		while( iPointer < ( m_XDimension - 1 ) && jPointer < ( m_YDimension - 1 ) ) 
		{
			// mark the equiv
			_equiv[ iPointer ] = jPointer;

			// get the next cell IDs
			// TempBufferIPointer :- needed as the jPointer-lookup needs the older iPointer
			int TempBufferIPointer = m_PathStoreMatrix[iPointer][jPointer][0];
			jPointer = m_PathStoreMatrix[iPointer][jPointer][1];
			iPointer = TempBufferIPointer;
		}

		// mark the equiv
		_equiv[ iPointer ] = jPointer;

		// all equivs are now set, we are done
	}

	void SimpleBestPath::FillScoreMatrix()
	{
		const BioSequence &seq1 = m_AlignDef->getSeq1();
		const BioSequence &seq2 = m_AlignDef->getSeq2();

		// make the scorematrix
		for( int i = 0; i < m_XDimension; i++ )
		{
			for( int j = 0; j < m_YDimension; j++ )
			{
				if( seq1[i] == seq2[j] )
				{
					m_ScoreMatrix[i][j] = 2.0;
				}
				else
				{
					m_ScoreMatrix[i][j] = 0.0;
				}
			}
		}
	}

	SimpleAligner::SimpleAligner() : AlignerBase()
	{
	}

	void SimpleAligner::AlignCore( AlignmentDef &_AlignDef ) const
	{
		// _AlignDef.ResetAlignment() is called by Align() in base class
		
		SimpleBestPath bp( 1, _AlignDef );
		bp.getBestPath(); // Class 'BestPathBase' calls 'FillScoreMatrix()' during getBestPath()
		bp.FillEquivList();
	}

	ManualAligner::ManualAligner() : AlignerBase()
	{
	}

	void ManualAligner::addMatch( size_t _seq1Index, size_t _seq2Index )
	{
		ManualMatch m;
		m.length = 1;
		m.seq1Index = _seq1Index;
		m.seq2Index = _seq2Index;
		m_matches.push_back( m );
	}

	void ManualAligner::addMatch( size_t _seq1StartIndex, size_t _seq2StartIndex, size_t _length )
	{
		ManualMatch m;
		m.length = _length;
		m.seq1Index = _seq1StartIndex;
		m.seq2Index = _seq2StartIndex;
		m_matches.push_back( m );
	}

	void ManualAligner::clearmatches()
	{
		m_matches.clear();
	}

	void ManualAligner::AlignCore( AlignmentDef &_AlignDef ) const
	{
		// _AlignDef.ResetAlignment() is called by Align() in base class

		int *equiv = getEquiv( _AlignDef );
		for( size_t i = 0; i < m_matches.size(); i++ )
		{
			const ManualMatch& mm = m_matches[i];

			if( mm.length == 0 ) 
				continue; // ignore
			
			ASSERT(mm.seq1Index+mm.length < _AlignDef.getSeq1().size(), ProcedureException, "Residue definition is out of range!");
			ASSERT(mm.seq2Index+mm.length < _AlignDef.getSeq2().size(), ProcedureException, "Residue definition is out of range!");
			bool error = false;
			for( int j = 0; j < mm.length; j++ )
			{
				if( equiv[j+mm.seq1Index] != -1 )
				{
					error = true;
				}
				equiv[j+mm.seq1Index] = j+mm.seq2Index;
			}
			if( error )
			{
				printf("Warning: Manual alignment contains overlapping sections! Overwriting has occured in stretch:\n");
				mm.info();
			}
		}
	}

	DirectAligner::DirectAligner()
	{
	}

	void DirectAligner::AlignCore( AlignmentDef &_AlignDef ) const
	{
		// _AlignDef.ResetAlignment() is called by Align() in base class

		int *equiv = getEquiv( _AlignDef );

		int len1 = _AlignDef.getSeq1().size();
		int len2 = _AlignDef.getSeq2().size();
		
		for( int i = 0; i < len1; i++ )
		{
			if( i < len2 )
			{
				equiv[i] = i;
			}
			else
			{
				break; // we are past the end of seq2
			}
		}
	}
}

