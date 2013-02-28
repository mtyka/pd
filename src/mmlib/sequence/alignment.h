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

#ifndef __ALIGNMENT_H
#define __ALIGNMENT_H

// Essential Headers
#include "maths/tntjama/tnt_array2d.h"
#include "maths/tntjama/tnt_array3d.h"
#include "maths/tntjama/tnt_array2d_utils.h"
#include "maths/tntjama/tnt_array3d_utils.h"
#include "sequence/sequence.h"

namespace Sequence
{
	class PD_API AlignmentDef; // Forward declaration






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
	class PD_API BestPathBase // Pure abstract base class
	{
	public:
		BestPathBase( int gapPenalty, size_t xDimension, size_t yDimension );
		virtual void FillScoreMatrix() = 0; // override for in derived classes for custom score behaviour
		void getBestPath(); // finds the sum route through the matrix and then fills the m_PathStoreMatrix to record the best path

	protected:
		int* getEquiv( AlignmentDef &_AlignDef ); /// Obtain a writable pointer to the current AlignDef

		int m_GapPenalty;
		int m_XDimension;
		int m_YDimension;
		int m_BestPathStartCellIDX;
		int m_BestPathStartCellIDY;

		TNT::Array2D<double> m_ScoreMatrix;
		TNT::Array3D<int> m_PathStoreMatrix;
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
	class PD_API AlignmentDef
	{
		friend class PD_API AlignerBase;
		friend class PD_API BestPathBase;
	public:
		AlignmentDef( const Sequence::BioSequence &_Seq1, const Sequence::BioSequence &_Seq2 ) : m_Seq1(&_Seq1), m_Seq2(&_Seq2)
		{
			m_AllocatedMemory = 0;
			m_Equiv = NULL;
			m_Aligned = false;
		}

		AlignmentDef( const AlignmentDef& _Clone ) : m_Seq1(_Clone.m_Seq1), m_Seq2(_Clone.m_Seq2)
		{	
			CloneAlignment(_Clone);
		}

#ifndef SWIG
		AlignmentDef& operator= (const AlignmentDef &_Clone)
		{
			m_Seq1 = _Clone.m_Seq1;
			m_Seq2 = _Clone.m_Seq2;
			CloneAlignment(_Clone);
			return *this;
		}
#endif

		virtual ~AlignmentDef()
		{
			delete[] m_Equiv;
		}

		void printToScreen() const;
		inline const int* getEquivelency() const { return m_Equiv; }
		inline const Sequence::BioSequence &getSeq1() const { return *m_Seq1; }
		inline const Sequence::BioSequence &getSeq2() const { return *m_Seq2; }

		bool isSeqAligned() const { return m_Aligned; } ///< Has this class been aligned?
		void ResetAlignment(); ///< Return this class to an unaligned state

	protected:
		inline int* getEquivelency() { return m_Equiv; }

		// Helper Functions
		void CloneAlignment( const AlignmentDef& _Clone );
		int nextEquivFromPoint( int lookFrom ) const;

		// Member Data
		const Sequence::BioSequence *m_Seq1;
		const Sequence::BioSequence *m_Seq2;
		size_t m_AllocatedMemory; // Record this to avoid a memory leak upon Copy Constructor if _Seq1's length changes
		int* m_Equiv;
	private:
		bool m_Aligned;
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
	class PD_API ExpSeqPair: public AlignmentDef
	{
	public:

		ExpSeqPair( 
			const Sequence::BioSequence &_BioSequence, // AlignmentDef m_Seq1 is the Sequence of the biological chain
			const Sequence::BioSequence &_StructuralSequence // AlignmentDef m_Seq2 is The Sequence of the biological chain for which atomic structure is known
			) 
			: AlignmentDef(_BioSequence,_StructuralSequence)
		{
		}

		virtual ~ExpSeqPair()
		{
		}

		/// check correct correspondence of the sequence:
		/// There should be no mis-matching residues in a biological/structural sequence pairing!
		/// Only completely missing residues are allowed, and then only in the structural sequence.
		bool isValid() const;
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
	class PD_API AlignerBase
	{
	public:
		AlignerBase()
		{
		}

		void Align( AlignmentDef &_AlignDef ) const;

	protected:
		int* getEquiv( AlignmentDef &_AlignDef ) const; /// Obtain a writable pointer to the current AlignDef
		virtual void AlignCore( AlignmentDef &_AlignDef ) const = 0; ///< Called by public Align() function - overridden in derived classes
		void PoxyAlign( AlignmentDef &_AlignDef ) const; ///< A ribbish direct alignment if the sequence is too short to allow BestPathBase to work
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
	class SimpleBestPath : public BestPathBase
	{
		friend class SimpleAligner;
	protected:
		SimpleBestPath( int gapPenalty, AlignmentDef &_AlignDef );
		void FillEquivList();	
		// Member Functions
		virtual void FillScoreMatrix();		
		AlignmentDef *m_AlignDef;
	};






//-------------------------------------------------
//
/// \brief  Generates a simple alignement between two sequences 
///
/// \details 
/// 
/*!
\code

void Test_Alignment_Code()
{
	FFParamSet ffps;
	ffps.readLib("lib/amber03aa.ff");

	// Forcefield resolution
	Sequence::BioSequence s1(ffps);
	char* seq1 = "*s-ASP-GLU-AlA-Gly-ASP-D-GLU-ALA-GLU-A*";
	printf("Resolving using 'FFParamSet':\n");
	printf("Seq1: '%s'\n",seq1);
	s1.setTo(seq1);
	printf("ResolvesTo: '");
	s1.printToScreen();
	printf("'\n\n");

	Library::NamingConventions* map2 = Library::NamingConventions::getSingleton();
	Sequence::BioSequence s2(*map2);
	char* seq2 = "GLU-AlA-G-ASP-D-AlA-GLU";
	printf("Resolving using Library::NamingConventions:\n");
	printf("Seq2: '%s'\n",seq2);
	s2.setTo(seq2);
	printf("ResolvesTo: '");
	s2.printToScreen();
	printf("'\n\n");

	Sequence::ExpSeqPair expPair(s1,s2);
	Sequence::SimpleAligner ali;
	ali.Align(expPair);
	expPair.printToScreen();

	return;
}
\endcode

*/
	///
	/// \author Jon Rea 
	///
	///
	class PD_API SimpleAligner : public AlignerBase
	{
	public:
		SimpleAligner();
	protected:
		virtual void AlignCore( AlignmentDef &_AlignDef ) const;
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
	class PD_API ManualAligner : public AlignerBase
	{
		struct ManualMatch
		{
			size_t seq1Index;
			size_t seq2Index;
			size_t length;
			void info() const
			{
				printf("Start Index S1:%d, S2:%d, Length:%d\n",seq1Index,seq2Index,length);
			}
		};
	public:
		ManualAligner();
		void addMatch( size_t _seq1Index, size_t _seq2Index );
		void addMatch( size_t _seq1StartIndex, size_t _seq2StartIndex, size_t _length );
		void clearmatches();
	protected:
		virtual void AlignCore( AlignmentDef &_AlignDef ) const;
		std::vector<ManualMatch> m_matches;
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
	class PD_API DirectAligner : public AlignerBase
	{
	public:
		DirectAligner();
	protected:
		virtual void AlignCore( AlignmentDef &_AlignDef ) const;
	};
}

#endif

