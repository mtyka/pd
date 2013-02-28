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

#ifndef __SEQUENCE_H
#define __SEQUENCE_H

// Essential system headers
#include <string>
#include <vector>
#include <iostream>

// Essential Headers
#include "library/residues.h"
#include "library/mapper.h"

// Forward Declarations
namespace Library
{
	class PD_API AliasMapper;
	class PD_API ClassMapper;
}

class PD_API FFParamSet;

namespace Sequence
{
	//-------------------------------------------------
	//
	/// \brief  Handles Biological Residue Sequences
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	///
	/// \author  Jon Rea 
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class BioSource
	{
	public:
		BioSource(
			int _BiologicalIndex,
			char _BiologicalICode,
			const std::string &_ResidueName
			);

		bool isEqual( const BioSource& residue ) const;
		void setTo( const BioSource& source );

		int getBiologicalIndex() const { return m_BiologicalIndex; }
		char getBiologicalICode() const { return m_BiologicalICode; }

		///  Returns the name of this residue, disragarding its chain placement prefix
		std::string getResName() const { return m_ResidueName; } 

	protected:
		int m_BiologicalIndex;

		/// Insertion code
		char m_BiologicalICode; 
		std::string m_ResidueName;
	};



	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	///
	/// \author  Jon Rea 
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class ResidueInfo : public BioSource
	{
		friend class BioSequence;

	public:
		ResidueInfo(
			char _SingleLetter,
			unsigned int _SequentialIndex,
			const std::string &_ResidueName );

		ResidueInfo(
			char _SingleLetter,
			unsigned int _SequentialIndex,
			int _BiologicalIndex,
			char _BiologicalICode,
			const std::string &_ResidueName
			);

		void setTo( const ResidueInfo& _Clone );

		/// // For 'N' and 'C' terminal residues
		void setPrefix( Library::PrePostFix _Prefix ) { m_Prefix = _Prefix; } 
		Library::PrePostFix getPrefix() const { return m_Prefix; }
		void ForceRename(const std::string &_ResidueName) { m_ResidueName = _ResidueName; }

		bool isEqual( const ResidueInfo& residue, bool useBioInfo ) const;

		// public accessors
		char getSingleResLetter() const { return m_SingleLetter; }
		Library::StandardResidues getResidueType() const { return m_ResidueType; }
		unsigned int getSequentialIndex() const { return m_SequentialIndex; }


		/// Returns the name of this residue, taking into account its placement within the chain
		std::string getFullName() const { return m_Prefix.build(m_ResidueName); } 

	protected:
		char m_SingleLetter;
		Library::StandardResidues m_ResidueType;
		unsigned int m_SequentialIndex;

		/// Determins this residues placement within the chain e.g. An alanine 
		/// residue at the N-terminus is termed NALA, not ALA, with the prefix 'N'
		Library::PrePostFix m_Prefix; 
	};


	//-------------------------------------------------
	//
	/// \brief Defines & Holds an arbitrary sequence of a polymer 
	///
	/// \details 
	/// Syntax of the parse Function
	///
	/// Strings are termed 'advanced sequence strings'
	///
	/// simple: NHIS-ALA-GLN-CYS-HID-TYR-NME  with residues deliminated by '-'s
	/// simplified: *H-(ANC)-HID-(Y)-NME      with single letter codes in brackets
	/// translation from A to ALA is via ALIAS defs in ff def file
	///
	/// \author  Jon Rea 
	///
	class PD_API BioSequence
	{
	public:
		// Constructor Logic

		/// This should be for normal use. The other constructors use the Library::NamingConventions class singleton
		BioSequence( const FFParamSet &_ParserMapper ); 
		BioSequence(); /// < Use the default singleton AliasMapper. This has no forcefield knowledge and will fail for all but the most standard residues.
		BioSequence( const Library::AliasMapper &_ParserNameMapper );
		BioSequence( const Library::ClassMapper &_ParserClassMapper );
		BioSequence( const Library::AliasMapper &_ParserNameMapper, const Library::ClassMapper &_ParserClassMapper );

		// Operators
#ifndef SWIG
		// Function defined below class bods allows use of 'BioSequence' with cout.
		// 'friend' required as this function needs to access private class member data
		friend std::ostream& operator<<(std::ostream &s, const BioSequence &_Print);

		BioSequence& operator= ( const BioSequence& _Seq );

		/// Yes this should be a 'char' - we want to be able to use the class as a string. 
		/// To get a ResidueInfo use: const ResidueInfo& getResidue(size_t _Index) const;
		char operator[](size_t index) const;
#endif
		bool compare( const BioSequence& _Seq ) const;

		// setup Functions
		void setNamingException( bool _enable ) { m_NamingExceptions = _enable; }
		void setPrefix(Library::PrePostFix _First, Library::PrePostFix _Last);
		void setPrefix( size_t index, Library::PrePostFix _Fix );
		void removeAllPrefixes();

		void setChainID( char chainID ) { m_ChainID = chainID; }
		char getChainID() const { return m_ChainID; }

		void setTo(const std::string &_FormattedSequence);
		void setTo( const BioSequence &_FullSequence );

		void append( const std::string &_FormattedSequence );
		void append( const BioSequence &_FullSequence );
		void append( const BioSequence &_Seq, size_t _Start, size_t _Length );
		void append( const ResidueInfo &_Res );

		/// Re-number the sequential index of all residues in this collection, beginning from 0.
		void reSequence(); 

		// Still to implement mirrors of all the above appends ...
		// void insert(  )

		void erase( size_t _Index );
		void erase( size_t _Index, size_t _Length );

		void clear();

		// Output Functions
		void printToScreen( char _Delimiter = '-', int _ResWidth = -1, int _MaxWidth = 70 ) const;
		void printToBuffer( char *buffer, char _Delimiter = '-' ) const;
		std::string printToString( char _Delimiter = '-' ) const;
		std::string printToStringSingle() const;
		std::vector<std::string> makeResidueStrings() const;

		// Acceessors
		BioSequence makeSubSequence( size_t _StartIndex, size_t _Length ) const;
		inline size_t size() const { return m_ResList.size(); }
		std::string getFullName(size_t _Index) const;
		std::string getResName(size_t _Index) const;
		char getSingleResLetter(size_t _Index) const;
		const ResidueInfo& getResidue(size_t _Index) const;

	protected:
		// Helper Functions
		int interpretMultiRes( const std::string &_FormattedSequence, int& pscan ) const;
		void addParsedRes( const StringBuilder& _name );

		// Member Data
		std::vector<ResidueInfo> m_ResList;
		const Library::AliasMapper *m_ParserNameMapper;
		const Library::ClassMapper *m_ParserClassMapper;
		bool m_NamingExceptions;
		char m_ChainID;
	};

	// Stream Operator Overloading
	std::ostream& operator<<(std::ostream &s, const BioSequence &_Print);
	bool operator==(const BioSequence& _Seq1,const BioSequence& _Seq2);



	//-------------------------------------------------
	//
	/// \brief  A container of sequences 
	///
	/// \details 
	///    
	///
	/// \author  Jon Rea 
	///
	template< class T >
	class BioSequenceCollection
	{
	public:
		BioSequenceCollection();

		Sequence::BioSequence &getSequence( size_t _index );
		const Sequence::BioSequence &getSequence( size_t _index ) const;

		Sequence::BioSequence &getSequence( const T &_indexer );
		const Sequence::BioSequence &getSequence( const T &_indexer ) const;

		Sequence::BioSequence &getLastSequence();

		bool HasIndex( const T &_indexer ) const;
		T &getIndexer( size_t _index );

		size_t size()const { return m_Sequence.size(); }

		void printScreen() const;

	protected:
		// Helper Functions
		void addSequence( T &_indexer, BioSequence &_s );
		void addSequence( T &_indexer ); /// add at '_indexer' with an empty sequence
		void clear();

	private:
		// Member Data
		std::vector<T> m_Indexer;
		std::vector<Sequence::BioSequence> m_Sequence;
	};

	// ---------------------------------------------------------------------
	// Explicit Templated function declarations in header file are required.
	// See 'Export' keyword in the following URL:
	// http://www.codeguru.com/Cpp/COM-Tech/atl/tutorials/article.php/c3617/
	// ---------------------------------------------------------------------

	template < class T >
	BioSequenceCollection<T>::BioSequenceCollection()
	{
	}

	template < class T >
	void BioSequenceCollection<T>::addSequence( T &_indexer, BioSequence &_s )
	{
		m_Indexer.push_back(_indexer);
		m_Sequence.push_back(_s);
	}

	template < class T >
	void BioSequenceCollection<T>::addSequence( T &_indexer )
	{
		m_Indexer.push_back(_indexer);
		m_Sequence.push_back(BioSequence());
	}

	template < class T >
	void BioSequenceCollection<T>::clear()
	{
		m_Indexer.clear();
		m_Sequence.clear();
	}

	template < class T >
	BioSequence &BioSequenceCollection<T>::getSequence( const T &_indexer )
	{
		for( size_t i = 0; i < m_Sequence.size(); i++ )
		{
			if( _indexer == m_Indexer[i] ) return m_Sequence[i];
		}
		THROW(ArgumentException,"Index not found for array!");
	}

	template < class T >
	const BioSequence &BioSequenceCollection<T>::getSequence( const T &_indexer ) const
	{
		for( size_t i = 0; i < m_Sequence.size(); i++ )
		{
			if( _indexer == m_Indexer[i] ) return m_Sequence[i];
		}
		THROW(ArgumentException,"Index not found for array!");
	}

	template < class T >
	BioSequence &BioSequenceCollection<T>::getSequence( size_t _index )
	{
		if( _index >= m_Sequence.size() )
			THROW(OutOfRangeException,"Index is outside array bounds");
		return m_Sequence[_index];
	}

	template < class T >
	const BioSequence &BioSequenceCollection<T>::getSequence( size_t _index ) const
	{
		if( _index >= m_Sequence.size() )
			THROW(OutOfRangeException,"Index is outside array bounds");
		return m_Sequence[_index];
	}

	template < class T >
	BioSequence &BioSequenceCollection<T>::getLastSequence()
	{
		if( m_Sequence.size() == 0 )
			THROW(OutOfRangeException,"Last element cannot be returned, array length is '0'");
		return m_Sequence[m_Sequence.size()-1];
	}

	template < class T >
	T &BioSequenceCollection<T>::getIndexer( size_t _index )
	{
		if( _index >= m_Sequence.size() )
			THROW(OutOfRangeException,"Index is outside array bounds");
		return m_Indexer[_index];
	}

	template < class T >
	bool BioSequenceCollection<T>::HasIndex( const T &_indexer ) const
	{
		for( size_t i = 0; i < m_Indexer.size(); i++ )
		{
			if( _indexer == m_Indexer[i] ) return true;
		}
		return false;
	}

	// 'inline' is apparently required ... oddness!
	// Quote: "Function Templates are an exempt of ODR (one definition rule) and may be more than one
	// definition of them in different TU's (Type Unknowns). Full function template
	// specialization is not a template, rather an ordinary function, so you
	// need to use inline keyword not to violate ODR if you want to put them
	// in a header file included into several transaction units."
	template<>
	inline void BioSequenceCollection<char>::printScreen() const
	{
		for( size_t i = 0; i < m_Sequence.size(); i++ )
		{
			std::cout << "ChainID: '" << m_Indexer[i] << '\'' << std::endl;
			std::cout << m_Sequence[i] << std::endl << std::endl;
		}
	}

	/// A template specialisation for the <char> based sequence collection
	/// Also see comment for 'printScreen()' above.
	template<>
	inline void BioSequenceCollection<char>::addSequence( char &_indexer )
	{
		m_Indexer.push_back(_indexer);
		m_Sequence.push_back(BioSequence());
		m_Sequence[m_Sequence.size()-1].setChainID(_indexer); // assign the chainID to the sequence itself
	}

	template < class T >
	void BioSequenceCollection<T>::printScreen()const
	{
		for( size_t i = 0; i < m_Sequence.size(); i++ )
		{
			std::cout << m_Indexer[i] << std::endl;
			std::cout << m_Sequence[i] << std::endl << std::endl;
		}
	}
}

#endif

