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

#include "sequence.h"

#include "tools/stringbuilder.h"
#include "library/nameconv.h"
#include "library/residues.h"
#include "forcefields/ffparam.h"

using namespace std;
using namespace Library;

namespace Sequence
{
	// --------------------------------------------------------------------------------
	// Handles Protein Sequences
	// --------------------------------------------------------------------------------

	BioSource::BioSource(
		int _BiologicalIndex,
		char _BiologicalICode,
		const std::string &_ResidueName
		)
	{
		m_BiologicalIndex = _BiologicalIndex;
		m_BiologicalICode = _BiologicalICode;
		m_ResidueName = _ResidueName;
	}

	bool BioSource::isEqual( const BioSource& residue ) const
	{
		if( m_BiologicalIndex != residue.m_BiologicalIndex ) return false;
		if( m_BiologicalICode != residue.m_BiologicalICode ) return false;
		if( m_ResidueName.length() == 1 && m_ResidueName[0] == '*' ) return true;
		if( 0 != m_ResidueName.compare(residue.m_ResidueName) ) return false;
		return true;
	}

	void BioSource::setTo( const BioSource& source )
	{
		m_BiologicalIndex = source.m_BiologicalIndex;
		m_BiologicalICode = source.m_BiologicalICode;
		m_ResidueName = source.m_ResidueName;
	}



	ResidueInfo::ResidueInfo( char _SingleLetter, unsigned int _SequentialIndex,
		int _BiologicalIndex, char _BiologicalICode, const std::string &_ResidueName )
		: BioSource( _BiologicalIndex, _BiologicalICode, _ResidueName )
	{
		m_SingleLetter = _SingleLetter;
		m_ResidueType = Library::getResidue(_SingleLetter);
		m_SequentialIndex = _SequentialIndex;
		m_Prefix = Library::PrePostFix::NoFix;
	}

	ResidueInfo::ResidueInfo( char _SingleLetter, unsigned int _SequentialIndex, const std::string &_ResidueName )
		: 	BioSource( _SequentialIndex, ' ', _ResidueName )
	{
		m_SingleLetter = _SingleLetter;
		m_ResidueType = Library::getResidue(_SingleLetter);
		m_SequentialIndex = _SequentialIndex;
		m_Prefix = Library::PrePostFix::NoFix;
	}

	void ResidueInfo::setTo( const ResidueInfo& _Copy )
	{
		m_SingleLetter = _Copy.m_SingleLetter;
		m_ResidueType = _Copy.m_ResidueType;
		m_SequentialIndex = _Copy.m_SequentialIndex;
		m_BiologicalIndex = _Copy.m_BiologicalIndex;
		m_BiologicalICode = _Copy.m_BiologicalICode;
		m_ResidueName = _Copy.m_ResidueName;
		m_Prefix = _Copy.m_Prefix;
	}

	bool ResidueInfo::isEqual( const ResidueInfo& residue, bool useBioInfo ) const
	{
		if( useBioInfo )
		{
			if( !BioSource::isEqual(residue) ) return false;
		}

		if( m_SingleLetter != residue.m_SingleLetter ) return false;
		if( m_ResidueType != residue.m_ResidueType ) return false;
		if( m_SequentialIndex != residue.m_SequentialIndex ) return false;

		return true;
	}


	BioSequence::BioSequence() 
		: m_NamingExceptions(true), m_ChainID(' ')
	{
		m_ParserNameMapper = Library::NamingConventions::getSingleton();
		m_ParserClassMapper = Library::NamingConventions::getSingleton();
	}

	BioSequence::BioSequence( const AliasMapper &_ParserNameMapper ) 
		: m_NamingExceptions(true), m_ChainID(' ')
	{
		m_ParserNameMapper = &_ParserNameMapper;
		m_ParserClassMapper = Library::NamingConventions::getSingleton();
	}

	BioSequence::BioSequence( const ClassMapper &_ParserClassMapper ) 
		: m_NamingExceptions(true), m_ChainID(' ')
	{
		m_ParserNameMapper = Library::NamingConventions::getSingleton();
		m_ParserClassMapper = &_ParserClassMapper;
	}

	BioSequence::BioSequence( const AliasMapper &_ParserNameMapper, const ClassMapper &_ParserClassMapper ) 
		: m_NamingExceptions(true), m_ChainID(' ')
	{
		m_ParserNameMapper = &_ParserNameMapper;
		m_ParserClassMapper = &_ParserClassMapper;
	}

	BioSequence::BioSequence( const FFParamSet &_ParserMapper ) 
		: m_NamingExceptions(true), m_ChainID(' ')
	{
		m_ParserNameMapper = &_ParserMapper;
		m_ParserClassMapper = &_ParserMapper;
	}

	BioSequence& BioSequence::operator= ( const BioSequence& _Seq )
	{
		m_NamingExceptions = _Seq.m_NamingExceptions;
		m_ChainID = _Seq.m_ChainID;
		m_ParserNameMapper = _Seq.m_ParserNameMapper;
		m_ParserClassMapper = _Seq.m_ParserClassMapper;
		clear();
		append(_Seq);
		return *this;
	}

	bool BioSequence::compare( const BioSequence& _Seq ) const
	{
		size_t count = size();
		if( count != _Seq.size() ) return false;
		for( size_t i = 0; i < count; i++ )
		{
			if( 0 != m_ResList[i].m_ResidueName.compare(_Seq.m_ResList[i].m_ResidueName) )
			{
				return false;
			}
		}
		return true;
	}

	char BioSequence::operator[](size_t index)const
	{
		if( index >= m_ResList.size() ) THROW(OutOfRangeException,"Index was outside array bounds!");
		return m_ResList[index].getSingleResLetter();
	}

	void BioSequence::addParsedRes( const StringBuilder& _name )
	{
		ASSERT(m_ParserNameMapper!=NULL, CodeException, "Internal member 'm_ParserNameMapper' should not be null, ever!");
		ASSERT( _name.size() != 0, ParseException, "Token length is 0!");

		bool polymerStart = _name[0] == '*';
		if( polymerStart && _name.size() == 1 ) throw ParseException("Residue name '*' is invalid");
		bool polymerEnd = false;
		if( _name.size() > 1 ) polymerEnd = _name[_name.size()-1] == '*';
		if( polymerStart && polymerEnd && _name.size() == 2 ) throw ParseException("Residue name '**' is invalid");
		if( polymerStart && polymerEnd ) throw ParseException("Residue name '*xxxx*' is invalid. A single residue cannot be the polymer start and end by definition.");
		if( polymerStart && m_ResList.size() != 0 ) throw ParseException("Polymer start residue indicator '*' specified for a non-terminal residue!");

		std::string lookupName = _name.toString();

		// We now just want the name itself with no operator
		lookupName = trim(lookupName,"*");

		// Perform an alias lookup first, it may return false, but this doesn't matter...		
		m_ParserNameMapper->lookupAlias(lookupName,lookupName);

		if ( lookupName.size() == 1 )
		{
			std::string longName("???");
			char singleID = lookupName[0];
			if(!m_ParserNameMapper->lookupLongName(singleID,longName) && m_NamingExceptions )
			{
				throw(ParseException("Mappers lookupLongName() failed for name '" + std::string(1,singleID) + "'!"));
			}
			m_ResList.push_back( ResidueInfo( singleID, m_ResList.size(),longName) );
		}
		else
		{
			if( !polymerStart && !polymerEnd )
			{
				// It could unfortunatly be something like NARG, which must be split 
				// into PrePostFix("N",false) and ARG
				m_ParserClassMapper->detectSplitRequirement( lookupName, polymerStart, polymerEnd );
			}

			char singleID = '?';
			if(!m_ParserNameMapper->lookupShortName(lookupName,singleID) && m_NamingExceptions )
			{
				throw(ParseException("Mappers lookupShortName() failed for name '" + lookupName + "'!"));
			}
			m_ResList.push_back(ResidueInfo( singleID, m_ResList.size(), lookupName ));
		}

		if( polymerStart )
		{
			ResidueInfo& r = m_ResList[m_ResList.size()-1]; // the one we just added
			Library::PrePostFix preFix;
			if(!m_ParserClassMapper->getStartTerminusFromResidue(r.getResName(),preFix) )
			{
				throw(ParseException("Internal class mapper is unable to determine the polymer start prefix. Did you include the definitions \"CLASS_PREFIX\" in the forcefield definition files you have read in ? (for example for peptides: CLASS_PREFIX polypeptide N C )"));
			}	
			r.setPrefix( preFix );
		}
		if( polymerEnd )
		{
			for( size_t i = 1; i < m_ResList.size(); i++ )
			{
				if( m_ResList[i].getPrefix().isDefined() != 0 )
				{
					throw(ProcedureException("Biosequence parsing has yielded multiple terminii definitions!"));
				}
			}
			ResidueInfo& r = m_ResList[m_ResList.size()-1]; // the one we just added
			Library::PrePostFix postFix;
			if(!m_ParserClassMapper->getEndTerminusFromResidue(r.getResName(),postFix) )
			{
				throw(ParseException("Internal class mapper is unable to determine the polymer end prefix"));
			}
			r.setPrefix( postFix );
		}
	}

	void BioSequence::setPrefix(Library::PrePostFix _First, Library::PrePostFix _Last)
	{
		removeAllPrefixes();
		if( m_ResList.size() < 2 ) return;
		m_ResList[0].setPrefix(_First);
		m_ResList[m_ResList.size()-1].setPrefix(_Last);
	}

	void BioSequence::setPrefix( size_t index, Library::PrePostFix _Fix )
	{
		m_ResList[index].setPrefix(_Fix);
	}

	void BioSequence::removeAllPrefixes()
	{
		for( size_t i = 0; i < m_ResList.size(); i++ )
		{
			m_ResList[i].setPrefix(Library::PrePostFix::NoFix);
		}
	}

	void BioSequence::setTo(const std::string &_FormattedSequence)
	{
		clear();
		append(_FormattedSequence.c_str());
	}

	void BioSequence::setTo( const BioSequence &_FullSequence )
	{
		clear();
		append(_FullSequence);
	}

	int BioSequence::interpretMultiRes( const std::string &_FormattedSequence, int& pscan ) const
	{
		StringBuilder sbNumber(4);
		if(_FormattedSequence[pscan] != '[') // we should be sitting on one of these
			THROW(CodeException,"Code path failure whilst interpreting the '[]' operator\n");
		while(true)
		{
			pscan++;
			if(_FormattedSequence[pscan] == '\0') // end of string
				THROW(ParseException,"SYNTAX ERROR: Unexpected end of string whilst interpreting the '[]' operator\n");
			if(_FormattedSequence[pscan] == ']') // end of string
			{
				if( sbNumber.size() == 0 )
				{
					THROW(ParseException,"SYNTAX ERROR: No number within the '[]' operator\n");
				}
				int count = sbNumber.parseInt();
				if( count <= 0 )
				{
					THROW(ParseException,"SYNTAX ERROR: Non-sensical count within the '[]' operator\n");
				}
				return count;
			}
			if( !isdigit( _FormattedSequence[pscan] ) )
				THROW(ParseException,"SYNTAX ERROR: Non-numeric character within the '[]' operator\n");
			sbNumber.append( _FormattedSequence[pscan] );
		}
	}

	void BioSequence::append(const std::string &_FormattedSequence)
	{
		if(_FormattedSequence.size()==0 )
			return;

		int pscan = -1; // asequence scanning position
		bool bracketon = false; // are we inside a bracket ?
		StringBuilder ident(8);
		BioSequence innerSequenceGroup(*m_ParserNameMapper,*m_ParserClassMapper); // provision for <> operator groups

		// Now take apart the sequence string...
		while(true)
		{
			pscan++;
			if(_FormattedSequence[pscan] == '\0') // end of string
			{
				if(!bracketon) // when outside the bracket
				{
					if(ident.size() != 0) // in case of an unfinished token
					{
						addParsedRes(ident);
					}
					return; // Finished!
				}
				else
				{
					THROW(ParseException,"SYNTAX ERROR: closing ')' expected\n");
				}
			}
			else if(_FormattedSequence[pscan] == '>') // the MultiRes operator :-D
			{
				THROW(ParseException,"SYNTAX ERROR: closing '>' expected\n");
			}
			else if(_FormattedSequence[pscan] == '<') // the MultiRes operator :-D
			{
				if( ident.size() != 0 )
					THROW(ParseException,"SYNTAX ERROR: missing '-' before operator '<'\n");
				StringBuilder sb;
				sb.setTo( _FormattedSequence, pscan+1, _FormattedSequence.size() - pscan - 1 );

				// Identify the ending '>' for the current '<' taking nesting into account
				size_t seekID = 0;
				size_t endRange = SIZE_T_FAIL;
				while(true)
				{
					size_t nextStart = sb.XthOf(seekID,'<');
					endRange = sb.XthOf(seekID,'>');
					if( endRange == SIZE_T_FAIL )
					{
						THROW(ParseException,"SYNTAX ERROR: Extected to find a closing '>' operator\n");
					}
					if( nextStart == SIZE_T_FAIL || nextStart > endRange )
					{
						break;
					}
					seekID++;
				}
				D_ASSERT(endRange!=SIZE_T_FAIL, CodeException, "Assumption failure");

				sb.TruncateRightTo(endRange);
				if( bracketon )
				{
					sb.insert(0,'(');
					sb.append(')');
				}
				innerSequenceGroup.setTo(sb.toString());
				pscan += endRange + 1;
				int count = 1;
				if( _FormattedSequence[pscan+1] == '[' )
				{
					pscan++;
					count = interpretMultiRes( _FormattedSequence, pscan );
				}
				for( int i = 0; i < count; i++ )
				{
					append(innerSequenceGroup);			
				}
			}
			else if(_FormattedSequence[pscan] == '[') // the MultiRes operator :-D
			{
				if( ident.size() == 0 )
				{
					THROW(ParseException,"SYNTAX ERROR: Expected identifier before '[]' operator\n");
				}
				int count = interpretMultiRes( _FormattedSequence, pscan );
				for( int i = 0; i < count; i++ )
				{
					addParsedRes(ident);					
				}
				ident.clear();
			}
			else if(bracketon) 
			{ 
				if(_FormattedSequence[pscan] == '(') 
					THROW(ParseException,"SYNTAX ERROR: opening '(' unexpected\n");
				if(_FormattedSequence[pscan] == '*') 
					THROW(ParseException,"SYNTAX ERROR: Terminii operator '*' unexpected in bracketed single-char section\n");
				// when inside the bracket
				if(_FormattedSequence[pscan] == ')') 
				{
					bracketon = false;
					continue;
				}
				ident.setTo(_FormattedSequence[pscan]); // else take character and make it into a string
				if( !isalnum(ident[0]) && (ident[0]!='_') )
					THROW( ParseException, "SYNTAX ERROR: illegal character\n");
				int count = 1;
				if( _FormattedSequence[pscan+1] == '[' )
				{
					pscan++;
					count = interpretMultiRes( _FormattedSequence, pscan );
				}
				for( int i = 0; i < count; i++ )
				{
					addParsedRes(ident);					
				}			
				ident.clear(); // reset
				continue;
			}
			else
			{
				// when outside the bracket
				if(_FormattedSequence[pscan] == ')') 
					THROW(ParseException,"SYNTAX ERROR: closing ')' unexpected\n");
				else if(_FormattedSequence[pscan] == '-')
				{ 
					// found delimiter
					if(ident.size() > 0)
					{
						addParsedRes(ident);
						ident.clear();
					}
					continue;
				}
				else if((ident.size() == 0) && // if at beginning of token
					(_FormattedSequence[pscan] == '('))
				{
					bracketon = true; // inside a bracket now different rules :-)
					continue;
				}

				// otherwise accumulate letter into identifier
				char c = toupper(_FormattedSequence[pscan]);
				if( !isalnum(c) && (c!='*') && (c!='_') )
					THROW( ParseException, "SYNTAX ERROR: illegal character\n");
				ident.append(c);
			}
		}
	}

	void BioSequence::append( const ResidueInfo &_Res )
	{
		m_ResList.push_back(_Res);
	}

	void BioSequence::append( const BioSequence &_Seq )
	{
		for( size_t i = 0; i < _Seq.size(); i++ )
		{
			m_ResList.push_back(_Seq.getResidue(i));
		}
	}

	void BioSequence::reSequence()
	{
		for( size_t i = 0; i < m_ResList.size(); i++ )
		{
			m_ResList[i].m_SequentialIndex = i;
		}
	}

	void BioSequence::append( const BioSequence &_Seq, size_t _Start, size_t _Length )
	{
		size_t toLength = _Start + _Length;
		for( size_t i = _Start; i < toLength; i++ )
		{
			m_ResList.push_back(_Seq.getResidue(i));
		}
	}

	void BioSequence::erase( size_t _Index )
	{
		m_ResList.erase( m_ResList.begin()+_Index );
	}

	void BioSequence::erase( size_t _Index, size_t _Length )
	{
		m_ResList.erase( m_ResList.begin()+_Index, m_ResList.begin()+_Index+_Length );
	}

	void BioSequence::clear()
	{
		m_ResList.clear();
	}

	BioSequence BioSequence::makeSubSequence( size_t _StartIndex, size_t _Length )const
	{
		ASSERT( _StartIndex + _Length <= m_ResList.size(), OutOfRangeException, "makeSubSequence is out of range" );
		BioSequence bio;
		for( size_t i = 0; i < _Length; i++ )
		{
			bio.append(m_ResList[_StartIndex+i]);
		}
		return bio;
	}

	void BioSequence::printToScreen( char _Delimiter, int _ResWidth, int _MaxWidth )const
	{
		if( m_ResList.size() == 0 )
		{
			return;
		}

		int printedChars = 0;
		int alloc = 0;

		size_t length = m_ResList.size();
		for( size_t i = 0; i < length; i++ )
		{
			std::string name =  m_ResList[i].getFullName();
			if( _ResWidth > 0 && name.size() < _ResWidth )
			{
				if( name.size() == 3 ) name = ' ' + name;
				name.resize(_ResWidth,' ');
			}
			alloc = name.size() + 1;
			if( alloc + printedChars >= _MaxWidth )
			{
				cout << endl;
				printedChars = 0;
			}
			printedChars += alloc;
			cout << name;
			if( i != length-1 )
			{
				cout << _Delimiter;
			}
		}
	}

	void BioSequence::printToBuffer( char *buffer, char _Delimiter )const
	{
		if( m_ResList.size() == 0 )
		{
			buffer[0] = 0;
			return;
		}
		size_t length = m_ResList.size() - 1;
		size_t i, strLen, buffPos = 0;
		for( i = 0; i < length; i++ )
		{
			const std::string& resName = m_ResList[i].getFullName();
			strLen = resName.length();
			memcpy(&buffer[buffPos],resName.c_str(),strLen);
			buffPos+=length;
			buffer[buffPos] = _Delimiter;
			buffPos++;
		}
		const std::string& resName = m_ResList[i].getFullName();
		strLen = resName.length();
		memcpy(&buffer[buffPos],resName.c_str(),strLen);
		buffPos+=length;
		buffer[buffPos] = 0;
	}

	std::string BioSequence::printToString( char _Delimiter )const
	{
		if( m_ResList.size() == 0 )
		{
			return string("");
		}
		size_t length = m_ResList.size() - 1;
		StringBuilder seqBuild( 4 * m_ResList.size() + 10 );
		size_t i;
		for( i = 0; i < length; i++ )
		{
			seqBuild.append(m_ResList[i].getFullName());
			seqBuild.append(_Delimiter);
		}
		seqBuild.append(m_ResList[i].getFullName());
		return seqBuild.toString();
	}

	std::string BioSequence::printToStringSingle()const
	{
		if( m_ResList.size() == 0 )
		{
			return string("");
		}
		size_t length = m_ResList.size() - 1;
		StringBuilder seqBuild( 4 * m_ResList.size() + 10 );
		size_t i;
		for( i = 0; i < length; i++ )
		{
			seqBuild.append(m_ResList[i].getSingleResLetter());
		}
		seqBuild.append(m_ResList[i].getSingleResLetter());
		return seqBuild.toString();
	}

	std::vector<std::string> BioSequence::makeResidueStrings()const
	{
		std::vector<std::string> resStrings;
		size_t length = m_ResList.size();
		for(size_t i = 0; i < length; i++ )
		{
			resStrings.push_back(m_ResList[i].getFullName());
		}
		return resStrings;
	}

	std::string BioSequence::getFullName(size_t _Index)const
	{
		if( _Index >= m_ResList.size() ) THROW(OutOfRangeException,"Residue index is not within the bounds of the sequence!");
		return m_ResList[_Index].getFullName();
	}

	std::string BioSequence::getResName(size_t _Index)const
	{
		if( _Index >= m_ResList.size() ) THROW(OutOfRangeException,"Residue index is not within the bounds of the sequence!");
		return m_ResList[_Index].getResName();
	}

	char BioSequence::getSingleResLetter(size_t _Index)const
	{
		if( _Index >= m_ResList.size() ) THROW(OutOfRangeException,"Residue index is not within the bounds of the sequence!");
		return m_ResList[_Index].getSingleResLetter();
	}

	const ResidueInfo& BioSequence::getResidue(size_t _Index)const
	{
		if( _Index >= m_ResList.size() ) THROW(OutOfRangeException,"Residue index is not within the bounds of the sequence!");
		return m_ResList[_Index];
	}

	std::ostream& operator<<(std::ostream &s, const BioSequence &_Print)
	{
		size_t length = _Print.m_ResList.size() - 1;
		size_t i;
		for( i = 0; i < length; i++ )
		{
			s << _Print.m_ResList[i].getFullName() << '-';
		}
		s << _Print.m_ResList[i].getFullName();
		return s;
	}

	bool operator==(const BioSequence& _Seq1,const BioSequence& _Seq2)
	{
		return _Seq1.compare(_Seq2);
	}
}


