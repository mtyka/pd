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

#include <fstream>
#include "tools/io.h"






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
class PD_API AlignFile
{
public:
	AlignFile( AlignmentDef &_AlignDef ) :
		m_AlignDef(_AlignDef),
		m_Equiv(_AlignDef.GetEquivelency()),
		m_EquivLength(_AlignDef.GetSeq1().GetResidueCount())
	{
	}

	void AlignFile::Save( const std::string &_FileName )
	{
		ofstream _filestream;

		try
		{
			_filestream.open(filename.c_str(), ios::out );
			if( !_filestream.is_open() )
			{
				THROW(IOException,"'AlignFile::Save()' could not open the file for writing!");
			}

			// Pickle a header and the required information for alignment recovery:
			_filestream << "MANUAL_SEQUENCE" << std::endl;
			_filestream << IO::Pickle<int>(m_AlignDef.GetSeq1().GetResidueCount());
			_filestream << IO::Pickle<int>(m_AlignDef.GetSeq2().GetResidueCount());
			_filestream << m_AlignDef.GetSeq1() << std::endl;
			_filestream << m_AlignDef.GetSeq2() << std::endl;
			_filestream << IO::PickleArray<int>(m_Equiv,m_EquivLength);
		}
		catch(ExceptionBase ex)
		{
			if( _filestream.is_open() ) _filestream.close(); // release the file handle!
			throw; // rethrow our exception - we just wanted to catch it to close the fileStream...
		}

		if( _filestream.is_open() ) _filestream.close(); // release the file handle!
	}

	void AlignFile::Load( const std::string &_FileName )
	{
		THROW(NotImplementedException,"");
		//ifstream _filestream;
	}

	void SetEquivelency( size_t _Seq1Index, size_t _Seq2Index )
	{
		m_Equiv[_Seq1Index] = _Seq2Index;
	}

	void ResetAt( size_t _Seq1Index )
	{
		m_Equiv[_Seq1Index] = -1;
	}

	void ResetAll()
	{
		for( size_t i = 0; i < m_EquivLength; i++ )
		{
			m_Equiv[i] = -1;
		}
	}

protected:
	// Member Data
	AlignmentDef &m_AlignDef;
	int* m_Equiv; // points to memory in m_AlignDef
	size_t m_EquivLength;
};

