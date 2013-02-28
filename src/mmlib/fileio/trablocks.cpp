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

#include "trablocks.h"

namespace IO 
{
	// ---------------------------------------------------------
	// BTF_Block_Comment: Begin class defintion
	// ---------------------------------------------------------

	BTF_Block_Comment::BTF_Block_Comment( size_t stringCapacity ) : 
		BTF_Block(), // Initialise the base class with no params
		StringBuilder(1) // Initialise with no capacity
	{		
		m_TraCapacity = stringCapacity;
		if( m_TraCapacity % 4 != 0 )
		{
			// then ensure correct binary alignment!
			m_TraCapacity += 4 - (m_TraCapacity % 4);
		}
		growCapacity( m_TraCapacity ); // Increase capacity by the amount we want to hold
		append((char)'\0',(size_t)m_TraCapacity); // initialise to all '/0's
		clear();
	}

	BTF_Block_Comment* BTF_Block_Comment::clone() const
	{
		BTF_Block_Comment* block = new BTF_Block_Comment(m_TraCapacity);
		block->setTo(*this);
		return block; 
	}

	void BTF_Block_Comment::appendHeader( std::ofstream &_traFile )const
	{
		_traFile.write("EXTDCMNT",8);
		_traFile.write((char*)&m_TraCapacity,sizeof(int));
	}

	void BTF_Block_Comment::appendData( std::ofstream &_traFile )
	{
		setBufferAt(size(),'\0'); // force a zero termination of the internal buffer, StringBuilder does NOT do this itself internally unless toString is called
		_traFile.write(buffer(),(std::streamsize)m_TraCapacity); // The buffer has been ensured to be the correct minimum size in the constructor - it can't shrink
		clear(); // zero length again
	}

	int BTF_Block_Comment::getHeaderSize()const
	{
		return (8 * sizeof(char)) + sizeof(int); // "EXTDCMNT" followed by an integer = 12 bytes
	}

	int BTF_Block_Comment::getBlockSize()const
	{
		return (int)m_TraCapacity * sizeof(char);
	}

	size_t BTF_Block_Comment::getTraCapactity()const
	{
		return m_TraCapacity;
	}

	// ---------------------------------------------------------
	// BTF_Block_Comment: End class definition
	// ---------------------------------------------------------






	// ---------------------------------------------------------
	// BTF_Block_Vector: Begin class defintion
	// ---------------------------------------------------------

	BTF_Block_Vector::BTF_Block_Vector( int vectorCount ):
		BTF_Block(), // call the base with no params
		m_VectorCount( vectorCount ),
		m_CurrentVectorIndex(0)
	{
		m_Vectors.resize(vectorCount);
	}

	BTF_Block_Vector* BTF_Block_Vector::clone() const
	{
		BTF_Block_Vector* block = new BTF_Block_Vector(m_VectorCount);
		block->m_Vectors = m_Vectors;
		return block; 
	}

	void BTF_Block_Vector::appendHeader( std::ofstream &_traFile )const
	{
		_traFile.write("EXTDVECT",8);
		_traFile.write((char*)&m_VectorCount,sizeof(int));
	}

	void BTF_Block_Vector::ForceClear()
	{
		// blank the ones that have been alocated in this frame
		for( int i = 0; i < m_CurrentVectorIndex; i++ )
		{
			m_Vectors[i].clear();
		}
		m_CurrentVectorIndex = 0; // reset this;
	}

	void BTF_Block_Vector::appendData( std::ofstream &_traFile )
	{
		// write the current data
		_traFile.write((char*)&m_Vectors[0],sizeof(DrawingVector)*m_VectorCount);
		ForceClear();
	}

	int BTF_Block_Vector::getHeaderSize()const
	{
		return (8 * sizeof(char)) + sizeof(int); // "EXTDVECT" followed by an integer = 12 bytes
	}

	int BTF_Block_Vector::getBlockSize()const
	{
		return m_VectorCount * sizeof(DrawingVector);
	}

	int BTF_Block_Vector::getVectorCount()const
	{
		return m_VectorCount;
	}

	DrawingVector* BTF_Block_Vector::request()
	{
		if( m_CurrentVectorIndex < m_VectorCount )
		{
			return &m_Vectors[m_CurrentVectorIndex++];
		}
		else
		{
			return NULL;
		}
	}

	// ---------------------------------------------------------
	// BTF_Block_Vector: End class defintion
	// ---------------------------------------------------------
}

