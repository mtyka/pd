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

#ifndef __BTF_BLOCKS_H
#define __BTF_BLOCKS_H

#include <fstream>
#include <string>
#include <vector>

#include "primitives.h"
#include "interface.h"
#include "tools/stringbuilder.h"

namespace IO 
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
	class PD_API BTF_Block // the base abstract data class PD_API for all tra-blocks
	{
		friend class PD_API OutTra_BTF;

	public:
		BTF_Block(){}
		virtual ~BTF_Block() {}
		virtual BTF_Block* clone() const = 0;

	protected:
		virtual void appendHeader( std::ofstream &_traFile ) const = 0; // pure virtual function uses = 0 ...
		virtual void appendData( std::ofstream &_traFile ) = 0;
		virtual int getHeaderSize() const = 0;
		virtual int getBlockSize() const = 0;
	};


	/// \brief  Allows lines to be drawn in a TrajectoryFormat file
	/// \details Allows lines to be drawn in a TrajectoryFormat file. Used by OutTra_BTF. 
	/// \author Jon Rea 
	class PD_API BTF_Block_Vector: public BTF_Block, public IDrawProvider
	{
	public:
		BTF_Block_Vector( int vectorCount = 20 );
		virtual ~BTF_Block_Vector(){}
		virtual BTF_Block_Vector* clone() const;

		int getVectorCount() const;
		virtual DrawingVector* request(); ///< returns NULL if we run out of vectors
		void ForceClear();

	protected:
		virtual void appendHeader( std::ofstream &_traFile ) const;
		virtual void appendData( std::ofstream &_traFile );
		virtual int getHeaderSize() const;
		virtual int getBlockSize() const;

		int m_VectorCount;
		int m_CurrentVectorIndex;
		std::vector<DrawingVector> m_Vectors;
	};


	//-------------------------------------------------
	//
	/// \brief Effectively a string holder, which can be told to append to the trajectory ...
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	///
	/// \author Jon Rea 
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class PD_API BTF_Block_Comment: public StringBuilder, public BTF_Block
	{
	public:
		BTF_Block_Comment( size_t traCapacity = 128 );
		virtual ~BTF_Block_Comment() {}
		virtual BTF_Block_Comment* clone() const;
		size_t getTraCapactity() const; ///< The underlying StringBuilder can be any size, but how much can we write to the TrajectoryFormat file?

	protected:
		virtual void appendHeader( std::ofstream &_traFile ) const;
		virtual void appendData( std::ofstream &_traFile );
		virtual int getHeaderSize() const;
		virtual int getBlockSize() const;

		size_t m_TraCapacity;
	};

	//-------------------------------------------------
	//
	/// \brief  derive a class PD_API from this to allow easy TrajectoryFormat-Comment support within that class
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
	class PD_API BTF_Block_Comment_User
	{
	protected:
		BTF_Block_Comment *m_ExtComments;
	public:
		BTF_Block_Comment_User(): m_ExtComments(NULL)
		{
		}
		inline bool hasTraComments() const { return m_ExtComments != NULL; }
		inline void enableComments( BTF_Block_Comment *comments ) { m_ExtComments = comments; }
		inline BTF_Block_Comment* getTraComments() { return m_ExtComments; } ///> Can return NULL, ask hasTraComments()
	};

	/// \brief  Defines a class that uses a BTF_Block_Vector.
	/// \details A client class can derive from BTF_Block_Vector_User to allow the user to supply 
	/// the class with a m_ExtVectors.
	/// \author Jon Rea 
	class PD_API BTF_Block_Vector_User
	{
	protected:
		BTF_Block_Vector *m_ExtVectors;
	public:
		BTF_Block_Vector_User() : m_ExtVectors(0)
		{
		}
		inline bool hasTraVectors() const { return m_ExtVectors != NULL; }
		inline void enableVectors( BTF_Block_Vector *vectors ) { m_ExtVectors = vectors; }
		inline BTF_Block_Vector* getTraVectors() { return m_ExtVectors; } ///> Can return NULL, ask hasTraVectors()
	};
}  // namespace IO

#endif

