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

#ifndef __BTF_H
#define __BTF_H

#include <vector>
#include <fstream>

#include "intra.h"
#include "outtra.h"
#include "tratypes.h"
#include "trablocks.h"
#include "workspace/workspace.fwd.h"

class PD_API System;
class PD_API MoleculeBase;

namespace Sequence
{
	class PD_API BioSequence;
}

namespace IO 
{
	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
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
	class PD_API BTF_Base
	{
	public:
		BTF_Base(){}
		virtual BTF_Base* clone() const = 0;

	protected:
		BTF_Header m_Header;
		BTF_Energy m_Ene;
	};


	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
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
	class PD_API OutTra_BTF : public BTF_Base, public OutputTrajectoryFile
	{
	public:
		OutTra_BTF(const std::string &_filestem, WorkSpace& _wspace);
		virtual OutTra_BTF* clone() const { return new OutTra_BTF(*this); }

		int create(BFT_StandardIncludes _includes);

		virtual int create();
		virtual int append();

		void addOwnedBlock( BTF_Block* _Block ); ///< Add a block that is memory managed by this class
		void addBlock( BTF_Block& _Block ); ///< Add a block that is memory managed elsewhere

	protected:
		std::string filename;
		ObjectContainer<BTF_Block> m_Blocks;
	};


	//-------------------------------------------------
	//
	/// \brief  A base class for files that read in trajectories.
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
	class PD_API BTF_ImportBase : public BTF_Base
	{
	public:
		enum AtomFilter
		{
			CA = 1,
			N = 2,
			O = 4,
			C = 8,
			Sidechain = 16,
			Hydrogens = 32,

			Backbone = CA | N | C | O,
			HeavyAtom = Backbone | Sidechain,
			All = Backbone | Sidechain | Hydrogens
		};

		BTF_ImportBase( const std::string &_fileName, bool _FullFileValidation);
		virtual BTF_ImportBase* clone() const = 0;

		void reOpen(const std::string &_fileName, bool _FullFileValidation = true);
		void validateFile(); ///< Performs a validation on all the Sanity-Tags that should be contained within the tar file
		void info( bool _verbose );

		int recountEntries(); ///< Just in case the file is being written to, we can do an internal recount - the header and system cannot change.
		int getEntryCount() const; ///< The total number of entries in the file
		int size() const;

		const Sequence::BioSequence& getSequence() const;

	protected:
		// Parsing
		void assertTag( std::ifstream &_stream, std::string _tag ); // parsing assistance function
		void load(); // Load the header and system definitions from the binary file.

		/// loadVectors() is applicable to any dvector  array stored in the tra file,
		/// that is the same length as that of m_Header.atoms.
		/// This includes the atom positions, the forcevectors, and velocities if these are added later ...
		void loadVectors( std::ifstream &_traFile, std::vector<Maths::dvector>& _storage, AtomFilter _filter, int _molNum );
		int passesFilter( int _Index, AtomFilter _filter, int _molNum );

		void seekToEntry( std::ifstream &_stream, int _entry );

		std::string m_Filename;
		int m_Entries;
		std::vector<BTF_SystemDefinitionEntry> sysDefs;
	};


	//-------------------------------------------------
	//
	/// \brief  A class to create a System from a given entry in a TrajectoryFormat file.
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
	class PD_API InTra_BTF: public BTF_ImportBase, public InputTrajectory_RandomAccess
	{
	public:
		InTra_BTF( const std::string &_fileName, bool _ValidateFile = true);
		virtual InTra_BTF* clone() const { return new InTra_BTF(*this); }

		/// The traditional InTra_BTF call to import into a system
		void loadIntoSystem( System &sysspec, int position );

		/// \brief Reads a structure and returns a SnapShot with the coordinates (and the box geometry)
		/// \return Returns true if end of file was reached during read, false otherwise.
		virtual bool readNext( SnapShot &ss );

		/// \brief Skips an entry in the file 
		/// \return Returns true if end of file was reached during read, false otherwise.
		virtual bool skip();

		/// \brief Returns true if file handle has reached end of file. 
		/// \return Returns true if end of file was reached during read, false otherwise.
		virtual bool isEndOfFile() const;

		/// \brief Go back to the start of the tra
		virtual void reset();

		virtual void readRandomAccess( SnapShot &ss, size_t entry );

		/// \brief works out how many entries are in the trajectory
		virtual size_t nEntries() const;

	protected:
		SnapShot makeSnapShot( size_t entry );

		size_t currentPos;
	};


	//-------------------------------------------------
	//
	/// \brief A class to allow analysis of a TrajectoryFormat file in the absence of the FFPS or FF classes - e.g. geometry analysis and distance calculations
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
	class PD_API BTF_Tools: public BTF_ImportBase
	{
	public:

		BTF_Tools(std::string &_fileName, bool _FullFileValidation = true);
		virtual BTF_Tools* clone() const { return new BTF_Tools(*this); }

		void getCoordinates( int _Entry );
		void getCoordinates( int _Entry, AtomFilter _filter );
		void getCoordinates( int _Entry, AtomFilter _filter, int _molNum );

		void savePDBFile( std::string &_FileName, int _Entry );
		void savePDBFile( std::string &_FileName, int _Entry, AtomFilter _filter );
		void savePDBFile( std::string &_FileName, int _Entry, AtomFilter _filter, int _molNum );

	protected:
		std::vector<Maths::dvector> positions;
		std::vector<Maths::dvector> forces;
		BTF_Energy ene;
	};

	int loadtra(System &sysspec, const std::string &filename, const std::string &entrystring);


	//-------------------------------------------------
	//
	/// \brief Saves a WorkSpace, System or Molecule to a BTF file
	///       
	/// \details This class is used to create BTF files from a WorkSpace, System or Molecule.
	///          Note that this only creates single files, not trajectories. For Trajectories use
	///          the class hierarchy deriving from OutputTrajectory.
	///          FileOutut_BTF can be called in two different ways. Either a system/workspace is passed as an argument
	///          or the "save" function of a system/workspace is used. Both lead to the same result:
	///          THe constructor requires a filename or a filestem (the extension is added automatically if its missing). 
	/*! 
	\code

	// examples
	OutputFile_BTF  mypdb("test.pdb");
	mypdb.save( mysystem );
	mypdb.save( myworkspace );
	mypdb.save( mymolecule );

	mysystem.save( OutputFile_BTF("test.pdb") );
	myworkspace.save( OutputFile_BTF("test.pdb") );
	mymolecule.save( OutputFile_BTF("test.pdb") );

	\endcode
	*/
	/// \author Mike Tyka 
	///
	class PD_API OutputFile_BTF: public OutputFile
	{
	public:
		OutputFile_BTF(const std::string& _filestem): OutputFile( _filestem ) {}
		virtual OutputFile* clone() const { return new OutputFile_BTF( *this ); } 

		/// Saves a work space to a BTF file
		virtual void save( WorkSpace &_wspace );
	};
}

#endif

