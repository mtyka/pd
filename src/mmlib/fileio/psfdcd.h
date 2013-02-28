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

#ifndef __PSFDCD_H
#define __PSFDCD_H

#include "outtra.h"
#include "workspace/workspace.fwd.h"


namespace IO 
{

	//-------------------------------------------------
	//
	/// \brief Saves a WorkSpace, System or Molecule to a PSF file
	///       
	/// \details This class is used to create PSF files from a WorkSpace, System or Molecule.
	///          Note that this only creates single files, not trajectories. For Trajectories use
	///          the class hierarchy deriving from OutputTrajectory.
	///          FileOutut_PSF can be called in two different ways. Either a system/workspace is passed as an argument
	///          or the "save" function of a system/workspace is used. Both lead to the same result:
	///          THe constructor requires a filename or a filestem (the extension is added automatically if its missing). 
	/*! 
	\code

	// examples
	OutputFile_PSF  mypdb("test.pdb");
	mypdb.save( mysystem );
	mypdb.save( myworkspace );
	mypdb.save( mymolecule );

	mysystem.save( OutputFile_PSF("test.pdb") );
	myworkspace.save( OutputFile_PSF("test.pdb") );
	mymolecule.save( OutputFile_PSF("test.pdb") );

	\endcode
	*/
	class PD_API OutputFile_PSF: public OutputFile
	{
	public:

		OutputFile_PSF(const std::string& _filestem): OutputFile( _filestem ) {}
		virtual OutputFile* clone() const { return new OutputFile_PSF( *this ); } 

		/// Saves a work space to a PSF file
		virtual void save( WorkSpace &_wspace);
	};






	//-------------------------------------------------
	//
	/// \brief Saves a WorkSpace, System or Molecule to a DCD file
	///       
	/// \details This class is used to create DCD files from a WorkSpace, System or Molecule.
	///          Note that this only creates single files, not trajectories. For Trajectories use
	///          the class hierarchy deriving from OutputTrajectory.
	///          FileOutut_DCD can be called in two different ways. Either a system/workspace is passed as an argument
	///          or the "save" function of a system/workspace is used. Both lead to the same result:
	///          THe constructor requires a filename or a filestem (the extension is added automatically if its missing). 
	/*! 
	\code

	// examples
	OutputFile_DCD  mypdb("test.pdb");
	mypdb.save( mysystem );
	mypdb.save( myworkspace );
	mypdb.save( mymolecule );

	mysystem.save( OutputFile_DCD("test.pdb") );
	myworkspace.save( OutputFile_DCD("test.pdb") );
	mymolecule.save( OutputFile_DCD("test.pdb") );

	\endcode
	*/
	class PD_API OutputFile_DCD: public OutputFile
	{
	public:
		OutputFile_DCD(const std::string& _filestem): OutputFile( _filestem ) {}
		virtual OutputFile* clone() const { return new OutputFile_DCD( *this ); } 

		/// Saves a work space to a DCD file
		virtual void save( WorkSpace &_wspace);
	};




	//-------------------------------------------------
	//
	/// \brief Implements the NAMD trajectory format, i.e. a .psf storing the topology + a .dcd file storing the coordinates 
	///
	/// \details Works like all OutputTrajectories, see base class OutputTrajectoryFile for usage.  
	///
	/// \author Mike Tyka  
	///
	///
	class PD_API OutTra_NAMD: public OutputTrajectoryFile {
	public:
		OutTra_NAMD(
			const std::string &_filestem,
			WorkSpace &_wspace)
		:	OutputTrajectoryFile(_filestem, _wspace)
		{
			filename_PSF = filestem + std::string(".psf");
			filename_PDB = filestem + std::string(".pdb");
			filename_DCD = filestem + std::string(".dcd");
			created = false;
		}
		
		virtual OutTra_NAMD* clone() const { return new OutTra_NAMD(*this); }

		virtual int create();
		virtual int append();

		friend class OutputFile_PSF;
		friend class OutputFile_DCD;

	private:

		// PSF Members
		void printPSFfile( const char *filename );
		int printPSFfile_structure(FILE *file);
		int printPSFfile_structurePDB(FILE *file);

		// DCD Members
		int write_DCD_header(const char *filename);
		int append_DCD_step(const char *filename);

		std::string filename_PSF;
		std::string filename_PDB;
		std::string filename_DCD;

		bool created;
	};



} //namespace IO




#endif
