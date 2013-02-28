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

#ifndef __GROMACS_H
#define __GROMACS_H

#include "global.h"

#include "mmlib/fileio/outtra.h"
#include "mmlib/fileio/intra.h"
#include "mmlib/workspace/workspace.h"
#include "mmlib/workspace/space.h"

#include "rpcxdr.h"

namespace IO 
{


	class PD_API OutTra_GRO: public OutputTrajectoryFile 
	{
	public:
		OutTra_GRO(
			std::string _filestem,
			WorkSpace &_wspace 
		):
			OutputTrajectoryFile(_filestem, _wspace )
		{
			filename_GRO = filestem + std::string(".gro");
			created = false;
		}

		virtual OutTra_GRO* clone() const { return new OutTra_GRO(*this); }

		virtual int create();
		virtual int append(){ if(!created) create(); return 0; };

	private:
		std::string filename_GRO;
	};




	class PD_API OutTra_XTC: public OutputTrajectoryFile
	{
	public:
		OutTra_XTC(
			std::string _filestem,
			WorkSpace &_wspace
		):
			OutputTrajectoryFile(_filestem, _wspace)
		{
			xdrs = NULL;
			filename_XTC = filestem + std::string(".xtc");
			m_StepCount = 0;
			created = false;
		}

		~OutTra_XTC();

		virtual OutTra_XTC* clone() const { return new OutTra_XTC(*this); }

		virtual int create();
		virtual int append();

	private:
		std::string filename_XTC;
		XDR *xdrs;

		int m_StepCount;
	};



		
	//-------------------------------------------------
	//
	/// \brief Saves a WorkSpace, System or Molecule to a GRO file
	///       
	/// \details This class is used to create GRO files from a WorkSpace, System or Molecule.
	///          Note that this only creates single files, not trajectories. For Trajectories use
	///          the class hierarchy deriving from OutputTrajectory.
	///          FileOutut_GRO can be called in two different ways. Either a system/workspace is passed as an argument
	///          or the "save" function of a system/workspace is used. Both lead to the same result:
	///          THe constructor requires a filename or a filestem (the extension is added automatically if its missing). 
/*! 
\code

// examples
OutputFile_GRO  mypdb("test.pdb");
mypdb.save( mysystem );
mypdb.save( myworkspace );
mypdb.save( mymolecule );

mysystem.save( OutputFile_GRO("test.pdb") );
myworkspace.save( OutputFile_GRO("test.pdb") );
mymolecule.save( OutputFile_GRO("test.pdb") );

\endcode
*/
	/// \author Mike Tyka 
	///
	class PD_API OutputFile_GRO: public OutputFile
	{
	public:
		OutputFile_GRO(const std::string& _filestem): OutputFile( _filestem ) {}
		virtual OutputFile* clone() const { return new OutputFile_GRO( *this ); } 

		/// Saves a work space to a GRO file
		virtual void save( WorkSpace &_wspace );
	};




		
	//-------------------------------------------------
	//
	/// \brief Saves a WorkSpace, System or Molecule to a XTC file
	///       
	/// \details This class is used to create XTC files from a WorkSpace, System or Molecule.
	///          Note that this only creates single files, not trajectories. For Trajectories use
	///          the class hierarchy deriving from OutputTrajectory.
	///          FileOutut_XTC can be called in two different ways. Either a system/workspace is passed as an argument
	///          or the "save" function of a system/workspace is used. Both lead to the same result:
	///          THe constructor requires a filename or a filestem (the extension is added automatically if its missing). 
/*! 
\code

// examples
OutputFile_XTC  mypdb("test.pdb");
mypdb.save( mysystem );
mypdb.save( myworkspace );
mypdb.save( mymolecule );

mysystem.save( OutputFile_XTC("test.pdb") );
myworkspace.save( OutputFile_XTC("test.pdb") );
mymolecule.save( OutputFile_XTC("test.pdb") );

\endcode
*/
	/// \author Mike Tyka 
	///
	class PD_API OutputFile_XTC: public OutputFile
	{
	public:
		OutputFile_XTC(const std::string& _filestem): OutputFile( _filestem ) {}
		virtual OutputFile* clone() const { return new OutputFile_XTC( *this ); } 

		/// Saves a work space to a XTC file
		virtual void save( WorkSpace &_wspace );
	};






	class PD_API InTra_XTC: public InputTrajectory
	{
	public:

		InTra_XTC(const std::string& _filename):
		InputTrajectory( _filename)
		{
		  xdrs = NULL;
		  open();
		}

		virtual ~InTra_XTC () 
		{ 
		  close();  
		}

		virtual InTra_XTC* clone() const { return new InTra_XTC(*this); }

		virtual bool readNext( SnapShot &ss );
		virtual bool skip();
		virtual bool isEndOfFile() const;
		virtual void reset(){ close(); open(); };
	protected:
		virtual int open();
		virtual int close();

	private:

		XDR *xdrs;
	};


}

#endif

