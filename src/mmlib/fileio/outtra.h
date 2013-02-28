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

#ifndef __OUTTRA_H
#define __OUTTRA_H

#include <string>
#include "object.h"

#include "system/system.fwd.h"
#include "system/molecule.fwd.h"
#include "workspace/workspace.fwd.h"
#include "workspace/componentbase.h"

namespace IO 
{

	//-------------------------------------------------
	//
	/// \brief  Base class for all Output Trajectory Formats. This class essentially 
	///         attaches to a workspace and one can use the functions create, attach 
	///         and save to save snapshots of trjectories into files or other memory structures.
	///          
	/// \details The difference between this class and OutputFile is that it supports multiple
	///         writes to the same file but it assumes that the number of atoms and the topology
	///         are unchanged between frames, as this is common to most trajectory formats.
	///         Hence this class "attaches" permanently to a WorkSpace (and not a System)
	///         A reference to a WorkSpace must be given at creation time.
	///
	/// \author Mike Tyka & Jon Rea 
	///

	class PD_API OutputTrajectory: public Object, 
		public WorkSpaceOperatorBase
	{
	public:
		OutputTrajectory( WorkSpace &_wspace );
		virtual ~OutputTrajectory () {}
		virtual OutputTrajectory* clone() const = 0;

		/// Create the trajectory
		virtual int create() = 0;

		/// Append a frame to the trajectory
		virtual int append() = 0;

	protected:

		/// flags if the trajectory has been created already (in most cases this means
		/// whether the header has already been written. append() uses this boolean to assess
		/// if create() has yet to be called.
		bool created;
	};






	//-------------------------------------------------
	/// \brief  Base class for all Output Trajectory Formats which write to files. This class essentially 
	///         attaches to a workspace and one can use the functions create, attach 
	///         and save to save snapshots of trjectories into files.
	///          
	/// \details The difference between this class and OutputFile is that it supports multiple
	///         writes to the same file but it assumes that the number of atoms and the topology
	///         are unchanged between frames, as this is common to most trajectory formats.
	///         Hence this class "attaches" permanently to a WorkSpace (and not a System)
	///         A reference to a WorkSpace must be given at creation time.
	///
	/// \author Mike Tyka & Jon Rea 
	///

	class PD_API OutputTrajectoryFile: public OutputTrajectory
	{
	public:
		OutputTrajectoryFile(const std::string& _filestem, WorkSpace &_wspace);
		virtual ~OutputTrajectoryFile () {}
		virtual OutputTrajectoryFile* clone() const = 0;

		/// Create the trajectory
		virtual int create() = 0;

		/// Append a frame to the trajectory
		virtual int append() = 0;

	protected:
		std::string filestem;
	};

}

#ifdef SWIG
%template(ObjectContainer_OutputTrajectory) ObjectContainer<IO::OutputTrajectory>;
#endif

namespace IO
{




	//-------------------------------------------------
	//
	/// \brief  This class ia a container for multiple OutputTrajectories. 
	/// \details Its main function 
	///         is inside Workspace where it hold all the trajetories currently loaded into that
	///         WorkSpace, but it can also be used outside to group trajectories. It itself derives from
	///         OutputTrajectory and thus can be used in all situations where an OutputTrajectory is required.
	///
	/// \author Mike Tyka & Jon Rea 
	///
	class PD_API OutputTrajectoryContainer: 
		public OutputTrajectory, 
		public ObjectContainer<OutputTrajectory> 
	{
	public:
		/// OutputTrajectoryContainer must take a WorkSpace at creation time
		OutputTrajectoryContainer(WorkSpace &_wspace);
		virtual ~OutputTrajectoryContainer () {}
		virtual OutputTrajectoryContainer* clone() const { return new OutputTrajectoryContainer(*this); }

		/// Create all trajectories currently stored
		virtual int create();

		/// Append a frame to all trajectories currently stored
		virtual int append();

	private:

		/// Calc some stuff before subclasses are asked to append
		void prepare(); 
	};



	//-------------------------------------------------
	//
	/// \brief  Base Class for saving the state of the WorkSpace, System or a Molecule. 
	/// \details This class is used to create files from a WorkSpace, System or Molecule.
	///          Note that this only creates single files, not trajectories. For Trajectories use
	///          the class hierarchy deriving from OutputTrajectory.
	///          OutputFile can be called in two different ways. Either a system/workspace is passed as an argument
	///          or the "save" function of a system/workspace is used. Both lead to the same result:
	///          The constructor requires a filename or a filestem (the extension is added automatically if its missing). 
	///
	///          Intended use of the derived class (here for a format XXX) is demonstrated below: 
				 /*! 
				 \code 

				 // examples
				 OutputFile_XXX  mypdb("test.pdb");
				 mypdb.save( mysystem );
				 mypdb.save( myworkspace );
				 mypdb.save( mymolecule );

 				 mysystem.save( OutputFile_XXX("test.pdb") );
				 myworkspace.save( OutputFile_XXX("test.pdb") );
				 mymolecule.save( OutputFile_XXX("test.pdb") );

				 \endcode
				 */
	/// \details 
	///
	/// \author Mike Tyka 
	///
	class PD_API OutputFile
	{
	public:
		OutputFile(const std::string& _filestem);
		virtual OutputFile* clone() const = 0;

		friend class ::System;
		friend class ::Molecule;
		friend class ::WorkSpace;

		virtual void save( WorkSpace &_wspace) = 0;
		virtual void save( System &_system);
		virtual void save( Molecule &_system);

	protected:
		std::string filestem;
	};

} // namespace IO

#endif

