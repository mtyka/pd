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
#include "outtra.h"

#include "system/system.h"
#include "system/molecule.h"
#include "workspace/workspace.h"
#include "workspace/componentbase.h"

namespace IO
{
	// ------------------
	//  OutputTrajectory
	// ------------------

	OutputTrajectory::OutputTrajectory(WorkSpace &_wspace) :
		WorkSpaceOperatorBase( _wspace ),
		created(false)
	{
	}

	// ----------------------
	//  OutputTrajectoryFile
	// ----------------------

	OutputTrajectoryFile::OutputTrajectoryFile(
		const std::string& _filestem, 
		WorkSpace & _wspace
	) :
		filestem(_filestem), 
		OutputTrajectory( _wspace )
	{
	}

	// ---------------------------
	//  OutputTrajectoryContainer
	// ---------------------------

	OutputTrajectoryContainer::OutputTrajectoryContainer( 
		WorkSpace & _wspace 
	) :
		OutputTrajectory( _wspace )
	{
	}

	void OutputTrajectoryContainer::prepare()
	{
		if( size() > 0 ) 
		{
			getWSpace().calcCRMS_HeavyAtom();
		}
	}

	int OutputTrajectoryContainer::create()
	{
		prepare();
		for(unsigned i=0;i<size();i++) 
			element(i).create();
		return 0;
	}

	int OutputTrajectoryContainer::append()
	{
		prepare();
		for(unsigned i=0;i<size();i++) 
		{
			element(i).append();
		}
		return 0;
	}

	// -------------
	//  OutputFile
	// -------------

	OutputFile::OutputFile(const std::string& _filestem)
	{
		filestem = _filestem;
	}

	void OutputFile::save( System &_system)
	{
		WorkSpace wspace( _system );
		save( wspace );
	}

	void OutputFile::save( Molecule &_molecule)
	{
		System mysystem( _molecule.ffps() );
		mysystem.add( _molecule );
		save( mysystem );
	}
} // namespace IO


