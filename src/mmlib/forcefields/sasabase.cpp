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
#include "pickers/basicpickers.h"
#include "sasabase.h"

namespace Physics
{
	FF_SurfaceArea::FF_SurfaceArea(WorkSpace &newwspace) 
		: ForcefieldBase( newwspace )
	{
		setToDefaults();
	}

	FF_SurfaceArea::~FF_SurfaceArea()
	{
	}

	void FF_SurfaceArea::setToDefaults()
	{
		// At the moment, we will default to choosing all atoms. This will of coarse yield a 
		// nonsense answer if the system contains explicit solvent. There has been a mild debate as to
		// wether or not a PickNotSolvent() class would be a good plan or wether it is better for the
		// user to set the for and against pickers to a custom selection of only a protein molecule or
		// a all molecules in a complex. If it is decided that a pick solvent is needed, then a new atom
		// flag isSolvent() may be the way to do - this would be set post-molecule-build using the
		// information contained within the current class system (file: default.class).
		m_ForPicker = PickAllAtoms();
		m_AgainstPicker = PickAllAtoms();
	}

	void FF_SurfaceArea::setForPicker( const PickBase& _ForPicker )
	{
		m_ForPicker = _ForPicker;
	}

	void FF_SurfaceArea::setAgainstPicker( const PickBase& _AgainstPicker )
	{
		m_AgainstPicker = _AgainstPicker;
	}
}

