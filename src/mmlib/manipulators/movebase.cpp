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
#include "manipulators/movebase.h"
#include "workspace/workspace.h"

namespace Manipulator
{
	MoveBase::MoveBase(WorkSpace& _wspace) 
		: wspace(&_wspace)
	{
	}

	MoveBase::~MoveBase()
	{
	}

	// Applies the move 'rounds' times for testing purposes
	void MoveBase::test(int rounds)
	{
		for(int i = 0; i < rounds; i++) 
		{
			apply();
			wspace->outtra.append();
		}
	}

	MoveSet::MoveSet(WorkSpace& newwspace)
		: MoveBase(newwspace) 
	{
		name = "MoveSet";
	}
		
	MoveSet::~MoveSet()
	{
	}

	MoveSet* MoveSet::clone() const 
	{ 
		return new MoveSet(*this); 
	}

	int MoveSet::apply()
	{
		for(size_t i = 0; i < size(); i++)
		{
			element(i).apply();
		}
		return 0;
	}
}

