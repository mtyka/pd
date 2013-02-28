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
#include "filterbase.h"

FilterBase::FilterBase()
{
}

FilterBase::~FilterBase()
{
}

void FilterBase::setInternalName()
{
	name = "FilterBase";
}

MolFilterBase::MolFilterBase() 
	: FilterBase(), 
	m_Mol(NULL) 
{
}

MolFilterBase::MolFilterBase( const MoleculeBase& _mol ) 
	: FilterBase(), 
	m_Mol(&_mol) 
{
}

void MolFilterBase::setInternalName()
{
	name = "MolFilterBase";
}

FilterContainer::FilterContainer()
{
}

void FilterContainer::setInternalName()
{
	name = "FilterContainer";
}

bool FilterContainer::passes()
{
	for( size_t i = 0; i < size(); i++ )
	{
		if( !element(i).passes() ) 
			return false;
	}
	return true;
}

std::string FilterContainer::reason()
{
	StringBuilder sb;
	bool foundBad = false;
	for( size_t i = 0; i < size(); i++ )
	{
		if( !element(i).passes() ) 
		{
			if( !foundBad ) 
			{
				sb.appendFormat("Filter Collection `%s' failed:\n")(name);
				foundBad = true;
			}
			sb.appendFormat("  (%d/%d) `%s' failed: %s\n")
				(i+1)
				(size())
				(element(i).name.c_str())
				(element(i).reason().c_str());						
		}
	}
	if( !foundBad ) sb.appendFormat("All %d Internal filters passed\n")(size());
	return sb.toString();
}

