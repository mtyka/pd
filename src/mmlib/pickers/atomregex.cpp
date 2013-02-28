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
#include "tools/stringbuilder.h"
#include "system/fundamentals.h"
#include "atomregex.h"

// ------------------------------------
// Begin Class PickAtomRegEx
// ------------------------------------

PickAtomRegEx::PickAtomRegEx( const std::string& _pattern ) 
{ 
	setPattern(_pattern); 
}

PickAtomRegEx* PickAtomRegEx::clone() const
{
	return new PickAtomRegEx(*this);
}

void PickAtomRegEx::setPattern( const std::string& _pattern )
{
	m_AtomNames.clear();
	m_Pattern = _pattern;
	StringBuilder sb( _pattern );
	size_t index;
	while( (index = sb.LastOf(PickAtomRegExDelimiter)) != SIZE_T_FAIL )
	{
		m_AtomNames.push_back(sb.toString(index,sb.size()-index));
		sb.erase(index,sb.size()-index);
	}
	m_AtomNames.push_back(sb.toString());
}

bool PickAtomRegEx::matches( const Particle& particle ) const
{
	for( size_t i = 0; i < m_AtomNames.size(); i++ )
	{
		if( 0 == particle.pdbname.compare( m_AtomNames[i] ) ) 
		{
			return true;
		}
	}
	return false;
}

