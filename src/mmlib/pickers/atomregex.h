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

#ifndef __ATOM_REGEX_H
#define __ATOM_REGEX_H

#include <string>
#include <vector>

#include "pickers/pickbase.h" // PickBase base class

const char PickAtomRegExDelimiter = '-';

//-------------------------------------------------
/// \brief  A simple RegEx atom picker
/// \details Currently this is the simplest of the simple REGEX for atoms, only individual names can be
/// expressed, separated by the 'PickAtomRegExDelimiter'.
/// This class should be extended to allow much more complex atom and parent residue expressions
/// but at the moment I don't have the time to code it...
/// \author Jon Rea 
/// \todo State of Development: Primordial Soup
/// \bug More than likely
class PickAtomRegEx : public PickBase
{
public:
	// Public constructor logic
	PickAtomRegEx();
	PickAtomRegEx( const std::string& _pattern );
	virtual PickAtomRegEx* clone() const;

	// Public function calls
	void setPattern( const std::string& _pattern ); ///< set the atom selection patten
	virtual bool matches( const Particle& particle ) const; ///< Returns true if the passed reference matches the internal pattern
protected:
	std::string m_Pattern;
	std::vector<std::string> m_AtomNames;
};

#endif

