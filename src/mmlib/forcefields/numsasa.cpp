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
#include <fstream>
#include "tools/stringtool.h"
#include "pickers/pickbase.h"
#include "workspace/workspace.h"
#include "numsasa.h"

NumSASA::NumSASA()
{
}

NumSASA::NumSASA( const std::string& _DatFileName )
{
	readDat(_DatFileName);
}

void NumSASA::assertData()
{
	//ASSERT( m_WSpace != NULL, ProcedureException, "The internal workspace pointer is null; setTo() has not been called");
	//ASSERT( m_Sphere.size() > 0, ProcedureException, "The Numeric SASA has not been initialised with its sphere data");
	//ASSERT( m_AtomIndexes.size() > 0, ProcedureException, "There is no sasa data; has calc been called? or does the  
}

void NumSASA::readDat(const std::string& _DatFileName)
{
	//m_Sphere.clear();
}

void NumSASA::setTo( const WorkSpace& _WSpace )
{
	//m_Picker = PickAllAtoms();
}

void NumSASA::detail() const
{
}

void NumSASA::info() const
{
}

void NumSASA::calc()
{
}


