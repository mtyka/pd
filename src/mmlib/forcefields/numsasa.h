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

#ifndef __NUMSASA_H
#define __NUMSASA_H

#include "workspace/workspace.fwd.h"






//-------------------------------------------------
//
/// \brief  Implements calculation of numerical SASAs
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API NumSASA
{
public:
	NumSASA();
	NumSASA( const std::string& _DatFileName );

	// read out data file
	void readDat(const std::string& _DatFileName);
	void setTo( const WorkSpace& _WSpace );

	// Info
	void detail() const;
	void info() const;

	/// Calculate the SASA
	void calc();
	
private:
	void assertData();

	//std::vector<dvector> m_Sphere;
};

#endif


