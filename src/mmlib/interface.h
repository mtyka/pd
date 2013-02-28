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

#ifndef __INTERFACE_H
#define __INTERFACE_H

// Defines a list of public ABC interfaces 
// that can be used throughout PD

// Forward declarations
struct DrawingVector; // defined in primitives.h

/// \brief  IDrawProvider
/// \details Defiens a class interface for those which can draw lines for user feedback
/// \author  Jon Rea 
class IDrawProvider
{
public:
	virtual DrawingVector* request() = 0;
};

/// \brief ICloneable - any clonable class should define itself as ICloneable
/// \details ICloneable - any clonable class should define itself as ICloneable
/// \author  Jon Rea 
class ICloneable
{
public:
	virtual ICloneable* clone() const = 0;
};

#endif

