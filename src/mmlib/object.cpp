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

long Object::object_nextid = 0;

Object::Object()
{
	object_id = object_nextid;
	object_nextid++;
	setInternalName();
}

Object::Object( const Object & _clone )
{
	(*this) = _clone;
}

Object & Object::operator= ( const Object & _clone )
{
	if( &_clone == this ) return *this; // Handle self-assignement
	name = _clone.name ;							 
	object_id = object_nextid;
	object_nextid++;
	return (*this);
}

bool operator== (const Object &a, const Object &b)
{
	return a.object_id == b.object_id;
}

void Object::setInternalName()
{
	name = "Object";
}

