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

#ifndef __TRAITS_H
#define __TRAITS_H

/// trait used in nonbonded_periodic.h to select which block 
/// of inline code to compile in
template< typename T_Type1, typename T_Type2 >
struct is_same_type{ static const bool value = false; };
/// This template specialisation essentially acts as a compile-time
/// if statment: if you access is_same_type<typea,typeb>::value it will
/// be true if typea and typeb are the same Type and false otherwise.
/// used in an if statement the compiler can compile away the if
/// because the result is already known at compile time (as if you
/// said if(true){ ... }
template< typename T >
struct is_same_type<T,T>{
	static const bool value = true;
}; 

#endif

