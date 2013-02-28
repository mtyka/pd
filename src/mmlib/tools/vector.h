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

#ifndef __VECTOR_TOOL_H
#define __VECTOR_TOOL_H

#include <vector>
#include <algorithm>
#include <functional>

///\brief Removes elements from a vector 
///\details 
/// Uses the std::remove_if algorithm (remove-erase-idiom)
/// Implements the link between a member function that returns bool for a const element where
/// element is the content type of a std::vector<element>. The element is copied in the 
/// procedure, but this is usually not limiting in any way.
/// Usage Example:
/// removematches( myVector, &MyClass::matches, this ); 
/// where 'this' is a 'this pointer' to the instance of MyClass and shouldRemove() is the member function of type:
/// bool MyClass::shouldRemove( const TArg bond ) const;
template< typename TParentClass, typename TArg >
typename std::vector<TArg>::iterator removematches( 
	std::vector<TArg>& vect, 																									
	bool(TParentClass::*matchPFunc)(const TArg) const, 
	TParentClass* matchFuncClassInstance )
{
	return vect.erase( 
		remove_if( 
			vect.begin(), 
			vect.end(), 
			std::bind1st(std::mem_fun(matchPFunc), 
			matchFuncClassInstance)
			),
		vect.end() 
		);
}

///\brief Decreases the capacity of a vector back down to that of its size 
///\details 
/// *-* TAKE CARE *-*: this function will invalidate all pointers and itterators refereing to 
/// this std::vector<T> bloatedVector!
/// This creates a new temporary unnamed vector from the bloated vector and then swaps 
/// the internal memory buffers.
/// The result is that the capacity of the the bloatedVector capactity() == size()
template < typename T >
inline void TrimVectorCapacity( std::vector<T>& bloatedVector )
{
	std::vector<T>(bloatedVector).swap(bloatedVector);
}

///\brief Returns true if the given array contains the given value
template < typename T >
inline bool VectorContains( const std::vector<T>& vect, T value )
{
	typename std::vector<T>::const_iterator iter;
	for( iter = vect.begin(); iter != vect.end(); iter++ )
	{
		if( *iter == value ) return true;
	}
	return false;
}

template < typename T, typename U >
inline size_t FindFirstInVector( const std::vector<T>& vect, U value )
{
	for( size_t i = 0; i < vect.size(); i++ )
	{
		if( vect[i] == value )
		{
			return i;
		}
	}
	return SIZE_T_FAIL;
}

template < typename T >
inline T Average( std::vector<T>& vect )
{
	size_t size = vect.size();
	T sum = (T)0.0;
	for( size_t i = 0; i < size; i++ )
	{
		sum += vect[i];
	}
	return (T)((double)sum / (double)size);
}

// I would template these, other than the fact that you need boost to 
// get _MIN and _MAX in templated form for a given type...

inline double Min( std::vector<double>& vect )
{
	double min = DBL_MAX;
	size_t size = vect.size();
	for( size_t i = 0; i < size; i++ )
	{
		min = std::min( vect[i], min );
	}
	return min;
}

inline double Max( std::vector<double>& vect )
{
	double max = DBL_MIN;
	size_t size = vect.size();
	for( size_t i = 0; i < size; i++ )
	{
		max = std::max( vect[i], max );
	}
	return max;
}

#endif

