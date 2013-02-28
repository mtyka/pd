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

#ifndef __STRING_BUILDER_H
#define __STRING_BUILDER_H

#include <string>
#include <vector>
#include <ostream>

#include "typedefs.h"

namespace LokiCore
{
	template <class Device, class Char> struct PrintfState;
}

// ----------------------------------------------------------------------------------
// --- A special custom string manipulation class
// ----------------------------------------------------------------------------------

// Class Declaration
class PD_API StringBuilder
{
public:
	// Allow use with cout, this function needs to access private member data

	static const size_t SB_DEFAULT_CAPACITY = 128; // seems like a normal size of string...

	StringBuilder( size_t _InitialCapacity = SB_DEFAULT_CAPACITY );
	StringBuilder( char _Fill, size_t _count );
	StringBuilder( const std::string &_Clone );
	StringBuilder( const StringBuilder &_Clone ); // Copy constructor
	virtual ~StringBuilder();

#ifndef SWIG
	friend std::ostream& operator<<(std::ostream &s, const StringBuilder &_Print);
	friend std::istream& operator<<(StringBuilder &_Print, std::istream &s);
	StringBuilder& operator=(const StringBuilder &_Clone);
	StringBuilder& operator=(const std::string &_Clone);
	const char& operator[](size_t i) const;
#endif

	void clear();
	void growCapacity( double _reallocFactor = 1.5 ) const;
	void growCapacity( size_t _capacityIncrease ) const;

	void replace( size_t _index, char _withChar );
	void replace( size_t _index, char _withChar, size_t _count );
	void replace( size_t _index, const char* _string, size_t _length );
	void replace( size_t _index, const std::string &_string );
	void replace( size_t _index, const StringBuilder &_string );
	
	void setTo( char _c );
	void setTo( char _c, size_t _count );
	void setTo( const char* const _string, size_t _length );
	void setTo( const std::string &_string );
	void setTo( const std::string &_string, size_t _index, size_t _length );
	void setTo( const StringBuilder &_string );
	void setTo( const StringBuilder &_string, size_t _index, size_t _length );

	// Print into the stringbuilder in a formatted way using the LokiCore
	LokiCore::PrintfState<StringBuilder&, char> setFormat(const std::string& format);
	LokiCore::PrintfState<StringBuilder&, char> appendFormat(const std::string& format);

	void append( char _c );
	void append( char _c, size_t _count );
	void append( const char* const _string, size_t _length );
	void append( const std::string &_string );
	void append( const std::string &_string, size_t _index, size_t _length );
	void append( const StringBuilder &_string );
	void append( const StringBuilder &_string, size_t _index, size_t _length );

	void append( int a, char *formatstring="%d" );
	void append( double a, char *formatstring="%8.3lf" );
	void append( long a, char *formatstring="%ld" );

	void endl();

	// 'appendingReplace' calls are the same as 'replace' except that you can
	// write past the end of the current string.
	void appendingReplace( size_t _index, char _withChar, size_t _count );
	void appendingReplace( size_t _index, const char* _string, size_t _length );
	void appendingReplace( size_t _index, const std::string &_string );
	void appendingReplace( size_t _index, const StringBuilder &_string );

	void insert( size_t _index, char _withChar );
	void insert( size_t _index, char _withChar, size_t _count );
	void insert( size_t _index, const char* _string, size_t _length );
	void insert( size_t _index, const std::string &_string );
	void insert( size_t _index, const StringBuilder &_string );

	void erase( size_t _index );
	void erase( size_t _index, size_t _count );

	size_t FirstOf( char _delimiter ) const;
	size_t FirstNotOf( char _delimiter ) const;
	size_t LastOf( char _delimiter ) const;
	size_t LastNotOf( char _delimiter ) const;

	size_t FirstOf( const std::string &_delimiters ) const;
	size_t FirstNotOf( const std::string &_delimiters ) const;
	size_t LastOf( const std::string &_delimiters ) const;
	size_t LastNotOf( const std::string &_delimiters ) const;

	size_t XthOf( size_t instance, char _delimiter ) const;
	size_t XthNotOf( size_t instance, char _delimiter ) const;
	size_t XthOf( size_t instance, const std::string &_delimiters ) const;
	size_t XthNotOf( size_t instance, const std::string &_delimiters ) const;

	void Trim( char _delimiter );
	void TrimLeft( char _delimiter );
	void TrimRight( char _delimiter );

	void Trim( const std::string &_delimiters = TOKEN_WHITESPACE );
	void TrimLeft( const std::string &_delimiters = TOKEN_WHITESPACE );
	void TrimRight( const std::string &_delimiters = TOKEN_WHITESPACE );

	void PadRight( size_t _toLength, char _withChar = ' ' );
	void PadLeft( size_t _toLength, char _withChar = ' ' );
	void TruncateRightTo( size_t _toLength );
	void TruncateLeftTo( size_t _toLength );
	void TruncateRightBy( size_t _removeCount );
	void TruncateLeftBy( size_t _removeCount );

	void removeAll( char _delimiter );
	void removeAll( const std::string &_delimiters );

	bool compare( char _c, size_t _index, bool _ignoreCase = false ) const;
	bool compare( const char *_stringB, size_t _indexA, size_t _length, size_t _indexB, bool _ignoreCase = false ) const;
	bool compare( const std::string &_stringB, size_t _indexA, size_t _length, size_t _indexB, bool _ignoreCase = false ) const;
	bool compare( const StringBuilder &_stringB, size_t _indexA, size_t _length, size_t _indexB, bool _ignoreCase = false ) const;

	int compareValue(const StringBuilder &_stringB, size_t _indexA, size_t _length, size_t _indexB, bool _ignoreCase = false ) const;

	void toScreen() const;
	void toScreen( size_t _index, size_t _length ) const;
	std::string toString() const;
	std::string toString( size_t _index, size_t _length ) const;
	std::vector<std::string> tokenise( const std::string &_delimiters = TOKEN_WHITESPACE ) const;

	int parseInt() const;
	int parseInt( size_t _index ) const;
	int parseInt( size_t _index, size_t _length ) const;

	double parseDouble() const;
	double parseDouble( size_t _index ) const;
	double parseDouble( size_t _index, size_t _length ) const;

	size_t capacity() const { return m_alloc; }
	size_t size() const { return m_pos; }

protected:
	// Helper Functions
	size_t newAlloc( size_t _Current, double _Factor ) const;
	const char* buffer() const;
	void setBufferAt(size_t _index, char _value );

	/// NOT *required* in the slightest, but makes debug in VisStudio easier as the visualiser shows 
	/// the actual string, not the internal buffer contents, which is often invalid.
	inline void zeroTerminateIfDebug() const
	{
#ifdef _DEBUG
		if( m_pos == m_alloc ) growCapacity((size_t)1);
		m_buffer[m_pos] = '\0'; 
#endif
	} 

private:// Derived classes must use the public access functions
	// Member Data
	// NOTE: functions like parseInt require '\0' termination of the internal string, and really should be const,
	// therefore we must make the internal data arrays mutable. They are never accessed from the outside 
	// directly anyway, and so to the user, the StringBuilder does indeed remain constant on Parse().
	mutable char* m_buffer; /// The char buffer - a pointer to heap memory. This memory is dynamically reallocated as required.
	mutable size_t m_alloc; /// The amount of memory that is currently allocated on the heap
	size_t m_pos; /// The current position within the internal 'char' buffer	
};

// Stream Operator Overloading
#ifndef SWIG
std::ostream& operator<<(std::ostream &s, const StringBuilder &_Print);
std::istream& operator<<(StringBuilder &_Print, std::istream &s);
#endif

#endif

