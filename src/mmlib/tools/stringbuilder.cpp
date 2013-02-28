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
#include "safeformat.h" // Get the LokiCore for formatted printing
#include "stringbuilder.h"

StringBuilder::StringBuilder( size_t _InitialCapacity )
{
	m_pos = 0;
	m_alloc = _InitialCapacity;
	m_buffer = new char[m_alloc];
	zeroTerminateIfDebug();
}

StringBuilder::StringBuilder( char _Fill, size_t _count )
{
	m_pos = _count;
	if( _count > SB_DEFAULT_CAPACITY )
	{
		m_alloc = newAlloc(_count,1.2);
	}
	else
	{
		m_alloc = SB_DEFAULT_CAPACITY;
	}
	m_buffer = new char[m_alloc];
	memset(m_buffer,_Fill,_count);
	zeroTerminateIfDebug();
}

StringBuilder::StringBuilder( const std::string &_Clone )
{
	size_t _count = _Clone.size();
	m_pos = _count;
	if( _count > SB_DEFAULT_CAPACITY )
	{
		m_alloc = newAlloc(_count,1.2);
	}
	else
	{
		m_alloc = SB_DEFAULT_CAPACITY;
	}
	m_buffer = new char[m_alloc];
	memcpy(m_buffer,_Clone.c_str(),_count);
	zeroTerminateIfDebug();
}

StringBuilder::StringBuilder( const StringBuilder &_Clone )
{
	m_pos = _Clone.m_pos;
	m_alloc = _Clone.m_alloc;
	m_buffer = new char[m_alloc];
	memcpy(m_buffer,_Clone.m_buffer,m_pos);
	zeroTerminateIfDebug();
}

StringBuilder::~StringBuilder()
{
	delete[] m_buffer;
}

StringBuilder& StringBuilder::operator=(const StringBuilder &_Clone)
{
	setTo(_Clone);
	return *this;
}

StringBuilder& StringBuilder::operator=(const std::string &_Clone)
{
	setTo( _Clone );
	return *this;
}

const char& StringBuilder::operator[](size_t i)const
{
	if( i >= m_pos )
	{
		THROW(OutOfRangeException,"StringBuilder: Index is outside the current array bounds!");
	}
	return m_buffer[i];
}

void StringBuilder::clear()
{
	m_pos = 0;
	zeroTerminateIfDebug();
}

size_t StringBuilder::newAlloc( size_t _Current, double _Factor ) const
{
	return (size_t)((double)_Current * _Factor);
}

const char* StringBuilder::buffer() const
{
	return m_buffer;
}

void StringBuilder::setBufferAt(size_t _index, char _value )
{
	if( _index >= m_alloc )
	{
		growCapacity( m_alloc - _index );
	}
	m_buffer[_index] = _value;
}

void StringBuilder::growCapacity( double _ReallocFactor ) const
{
	char* oldArray = m_buffer;
	m_alloc = newAlloc(m_alloc, _ReallocFactor);
	m_buffer = new char[m_alloc];
	if( m_pos > 0 )
		memcpy(m_buffer,oldArray,m_pos);
	delete[] oldArray;
	zeroTerminateIfDebug();
}

void StringBuilder::growCapacity( size_t _capacityIncrease ) const
{
	char* oldArray = m_buffer;
	m_alloc += _capacityIncrease;
	m_buffer = new char[m_alloc];
	if( m_pos > 0 )
		memcpy(m_buffer,oldArray,m_pos);
	delete[] oldArray;
	zeroTerminateIfDebug();
}

void StringBuilder::replace( size_t _index, char _withChar )
{
	if( _index >= m_pos ) throw OutOfRangeException("replace() call is outside of current array bounds!");
	m_buffer[_index] = _withChar;
}

void StringBuilder::replace( size_t _index, char _withChar, size_t _count )
{
	if( (_index + _count) > m_pos ) throw OutOfRangeException("replace() call is outside of current array bounds!");
	memset(&m_buffer[_index],_withChar,_count);
}

void StringBuilder::replace( size_t _index, const char* _string, size_t _length )
{
	if( (_index + _length) > m_pos ) throw OutOfRangeException("replace() call is outside of current array bounds!");
	memcpy(&m_buffer[_index],_string,_length);
}

void StringBuilder::replace( size_t _index, const std::string &_string )
{
	size_t length = _string.length();
	if( (_index + length) > m_pos ) throw OutOfRangeException("replace() call is outside of current array bounds!");
	memcpy(&m_buffer[_index],_string.c_str(),length);
}

void StringBuilder::replace( size_t _index, const StringBuilder &_string )
{
	if( (_index + _string.m_pos) > m_pos ) throw OutOfRangeException("replace() call is outside of current array bounds!");
	memcpy(&m_buffer[_index],_string.m_buffer,_string.m_pos);
}

void StringBuilder::setTo( char _c )
{
	clear();
	append(_c);
}

void StringBuilder::setTo( char _c, size_t _count )
{
	clear();
	append(_c,_count);
}

void StringBuilder::setTo( const char* const _string, size_t _length )
{
	clear();
	append(_string,_length);
}

void StringBuilder::setTo( const std::string &_string )
{
	clear();
	append(_string);
}

void StringBuilder::setTo( const std::string &_string, size_t _index, size_t _length )
{
	clear();
	append(_string,_index,_length);
}

void StringBuilder::setTo( const StringBuilder &_string )
{
	// This function is used by the copy constructor
	if( m_alloc < _string.m_pos ) // only need enough room to host the clones current string, not the whole allocation
	{
		m_alloc = _string.m_alloc;
		delete[] m_buffer;
		m_buffer = new char[m_alloc];
	}
	m_pos = _string.m_pos;
	memcpy(m_buffer,_string.m_buffer,m_pos);
	zeroTerminateIfDebug();
}

void StringBuilder::setTo( const StringBuilder &_string, size_t _index, size_t _length )
{
	clear();
	append(_string,_index,_length);
}

LokiCore::PrintfState<StringBuilder&, char> StringBuilder::setFormat(const std::string& format) 
{
	clear();
	return LokiCore::PrintfState<StringBuilder&, char>(*this, format.c_str());
}	

LokiCore::PrintfState<StringBuilder&, char> StringBuilder::appendFormat(const std::string& format) 
{
	return LokiCore::PrintfState<StringBuilder&, char>(*this, format.c_str());
}

void StringBuilder::append( char _c )
{
	if( m_pos == m_alloc ) growCapacity();
	m_buffer[m_pos++] = _c;
	zeroTerminateIfDebug();
}

void StringBuilder::append( char _c, size_t _count )
{
	while( (m_pos + _count) > m_alloc ) growCapacity();
	memset(&m_buffer[m_pos],_c,_count);
	m_pos += _count;
	zeroTerminateIfDebug();
}

void StringBuilder::append( const char* const _string, size_t _length )
{
	while( (m_pos + _length) > m_alloc ) growCapacity();
	memcpy(&m_buffer[m_pos],_string,_length);
	m_pos += _length;
	zeroTerminateIfDebug();
}

void StringBuilder::append( const std::string &_string )
{
	size_t length = _string.size();
	while( (m_pos + length) > m_alloc ) growCapacity();
	memcpy(&m_buffer[m_pos],_string.c_str(),length);
	m_pos += length;
	zeroTerminateIfDebug();
}

void StringBuilder::append( const std::string &_string, size_t _index, size_t _length )
{
	size_t strLen = _string.length();
	if( strLen < (_index+_length) ) throw OutOfRangeException("Range given is outside of the range of the given string!");
	while( (m_pos + _length) > m_alloc ) growCapacity();
	memcpy(&m_buffer[m_pos],&_string.c_str()[_index],_length);
	m_pos += _length;
	zeroTerminateIfDebug();
}

void StringBuilder::append( const StringBuilder &_string )
{
	while( (m_pos + _string.m_pos) > m_alloc ) growCapacity();
	memcpy(&m_buffer[m_pos],_string.m_buffer,_string.m_pos);
	m_pos += _string.m_pos;
	zeroTerminateIfDebug();
}

void StringBuilder::append( const StringBuilder &_string, size_t _index, size_t _length )
{
	size_t strLen = _string.size();
	if( strLen < (_index+_length) ) throw OutOfRangeException("Range given is outside of the range of the given string!");
	while( (m_pos + _length) > m_alloc ) growCapacity();
	memcpy(&m_buffer[m_pos],&_string.m_buffer[_index],_length);
	m_pos += _length;
	zeroTerminateIfDebug();
}

void StringBuilder::append( int a, char *formatstring )
{
	// System defines: DBL_MAX 1.7976931348623158e+308
	const int ensureCapacity = 23; // make sure we can store 23 characters in our buffer
	while( (m_pos + ensureCapacity) > m_alloc ) growCapacity();
	int addedCharacters = sprintf(&m_buffer[m_pos],formatstring,a);
	m_pos += addedCharacters;
	zeroTerminateIfDebug();
}

void StringBuilder::append( double a, char *formatstring )
{
	// System defines: INT_MAX 2147483647
	const int ensureCapacity = 10; // make sure we can store 10 characters in our buffer
	while( (m_pos + ensureCapacity) > m_alloc ) growCapacity();
	int addedCharacters = sprintf(&m_buffer[m_pos],formatstring,a);
	m_pos += addedCharacters;
	zeroTerminateIfDebug();
}

void StringBuilder::append( long a, char *formatstring )
{
	// System defines: LLONG_MAX 9223372036854775807
	const int ensureCapacity = 20; // make sure we can store 20 characters in our buffer
	while( (m_pos + ensureCapacity) > m_alloc ) growCapacity();
	int addedCharacters = sprintf(&m_buffer[m_pos],formatstring,a);
	m_pos += addedCharacters;
	zeroTerminateIfDebug();
}

void StringBuilder::endl()
{
	append( TOKEN_LINEFEED );
}

void StringBuilder::appendingReplace( size_t _index, char _withChar, size_t _count )
{
	if( _index > m_pos ) throw OutOfRangeException("appendingReplace() call index is outside the current array bounds!");
	size_t endIndex = _index + _count;
	while( endIndex > m_alloc ) growCapacity();
	if( endIndex > m_pos ) m_pos = endIndex;
	memset(&m_buffer[_index],_withChar,_count);
	zeroTerminateIfDebug();
}

void StringBuilder::appendingReplace( size_t _index, const char* _string, size_t _length )
{
	if( _index > m_pos ) throw OutOfRangeException("appendingReplace() call index is outside the current array bounds!");
	size_t endIndex = _index + _length;
	while( endIndex > m_alloc ) growCapacity();
	if( endIndex > m_pos ) m_pos = endIndex;
	memcpy(&m_buffer[_index],_string,_length);
	zeroTerminateIfDebug();
}

void StringBuilder::appendingReplace( size_t _index, const std::string &_string )
{
	if( _index > m_pos ) throw OutOfRangeException("appendingReplace() call index is outside the current array bounds!");
	size_t length = _string.length();
	size_t endIndex = _index + length;
	while( endIndex > m_alloc ) growCapacity();
	if( endIndex > m_pos ) m_pos = endIndex;
	memcpy(&m_buffer[_index],_string.c_str(),length);
	zeroTerminateIfDebug();
}

void StringBuilder::appendingReplace( size_t _index, const StringBuilder &_string )
{
	if( _index > m_pos ) throw OutOfRangeException("appendingReplace() call index is outside the current array bounds!");
	size_t endIndex = _index + _string.m_pos;
	while( endIndex > m_alloc ) growCapacity();
	if( endIndex > m_pos ) m_pos = endIndex;
	memcpy(&m_buffer[_index],_string.m_buffer,_string.m_pos);
	zeroTerminateIfDebug();
}

void StringBuilder::insert( size_t _index, char _Char )
{
	if( m_pos == m_alloc ) growCapacity();
	for( int i = (int)m_pos; i > _index; i-- ) // This must be an 'int' if _index == 0!!
	{
		m_buffer[i] = m_buffer[i-1];
	}
	m_buffer[_index] = _Char;
	m_pos++;
	zeroTerminateIfDebug();
}

void StringBuilder::insert( size_t _index, char _Char, size_t _count )
{
	while( (m_pos + _count) > m_alloc ) growCapacity();
	size_t j;
	for( int i = (int)m_pos; i > _index; i-- ) // This must be an 'int' if _index == 0!!
	{
		j = i - 1;
		m_buffer[_count+j] = m_buffer[j];
	}
	memset(&m_buffer[_index],_Char,_count);
	m_pos += _count;
	zeroTerminateIfDebug();
}

void StringBuilder::insert( size_t _index, const char* _string, size_t _length )
{
	while( (m_pos + _length) > m_alloc ) growCapacity();
	int j;
	for( int i = (int)m_pos; i > _index; i-- ) // This must be an 'int' if _index == 0!!
	{
		j = i - 1;
		m_buffer[_length+j] = m_buffer[j];
	}
	memcpy(&m_buffer[_index],_string,_length);
	m_pos += _length;
	zeroTerminateIfDebug();
}

void StringBuilder::insert( size_t _index, const std::string &_string )
{
	int j;
	size_t _length = _string.length();
	while( (m_pos + _length) > m_alloc ) growCapacity();
	for( int i = (int)m_pos; i > _index; i-- ) // This must be an 'int' if _index == 0!!
	{
		j = i - 1;
		m_buffer[_length+j] = m_buffer[j];
	}
	memcpy(&m_buffer[_index],_string.c_str(),_length);
	m_pos += _length;
	zeroTerminateIfDebug();
}

void StringBuilder::insert( size_t _index, const StringBuilder &_string )
{
	int j;
	while( (m_pos + _string.m_pos) > m_alloc ) growCapacity();
	for( int i = (int)m_pos; i > _index; i-- ) // This must be an 'int' if _index == 0!!
	{
		j = i - 1;
		m_buffer[_string.m_pos+j] = m_buffer[j];
	}
	memcpy(&m_buffer[_index],_string.m_buffer,_string.m_pos);
	m_pos += _string.m_pos;
	zeroTerminateIfDebug();
}

void StringBuilder::erase( size_t _index )
{
	if( m_pos == 0 || _index >= m_pos ) throw OutOfRangeException("erase() call is outside of current array bounds!");
	m_pos--;
	for( size_t i = _index; i < m_pos; i++ )
	{
		m_buffer[i] = m_buffer[i+1];
	}
	zeroTerminateIfDebug();
}

void StringBuilder::erase( size_t _index, size_t _count )
{
	if( (_index + _count) > m_pos ) throw OutOfRangeException("erase() call is outside of current array bounds!");
	m_pos -= _count;
	for( size_t i = _index; i < m_pos; i++ )
	{
		m_buffer[i] = m_buffer[i+_count];
	}
	zeroTerminateIfDebug();
}

size_t StringBuilder::FirstOf( char _delimiter ) const
{
	for( size_t i = 0; i < m_pos; i++ )
	{
		if( m_buffer[i] == _delimiter ) return i;
	}
	return SIZE_T_FAIL;
}

size_t StringBuilder::FirstNotOf( char _delimiter ) const
{
	for( size_t i = 0; i < m_pos; i++ )
	{
		if( m_buffer[i] != _delimiter ) return i;
	}
	return SIZE_T_FAIL;
}

size_t StringBuilder::LastOf( char _delimiter ) const
{
	for( int i = (((int)m_pos)-1); i >= 0; i-- ) // This must be an 'int' i can be zero, otherwise an infinate loop is created (size_t cant be <0, int can)
	{
		if( m_buffer[i] == _delimiter ) 
			return i;
	}
	return SIZE_T_FAIL;
}

size_t StringBuilder::LastNotOf( char _delimiter ) const
{
	if( m_pos == 0 ) return SIZE_T_FAIL;
	for( int i = (((int)m_pos)-1); i >= 0; i-- ) // This must be an 'int' i can be zero, otherwise an infinate loop is created (size_t cant be <0, int can)
	{
		if( m_buffer[i] != _delimiter ) 
			return i;
	}
	return SIZE_T_FAIL;
}

inline bool isDelimiter( char c, const std::string& _delimiters )
{
	for( int j = 0; j < _delimiters.length(); j++ )
	{
		if( c == _delimiters[j] )
		{
			return true;
		}
	}
	return false;
}

size_t StringBuilder::FirstOf( const std::string &_delimiters ) const
{
	size_t index, indexBest = SIZE_T_FAIL;
	for( size_t i = 0; i < _delimiters.length(); i++ )
	{
		index = FirstOf(_delimiters[i]);
		if( index < indexBest ) 
			indexBest = index;
	}
	return indexBest;
}

size_t StringBuilder::FirstNotOf( const std::string &_delimiters ) const
{
	for( size_t i = 0; i < m_pos; i++ )
	{
		if(!isDelimiter(m_buffer[i],_delimiters))
			return i;
	}
	return SIZE_T_FAIL;
}

size_t StringBuilder::LastOf( const std::string &_delimiters ) const
{
	size_t index, indexBest = SIZE_T_FAIL;
	for( size_t i = 0; i < _delimiters.length(); i++ )
	{
		index = LastOf(_delimiters[i]);
		if( index != SIZE_T_FAIL && ((index > indexBest) || (indexBest==SIZE_T_FAIL)) ) 
		{
			indexBest = index;
		}
	}
	return indexBest;
}

size_t StringBuilder::LastNotOf( const std::string &_delimiters ) const
{
	for( int i = (((int)m_pos)-1); i >= 0; i-- ) // This must be an 'int' i can be zero, otherwise an infinate loop is created (size_t cant be <0, int can)
	{
		if(!isDelimiter(m_buffer[i],_delimiters))
			return i;
	}
	return SIZE_T_FAIL;
}

size_t StringBuilder::XthOf( size_t instance, char _delimiter ) const
{
	size_t count = 0;
	for( size_t i = 0; i < m_pos; i++ )
	{
		if( m_buffer[i] == _delimiter ) 
		{
			if( count < instance ) 
			{
				count++;
			}
			else
			{
				return i;
			}
		}
	}
	return SIZE_T_FAIL;
}

size_t StringBuilder::XthNotOf( size_t instance, char _delimiter ) const
{	
	size_t count = 0;
	for( size_t i = 0; i < m_pos; i++ )
	{
		if( m_buffer[i] != _delimiter ) 
		{
			if( count < instance ) 
			{
				count++;
			}
			else
			{
				return i;
			}
		}
	}
	return SIZE_T_FAIL;
}

size_t StringBuilder::XthOf( size_t instance, const std::string &_delimiters ) const
{
	size_t count = 0;
	for( size_t i = 0; i < m_pos; i++ )
	{
		if(isDelimiter(m_buffer[i],_delimiters)) 
		{
			if( count < instance ) 
			{
				count++;
			}
			else
			{
				return i;
			}
		}
	}
	return SIZE_T_FAIL;
}

size_t StringBuilder::XthNotOf( size_t instance, const std::string &_delimiters ) const
{
	size_t count = 0;
	for( size_t i = 0; i < m_pos; i++ )
	{
		if(!isDelimiter(m_buffer[i],_delimiters)) 
		{
			if( count < instance ) 
			{
				count++;
			}
			else
			{
				return i;
			}
		}
	}
	return SIZE_T_FAIL;
}

void StringBuilder::Trim( char _delimiter )
{
	TrimLeft(_delimiter);
	TrimRight(_delimiter);
}

void StringBuilder::TrimLeft( char _delimiter )
{
	size_t index = FirstNotOf(_delimiter);
	if( index != SIZE_T_FAIL ) erase(0,index);
}

void StringBuilder::TrimRight( char _delimiter )
{
	size_t index = LastNotOf(_delimiter);
	if( index != SIZE_T_FAIL ) erase(index,m_pos-index);
}

void StringBuilder::Trim( const std::string &_delimiters )
{
	TrimRight(_delimiters);
	TrimLeft(_delimiters);	
}

void StringBuilder::TrimLeft( const std::string &_delimiters )
{
	size_t index = FirstNotOf(_delimiters);
	if( index != SIZE_T_FAIL ) erase(0,index);
	else clear();
}

void StringBuilder::TrimRight( const std::string &_delimiters )
{
	size_t index = LastNotOf(_delimiters);
	if( index != SIZE_T_FAIL ) erase(index+1,m_pos-index-1);
}

void StringBuilder::PadRight( size_t _toLength, char _withChar )
{
	if( m_pos < _toLength ) append(_withChar,_toLength-m_pos);
}

void StringBuilder::PadLeft( size_t _toLength, char _withChar )
{
	if( m_pos < _toLength ) insert(0,_withChar,_toLength-m_pos);
}

void StringBuilder::TruncateRightTo( size_t _toLength )
{
	if( m_pos > _toLength ) m_pos = _toLength;
	zeroTerminateIfDebug();
}

void StringBuilder::TruncateLeftTo( size_t _toLength )
{
	if( m_pos > _toLength ) erase(0,m_pos-_toLength);
}

void StringBuilder::TruncateRightBy( size_t _removeCount )
{
	if( (int)m_pos - (int)_removeCount < 0 ) 
		throw OutOfRangeException("StringBuilder::TruncateRightBy() is longer than the string!");
	m_pos -= _removeCount;
	zeroTerminateIfDebug();
}

void StringBuilder::TruncateLeftBy( size_t _removeCount )
{
	if( (int)m_pos - (int)_removeCount < 0 ) 
		throw OutOfRangeException("StringBuilder::TruncateLeftBy() is longer than the string!"); 
	erase(0,_removeCount);
}

void StringBuilder::removeAll( char _delimiter )
{
	size_t index;
	while(SIZE_T_FAIL != (index = FirstOf(_delimiter))) erase(index);
}

void StringBuilder::removeAll( const std::string &_delimiters )
{
	size_t index;
	for( size_t i = 0; i < _delimiters.length(); i++ )
	{
		while(SIZE_T_FAIL != (index = FirstOf(_delimiters[i]))) erase(index);
	}
}

bool StringBuilder::compare( char _c, size_t _index, bool _ignoreCase )const
{
	if( _index >= m_pos ) throw OutOfRangeException("StringBuilder::Compare()");
	if( _ignoreCase )
	{
		return isSameCharIgnoringCase(_c,m_buffer[_index]);
	}
	else
	{
		return _c == m_buffer[_index];
	}
}

bool StringBuilder::compare( const char *_stringB, size_t _indexA, size_t _length, size_t _indexB, bool _ignoreCase ) const
{
	if( ((_indexA + _length) > m_pos) || (_indexB+_length) > strlen(_stringB) ) 
		throw OutOfRangeException("StringBuilder::Compare() bad range");
	if( _ignoreCase )
	{
		for( size_t i = 0; i < _length; i++ )
			if( !isSameCharIgnoringCase(_stringB[_indexB+i],m_buffer[_indexA+i]) ) 
				return false;
	}
	else
	{
		for( size_t i = 0; i < _length; i++ )
			if( _stringB[_indexB+i] != m_buffer[_indexA+i] ) 
				return false;
	}
	return true;
}

bool StringBuilder::compare( const std::string &_stringB, size_t _indexA, size_t _length, size_t _indexB, bool _ignoreCase )const
{
	if( ((_indexA+_length) > m_pos) || (_indexB+_length) > _stringB.size() ) 
		throw OutOfRangeException("StringBuilder::Compare() bad range");
	if( _ignoreCase )
	{
		for( size_t i = 0; i < _length; i++ )
			if( !isSameCharIgnoringCase(_stringB[_indexB+i],m_buffer[_indexA+i]) ) 
				return false;
	}
	else
	{
		for( size_t i = 0; i < _length; i++ )
			if( _stringB[_indexB+i] != m_buffer[_indexA+i] ) 
				return false;
	}
	return true;
}

bool StringBuilder::compare( const StringBuilder &_stringB, size_t _indexA, size_t _length, size_t _indexB, bool _ignoreCase )const
{
	if( ((_indexA+_length) > m_pos) || (_indexB+_length) > _stringB.m_pos ) 
		throw OutOfRangeException("StringBuilder::Compare() bad range");
	if( _ignoreCase )
	{
		for( size_t i = 0; i < _length; i++ )
			if( !isSameCharIgnoringCase(_stringB[_indexB+i],m_buffer[_indexA+i]) ) 
				return false;
	}
	else
	{
		for( size_t i = 0; i < _length; i++ )
			if( _stringB[_indexB+i] != m_buffer[_indexA+i] ) 
				return false;
	}
	return true;
}

int StringBuilder::compareValue(const StringBuilder &_stringB, size_t _indexA, size_t _length, size_t _indexB, bool _ignoreCase ) const
{
	if( ((_indexA+_length) > m_pos) || (_indexB+_length) > _stringB.m_pos ) 
		throw OutOfRangeException("StringBuilder::Compare() bad range");
	if( _ignoreCase )
	{
		for( size_t i = 0; i < _length; i++ )
			if( !isSameCharIgnoringCase(_stringB[_indexB+i],m_buffer[_indexA+i]) ) 
				return _stringB[_indexB+i] > m_buffer[_indexA+i] ? -1 : 1;
	}
	else
	{
		for( size_t i = 0; i < _length; i++ )
			if( _stringB[_indexB+i] != m_buffer[_indexA+i] )
				return _stringB[_indexB+i] > m_buffer[_indexA+i] ? -1 : 1;
	}
	return 0;
}

void StringBuilder::toScreen() const
{
	if( m_pos == 0 ) return;
	if( m_pos == m_alloc ) growCapacity((size_t)1); // make sure we have enough room for the following temproary '/0' termination
	m_buffer[m_pos] = '\0'; // *must* zero terminate
	printf("%s",m_buffer);
}

void StringBuilder::toScreen( size_t _index, size_t _length ) const
{
	if( (_index + _length) > m_pos ) throw OutOfRangeException("toString( size_t _index, size_t _length ) call is outside of current array bounds!");
	if( _length == 0 ) return;
	if( _index+_length == m_alloc ) growCapacity((size_t)1); // make sure we have enough room for the following temproary '/0' termination
	char store = m_buffer[_index+_length]; // store the value
	m_buffer[_index+_length] = '\0'; // *must* zero terminate
	printf("%s",&m_buffer[_index]); // perform the print operation on the now terminated string
	m_buffer[_index+_length] = store;
}

std::string StringBuilder::toString() const
{
	if( m_pos == 0 )
	{
		return std::string("");
	}
	else
	{
		return std::string(m_buffer,m_pos);
	}
}

std::string StringBuilder::toString( size_t _index, size_t _length ) const
{
	if( (_index + _length) > m_pos ) throw OutOfRangeException("toString( size_t _index, size_t _length ) call is outside of current array bounds!");
	if( m_pos == 0 || _length == 0 )
	{
		return std::string("");
	}
	else
	{
		return std::string(&m_buffer[_index],_length);
	}
}

std::vector<std::string> StringBuilder::tokenise( const std::string &_delimiters ) const
{
	std::vector<std::string> parts;
	size_t i;
	size_t j = 0;
	size_t x = 0;
	StringBuilder sub;
	while( SIZE_T_FAIL != (i = XthOf(x++,_delimiters) ) )
	{
		sub.setTo(*this,j,i-j);
		j = i;
		sub.Trim(_delimiters);
		if( sub.size() > 0 )
		{
			parts.push_back(sub.toString());
		}
	}
	if( j != m_pos )
	{
		sub.setTo(*this,j,m_pos-j);
		sub.Trim(_delimiters);
		if( sub.size() > 0 )
		{
			parts.push_back(sub.toString());
		}
	}
	return parts;
}

int StringBuilder::parseInt() const
{
	if( m_pos == 0 ) throw ParseException("StringBuilder::parseInt(): Cannot parse a '0' length string!");
	if( m_pos == m_alloc ) growCapacity((size_t)1); // make sure we have enough room for the following temproary '/0' termination
	int lvalue;
	if(sscanf(m_buffer,"%d",&lvalue)!=1) throw ParseException("StringBuilder::parseInt()");
	return lvalue;
}

int StringBuilder::parseInt( size_t _index ) const
{
	if( m_pos == 0 || _index >= m_pos-1 ) throw ParseException("StringBuilder::parseInt(): Cannot parse a '0' length string!");
	if( m_pos == m_alloc ) growCapacity((size_t)1); // make sure we have enough room for the following temproary '/0' termination
	m_buffer[m_pos] = 0; // ensure a '0' terminated string for sscanf() but do not increase the size of the string
	int lvalue;
	if(sscanf(&m_buffer[_index],"%d",&lvalue)!=1) throw ParseException("StringBuilder::parseInt()");
	return lvalue;
}

int StringBuilder::parseInt( size_t _index, size_t _length ) const
{
	size_t tempTerminatorIndex = _index + _length;
	if( _length == 0 || tempTerminatorIndex > m_pos ) throw ParseException("StringBuilder::parseInt(): Cannot parse a '0' length string!");
	if( tempTerminatorIndex+1 > m_alloc ) growCapacity((size_t)(tempTerminatorIndex+1 - m_alloc)); // make sure we have enough room for the following temproary '/0' termination
	char cache = m_buffer[tempTerminatorIndex]; // store what was here before
	m_buffer[tempTerminatorIndex] = 0; // ensure a '0' terminated string for sscanf() but do not increase the size of the string
	int lvalue;
	int result = sscanf(&m_buffer[_index],"%d",&lvalue);
	m_buffer[tempTerminatorIndex] = cache; // put it back!
	if(result!=1) throw ParseException("StringBuilder::parseInt()");
	else return lvalue;
}

double StringBuilder::parseDouble() const
{
	if( m_pos == 0 ) throw ParseException("StringBuilder::parseDouble(): Cannot parse a '0' length string!");
	if( m_pos == m_alloc ) growCapacity((size_t)1); // make sure we have enough room for the following temproary '/0' termination
	m_buffer[m_pos] = 0; // ensure a '0' terminated string for sscanf() but do not increase the size of the string
	double lvalue;
	if(sscanf(m_buffer,"%lf",&lvalue)!=1) throw ParseException("StringBuilder::parseDouble()");
	return lvalue;
}

double StringBuilder::parseDouble( size_t _index ) const
{
	if( m_pos == 0 || _index >= m_pos-1 ) throw ParseException("StringBuilder::parseDouble( size_t _index )");
	if( m_pos == m_alloc ) growCapacity((size_t)1); // make sure we have enough room for the following temproary '/0' termination
	m_buffer[m_pos] = 0; // ensure a '0' terminated string for sscanf() but do not increase the size of the string
	double lvalue;
	if(sscanf(&m_buffer[_index],"%lf",&lvalue)!=1) throw ParseException("StringBuilder::parseDouble()");
	return lvalue;
}

double StringBuilder::parseDouble( size_t _index, size_t _length ) const
{
	size_t tempTerminatorIndex = _index + _length;
	if( _length == 0 || tempTerminatorIndex > m_pos ) throw ParseException("StringBuilder::parseDouble(): Cannot parse a '0' length string!");
	if( tempTerminatorIndex+1 > m_alloc ) growCapacity((size_t)(tempTerminatorIndex+1 - m_alloc)); // make sure we have enough room for the following temproary '/0' termination
	char cache = m_buffer[tempTerminatorIndex]; // store what was here before
	m_buffer[tempTerminatorIndex] = 0; // ensure a '0' terminated string for sscanf() but do not increase the size of the string
	double lvalue;
	int result = sscanf(&m_buffer[_index],"%lf",&lvalue);
	m_buffer[tempTerminatorIndex] = cache; // put it back!
	if(result!=1) throw ParseException("StringBuilder::parseInt()");
	else return lvalue;
}

std::ostream& operator<<(std::ostream &s, const StringBuilder &_Print)
{
 
	if( _Print.m_pos > 0 )
	{
		if( _Print.m_pos == _Print.m_alloc ) _Print.growCapacity((size_t)1); // make sure we have enough room for the following temproary '/0' termination
		_Print.m_buffer[_Print.m_pos] = 0; // zero terminate the current string, but no need to increment m_pos
		s << _Print.m_buffer;
	}
	
	return s;
}

std::istream& operator<<(StringBuilder &_Print, std::istream &s)
{
	_Print.clear(); // Reset the buffer
	while(true)
	{
		if( s.eof() )
			break;

		char c = s.get();

		if( c == '\n' || c == '\0' ) 
			break;
		else if( c == '\r' ) 
			continue;

		if( _Print.m_alloc - _Print.m_pos == 0 ) 
			_Print.growCapacity(1.2);

		_Print.m_buffer[_Print.m_pos++] = c;

		_Print.zeroTerminateIfDebug();
	}
	return s;
}

