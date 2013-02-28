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

#ifndef __EXCEPTION_H
#define __EXCEPTION_H

#include <exception>
#include <string>

// Not for end-user usage:
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define FILE_POS "\nCode file: '" __FILE__ "' Line: '" TOSTRING(__LINE__) "'"

// We can automatically add the '__FILE__' and '__LINE__'
// preprocessor macros to a given excpetion with a set of defines
// In-code usage:
// THROW(AnyYouLikeException,"Needs Implementation")
// THROW2(AnyYouLikeException,"Needs Implementation",arg1,arg2)
#define THROW(Exception_Type,Message) throw Exception_Type(std::string(Message)+std::string(FILE_POS))
#define THROW1(Exception_Type,Message,a) throw Exception_Type(std::string(Message)+std::string(FILE_POS),(a))
#define THROW2(Exception_Type,Message,a,b) throw Exception_Type(std::string(Message)+std::string(FILE_POS),(a),(b))
#define THROW3(Exception_Type,Message,a,b,c) throw Exception_Type(std::string(Message)+std::string(FILE_POS),(a),(b),(c))

// Assertions, same as THROW, but has a built in condition. Exception is thrown if the condtion evaluates to 'false'.
// These should NOT be used in inner loops due to performance impact, only in setup-code
#define ASSERT(CONDITION,Exception_Type,Message) if(!(CONDITION)) throw Exception_Type(std::string(Message)+std::string(FILE_POS))
#define ASSERT1(CONDITION,Exception_Type,Message,a) if(!(CONDITION)) throw Exception_Type(std::string(Message)+std::string(FILE_POS),(a))
#define ASSERT2(CONDITION,Exception_Type,Message,a,b) if(!(CONDITION)) throw Exception_Type(std::string(Message)+std::string(FILE_POS),(a),(b))
#define ASSERT3(CONDITION,Exception_Type,Message,a,b,c) if(!(CONDITION)) throw Exception_Type(std::string(Message)+std::string(FILE_POS),(a),(b),(c))

// Debug-only assertions - debug mode error checking with no runtime performance hit
// the exceptions just compile away in release mode ...
// It is fine to use these in inner loops to make sure that any assumption made are valid.
#ifdef _DEBUG
	#define D_ASSERT(CONDITION,Exception_Type,Message) if(!(CONDITION)) throw Exception_Type(std::string(Message)+std::string(FILE_POS))
	#define D_ASSERT1(CONDITION,Exception_Type,Message,a) if(!(CONDITION)) throw Exception_Type(std::string(Message)+std::string(FILE_POS),(a))
	#define D_ASSERT2(CONDITION,Exception_Type,Message,a,b) if(!(CONDITION)) throw Exception_Type(std::string(Message)+std::string(FILE_POS),(a),(b))
	#define D_ASSERT3(CONDITION,Exception_Type,Message,a,b,c) if(!(CONDITION)) throw Exception_Type(std::string(Message)+std::string(FILE_POS),(a),(b),(c))
#else
	#define D_ASSERT(CONDITION,Exception_Type,Message)
	#define D_ASSERT1(CONDITION,Exception_Type,Message,a)
	#define D_ASSERT2(CONDITION,Exception_Type,Message,a,b)
	#define D_ASSERT3(CONDITION,Exception_Type,Message,a,b,c)
#endif

/// \brief Base class PD_API for all pd exceptions...
class PD_API ExceptionBase: public std::exception
{
public:
	ExceptionBase(const std::string &_message = "");
	virtual ~ExceptionBase() throw();
	void Details();
protected:
	std::string message;
};

/// \brief The code for this action is not yet implemented...
class PD_API NotImplementedException: public ExceptionBase
{
public:
	NotImplementedException(const std::string &_message = "");
	virtual ~NotImplementedException() throw();
};

/// \brief There is not enough free memory to allow this allocation...
class PD_API OutOfMemoryException: public ExceptionBase
{
public:
	OutOfMemoryException(const std::string &_message = "");
	virtual ~OutOfMemoryException() throw();
};

/// \brief The code itself is procedurally is broken, i.e. an
/// assumption just fell over and was caught, debugging required!!
class PD_API CodeException: public ExceptionBase
{
public:
	CodeException(const std::string &message );
	virtual ~CodeException() throw();
};


/// \brief Thrown when a function argument is invalid or out of range
class PD_API ArgumentException: public ExceptionBase
{
public:
	ArgumentException(const std::string &_message = "");
	virtual ~ArgumentException() throw();
};

/// \brief Thrown when a function argument is NULL
class PD_API NullArgumentException: public ArgumentException
{
public:
	NullArgumentException(const std::string &_message = "");
	virtual ~NullArgumentException() throw();
};

/// \brief Thrown when in internal pointer is unexpectedly NULL
class PD_API NullInternalException: public ArgumentException
{
public:
	NullInternalException(const std::string &_message = "");
	virtual ~NullInternalException() throw();
};

/// \brief Thrown when a value is not within a given
/// set of bounds / restrained values
class PD_API OutOfRangeException: public ExceptionBase
{
public:
	OutOfRangeException(const std::string &_message = "");
	virtual ~OutOfRangeException() throw();
};

// \brief Helper functions for OutOfRangeException
// Note: printing directly to cerr rather than handing message to
// the exception. maybe it'd be better to make a stringstream & convert the message to
// a string and then pass it ?
/// Non-debug Assert that throws OutOfRangeException if out of range
template <class T>
void argcheck_range(const std::string &name,const T value,const T min,const T max){
	if((value<min)||(value>max)){
		std::cerr << "Parameter " << name << " must be within " << min << " and " << max << std::endl;
		std::cerr.flush();
		throw OutOfRangeException();
	}
}

/// \brief Non-debug Assert that throws OutOfRangeException if value is greater than min
template <class T>
void argcheck_gt(const std::string &name,const T value,const T min){
	if((value<min)){
		std::cerr << "Parameter " << name << " must be greater than " << min << std::endl;
		std::cerr.flush();
		throw OutOfRangeException();
	}
}

/// \brief Non-debug Assert that throws OutOfRangeException if value is less than max
template <class T>
void argcheck_lt(const std::string &name,const T value,const T max){
	if((value>max)){
		std::cerr << "Parameter " << name << " must be less than " << max << std::endl;
		std::cerr.flush();
		throw OutOfRangeException();
	}
}

/// \brief Non-debug Assert that throws OutOfRangeException if value1 is equal to value2
template <class T>
void argcheck_eq(const std::string &name1,const T value1,
								 const std::string &name2,const T value2){
	if((value1!=value2)){
		std::cerr << "Parameter " << name1 << " (" << value1 << ") " << " must be equal to" << name2 << " (" << value2 << ") " << std::endl;
		std::cerr.flush();
		throw OutOfRangeException();
	}
}

/// Non-debug Assert that throws OutOfRangeException if value1 is equal to value2
template <class T>
void argcheck_ne(const std::string &name1,const T value1,
								 const std::string &name2,const T value2){
	if((value1==value2)){
		std::cerr << "Parameter " << name1 << " (" << value1 << ") " << " must not be equal to" << name2 << " (" << value2 << ") " << std::endl;
		std::cerr.flush();
		throw OutOfRangeException();
	}
}


/// \brief ProcedureExceptions are used where a function fails due
/// to invalid data that is not an argument and there is no
/// code error, i.e. some internal data is invalid
class PD_API ProcedureException: public ExceptionBase
{
public:
	ProcedureException(const std::string &_message = "");
	virtual ~ProcedureException() throw();
};

/// \brief Thrown when a file cannot be found or when reading/writing is blocked
class PD_API IOException: public ExceptionBase
{
public:
	IOException(const std::string &_message = "");
	virtual ~IOException() throw();
};


/// \brief Thrown when a parse is requested on an invalid string / input
/// or some syntax error occured.
class PD_API ParseException: public ExceptionBase
{
public:
	ParseException(const std::string &_message = "");
	ParseException(const std::string &_message, 
								 const std::string &_filename,
								 const int          _linenumber = -1);
	virtual ~ParseException() throw();
};



void	PrintFileLineError(const std::string &_message, 
								 				 const std::string &_filename,
								 				 const int  _linenumber = -1);


void	PrintFileLineWarning(const std::string &_message, 
								 				 const std::string &_filename,
								 				 const int  _linenumber = -1);


class PD_API ParseErrorLogger 
{
public:
	ParseErrorLogger( const std::string &_Filename ):m_Filename(_Filename){
		m_LineNumber   = 0;	
		m_WarningCount = 0;
		m_ErrorCount   = 0;
	}

	void setLineNumber( int _LineNumber ) { m_LineNumber = _LineNumber; }
	void incLineNumber() { m_LineNumber++; };
	int getLineNumber() const { return m_LineNumber; };
	
	void setFileName( const std::string &_Filename) { m_Filename = _Filename; }

	void logWarning(  const std::string& _message ){
		PrintFileLineWarning( _message, m_Filename, m_LineNumber );
		m_WarningCount ++;
	}
	void logError(  const std::string& _message ){
		PrintFileLineError( _message, m_Filename, m_LineNumber );
		m_ErrorCount ++;
	}

	size_t getErrorCount(){ return m_ErrorCount; }
private:
	std::string m_Filename;
	int         m_LineNumber;

	size_t m_WarningCount;
	size_t m_ErrorCount;

};




/// \brief Thrown when a class is not properly initialised
/// before its 'internal bits' are used
class PD_API UninitialisedException: public ExceptionBase
{
public:
	UninitialisedException(const std::string &_message = "");
	virtual ~UninitialisedException() throw();
};


void testpdException();

#endif

