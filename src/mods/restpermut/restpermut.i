%include "std_string.i"
%include "exception.i"

#define PD_API
// Propagate Exceptions thrown in C++ to Python
%exception{
		try {
		$action
		} catch ( ExceptionBase ){
			SWIG_exception(SWIG_RuntimeError, "pd Exception occured" );
		}
}


%module restpermut 
%{

// MMlib.h includes all the headers that mmlib provides
#include "mmlib.h"

// Module specific headers to be included in the wrap file
#include "restpermut.h"
%}

// This import allows type information to be shared between the module and 
// the pd code. 
%import ../../pd/pd.i


// Module specific headers to be SWIG-parsed
%include restpermut.h


