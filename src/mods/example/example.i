%include "std_string.i"
%include "exception.i"

#define PD_API

// Propagate Exceptions thrown in C++ to Python
%exception{
		try {
		$action
		} catch ( ExceptionBase ){
			SWIG_exception(SWIG_RuntimeError, "PD Exception occured" );
		}
}


%module example
%{
// Include the PD/MMLIB headers
#include "mmlib.h"

// Include the module headers
#include "example.h"
%}

// This import allows type information to be shared between the module and 
// the pd code. 
%import ../../pd/pd.i

%include example.h


