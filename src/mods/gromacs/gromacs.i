%include "std_string.i"
%include "exception.i"
#define PD_API
// Propagete Exceptions thrown in C++ to Python
%exception{
		try {
		$action
		} catch ( ExceptionBase ){
			SWIG_exception(SWIG_RuntimeError, "PD Exception occured in Gromacs module" );
		}
}


%module gromacs
%{
// Include the PD/MMLIB headers
#include "mmlib.h"

// Include the module headers
#include "gromacs.h"
%}


// This import allows type information to be shared between the module and 
// the pd code. 
%import ../../pd/pd.i

// Include the module header files for parsing
%include gromacs.h


