%include "std_string.i"
%include "std_vector.i"

namespace std 
{
   %template(IntVector) vector<int>;
   %template(FloatVector) vector<float>;
   %template(DoubleVector) vector<double>;
}

%include "exception.i"


// Propagate Exceptions thrown in C++ to Python
%exception{
		try {$action} 
		catch ( ExceptionBase ){SWIG_exception(SWIG_RuntimeError, "PD Exception" );}
}

%feature("compactdefaultargs");

// Prevent the user dynamically add properties to wrapped objects
// as that is extremely error prone.
%pythonnondynamic;

// Set special properties for certain functions to keep the python object
// alive
%pythonprepend ObjectContainer::add %{if not hasattr(args[0],"store"): _swig_setattr(args[0], args[0] , "store", []);
        args[0].store.append(args[1])    %}
%pythonprepend ObjectContainer::lend %{print "lend() is DEPRECATED! Use add() instead!";if not hasattr(args[0],"store"): _swig_setattr(args[0], ObjectContainer_ForcefieldBase, "store", []);args[0].store.append(args[1])    %}

//  - Add code before the C++ function is called
//  %pythonappend    - Add code after the C++ function is called

// Define this just so SWIG doesnt complain with silly warnings
#define PD_API

%define DOCSTRING
"PD - Molecular Mechanics Library for Python ver beta0.01 
 (c) 2003-2006 M.Tyka, J.Rea"
%enddef

%module(docstring=DOCSTRING) pd
%{
#include "mmlib/mmlib.h"
#include "accessory.h"
%}

// Now include all the headers that SWIG will auto-analyse
// and generate python interfaces for.
// You MUST us the path relative to src, not relative to src/pd ! otherwise
// the makefile wont work properly

%include "mmlib/verbosity.h"
%include "mmlib/interface.h"
%include "mmlib/object.h"
%include "mmlib/primitives.h"

%include "mmlib/maths/maths.fwd.h"
%include "mmlib/maths/maths.h"
%include "mmlib/maths/maths_matrix3x3.h"
%include "mmlib/maths/maths_vector.h"
%include "mmlib/maths/graphtheory.h"
%include "mmlib/maths/fastrandom.h"

%template(dvector) Maths::Tvector<double>;
%template(fvector) Maths::Tvector<float>;

%include "mmlib/maths/histogram.h"
%include "mmlib/tools/statclock.h"
%include "mmlib/tools/streamwriter.h"
%include "mmlib/tools/stringbuilder.h"

%include "mmlib/system/fundamentals.h"
%include "mmlib/system/molecule.h"
%include "mmlib/system/system.h"
%include "mmlib/system/genpolymer.h"
%include "mmlib/pickers/pickbase.h"
%include "mmlib/pickers/basicpickers.h"
%include "mmlib/pickers/pickfromfile.h"
%include "mmlib/system/rebuilder.h"
%include "mmlib/protocols/temperature.h"

%include "mmlib/workspace/componentbase.h"
%include "mmlib/forcefields/forcefield.h"
%include "mmlib/monitors/monitorbase.h"
%include "mmlib/protocols/protocolbase.h"
%include "mmlib/manipulators/movebase.h"
%include "mmlib/filters/filterbase.h"

%include "mmlib/library/backbonetorsions.h"
%include "mmlib/library/angleset.h"
%include "mmlib/library/mapper.h"
%include "mmlib/library/nameconv.h"
%include "mmlib/library/rotamerlib.h"
%include "mmlib/library/rotamer_dunbrack.h"
%include "mmlib/library/rotamer_shetty.h"

%include "mmlib/sequence/sequence.h"
%include "mmlib/sequence/alignment.h"

%include "mmlib/fileio/outtra.h"
%include "mmlib/fileio/intra.h"
%include "mmlib/fileio/psfdcd.h"
%include "mmlib/fileio/infile.h"
%include "mmlib/fileio/pdb.h"
%include "mmlib/fileio/tratypes.h"
%include "mmlib/fileio/trablocks.h"
%include "mmlib/fileio/tra.h"

%include "mmlib/workspace/space.h"
%include "mmlib/workspace/bondorder.h"      
%include "mmlib/workspace/neighbourlist.h"
%include "mmlib/workspace/rotbond.h"
%include "mmlib/workspace/workspace.h"
%include "mmlib/workspace/snapshot.h"
%include "mmlib/workspace/pospointer.h"

%include "mmlib/sequence/sequence.h"
%include "mmlib/sequence/alignment.h"

%include "mmlib/library/residues.h"
%include "mmlib/library/nameconv.h"

%include "mmlib/forcefields/ffparam.h"
%include "mmlib/forcefields/nonbonded.h"
%include "mmlib/forcefields/nonbonded_ti.h"
%include "mmlib/forcefields/nonbonded_ti_linear.h"
%include "mmlib/forcefields/ffbonded.h"
%include "mmlib/forcefields/gbff.h"
%include "mmlib/forcefields/lcpo.h"

%include "mmlib/forcefields/restraintbase.h"
%include "mmlib/forcefields/restraint_positional.h"
%include "mmlib/forcefields/restraint_internal.h"
%include "mmlib/forcefields/restraint_torsional.h"
%include "mmlib/forcefields/restraint_atomdist.h"
%include "mmlib/forcefields/restraint_native_contact.h"

%include "mmlib/forcefields/ffcustom.h"
%include "mmlib/forcefields/ffsoftvdw.h"
%include "mmlib/forcefields/breakablebonded.h"

%include "mmlib/manipulators/rotamer_applicatorbase.h"
%include "mmlib/manipulators/basicmoves.h"
%include "mmlib/manipulators/rotamer_scwrl.h"

%include "mmlib/filters/basicfilters.h"

%include "mmlib/protocols/energy.h"
%include "mmlib/protocols/minimise.h"
%include "mmlib/protocols/md.h"
%include "mmlib/protocols/rerun.h"
%include "mmlib/protocols/remd.h"
%include "mmlib/protocols/torsionalminimisation.h"
%include "mmlib/protocols/montecarlo.h"
%include "mmlib/protocols/scpack.h"
%include "mmlib/protocols/dualffminimiser.h"

%include "mmlib/monitors/basicmonitors.h"
%include "mmlib/monitors/umbrella.h"
%include "mmlib/monitors/monexp.h"
%include "mmlib/protocols/nmode.h"

%include "pd/accessory.h"

%pythoncode %{
import sys

def after(optionname, default):
  result = default
  next = 0
  for arg in sys.argv:
    if next == 1:
      result = arg
      break
    if arg == optionname:
      next = 1
  return result

def isarg(optionname):
  for arg in sys.argv:
    if arg == optionname:
      return True
  return False
  
def flushed_info(line):
	info(line)
	sys.stdout.flush()
	
def usage(Object): Object.usage()	

import sys

def pd_excepthook(etype, value, tb):
  if str(value) != "PD Exception": cprint( "\n-----! Python Exception !-------------------------")
  if hasattr(etype,'__name__'): stype = etype.__name__
  else:                       stype = etype
  cprint("%s %s %s %s"%("Fatal", str(stype),"occured:", str(value)))
  n = 0
  stack_filepos = []
  stack_funcname = []
  stack_line = []
  max_filepos = 0
  max_funcname = 0
  while tb is not None:
      co = tb.tb_frame.f_code
      filename = co.co_filename

      try:
        fp = open(co.co_filename, 'rU')
        lines = fp.readlines()
        line = lines[tb.tb_lineno-1]
        fp.close()
        if line: line = line.strip()
        else: line = None
      except: line = "cant open file"
      if filename[-5:] != 'pd.py':       ## ignore internal call stack
        filepos =  '%s:%d'%(co.co_filename, tb.tb_lineno)
        if len(filepos) > max_filepos: max_filepos = len(filepos)
        if len(co.co_name) > max_funcname: max_funcname = len(co.co_name)
        stack_filepos.append(filepos)
        stack_funcname.append(co.co_name)
        stack_line.append(line)
      tb = tb.tb_next
      n=n+1
  n = len(stack_filepos)-1
  fstr1 = '%%%ds'%(max_filepos+0)
  fstr2 = '%%%ds'%(max_funcname+0)
  prefix = "         at:"
  while n>=0:
    if stack_funcname[n] != "?":
      cprint( "%s %s %s %s %s %s"%( prefix, fstr1%stack_filepos[n],"in", fstr2%stack_funcname[n], "(...) :", stack_line[n] ))
    else:
      cprint( "%s %s %s %s %s %s"%( prefix, fstr1%stack_filepos[n],"  ", fstr2%" ",               "      :", stack_line[n] ))
    prefix = "called from:"
    n=n-1
 
## replace the standard exception handler with our own pd-aware version
sys.excepthook = pd_excepthook

%}

