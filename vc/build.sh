#!/bin/sh
## Build from cammand line (sh script)

## More info here http://msdn2.microsoft.com/en-us/library/hw9dzw3c.aspx
## Paths may need to be adjusted

swig -I"../src" -c++ -fvirtual -python "../src/pd/pd.i" 
/c/Program\ Files/Microsoft\ Visual\ Studio\ 8/VC/vcpackages/vcbuild.exe pd.sln 'Release|Win32'
