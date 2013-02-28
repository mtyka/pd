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

/*_________________________________________________________________
 |
 | xdrf.h - include file for C routines that want to use the 
 |	    functions below.
*/

int xdropen(XDR *xdrs, const char *filename, const char *type);
int xdrclose(XDR *xdrs) ;
int xdr3dfcoord(XDR *xdrs, float *fp, int *size, float *precision) ;
