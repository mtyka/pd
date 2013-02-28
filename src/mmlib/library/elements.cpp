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

#include "elements.h"

namespace Library
{
	// ----------------------------------------------------------------------------------
	// --- The Elements
	// ----------------------------------------------------------------------------------

	// incomplete, stops before Lanthanoids, but hey .. :)
	const char *element[] = { "-",
		"H", "He",
		"Li","Be", "B" ,"C" , "N", "O", "F","Ne",
		"Na","Mg", "Al","Si", "P", "S","Cl","Ar",
		"K" ,"Ca","Sc","Ti","V", "Cr","Mn","Fe","Co","Ni","Cu","Zn", "Ga","Ge","As","Se","Br","Kr",
		"Rb","Sr","Y", "Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd", "In","Sn","Sb","Te","I", "Xe",
		"Cs","Ba" };

	int ElementSymbol2Z(const std::string &Symbol){
		char buffer1[15];
		char buffer2[15];
		strncpy(&buffer1[0],Symbol.c_str(),15);
		strlc(&buffer1[0]);
		for(int i=0;i<56;i++){
			strcpy(&buffer2[0],element[i]);
			strlc(&buffer2[0]);
			if(strcmp(&buffer1[0],&buffer2[0])==0){
				return i;
			}
		}
		return -1;
	}
}

