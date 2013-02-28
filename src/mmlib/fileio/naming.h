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

#ifndef __NAMING_H
#define __NAMING_H

// ----------------------------------------------------------------------------------------
// Naming.h
// DEPRECATED: Note this functionality should no longer be used!
// It has been replaced by functionality in sequence.h
// ----------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------
// --- Naming constants ---
// ----------------------------------------------------------------------------------

//extern const char Nameconvention_AAOrder[22];
//extern char *PDB_Nameconvention[22][50];
//extern char *XPLOR_Nameconvention[22][50];
//extern char *PDAMBER_Nameconvention[22][50];

// --------------------------------------------------------------------------------
// Convertion from aminoacid numbers to letter and vice versa
// --------------------------------------------------------------------------------

char getAALetterFromFullName(char *resname);
char getAANumberFromFullName(char *resname);
char getAALetter(int aanum);
// The following are unused:
//const char *getAANameFull(int aanum);
//int getAANumber(char aachar);

#endif


