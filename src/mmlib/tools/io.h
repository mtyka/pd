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

#ifndef __IO_H
#define __IO_H

#include "maths/maths.h"
#include <stdlib.h>
#include <string>

namespace IO // tools for file manipulation and information import
{
	// ----------------------------------------------------------------------------------
	// --- Handles conversion of binary data to a text-based storage container ----------
	// ----------------------------------------------------------------------------------

	// [todo] deprecated i believe ... (mook)
	void PD_API picklemem(const unsigned char* memory, unsigned length, std::string &jar);
	int PD_API depicklemem(unsigned char* memory, unsigned length, const std::string &jar);

	/// MIME Type 6bit encoding (bundles 3 resl bytes into 4 encoded bytes)
	void PD_API encode6bit(const unsigned char* memory, unsigned length, std::string &jar);
	/// MIME Type 6bit decoding (debundles 4 encoded bytes into 3 real bytes)
	int PD_API decode6bit(unsigned char** memory, const std::string &jar);

	// ----------------------------------------------------------------------------------
	// --- File Handing Functions & Classes ---------------------------------------------
	// ----------------------------------------------------------------------------------

	long PD_API getFileSize( const std::string &_filename );
	bool PD_API fileExists(const std::string &_filename);

	// [todo] need updating to std::string (mook)
	int PD_API fcopy(char *, char *); // copies a file from one to other
	int PD_API readFlatVectorFile(char *filename, Maths::dvector ** point, int *npoints);
	int PD_API readFlatFloatFile(char *filename, double **value, int *nvalues);

	// Generic File Handle
	// This is directly lifted from bjarne stroustrup's C++ Classic





//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API FilePtr{
	protected:
		FILE * p;
	public:
		FilePtr(const char *name, const char *mode){
			p = fopen(name, mode);
		}
		~FilePtr(){
			if(p != NULL)
				fclose(p);
		}

		bool loaded(){
			return p == NULL ? false : true;
		}
		operator FILE *(){
			return p;
		}
	};

	// Elaboration of the above that handles line by line text files





//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API TextFilePtr: public FilePtr{
	public:
		TextFilePtr(const char *name, const char *mode):FilePtr(name, mode) {
		};

		// this one just masks fgets
		int readline_raw(char *buffer, int max){
			fgets(buffer, max, p);
			return 1;
		}
		// this one will skip comment lines (default comment sign is '#'
		// returns the number of lines read
		int readline(char *buffer, int max){
			int i, lines;
			for(lines = 0; lines > (-1); lines++) {
				if(feof(p))
					return -1; // quit if EOF
				fgets(buffer, max, p);
				for(i = 0; i < max; i++) {
					if((buffer[i] == 0))
						break;
					if((buffer[i] == ' ') || (buffer[i] == '\t') || (buffer[i] < 32))
						continue; // ignore whitespaces
					break;
				}

				if((buffer[i] == '#'))
					continue; // ignore comment lines
				else if((buffer[i] == 0))
					continue; // empty line, ignore
				else
					break;
			}
			return lines + 1;
		}
	};

} // namespace IO
#endif


