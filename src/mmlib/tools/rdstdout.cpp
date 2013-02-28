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

#include "tools/rdstdout.h"

char *rdstdout_file = NULL;

void set_rdstdout_file(char *filename){
	rdstdout_file = filename;
}

void clearrdstdout(){
	if(rdstdout_file != NULL) {
		FILE *rdfile = fopen(rdstdout_file, "w");
		if(rdfile == NULL)
			return;
		fclose(rdfile);
	}
}

void rdprintf(const char *a){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}
