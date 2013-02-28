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

// tools/rdstdout.h
//
// redirects all printf statements into a file named rdstdout_file


#ifndef __RDSTDOUT_H
#define __RDSTDOUT_H

#include <stdio.h>


extern char *rdstdout_file;

void set_rdstdout_file(char *filename);
void clearrdstdout();
void rdprintf(const char *a);

template < class  B >
void rdprintf(const char *a, B b){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a, (B) b);
		fflush(rdfile);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}
template < class  B, class PD_API C > void rdprintf(const char *a, B b, C c){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a, (B) b, (C) c);
		fflush(rdfile);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}
template < class  B, class PD_API C, class PD_API D > void rdprintf(const char *a, B b, C c, D d){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a, (B) b, (C) c, (D) d);
		fflush(rdfile);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}
template < class  B, class PD_API C, class PD_API D, class PD_API E > void rdprintf(const char *a, B b, C c, D d, E e){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a, (B) b, (C) c, (D) d, (E) e);
		fflush(rdfile);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}
template < class  B, class PD_API C, class PD_API D, class PD_API E, class PD_API F > void rdprintf(const char *a, B b, C c, D d, E e, F f){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a, (B) b, (C) c, (D) d, (E) e, (F) f);
		fflush(rdfile);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}
template < class  B, class PD_API C, class PD_API D, class PD_API E, class PD_API F, class PD_API G > void rdprintf(const char *a, B b, C c, D d, E e, F f, G g){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a, (B) b, (C) c, (D) d, (E) e, (F) f, (G) g);
		fflush(rdfile);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}
template < class  B, class PD_API C, class PD_API D, class PD_API E, class PD_API F, class PD_API G, class PD_API H >
void rdprintf(const char *a, B b, C c, D d, E e, F f, G g, H h){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a, (B) b, (C) c, (D) d, (E) e, (F) f, (G) g, (H) h);
		fflush(rdfile);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}
template < class  B, class PD_API C, class PD_API D, class PD_API E, class PD_API F, class PD_API G, class PD_API H, class PD_API I >
void rdprintf(const char *a, B b, C c, D d, E e, F f, G g, H h, I i){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a, (B) b, (C) c, (D) d, (E) e, (F) f, (G) g, (H) h, (I) i);
		fflush(rdfile);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}
template < class  B, class PD_API C, class PD_API D, class PD_API E, class PD_API F, class PD_API G, class PD_API H, class PD_API I, class PD_API J >
void rdprintf(const char *a, B b, C c, D d, E e, F f, G g, H h, I i, J j){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a, (B) b, (C) c, (D) d, (E) e, (F) f, (G) g, (H) h, (I) i, (J) j);
		fflush(rdfile);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}
template < class  B, class PD_API C, class PD_API D, class PD_API E, class PD_API F, class PD_API G, class PD_API H, class PD_API I, class PD_API J, class PD_API K >
void rdprintf(const char *a, B b, C c, D d, E e, F f, G g, H h, I i, J j, K k){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a, (B) b, (C) c, (D) d, (E) e, (F) f, (G) g, (H) h, (I) i, (J) j, (K) k);
		fflush(rdfile);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}
template < class  B, class PD_API C, class PD_API D, class PD_API E, class PD_API F, class PD_API G, class PD_API H, class PD_API I, class PD_API J, class PD_API K, class PD_API L >
void rdprintf(const char *a, B b, C c, D d, E e, F f, G g, H h, I i, J j, K k, L l){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a, (B) b, (C) c, (D) d, (E) e, (F) f, (G) g, (H) h, (I) i, (J) j, (K) k, (L) l);
		fflush(rdfile);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}
template < class  B, class PD_API C, class PD_API D, class PD_API E, class PD_API F, class PD_API G, class PD_API H, class PD_API I, class PD_API J, class PD_API K, class PD_API L, class PD_API M >
void rdprintf(const char *a, B b, C c, D d, E e, F f, G g, H h, I i, J j, K k, L l, M m){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a, (B) b, (C) c, (D) d, (E) e, (F) f, (G) g, (H) h, (I) i, (J) j, (K) k, (L) l, (M) m);
		fflush(rdfile);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}
template < class  B, class PD_API C, class PD_API D, class PD_API E, class PD_API F, class PD_API G, class PD_API H, class PD_API I, class PD_API J, class PD_API K, class PD_API L, class PD_API M,
class PD_API N > void rdprintf(const char *a, B b, C c, D d, E e, F f, G g, H h, I i, J j, K k, L l, M m, N n){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a, (B) b, (C) c, (D) d, (E) e, (F) f, (G) g, (H) h, (I) i, (J) j, (K) k, (L) l, (M) m, (N) n);
		fflush(rdfile);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}
template < class  B, class PD_API C, class PD_API D, class PD_API E, class PD_API F, class PD_API G, class PD_API H, class PD_API I, class PD_API J, class PD_API K, class PD_API L, class PD_API M,
class PD_API N, class PD_API O > void rdprintf(const char *a, B b, C c, D d, E e, F f, G g, H h, I i, J j, K k, L l, M m, N n, O o){
	FILE *rdfile;
	if(rdstdout_file != NULL) {
		rdfile = fopen(rdstdout_file, "at");
	} else {
		rdfile = stdout;
	}
	if(rdfile != NULL) {
		fprintf(rdfile, a, (B) b, (C) c, (D) d, (E) e, (F) f, (G) g, (H) h, (I) i, (J) j, (K) k, (L) l, (M) m, (N) n,
			(O) o);
		fflush(rdfile);
		if(rdstdout_file != NULL)
			fclose(rdfile);
	}
}

// One can switch off redirection at compile time by defining NORDPRINTF
// Use fprintf(stdout, ... ) or frpintf(stdout, ... ) to force screen output independent
// of redirection

#ifndef NORDPRINTF
#define printf rdprintf
#endif

#endif
