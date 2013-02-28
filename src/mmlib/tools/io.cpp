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

#include "io.h"

using namespace Maths;

namespace IO // tools for file manipulation and information import
{
	void c_flush_stdout()
	{
		fflush(stdout);
	}

	void picklemem(const unsigned char* memory, unsigned length, std::string &jar)
	{
		jar = "";
		unsigned char directlimit = 34;

		unsigned i;
		unsigned ijar=0;
		unsigned linewidth=100;

		for(i=0;i<length;i++){
			if(memory[i] > directlimit){
				//printf("add %d \n",memory[i]);
				jar += memory[i];
				ijar ++;
			}else{
				//printf("%c",directlimit);
				jar += directlimit;
				ijar ++;
				if((ijar>0)&&((ijar%linewidth)==0)){
					jar += "\n";
				}
				unsigned char mod = (directlimit + 2 + memory[i]);
				//printf("add %d %d \n",directlimit,mod);
				ijar ++;
				jar += mod;
			}
			if((ijar>0)&&((ijar%linewidth)==0)){
				jar += "\n";
			}
		}
	}

	// returns -1 on error or the number of bytes written
	int depicklemem(unsigned char* memory, unsigned length, const std::string &jar){
		unsigned char directlimit = 34;

		unsigned char minusnext=0;

		size_t i;

		const unsigned char *jarmemory = (unsigned char *)jar.c_str();
		unsigned int pos=0;

		for(i=0;i<jar.size();i++){
			//printf("get %d \n", jarmemory[i]);
			if(pos>=length) break;
			if(jarmemory[i] < directlimit){
				//printf("unexpected char: %d ", (int) jarmemory[i] );
				continue;
			}
			if(jarmemory[i] > directlimit){
				memory[pos] = jarmemory[i]-minusnext;
				//printf("->%d\n",memory[pos]);
				minusnext=0;
				pos++;
				continue;
			}
			if(jarmemory[i] == directlimit){
				minusnext = directlimit+2;
			}
		}
		//memory[pos]=0;
		return pos;
	}
	inline unsigned char code_to_6bit(unsigned char _6bit)
	{
		//return _6bit;
		const unsigned char conversion [] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
		return conversion[ _6bit & 63 ];
	}

	inline unsigned char code_from_6bit(unsigned char _8bit)
	{
		//return _8bit;
		if( ( _8bit >= 'A') && (_8bit <= 'Z') ) return _8bit - 'A';
		if( ( _8bit >= 'a') && (_8bit <= 'z') ) return _8bit - 'a' + 26;
		if( ( _8bit >= '0') && (_8bit <= '9') ) return _8bit - '0' + 52;
		if( ( _8bit == '+') ) return 62;
		return 63;
	}

	inline void encode_24_to_32(
		unsigned char i0,
		unsigned char i1,
		unsigned char i2,
		unsigned char &o0,
		unsigned char &o1,
		unsigned char &o2,
		unsigned char &o3
	)
	{
		// converts 3 8 bit chars into 
	  o0 = code_to_6bit(i0 & 63);  // delete last two bits
	  o1 = code_to_6bit(((i1 << 2) | (i0 >> 6)) & 63);
	  o2 = code_to_6bit(((i1 >> 4) | ((i2 << 4) & 63)) & 63);
	  o3 = code_to_6bit(i2 >> 2);  // right shift by two 
	}

	inline void decode_32_to_24(
		unsigned char i0,
		unsigned char i1,
		unsigned char i2,
		unsigned char i3,
		unsigned char &o0,
		unsigned char &o1,
		unsigned char &o2
	)
	{
		i0 = code_from_6bit(i0);
		i1 = code_from_6bit(i1);
		i2 = code_from_6bit(i2);
		i3 = code_from_6bit(i3);
		o0 = i0 | (i1 << 6); 
		o1 = (i1 >> 2) | (i2 << 4);
		o2 = (i3 << 2) | (i2 >> 4);
	}



	void encode6bit(const unsigned char* memory, unsigned length, std::string &jar){
		jar = "";
		unsigned i;
		unsigned allfourcount=0;
		unsigned fourcount=0;
		unsigned linewidth=15;
		for(i=0;i<length;){
			unsigned char buffer[3] = {0,0,0};
			unsigned char outbuffer[4] = {0,0,0,0};
			int ibuf=0;
			for(;((i<length)&&(ibuf<3));i++){
				buffer[ibuf] = memory[i];
				ibuf++;
			}
			encode_24_to_32(buffer[0],buffer[1],buffer[2], 
				              outbuffer[0],outbuffer[1],outbuffer[2],outbuffer[3]);
			jar += outbuffer[0];
			jar += outbuffer[1];
			jar += outbuffer[2];
			jar += outbuffer[3];
			fourcount +=1;
			allfourcount +=1;	
			if(fourcount > linewidth){
				fourcount = 0;
				jar += '\n';
			}
		}
		jar += '\n';
	}

	int decode6bit(unsigned char** memory, const std::string &jar){
		//printf("-->%s\n",jar.c_str());
		const unsigned char *jarmemory = (unsigned char *)jar.c_str();
		unsigned memlength = (unsigned)jar.length() * 3 / 4 + 1;
		unsigned mempos = 0;
		*memory = new unsigned char [memlength];  
		for(size_t i=0;i<jar.size();){
			unsigned char inbuf[4] = {0,0,0,0};
			unsigned char ibuf=0;
			for(;((i<jar.size())&&(ibuf<4));i++){
				inbuf[ibuf] = jarmemory[i];
				if(inbuf[ibuf] < 32) continue;
				ibuf++;
			}

			decode_32_to_24(inbuf[0],inbuf[1],inbuf[2],inbuf[3],
											(*memory)[mempos+0],(*memory)[mempos+1],(*memory)[mempos+2]);

			mempos += 3;			
		}
		return mempos;
	}



	long getFileSize( const std::string &_filename )
	{
		std::ifstream in( _filename.c_str(), std::ios::in );
		if( !in.is_open() ) throw new IOException("getFileSize() failed. Filename is not valid or file is locked");
		in.seekg( 0, std::ios::end );
		std::ios::pos_type pos = in.tellg();
		in.close();
		return (long)pos;
	}

	bool PD_API fileExists(const std::string &_filename)
	{
		FILE *file;
		file = fopen(_filename.c_str(), "r");
		if(file == NULL)
			return false;
		fclose(file);
		return true;
	}

	int fcopy(char *dstfilename, char *srcfilename)
	{
		FILE *dstfile, *srcfile;
		char *memory;
		int bsize = 0xFFFF;
		int n;

		if(dstfilename == NULL)
			return -1;
		dstfile = fopen(dstfilename, "wb");
		if(dstfile == NULL) {
			printf("ERROR: Cannot open %s for writing \n", dstfilename);
			return -1;
		}

		if(srcfilename == NULL)
			return -1;
		srcfile = fopen(srcfilename, "rb");
		if(srcfile == NULL) {
			printf("ERROR: Cannot open %s for writing \n", srcfilename);
			return -1;
		}

		memory = new char[bsize];

		while(!feof(srcfile)) {
			n = (int) fread((void *) memory, 1, bsize, srcfile);
			fwrite((void *) memory, n, 1, dstfile);
			if(n < bsize)
				break; // probably end of file - end copying
		}

		fclose(srcfile);
		fclose(dstfile);
		delete[]memory;

		return 0;
	}

	int readFlatVectorFile(char *filename, dvector ** point, int *npoints){
		FILE *file;
		char buffer[256];
		int n, i;
		double x, y, z;

		printf("Opening %s .. \n", filename);
		file = fopen(filename, "r");
		if(file == NULL)
			return -1;


		// first line must contain integer with number of entries
		fgets(&buffer[0], 256, file);
		n = sscanf(&buffer[0], "%d", npoints);

		if((n != 1) || (*npoints > 100000)) {
			printf("ERROR: couldn't read number of entries in line 1\n");
			fclose(file);
			return -1;
		}


		*point = new dvector[*npoints];

		for(i = 0; i < *npoints; i++) {
			fgets(&buffer[0], 256, file);
			n = sscanf(&buffer[0], "%lf %lf %lf", &x, &y, &z);
			if(n != 3) {
				printf("SYNTAX ERROR: line %d error \n", i + 2);
				fclose(file);
				return -1;
			}
			(*point)[i].setTo(x, y, z);
		}

		fclose(file);

		return 0;
	}

	int readFlatFloatFile(char *filename, double **value, int *nvalues){
		FILE *file;
		char buffer[256];
		int n, i;
		double x;

		//printf("Opening %s .. \n",filename);
		file = fopen(filename, "r");
		if(file == NULL)
			return -1;


		// first line must contain integer with number of entries
		fgets(&buffer[0], 256, file);
		n = sscanf(&buffer[0], "%d", nvalues);

		if((n != 1) || (*nvalues > 100000)) {
			printf("ERROR: couldn't read number of entries in line 1\n");
			fclose(file);
			return -1;
		}


		*value = new double[*nvalues];

		for(i = 0; i < *nvalues; i++) {
			fgets(&buffer[0], 256, file);
			n = sscanf(&buffer[0], "%lf", &x);
			if(n != 1) {
				printf("SYNTAX ERROR: line %d error \n", i + 2);
				fclose(file);
				return -1;
			}
			(*value)[i] = x;
		}

		fclose(file);

		return 0;
	}

	/*
	int findCmdLineParam(int argc, char **argv, char *name){
	for(int i = 0; i < argc; i++)
	if(strcmp(name, argv[i]) == 0)
	return i;
	return -1;
	}

	int findCmdLineParam(int argc, char **argv, char *name, char **following){
	for(int i = 0; i < argc; i++)
	if(strcmp(name, argv[i]) == 0) {
	if(i < (argc - 1)) {
	if(*argv[i + 1] != '-')
	*following = argv[i + 1];
	else
	*following = NULL;
	return i;
	}
	}
	*following = NULL;
	return -1;
	}

	*/

} // namespace IO


