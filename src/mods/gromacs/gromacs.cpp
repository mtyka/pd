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

#include "rpcxdr.h"

extern "C" {
#include "xdrf.h"
}
#include "gromacs.h"

namespace IO 
{


	int OutTra_GRO::create()
	{
		FILE *file;
		file = fopen(filename_GRO.c_str(),"w");

		fprintf(file, "%s\n", "PD Trajectory");
		fprintf(file, "%5d\n", getWSpace().nAtoms() );
		for(size_t i=0; i<getWSpace().nAtoms(); i++) {

			fprintf(file, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
				getWSpace().atom[i].ir,
				trim(getWSpace().atom[i].parentname," ").c_str(),
				trim(getWSpace().atom[i].pdbname," ").c_str(),
				i,
				getWSpace().cur.atom[i].p.x/10.0,   // gromacs uses nanometers instead of angstroms
				getWSpace().cur.atom[i].p.y/10.0,
				getWSpace().cur.atom[i].p.z/10.0,
				getWSpace().cur.atom[i].v.x/1000.0, // from m/s to units of nm/ps 
				getWSpace().cur.atom[i].v.y/1000.0,
				getWSpace().cur.atom[i].v.z/1000.0);
		}

		ClosedSpace *space = castTo_ClosedSpace(&getWSpace().boundary());
		Maths::dvector A(0,0,0),B(0,0,0),C(0,0,0);

		if( space == NULL){
		}else{
			space->getBoxVectors(A,B,C);
		}
		fprintf(file, "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n",
			A.x, B.y, C.z,
			A.y, A.z, 
			B.x, B.z,
			C.y, C.z);
		fclose(file);

		created = true;

		return 0;
	}














	// Multiplication factor before integerisation of coordinates.
	// i.e. the resolution is 1/100 Angstrom = 1pico meter
	const float  XTC_Precision = 100.0;

	int OutTra_XTC::create()
	{
		delete xdrs; // just in case - to prevent a memory leak;
		xdrs = new XDR;
		if (xdropen(xdrs, filename_XTC.c_str() ,"w") == 0) {
			throw( IOException("ERROR opening " + filename_XTC + " for writing") );	
		}
		created = true;

		return 0;
	}


	OutTra_XTC::~OutTra_XTC(){
		if(!xdrclose(xdrs) ){ 
			throw( IOException("ERROR closing " + filename_XTC) );
		}
	}


	int OutTra_XTC::append()
	{
		if(!OutTra_XTC::created){
			create();
		}

		if(xdrs == NULL){
			throw( IOException("ERROR writing to " + filename_XTC + " for XTC trajectory") );	
		}

		int XTC_magic = 1995;                 // xtc format magic number
		int XTC_natoms = getWSpace().nAtoms();    
		int XTC_step = m_StepCount;
		float XTC_time = 0.0;
		m_StepCount++;

		// write magic number, natoms (constant), step number and time into trajectory
		if (!xdr_int(xdrs, &XTC_magic))    throw( IOException("ERROR writing magic number to " + filename_XTC) );	
		if (!xdr_int(xdrs, &XTC_natoms))   throw( IOException("ERROR writing number of atoms to " + filename_XTC) );	
		if (!xdr_int(xdrs, &XTC_step))     throw( IOException("ERROR writing step number to " + filename_XTC) );	
		if (!xdr_float(xdrs, &XTC_time))   throw( IOException("ERROR writing time stamp to " + filename_XTC) );	

		// Now write the box
		ClosedSpace *space = castTo_ClosedSpace(&getWSpace().boundary());
		Maths::dvector A(0,0,0),B(0,0,0),C(0,0,0);

		if( space == NULL){
		}else{
			space->getBoxVectors(A,B,C);
		}
		float XTC_box[3][3];
		XTC_box[0][0] = A.x;
		XTC_box[0][1] = A.y;
		XTC_box[0][2] = A.z;
		XTC_box[1][0] = B.x;
		XTC_box[1][1] = B.y;
		XTC_box[1][2] = B.z;
		XTC_box[2][0] = C.x;
		XTC_box[2][1] = C.y;
		XTC_box[2][2] = C.z;

		int i,j;
		for(i=0; i<3; i++)
		{
			for(j=0; j<3; j++)
			{
				if (!xdr_float(xdrs, &XTC_box[i][j]))
				{ 
					throw( IOException("ERROR writing box vectors to " + filename_XTC) );	
				}
			}
		}

		// put coordinates into a linear array in the form of xyzxyzxyz..xyz

		float *coordinate_array = new float [XTC_natoms * 3];
		for(i=0; i<XTC_natoms; i++)
		{
			coordinate_array[i*3 + 0] = getWSpace().cur.atom[i].p.x*0.1;  // convert to nano meters
			coordinate_array[i*3 + 1] = getWSpace().cur.atom[i].p.y*0.1;
			coordinate_array[i*3 + 2] = getWSpace().cur.atom[i].p.z*0.1;
		}

		// now write the coordinates to the file using 
		// Frans van Hoesel's compression algorithm
		float prec = XTC_Precision;
		if (xdr3dfcoord(xdrs, coordinate_array, &XTC_natoms, &prec) == 0) {
			throw( IOException("ERROR writing coordinates to " + filename_XTC) );	
		}

		delete [] coordinate_array;

		// make sure everything is written to the file by flushing
		fflush( (FILE *) xdrs->x_private );

		return 0;
	}








void OutputFile_GRO::save( WorkSpace &_wspace){
	OutTra_GRO outtra( filestem, _wspace );
	outtra.append( );
}


void OutputFile_XTC::save( WorkSpace &_wspace){
	OutTra_XTC outtra( filestem, _wspace );
	outtra.append( );
}
















	////////////////////////////////////////////////////////////////////


	int InTra_XTC::open(){
		delete xdrs; // just in case - to prevent a memory leak;
		xdrs = new XDR;
		if (xdropen(xdrs, filename.c_str() ,"r") == 0) {
			throw( IOException("ERROR opening " + filename + " for reading") );	
		}
		return 0;
	}

	int InTra_XTC::close(){
		if(!xdrclose(xdrs) ){ 
			throw( IOException("ERROR closing " + filename) );
		}
		return 0;
	}

	bool InTra_XTC::readNext( SnapShot &ss ){
		if(xdrs == NULL){
			throw( IOException("ERROR reading  " + filename + " for XTC trajectory") );	
		}

		int XTC_magic = 0;                 // xtc format magic number - expecting 1995 here
		int XTC_natoms = 0; 
		int XTC_step = 0; 
		float XTC_time = 0.0;

		// write magic number, natoms (constant), step number and time into trajectory
		if (!xdr_int(xdrs, &XTC_magic))    return true; 
		if (!xdr_int(xdrs, &XTC_natoms))   throw( IOException("ERROR reading number of atoms from " + filename) );	
		if (!xdr_int(xdrs, &XTC_step))     throw( IOException("ERROR reading step number from " + filename) );	
		if (!xdr_float(xdrs, &XTC_time))   throw( IOException("ERROR reading time stamp from " + filename) );	

		// test the data:
		if (XTC_magic != 1995)             throw( IOException("ERROR reading " + filename + ". Magic number is not 1995 as expected - are you using a newer/older format of XTC ?")  );
		if (XTC_natoms < 0)                throw( IOException("ERROR reading " + filename + ". THe number of atoms is negative. Is the file corrupt ?")  ); 

		// Now read the box vectors
		Maths::dvector A(0,0,0),B(0,0,0),C(0,0,0);
		float XTC_box[3][3];

		ss = SnapShot( XTC_natoms );

		int i,j;
		for(i=0; i<3; i++)
		{
			for(j=0; j<3; j++)
			{
				if (!xdr_float(xdrs, &XTC_box[i][j]))
				{ 
					throw( IOException("ERROR reading from box vectors to " + filename) );	
				}
			}
		}

		A.x = XTC_box[0][0];
		A.y = XTC_box[0][1];
		A.z = XTC_box[0][2];
		B.x = XTC_box[1][0];
		B.y = XTC_box[1][1];
		B.z = XTC_box[1][2];
		C.x = XTC_box[2][0];
		C.y = XTC_box[2][1];
		C.z = XTC_box[2][2];
		ss.A = A;
		ss.B = B;
		ss.C = C;

		// put coordinates into a linear array in the form of xyzxyzxyz..xyz

		float *coordinate_array = new float [XTC_natoms * 3];
		for(i=0; i<XTC_natoms; i++)
		{
			coordinate_array[i*3 + 0] = 0;  // convert from nano meters
			coordinate_array[i*3 + 1] = 0;
			coordinate_array[i*3 + 2] = 0;
		}

		// now read the coordinates from the file using 
		// Frans van Hoesel's compression algorithm
		float prec = XTC_Precision;
		if (xdr3dfcoord(xdrs, coordinate_array, &XTC_natoms, &prec) == 0) {
			fprintf(stderr,"error while writing coordinates\n");
		}
		for(i=0; i<XTC_natoms; i++)
		{
			ss.atom[i].p.x = coordinate_array[i*3 + 0]*10;  // convert from nano meters
			ss.atom[i].p.y = coordinate_array[i*3 + 1]*10;
			ss.atom[i].p.z = coordinate_array[i*3 + 2]*10;
		}

		delete [] coordinate_array;

		return false;
	}

	bool InTra_XTC::skip()
	{
		// naiive implementation of skip. just read and discard result.
		SnapShot ss;
		readNext(ss);	
		return false;
	}

	bool InTra_XTC::isEndOfFile() const
	{
		return (bool)feof( (FILE *) xdrs->x_private );
	}


}  //namespace IO 

