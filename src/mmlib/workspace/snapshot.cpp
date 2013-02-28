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

#include "snapshot.h"

#include "maths/maths.h"
#include "forcefields/ffparam.h"
#include "forcefields/forcefield.h"
#include "workspace/workspace.h"
#include "tools/io.h"

using namespace Maths;

SnapShot::SnapShot(size_t _natoms)
{
	ASSERT(_natoms > 0,ArgumentException,"Cannot create PhaseSpacePoint - natoms <= 0");

	natoms = _natoms;
	atom = new SnapShotAtom[ natoms ];
	epot=0;
	status1=0;
	status2=0;
	status3=0;
	status4=0;
	A.setTo(0,0,0);
	B.setTo(0,0,0);
	C.setTo(0,0,0);
}

SnapShot::~SnapShot()
{
	delete [] atom;
}

SnapShot::SnapShot(const SnapShot &newpsp)
{
	natoms = 0;
	atom = NULL;
	(*this) = newpsp;
}

void SnapShot::setTo(const SnapShot &newpsp)
{
	(*this) = newpsp;
}

void SnapShot::setToPartial(const SnapShot &newpsp)
{
	for(int i = 0;i < min(natoms,newpsp.natoms); i++)
	{
		atom[i] = newpsp.atom[i];
	}
}

SnapShot &SnapShot::operator= (const SnapShot &psp)
{
	// check for self assignement
	if(&psp == this)	return (*this);

	// reallocate only if the number of atoms changes
	if( natoms != psp.natoms )
	{
		delete[] atom;
		natoms = psp.natoms;	
		atom = new SnapShotAtom[ natoms ];
	}

	for( int i = 0; i < natoms; i++ )
	{
		atom[i] = psp.atom[i];
	}

	epot    = psp.epot;
	status1	= psp.status1;
	status2	= psp.status2;
	status3	= psp.status3;
	status4 = psp.status4;

	return (*this);
}



double SnapShot::cRMSFrom(const SnapShot &psp2) const
{	
	if( natoms != psp2.natoms ) throw("PhaseSpacePoints are incompatible (Total number of atoms doesnt match)");
	
	int i;
	double crms;
	int nrmsatoms = 0;

	dvector *atomarray1 = new dvector[natoms];
	dvector *atomarray2 = new dvector[natoms];

	for(i = 0; i < natoms; i++) 
	{
		atomarray1[nrmsatoms].setTo(atom[i].p.x, atom[i].p.y, atom[i].p.z);
		atomarray2[nrmsatoms].setTo(psp2.atom[i].p.x, psp2.atom[i].p.y, psp2.atom[i].p.z);
		nrmsatoms++;
	}

	crms = calcVectorCRMS(&atomarray1[0], &atomarray2[0], nrmsatoms);

	delete[]atomarray1;
	delete[]atomarray2;

	return crms;
}

double SnapShot::cRMSFrom(
						  const SnapShot &psp2,
						  const MoleculeBase& sys,
						  const std::string& name
						  ) const
{
	if( natoms != psp2.natoms ) throw("PhaseSpacePoints are incompatible (atom number doesnt match up)");
	if( natoms != sys.atom.size() ) throw("PhaseSpacePoints are incompatible with System (sys) (atom number doesnt match up)");

	int i;
	double crms;
	int nrmsatoms = 0;

	dvector *atomarray1 = new dvector[natoms];
	dvector *atomarray2 = new dvector[natoms];

	for(i = 0; i < natoms; i++) 
	{
		if(sys.atom[i].pdbname != name) continue;
		atomarray1[nrmsatoms].setTo(atom[i].p.x, atom[i].p.y, atom[i].p.z);
		atomarray2[nrmsatoms].setTo(psp2.atom[i].p.x, psp2.atom[i].p.y, psp2.atom[i].p.z);
		nrmsatoms++;
	}

	crms = calcVectorCRMS(&atomarray1[0], &atomarray2[0], nrmsatoms);

	delete[]atomarray1;
	delete[]atomarray2;

	return crms;
}



double SnapShot::dRMSFrom(const SnapShot &psp2) const{
	if( natoms != psp2.natoms ) throw("PhaseSpacePoints are incompatible (atom number doesnt match up)");
	int i, j;
	double distij1, distij2;
	double drms = 0.0;
	int count = 0;

	for(i = 0; i < natoms; i++) {
		for(j = i + 1; j < natoms; j++) {
			distij1 =  sqrt(sqr(atom[i].p.x - atom[j].p.x) +
											sqr(atom[i].p.y - atom[j].p.y) +
											sqr(atom[i].p.z - atom[j].p.z));
			distij2 =  sqrt(sqr(psp2.atom[i].p.x - psp2.atom[j].p.x) +
											sqr(psp2.atom[i].p.y - psp2.atom[j].p.y) +
											sqr(psp2.atom[i].p.z - psp2.atom[j].p.z));
			drms += sqr(distij1 - distij2);
			count++;
		}
	}

	drms /= (double) count;
	drms = sqrt(drms);
	return drms;
}





int SnapShot::insertPSPfragment(const SnapShot &psp2,
																						int startir, int endir,
																						const MoleculeBase &wspace){
	int i;
	dvector vsys_sten_vec;
	dvector vsys_sten_CCA;
	dvector vsys_sten_CN;
	matrix3x3 stdrot;

	// find end and start particles ...
	int sysCpre0 =	wspace.findParticle(startir - 1, "C");
	int sysN0 =			wspace.findParticle(startir, "N");
	int sysCA0 =		wspace.findParticle(startir, "CA");

	int sysNpostn = wspace.findParticle(endir, "N");
	int sysCn =			wspace.findParticle(endir - 1, "C");
	int sysCAn =		wspace.findParticle(endir - 1, "CA");

	if((sysCpre0 < 0)) {
		// printf("Alternative third ref. atom at beginning \n");
		sysCpre0 = wspace.findParticle(startir, "C");
	}
	if((sysNpostn < 0)) {
		// printf("Alternative third ref. atom at end\n");
		sysNpostn = wspace.findParticle(endir - 1, "N");
	}
	// printf("Segment to be replaced: %d %d %d --> %d %d %d \n",
	// sysCpre0,sysN0,sysCA0,sysNpostn,sysCn ,sysCAn );

	if((sysCpre0 < 0) ||
		(sysNpostn < 0) ||
		(sysN0 < 0) ||
		(sysCn < 0) ||
		(sysCA0 < 0) ||
		(sysCAn < 0))
		return -1;

	dvector psp1_sysCpre0(atom[sysCpre0].p);
	dvector psp1_sysN0(atom[sysN0].p);
	dvector psp1_sysCA0(atom[sysCA0].p);

	dvector psp1_sysNpostn(atom[sysNpostn].p);
	dvector psp1_sysCn(atom[sysCn].p);
	dvector psp1_sysCAn(atom[sysCAn].p);

	dvector psp2_sysCpre0(psp2.atom[sysCpre0].p);
	dvector psp2_sysN0(psp2.atom[sysN0].p);
	dvector psp2_sysCA0(psp2.atom[sysCA0].p);

	dvector psp2_sysNpostn(psp2.atom[sysNpostn].p);
	dvector psp2_sysCn(psp2.atom[sysCn].p);
	dvector psp2_sysCAn(psp2.atom[sysCAn].p);

	dvector trans; // translation dvector  and rotation matrix for fragment
	matrix3x3 rmat;
	dvector trans2; // ditto for rest of protein (after fragment)
	matrix3x3 rmat2;

	// superimpose first three atoms of fragment with last three atoms of
	// preceding protein segment
	superimpose(psp1_sysCA0,
		psp1_sysN0,
		psp1_sysCpre0,
		psp2_sysCA0,
		psp2_sysN0,
		psp2_sysCpre0,
		rmat);


	// work out the translation dvector  (after rotation by rmat) such that the
	// N atoms of the first fragment residue of protein and fragment are superimposed
	trans.setTo(psp2_sysN0);
	trans.mulmat(rmat);
	trans.sub(psp1_sysN0);

	// printf("Fragatom Transform: vec & matrix \n");
	// trans.info();
	// rmat.info();
	// rotate the fragment atoms into place

	dvector Temp1;
	for(i = wspace.res[startir].ifirst; i <= wspace.res[endir - 1].ilast; i++) {
		Temp1.setTo(psp2.atom[i].p);
		Temp1.mulmat(rmat);
		Temp1.sub(trans);
		atom[i].p.setTo(Temp1);
	}

	psp2_sysCAn.mulmat(rmat);
	psp2_sysCAn.sub(trans);
	psp2_sysCn.mulmat(rmat);
	psp2_sysCn.sub(trans);
	psp2_sysNpostn.mulmat(rmat);
	psp2_sysNpostn.sub(trans);

	// now superimpose first last atoms of the fragment with the first atoms of the
	// succeeding protein part
	superimpose(psp2_sysCAn, psp2_sysCn, psp2_sysNpostn, psp1_sysCAn, psp1_sysCn, psp1_sysNpostn, rmat2);

	// work out the translation dvector  (after rotation by rmat) such that the
	// C atoms of the last fragment residue of protein and fragment are superimposed
	trans2.setTo(psp1_sysCn);
	trans2.mulmat(rmat2);
	trans2.sub(psp2_sysCn);

	// rotate the fragment succeeding atoms into place
	for(i = wspace.res[min((int)wspace.res.size(), (int)endir)].ifirst; i < wspace.atom.size(); i++) {
		Temp1.setTo(atom[i].p);
		Temp1.mulmat(rmat2);
		Temp1.sub(trans2);
		atom[i].p.setTo(Temp1);
	}

	return 0;
}

























int SnapShot::writeMIME(const std::string &filename,const std::string &name, const char *mode){ // append by default
	FILE *file = fopen(filename.c_str(),mode);
	if(file == NULL) return -1;

	writeMIME(file,name);

	fclose(file);
	return 0;
}
int SnapShot::appendMIME(const std::string &filename,const std::string &name){ // append by default
	return writeMIME(filename,name,"a");
}
void SnapShot::printMIME(const std::string &name){
	writeMIME(stdout,name);
}

void SnapShot::writeMIME(FILE *file,const std::string &name){
	fprintf(file,"!PSP %s\n",name.c_str());

	// make some memory
	unsigned int length = 128 +            /// 100 bytes for the header info
												natoms*9*8;     /// rest for 16 doubles * natoms
	unsigned char *tempMemory = new unsigned char [length];

	// copy over the contents of 'this'
	*((int*)(&tempMemory[0])) = natoms;
	*((double*)(&tempMemory[4])) = epot;
	*((int*)(&tempMemory[12])) = status1;
	*((int*)(&tempMemory[16])) = status2;
	*((int*)(&tempMemory[20])) = status3;
	*((int*)(&tempMemory[24])) = status4;

	double *mempointer = (double*)(&tempMemory[128]);
	for(int i=0;i<natoms;i++){
		*mempointer = atom[i].p.x; mempointer++;
		*mempointer = atom[i].p.y; mempointer++;
		*mempointer = atom[i].p.z; mempointer++;
		*mempointer = atom[i].v.x; mempointer++;
		*mempointer = atom[i].v.y; mempointer++;
		*mempointer = atom[i].v.z; mempointer++;
		*mempointer = atom[i].f.x; mempointer++;
		*mempointer = atom[i].f.y; mempointer++;
		*mempointer = atom[i].f.z; mempointer++;
	}

	std::string jar;
	IO::encode6bit(tempMemory,length,jar);

	fprintf(file,"%s",jar.c_str());
	fprintf(file,"\n!ENDPSP\n");
}

int SnapShot::readMIME(const std::string &filename,const std::string &name){
	FILE *file;
	int returnval=0;
	file = fopen(filename.c_str(),"r");
	if(file==NULL){
		throw(IOException(std::string("ERROR: Error reading pickled structure from '") + filename + std::string("'")));
	}
	returnval = readMIME(file,name);
	if(returnval != 0){
		throw(ParseException(std::string("ERROR: Cannot find MIME structure '") + name + std::string("' in ") + filename));
	}
	fclose(file);
	printf("Loaded PSP from %s: %d atoms \n",filename.c_str(), natoms );
	return 0;
}

int SnapShot::readMIME(FILE *file,const std::string &name){
	while(!feof(file)){
		
		std::string line;
		if(getSingleLine(file, line)!=0) return -1;
		//printf(">%s",line.c_str());
		if(cmpstring(line.substr(0,4),"!PSP")){
 			std::vector <std::string> token;
			token = chopstr(line," \t\12\15");
			bool foundtag=false;
			if(cmpstring(name,""))			 foundtag = true;
			else {
				if(token.size() >= 2){
					if(cmpstring(token[1],name)) foundtag = true;
				}
			}

			if(foundtag){ // found the correct pickle
				std::string pickledata;
				// now gather the pickled data
				while(!feof(file)){
					std::string line;
					if(getSingleLine(file, line)!=0) break;
					//printf("=%s",line.c_str());
					if(cmpstring(line.substr(0,7),"!ENDPSP")) break;
					pickledata += line; // add
				}

				// now depickle the data
				// make some memory

				delete [] atom;

				unsigned int length;
				unsigned char *tempMemory;
				length = IO::decode6bit(&tempMemory,pickledata);

				// copy over the contents of 'this'
				natoms = *((int*)(&tempMemory[0]));
				epot = *((double*)(&tempMemory[4]));
				status1 = *((int*)(&tempMemory[12]));
				status2 = *((int*)(&tempMemory[16]));
				status3 = *((int*)(&tempMemory[20]));
				status4 = *((int*)(&tempMemory[24]));
				printf("Read Snapshot: Atoms: %d  Energy: %7.2f  \n",natoms,epot*Physics::PhysicsConst::Na *Physics::PhysicsConst::J2kcal );
				if(natoms <= 0){
					throw(ParseException("Cannot read MIME structure due to data corruption - natoms appears to be <= 0"));
				}
				atom = new SnapShotAtom [ natoms ];

				double *mempointer = (double*)(&tempMemory[128]);
				for(int i=0;i<natoms;i++){
					atom[i].p.x = *mempointer; mempointer++;
					atom[i].p.y = *mempointer; mempointer++;
					atom[i].p.z = *mempointer; mempointer++;
					atom[i].v.x = *mempointer; mempointer++;
					atom[i].v.y = *mempointer; mempointer++;
					atom[i].v.z = *mempointer; mempointer++;
					atom[i].f.x = *mempointer; mempointer++;
					atom[i].f.y = *mempointer; mempointer++;
					atom[i].f.z = *mempointer; mempointer++;
				}
				delete [] tempMemory;

				return 0;
			}
		}
	}
	return -1;
}


SnapShot readMIME(const std::string &filename,const std::string &name){
	SnapShot psp;
	psp.readMIME(filename,name);
	return psp;
}

int SnapShot::readOldPickle(const std::string &filename,const std::string &name){
	FILE *file;
	int returnval=0;
	file = fopen(filename.c_str(),"r");
	if(file==NULL){
		fprintf(stderr,"ERROR: Error reading file '%s' \n",filename.c_str());
		throw(IOException(""));
		return -1;
	}
	returnval = readOldPickle(file,name);
	if(returnval != 0){
		fprintf(stderr,"ERROR: Error reading pickled structure '%s'\n",name.c_str());
		throw(ArgumentException(""));
		return returnval;
	}
	fclose(file);
	printf("Loaded PSP from %s: %d atoms \n",filename.c_str(), natoms );
	return 0;
}

int SnapShot::readOldPickle(FILE *file,const std::string &name){
	while(!feof(file)){
		std::string line;
		if(::getSingleLine(file, line)!=0) return -1;

		if(cmpstring(line.substr(0,4),"!PSP")){
			std::vector <std::string> token;
			token = chopstr(line," \t\12\15");
			if(token.size() >= 2){
				bool foundtag=false;
				if(name == "")							 foundtag = true;
				if(cmpstring(token[1],name)) foundtag = true;

				if(foundtag){ // found the correct pickle
					std::string pickledata;
					// now gather the pickled data
					while(!feof(file)){
						std::string line;
						if(getSingleLine(file, line)!=0) break;
						if(cmpstring(line.substr(0,4),"!ENDPSP")) break;

						pickledata += line; // add
					}

					delete [] atom;

					unsigned int length = 280;
					unsigned char *tempMemory;
					tempMemory = new unsigned char [length];
					IO::depicklemem(tempMemory,length,pickledata);
					natoms = *((int*)(&tempMemory[0]));
					epot = *((double*)(&tempMemory[12]));
					if(natoms <= 0){
						throw(ParseException("Cannot read old pickle due to data corruption - natoms appears to be <= 0"));
					}
					delete [] tempMemory;
					length += natoms*9*8;  // enlarge the temp memory
					tempMemory = new unsigned char [length + 10];

					IO::depicklemem(tempMemory,length,pickledata);
					atom = new SnapShotAtom [ natoms ];

					double *mempointer = (double*)(&tempMemory[280]);
					for(int i=0;i<natoms;i++){
						atom[i].p.x = *mempointer; mempointer++;
						atom[i].p.y = *mempointer; mempointer++;
						atom[i].p.z = *mempointer; mempointer++;
						atom[i].v.x = *mempointer; mempointer++;
						atom[i].v.y = *mempointer; mempointer++;
						atom[i].v.z = *mempointer; mempointer++;
						atom[i].f.x = *mempointer; mempointer++;
						atom[i].f.y = *mempointer; mempointer++;
						atom[i].f.z = *mempointer; mempointer++;
					}
					delete [] tempMemory;
					return 0;
				}
			}
		}
	}
	return -1;
}




SnapShot readOldPickle(
	const std::string &filename,
	const std::string &name
)
{
	SnapShot psp;
	psp.readOldPickle(filename,name);
	return psp;
}


SnapShotLibrary::SnapShotLibrary()
{
}

SnapShotLibrary::~SnapShotLibrary()
{
}

void SnapShotLibrary::push_back(const SnapShot& s)
{
	data.push_back(s);
}

size_t SnapShotLibrary::dataSize() const
{
	return data.size();
}

const SnapShot &SnapShotLibrary::getData(size_t snapShotNumber)
{
	return data[snapShotNumber];
}


SnapShotOutTra::SnapShotOutTra(SnapShotLibrary& lib, WorkSpace& wspace):
	OutputTrajectory( wspace )
{
	m_lib = &lib;
	m_wspace = &wspace;
}

SnapShotOutTra::~SnapShotOutTra()
{
}

SnapShotOutTra* SnapShotOutTra::clone() const
{
	return new SnapShotOutTra(*this);
}

int SnapShotOutTra::create()
{
	return 0;
}

int SnapShotOutTra::append()
{
	m_lib->push_back ( m_wspace->save() );
	return 0;
}

