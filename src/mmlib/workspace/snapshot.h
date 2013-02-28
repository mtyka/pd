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

#ifndef __SNAPSHOT_H
#define __SNAPSHOT_H

// Essential Headers
#include <string>
#include <vector>
#include "stdio.h"
#include "fileio/outtra.h" // OutputTrajectory base class
#include "workspace/workspace.fwd.h"
#include "system/molecule.fwd.h"

//-------------------------------------------------
/// \brief Passive class that holds position velocity and forces of a particle.
/// \details Class used exclusively by SnapShot
/// \author Mike Tyka
class PD_API SnapShotAtom
{
 public:
	Maths::dvector p; ///< position
	Maths::dvector f; ///< force
	Maths::dvector v; ///< velocity
};


//-------------------------------------------------
//
/// \brief  Class that works in conjunction with WorkSpace and is used to hold a SnapShot of a simualation
///
/// \details  SnapShot represents a set of coordinate with a few associated pieces of data such as 
///           potential energy, and some custom data slots and most importantly periodic space vectors. 
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo -Do we keep the custom data ? Do we keep the epot ?
///
class PD_API SnapShot
{
public:
	/// constructor - PhaseSpacePoint has a minimum of 1 atom.
	SnapShot(size_t _natoms = 1);

	/// destructor - deletes atom array
	~SnapShot();

	/// copy constructor
	SnapShot(const SnapShot &newpsp);

	/// same as assignment operator - makes a deep copy
	void setTo(const SnapShot &newpsp);

	/// resets as many atom positions as possible - i.e.
	/// it allows assignement of technically incompatible PhaseSpacePoints
	void setToPartial(const SnapShot &newpsp);

#ifndef SWIG
	/// assignment operator
	SnapShot &operator= (const SnapShot &psp);
#endif

	size_t memuse(int level){
		level--;
		return sizeof(*this)+natoms*sizeof(SnapShotAtom);
	}

	// maybe make these external friend functions a la operators ?
	double dRMSFrom(const SnapShot &psp2) const; 
	double cRMSFrom(const SnapShot &psp2) const; 
	double cRMSFrom(
		const SnapShot &psp2, 
		const MoleculeBase& sys, 
		const std::string& name
		) const; 

	int insertPSPfragment(
		const SnapShot &psp2, 
		int startir, 
		int endir, 
		const MoleculeBase &wspace
		);
	
#ifndef SWIG
	int  readMIME(FILE *file,const std::string &name="");
	void writeMIME(FILE *file,const std::string &name="");
	int  readOldPickle(FILE *file,const std::string &name="");
#endif

	int  readMIME(const std::string &filename,const std::string &name="");
	void printMIME(const std::string &name="");
	int  writeMIME(const std::string &filename,const std::string &name="", const char *mode="w");
	int  appendMIME(const std::string &filename,const std::string &name="");
	int  readOldPickle(const std::string &filename,const std::string &name="");

	size_t nAtoms() const { return natoms; }

public:
	// Member Data
	/// atom positions, forces and velocities
	SnapShotAtom *atom;						
	
	/// number of atoms
	size_t natoms;

	// Extra information about the coordinate set saved

	/// potential energy
	double epot;														
	
	/// some data slots to attach custom information
	int status1, status2, status3, status4;    

	/// Periodic cell vectors. 
	Maths::dvector A,B,C; 
};

SnapShot readMIME(
	const std::string &filename,
	const std::string &name = ""
);

SnapShot readOldPickle(
	const std::string &filename,
	const std::string &name = ""
);







//-------------------------------------------------
//
///\brief The SnapShotLibrary base class creates a std::vector of SnapShots.
/// It can then add SnapShots to this array, get the size of the array and
/// get the SnapShot held at any position in the array.
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
class SnapShotLibrary
{
public:

	/// constructor
	SnapShotLibrary();

	///destructor
	virtual ~SnapShotLibrary();

	/// add a SnapShot to the std::vector data, a protected member variable
	void push_back(const SnapShot& s);
	
	/// get the size of the std::vector data, a protected member variable
	size_t dataSize() const;

	/// get the SnapShot held at data[snapShotNumber]
	const SnapShot &getData( size_t snapShotNumber );

protected:

	std::vector <SnapShot> data;
};






//-------------------------------------------------
//
///\brief The SnapShotOutTra class allows the user to append a SnapShot of
/// the current WorkSpace to a SnapShotLibrary.
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
class SnapShotOutTra: public IO::OutputTrajectory
{
public:
	/// constructor: takes a SnapShotLibrary and a workspace as arguments
	SnapShotOutTra(SnapShotLibrary& lib, WorkSpace& wspace);
	
	/// destructor
	virtual ~SnapShotOutTra();
	
	/// copy constructor
	virtual SnapShotOutTra *clone() const;

	/// empty create function: returns 0
	virtual int create();

	///\brief append function: take a SnapShot of the current WorkSpace and
	/// add it to the SnapShotLibrary
	virtual int append();

private:
	/// make a pointer to a SnapShotLibrary to be used in the append function
	SnapShotLibrary *m_lib;

	/// make a pointer to a WorkSpace to be used in the append function
	const WorkSpace *m_wspace;
};

#endif

