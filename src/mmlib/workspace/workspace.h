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

// (c) Michael Tyka & Jon Rea
/// \file workspace.h 
/// \brief Contains all the workspace main header definitions
/// \details Details on the workings and general properties of this file
/// \author Michael Tyka & Jon Rea
///

#ifndef __WORKSPACE_H
#define __WORKSPACE_H

// System includes
#include <vector>

// Forward declarations
#include "system/system.fwd.h"
#include "system/workspacecreator.fwd.h"
#include "workspace/space.fwd.h"
#include "workspace/neighbourlist.fwd.h"
#include "workspace/bondorder.fwd.h"
#include "workspace/rotbond.fwd.h"
#include "workspace/workspace.fwd.h"
#include "fileio/tra.fwd.h" 

// Used Types
#include "system/molecule.h" // Base class
#include "workspace/hamiltonian.h" // Member data
#include "workspace/snapshot.h" // Member data
#include "fileio/outtra.fwd.h" // Member data


//-------------------------------------------------
//
/// \brief Central simulation class 
///
/// \details 
/// Workspace is responsible for holding a simulation state. It
/// provides a linear atom array where each atom has properties such as
/// position, previous position, force, old force, velocity etc.. which the
/// Protocol s and Forcefield s can modify.
/// Workspace  thus acts as a scratch space for simulations.
/// Workspace  also
/// provides functions to calculate neighbour lists which the
/// Forcefield s can use to calculate non-bonded forces. Further a
/// large variety of structure modifying functions are available such as
/// atom movements, bond rotations, etc..
/// and other transformations as well as basic measuring functions (e.g. to
/// calculate inter-atomic distances or torsion angles). Finally,
/// Workspace contains a Trajectory container class
/// which allows various output trajectory formats (represented by
/// respective classes) to be 'loaded' and thus capture structural snapshots
/// of the simulation.
/// A Workspace is always created from a System where the
/// extracted atoms remain linked to the original atoms in the
/// System. Thus, following a simulation,
/// the original positions of the atoms in the System  can be
/// updated and other, different Workspaces created. For example a
/// simulation set up is conceivable in which firstly only the loops in a
/// protein are simulated separately (each with only its immediately
/// surrounding atoms), after which the entire protein is to be minimised or
/// relaxed. This is easily implementable in our setup because of the
/// two-way interaction between System and Workspace.
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo Many things, but the core outline of this class is set 
///
/// \bug BUGS? probably 
///
class PD_API WorkSpace : public MoleculeBase
{
	friend class PD_API WorkspaceCreatorBase;

public:
	// -----------------------
	//  Constructor Logic
	// -----------------------
	WorkSpace(const WorkspaceCreatorBase &wc);
	WorkSpace(System &sysspec);
	WorkSpace(const WorkSpace &cloneMe);
	virtual ~WorkSpace();

private:
	System *parentspec;

	void CommonInit();
	void reinitAll();
	void Allocate(int _natom);
	void parseGroups(); 

	RotBond           *default_rotbond;
	NeighbourListBase *default_nlist;
	BondOrder         *default_bondorder;

	/// the default boundary is InfiniteSpace
	Space             *default_boundary;    

	RotBond           *ptr_rotbondlist;
	NeighbourListBase *ptr_nlist;
	BondOrder         *ptr_bondorder;
	Space             *ptr_boundary;

	long m_CheckSum;

public:
	long getCheckSum() { return m_CheckSum; }

	virtual char getChainID( size_t iAtom ) const;

	// Non_constant Wrappers for the Operators
	inline       RotBond           &rotbond()       { return *ptr_rotbondlist; }
	inline       NeighbourListBase &nlist()         { return *ptr_nlist; }
	inline       Space             &boundary()      { return *ptr_boundary; }

	// Constant wrappers for the Operators
	inline const RotBond           &rotbond()   const { return *ptr_rotbondlist; } 
	inline const NeighbourListBase &nlist()     const { return *ptr_nlist; }
	inline const BondOrder         &bondorder() const { return *ptr_bondorder; }
	inline const Space             &boundary()  const { return *ptr_boundary; }

	// setting operators


	/// \brief  BRIEF DESCRIPTION
	///
	/// \details DETAILED DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	/// \param[in]  Name   INPUT PARAMETER DESCRIPTION
	/// \param[out] Name   OUTPUT PARAMETER DESCRIPTION
	/// \return            REtURN VALUE DESCRIPTION 
	///
	/// \author AUTHOR NAME(S)
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	void setRotatableBondList(RotBond *new_rotbond);
	void setNeighbourList(NeighbourListBase *new_nlist);
	void setSpace(Space *new_boundary);


	/// \brief Accessor for atom coordinates
	/// \details These functions from MoleculeBase allow derived classes to use separate atom coordinate
	/// arrays without loosing the functionality of the function in that class.
	/// In WorkSpace these will point to cur.atom.p rather than atom.pos such that the working
	/// Coordinates are modified instead.
	/// 
	/// \param[in]  index  The serial atom number of the atom
	/// \return     A reference to the atom position

	virtual Maths::dvector &atomxyz(
		/// The atom index requested
		const size_t index  
	) 
	{ 
		return cur.atom[index].p; 
	}

	/// \brief Const accessor for atom coordinates
	/// \details These functions from MoleculeBase allow derived classes to use separate atom coordinate
	/// arrays without loosing the functionality of the function in that class.
	/// In WorkSpace these will point to cur.atom.p rather than atom.pos such that the working
	/// Coordinates are modified instead.
	/// 
	/// \param[in]  index  The serial atom number of the atom
	/// \return     A reference to the atom position
	virtual const Maths::dvector &atomxyz(const size_t index) const 
	{ 
		return cur.atom[index].p; 
	}

	// Needs a decision on what to do with it. 
	// Hamiltonian kinda needs to stay - maybe in a slimmed down form.

	/// Energy of the current conformer
	mutable Hamiltonian ene; 

	// Base class overrides to set cRMS in the Hamiltonian

	/// calc the cRMS - the result is also places in Hamiltonian
	virtual double calcCRMS_AllAtom( bool useOnlyKnownAtom = false ) const; 

	/// calc the cRMS - the result is also places in Hamiltonian
	virtual double calcCRMS_HeavyAtom( bool useOnlyKnownAtom = false ) const; 

	/// calc the cRMS - the result is also places in Hamiltonian
	virtual double calcCRMS_CA( bool useOnlyKnownAtom = false ) const; 

	/// \brief This variable is set by some protocols to indiceted their step number
	/// \details
	/// This variable can be set by protocols to broadcast to forcefields and
	/// Other workspace operators which support multistepping which step
	/// of a calculation the system is at. Setting Step to 0 (it's default
	/// value) will cause everyone the recalculate properties fully while
	/// higehr values may mean that the components use less expensive, update-like
	/// operations (e.g. instead of recalculating the neighbor list, it'd merely
	/// recalculate the distnaces in the neighbour list)
	int Step;
	
	/// \brief  sets all energies to 0 and all force vectors to (0,0,0) 
	void zeroForces();

	/// \brief  prints a list of forces - mainly for debugging purposes
	void printForces(); 

	/// \brief  multiply velocities by some factor
	void scaleVelocities(double factor);

	/// \brief  returns volume of boundary or -1 if the boundary is infinite [Angstrom^3]
	double getVolume() const;

	/// \brief  returns density inside boundary or -1 if the boundary is infinite [g/cc]
	double getDensity() const;

	/// \brief scale the system by moving molecules, but don't scale the molcules intra molecularly
	void scaleSystemMolecular(double factor);

	/// scale the system by multiplying every atom position by factor
	void scaleSystem(double factor);

	/// moves any stray molecules back into the box 
	void cleanSpace();
	
// Trajectory output functions ---------------------------------------
public:
	/// adds a OutputTrajectoryFile to the outtra container
	/// - this function is merely a shortcut for use from python
	/// (since SWIG doesnt seem to support double dereferencing
	/// of the wspace.outtra.add(...) Type
	void addTra(IO::OutputTrajectory &newtra);

#ifndef SWIG
	/// For use in the c++ interface
	void addTraWithOwnership(IO::OutputTrajectory* newtra);
#endif

	/// Another shortcut function - adds a pd trajectory to this workspace
	/// - in most situations all one wants to do is to save trajectory entries in a file and that's it - this
	/// function makes that very simple:
	/// myworkspace.addStdTra("myfilestem");
	IO::OutTra_BTF& addStdTra(const std::string &_filestem);
	IO::OutTra_BTF& addExtdTra(

		/// The filestem to be written to
		const std::string &_filestem, 

		/// Should forcevectors and Phi/Psis be included too?
		bool _MoreModes = false 
		);
	
	/// Output trajectory container - trajectories can be added&removed
	/// Protocols can call the functions below to append trajectory snapshots
	IO::OutputTrajectoryContainer outtra;
	
	/// This saves the current state of the WorkSpace in a file format of the users choice
	///  For example  myworkspace.save( OutputFile_PDB("mypdbfile") ) 
	///  saves a PDB file.
	void save( IO::OutputFile &_output );

// -------------------------- Member Data ---------------------------------------
public:

	/// \name Member Data
	/// @{
	/// CURRENT Atom positions, forces and velocities
	SnapShot	cur; 
	
	/// OLD Atom positions, forces etc. (from previous Step)
	SnapShot	old; 

	/// Save snapshot of workspace - returns a copy of cur
	SnapShot save() const;

	/// Load a snapshot into the cur array (only if the atom numbers match)
	void load(const SnapShot &psp); 

	/// Load a snapshot into the cur array even if the atom numbers dont match
	void load_forced(const SnapShot &psp); 
	/// @}

	// This should be replaced with an array of atom IDENTIFICATION Numbers -
	// infact these would aready be in the atom array;
	
	/// Loads the workspacepositions positions back into the parent system description
	int updateSystemPositions();
	
	size_t memuse(int level);

// should be private:
	// Variables for use with reassigning the workspace back into the system

	/// temporary structure to hold the back indices
	std::vector<int> isysmol;   
	std::vector<int> isysatom;

  std::vector<unsigned>  group_index;
};

#endif
