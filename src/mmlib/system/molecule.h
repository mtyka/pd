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

#ifndef __MOLECULE_H
#define __MOLECULE_H

// System Includes
#include <stdlib.h>
#include <string>
#include <vector>

// MMLib Includes
#include "forcefields/ffparam.h"
#include "system/fundamentals.h"
#include "sequence/sequence.h"
#include "workspace/workspace.fwd.h"
#include "fileio/outtra.fwd.h"

// Class Forward-Declarations
class PD_API FFParamSet;
class PD_API System;
class PD_API MoleculeBase;
class PD_API Molecule;
class PD_API PickBase;

#include "molecule_base_array.h"

//-------------------------------------------------
//
/// \brief Defines the base class for all molecule types.
///
/// \details Derived from this class will be Molecule, representing class Systems internal container
/// and also the imutable workspace. Once allocated by a workspace creator the workspace cannot change size. 
/// This means that cached (performance enhancing) itterators that point to the workspace will not become invalid
/// following their allocation. To ensure that the WorkSpace is not changed following allocation, the public interface
/// to allocation functions must be restricted to friend classes only.
///
/// \author Mike Tyka & Jon Rea 
///


class PD_API MoleculeBase
{

#ifndef SWIG
	friend bool operator==(const MoleculeBase &mol1, const MoleculeBase &mol2 ); 
#endif

public:
	////////////////////////////////////////////////////////////////////////////////
	// Constructor logic

	/// Default Constructor needed by standard containers.
	MoleculeBase(): ptr_ffps(&nullffps)
	{
		loadedparams = false;
		m_BondingMode = false;
		name = "Molecule";
	}

	/// Constructor must take a reference to an ffps
	MoleculeBase(const FFParamSet &_ptr_ffps)
		: ptr_ffps(&_ptr_ffps)
	{
		loadedparams = false;
		m_BondingMode = false;
		name = "Molecule";
	}

	virtual ~MoleculeBase()
	{
	}

	virtual size_t memuse(int level);

	////////////////////////////////////////////////////////////////////////////////
	// Info and printing functions

	void info() const;
	void detail() const;
	void detail( const PickBase &picker ) const;

	void printBondConnectivity() const;


	////////////////////////////////////////////////////////////////////////////////
	// File IO 

	/// writes a PDB to a file. 
	void printPDB( const std::string& targetFile);
	
	/// prints a PDB to the screen !
	void printPDB();

	////////////////////////////////////////////////////////////////////////////////
	// Inspectors / Accessors


	/// this can be used to access ffps through the molecule
	const FFParamSet &ffps() const { return *ptr_ffps; } 

	/// Chain info
	virtual char getChainID( size_t iAtom ) const = 0;

	/// Sequence info
	const Sequence::BioSequence& getSequence() const;

	/// The true mappings for this are sent through the internal BioSequence
	std::string getResName( size_t _ir ) const; 

	/// The true mappings for this are sent through the internal BioSequence
	std::string getFullResName( size_t _ir ) const; 

	/// Returns the number of atoms in the Molecule
	inline size_t nAtoms() const { return atom.size(); }

	/// Returns the number of residue in the Molecule
	inline size_t nResidues() const { return res.size(); }

	/// Returns the number of (sub)molecule in the molecule. Note that the class MoleculeBase can represent multiple chemical molecules
	/// i.e. for example all the waters of a system could be represented by one Molecule. 
	inline size_t nMolecules() const { return mol.size(); }



	////////////////////////////////////////////////////////////////////////////////
	/// \name  Coordinate Accessors
	/// These two functions allow derived classes to use separate atom coordinate
	/// arrays without loosing the functionality of the function in this class.
	/// Thus all the functions in this class PD_API that modify atom coordinates should use
	/// these functions rather than accessing atom directly. Efficientcy *should*
	/// be comaprable. Custom overloaded functions can be used in derived classes
	/// if efficiency is critical.
	/// \{

	virtual Maths::dvector &atomxyz(const size_t index) 
	{ 
		return atom.at(index).pos(); 
	}

	virtual const Maths::dvector &atomxyz(const size_t index) const 
	{ 
		return atom.at(index).pos();
	}

	/// These two functions expose the internal atoms publically for position modification, but do NOT
	/// allow the size of the atom array to change. This is important in the derived workspace class.
	inline Particle &getAtom(const size_t index) 
	{
		return atom[index];
	}
	inline const Particle &getAtom(const size_t index) const 
	{
		return atom[index];
	}
	inline Residue &getRes(const size_t index) 
	{
		return res[index];
	}
	inline const Residue &getRes(const size_t index) const 
	{
		return res[index];
	}
	inline MoleculeRange &getMol(const size_t index) 
	{
		return mol[index];
	}
	inline const MoleculeRange &getMol(const size_t index) const 
	{
		return mol[index];
	}

	
	// For doxygen:
	/// \}





	////////////////////////////////////////////////////////////////////////////////
	/// \name  System Property Accessors
	/// \{
	double getTotalMass() const;
	double getTotalCharge() const;

	void resetAllAtomMovedFlags();

	/// calculates the inertia tensor of the molecule
	void calcInertiaTensor(Maths::matrix3x3 & I) const;
	void printInertiaInfo() const;
	double calcRotationalPartition(double temp, unsigned SymNumber) const;
	int findParticleBy_ffname(int ir, const std::string &name) const;
	int findParticle(int ir, const std::string &name) const;
	// For doxygen:
	/// \}




	////////////////////////////////////////////////////////////////////////////////
	/// \name Particle Property Accessors 
	/// \{
	double calcAtomDistance(int i, int j) const;
	double calcAtomSqrDistance(int i, int j) const;
	double calcAtomAngle(int i, int a, int j) const;

	bool isMolStartResIndex( int _ResIndex ) const;
	bool isMolEndResIndex( int _ResIndex ) const;


	// For doxygen:
	/// \}

	////////////////////////////////////////////////////////////////////////////////
	/// \name Coordinate Analysis / Modification Functions
	/// \{

	/// \brief returns the largest absolute vector that encompasses all the atoms - i.e. it returns the
	/// largets value of fabs(x), fabs(y) and fabs(z) for the set of atoms in this molecule
	Maths::dvector getEncompassingVector() const;

	/// gets maximum coordinate in system (i.e. largest value of |x|,|y| and |z|
	Maths::dvector getMaximumVector() const;

	/// gets minimum coordinate in system (i.e. smallest value of |x|,|y| and |z|
	Maths::dvector getMinimumVector() const;

	// Ridid-Body centre information
	Maths::dvector getCentreOfMass() const;
	Maths::dvector getCentreOfGeometry() const;

	/// Returns the centre of geometry for atoms flagging not isRebuildRequired, isValid and not isDummy
	Maths::dvector getCentreOfAvailableGeometry() const; 
	Maths::dvector getCentreOfGeometry(unsigned i) const;
	Maths::dvector getCentreOfGeometry(unsigned ifirst, unsigned ilast) const;
	Maths::dvector getCentreOfGeometry(const PickBase& _picker) const;

	// Reference coordinate sets
	/// setToGeometry: This function sets the atom[].p positions to the atom[].posGeom()         
	void setToGeometry();  

	/// setToReference: This function sets the atom[].p positions to the atom[].posRef()
	void setToReference(); 

	/// setAsReference: This function sets the atom[].posRef() positions to the atom[].p
	void setAsReference();

	// RMS calculations
	virtual double calcCRMS_AllAtom( bool useOnlyKnownAtom = false ) const;
	virtual double calcCRMS_HeavyAtom( bool useOnlyKnownAtom = false ) const;
	virtual double calcCRMS_CA( bool useOnlyKnownAtom = false ) const;

	// Cartesian movement
	//void move(double  x, double y, double z){ moveParticles(0,atom.size(), Maths::dvector( x, y, z) ); }
	void move(const Maths::dvector & disp){ moveParticles(0,atom.size(), disp); }
	void moveParticles(size_t start, size_t end, const Maths::dvector & disp);
	void moveMolecule(int imol, const Maths::dvector & disp);
	void zeroCentreOfAvailableGeometry();
	void zeroCentreOfGeometry();
	void zeroCentreOfMass();
	void alignAlongPrincipalAxes();

	inline void rotate(const Maths::matrix3x3 & rmat)
	{ 
		rotateParticles(0,atom.size(),getCentreOfGeometry(), rmat); 
	}
	inline void rotate(const Maths::dvector & centre, const Maths::matrix3x3 & rmat)
	{ 
		rotateParticles(0,atom.size(),centre, rmat); 
	}
	void rotateParticles(size_t start, size_t end, const Maths::dvector & centre, const Maths::matrix3x3 & rmat);
	void rotateMolecule(int imol, const Maths::dvector & centre, const Maths::matrix3x3 & rmat);

	/// Bond rotation
	int rotateBond(int i, int j, double angle);

	// Backbone Torsion Angle calculations (Protein dependent)
	void calcResiduePhiPsi(int ir, double &fphi, double &fpsi) const;
	void calcResidueOmega(int ir, double &fomega) const;
	void setAllResiduePhiPsi(double nphi, double npsi);
	void setAllResidueOmega(double nomega);
	void setProteinAlphaHelix();
	void setProteinReverseHelix();

	void setResiduePhiPsi(int ir, double nphi, double npsi);
	void setResidueOmega(int ir, double nomega);

	// For doxygen:
	/// \}

public:
	std::string name;

	// This is termporarily public, needed by InTra_BTF prior to its full revamp to use genPolymer
	// Setters
	int setBioSequence( const Sequence::BioSequence& sequenceInformation );

protected:

	////////////////////////////////////////////////////////////////////////////////
	// Protected Setup Functions 

	// Configuration
	int checkParameters();
	int loadParameters(const FFParamSet &ffps);

	/// calls detectResidueBoundaries() and detectMoleculeBoundaries()  
	void detectSubBoundaries();     

	/// Detects residue boundaries and fills the 'res' array with data.
	void detectResidueBoundaries();  

	/// Detects (sub)molecule boundaries and fills the 'mol' array with data.
	void detectMoleculeBoundaries(); 

	// Covalency:
	/// 1) Uses the atom parameters to create bonds (cov12atom)
	int createCovalentStructure();  
	/// 2) uses the information in cov12atom to fill cov13atom and cov14atom
	int calc1314BondOrders();     
	/// 3) Checks the covalent structure defined in cov12atom is consistent	
	int checkCovalentStructure(); 	 

	////////////////////////////////////////////////////////////////////////////////
	// Functions that change MolculeBase contents
	//

	/// Only full molecules are allowed to append to the base class internal containers
	void append(const MoleculeBase &appendmol); 
	
	/// Only full molecules are allowed to append to the base class internal containers
	void append(const System &appendsys); 
	
	void append(const Particle &newatom) { atom.push_back(newatom); }

	/// Begin bonding allows you to add multiple bonds without triggering internal bond order recalculation
	/// until all the bonds have been added. This is important if many bonds are being added to the structure.
	void beginBonding();

	/// End bonding mode initiated when beginBonding() was called. The internal bond order arrays will now be calculated
	void endBonding();

	void addBond(size_t indexi, size_t indexj );
	void addDisulphide(size_t indexi, size_t indexj );

	////////////////////////////////////////////////////////////////////////////////
	// protected Data

	Sequence::BioSequence m_Sequence;

private:

	////////////////////////////////////////////////////////////////////////////////
	// private Data

	bool loadedparams;
	bool m_BondingMode;

	// A MoleculeBase *must* be connected to a parameter set 
	// even if it's a default, empty one. 
	const FFParamSet *ptr_ffps;

public: 

#ifndef SWIG
	ParticleStore      atom; 
	ResidueStore       res;
	MoleculeRangeStore mol;
#endif
};


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
class Molecule : public MoleculeBase
{
	// Generate polymer makes new molecules - let it see our private parts.
	friend int GeneratorCore( 
		const Sequence::BioSequence &bioSeq, 
		Molecule &newmol,												
		const FFParamSet &ffps,	
		bool _Polymerise,
		Verbosity::Type OutputLevel);														

	// Generare polymer makes new molecules - let it see our private parts.
	friend int loadMolecule( 
		const std::string &MolName, 
		Molecule &newmol, 
		const FFParamSet &ffps);

	// Generare polymer makes new molecules - let it see our private parts.
	friend int loadMolecule(
		const Sequence::ResidueInfo& _Res,
		Molecule &newmol, 
		const FFParamSet &ffps);

	friend class System;

public:
	Molecule(const FFParamSet &_ptr_ffps): MoleculeBase(_ptr_ffps), chainID(' ')
	{
	}

	Molecule(const FFParamSet &_ptr_ffps, char _chainID): MoleculeBase(_ptr_ffps), chainID(_chainID)
	{
	}

	// Public append layer for full molecules *only* 
	// (WorkSpace is not allowed to append following allocation - its imutable)
	inline void append(const MoleculeBase &appendmol) { MoleculeBase::append(appendmol); }
	inline void append(const System &appendsys) { MoleculeBase::append(appendsys); }
	inline void append(const Particle &newatom) { MoleculeBase::append(newatom); }
	inline void beginBonding() { MoleculeBase::beginBonding(); }
	inline void endBonding() { MoleculeBase::endBonding(); }
	inline void addBond(size_t indexi, size_t indexj) { MoleculeBase::addBond(indexi, indexj); }
	inline void addDisulphide(size_t indexi, size_t indexj) { MoleculeBase::addDisulphide(indexi, indexj); }
	inline void addParticle(const Particle &newatom){ atom.push_back(newatom); }

	////////////////////////////////////////////////////////////////////////////////
	// File IO 

	/// This saves the current state of the system in a file format of the users choice
	///
	///  For example  mysystem.save( OutputFile_PDB("mypdbfile") ) 
	///  saves a PDB file.
	void save( IO::OutputFile &_output );

	// Chain info
	virtual char getChainID( size_t iAtom ) const { return chainID; }
	char getChainID() const { return chainID; }

private:
	char chainID;
};

#ifndef SWIG
bool operator==(const MoleculeBase &mol1, const MoleculeBase &mol2 ); 
#endif



#endif

