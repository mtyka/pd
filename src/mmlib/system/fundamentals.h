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

#ifndef __FUNDAMENTALS_H
#define __FUNDAMENTALS_H

#include <vector>
#include "primitives.h" // IndexPair base class

class MoleculeBase;
class MoleculeDefinition;
class WorkspaceCreatorBase;






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
class PD_API IndexPairDistance: public IndexPair
{
 public:
	IndexPairDistance(){};
	IndexPairDistance(int _i, int _j, double _l): IndexPair(_i,_j),l(_l) {}
	double l; ///< length/distance parameter
};







//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details 
/// The structure CovalentAtom is used for all bonded interaction and
/// stores information about any partner atoms within a molecule that
/// are restraint in any way comprising 1,2 (bonds), 1,3 (angles) and
/// 1,4 (torsions) atom pairs
/// for any atom x:
/// bond: x-i // sorry for awkward notation ...
/// angle x-i2-i
/// torsion x-i3-i2-i
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API CovalentAtom
{
public:
	CovalentAtom() : i(-1), i2(-1), i3(-1), d(0.0)
	{
	}
	int i; ///< index of partner atom (in 1,2 1,3 and 1,4 atoms)
	int i2; ///< index of pre-partner atom (only 1,3 and 1,4)
	int i3; ///< index of pre-pre-partner atom (only for 1,4 atoms);
	double d; ///< equilibrium distance of x to i
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
class PD_API Particle
{
public:
	friend class MoleculeBase;
	friend class WorkspaceCreatorBase;

	Particle();
	Particle(
		const std::string &_rawname,
		const std::string &_parentname
	);

	Particle(const Particle& _newpart){
		(*this) = _newpart;
	}
#ifndef SWIG
	const Particle& operator= ( const Particle& _newpart );
#endif
	// ---------------------------------------------
	// Parameters that should be readonly to the
	// vast majority of client classes.
	// Particle should be passed as const
	// ---------------------------------------------

	const MoleculeDefinition *parent; // links back to a parent - or a nullparameter.

	// Organisational properties
	int i;                      ///< index in array or whatever this particle is held in
	int ir;                     ///< index of residue in parent array
	int imol;                   ///< index of molecule in parent array
	int igroup;                 ///< index of charge group

	// --------------------------
	// Coordinate stuff:
	// --------------------------
private:

	/// pointer to coordinate. If Particle exists by itself or as part of a Molecule,
	/// the this pointer will point to the internal position store _pos. If part of
	/// a Workspace, then this position will point to wspace.cur[i].p - a linear array
	/// of atom positions specific to Workspace and there for fast access.
	Maths::dvector* _posPointer;         

	void setPosPointer(Maths::dvector& newpos){ _posPointer = &newpos; };
	void resetPosPointer(){ _posPointer = &_pos; };

	/// internal coordinate - used for independent particles and those which are part of a Molecule
	Maths::dvector _pos;         
	
	/// reference coordinate - Reference coordinates define the state that is used for RMS calculations
	Maths::dvector _posRef;      

	/// Geometric coordinate - position definitions with correct relative bond angles, lengths and connectivity
	Maths::dvector _posGeom;     

public:
	
	inline       Maths::dvector& pos()           { return *_posPointer; };
	inline const Maths::dvector& pos()     const { return *_posPointer; };
	inline       Maths::dvector& posRef()        { return _posRef;      };
	inline const Maths::dvector& posRef()  const { return _posRef;      };
	inline       Maths::dvector& posGeom()       { return _posGeom;     };
	inline const Maths::dvector& posGeom() const { return _posGeom;     };


public:
	// --------------------------
	// Build-in properties
	// --------------------------

	// Covalency
	std::vector<CovalentAtom> cov12atom;
	std::vector<CovalentAtom> cov13atom;
	std::vector<CovalentAtom> cov14atom;

	// Names
	std::string rawname;        ///< name of individual particle
	std::string pdbname;        ///< PDB name
	std::string type_name;      ///< name of fundamental Type (optional)
	std::string parentname;     ///< Full name of Parent (residue or molecule)
	std::string parentl3name;   ///< Three letter name of Parent (residue or molecule)
	char parentletter;          ///< Single letter name of Parent (residue or molecule)

	// Fundamental physical properties
	                            
	int Z;          ///< atomic number
	double mass;    ///< mass in kg
	double radius;  ///< VdW radius in Angstrom
	double epsilon; ///< Epsilon (VdW well depth) in J^0.5
	double charge;  ///< Formal/Partial charge	
	double epot;    ///< The current potential energy of this atom
	int FFType;     ///< forcefield Type number	

	// --------------------------
	// Custom property calls
	// --------------------------

	/// Custom Property - add one to this particle
	void addProperty(const std::string &_name,const std::string &_data);

	/// Custom Property - obtain the data string
	std::string getCustomProperty_String(const std::string &_name) const;

	/// Custom Property - obtain an integer
	int getCustomProperty_Int(const std::string &_name) const;

	/// Custom Property - obtain a double
	double getCustomProperty_Double(const std::string &_name) const;

	// --------------------------
	// Helper info function calls
	// --------------------------

	void detail() const;
	void info( 
		Verbosity::Type type = Verbosity::Normal, 
		int padNumbers = 4, 
		bool endLine = true ) const;

	size_t memuse(int level) const;

	// -------------
	// Flagset calls
	// -------------

	// Flag gets: General Properties
	inline bool isValid()           const {return (flags & flag_Valid)           == 0 ? false : true;}
	inline bool isUsed()            const {return (flags & flag_Used)            == 0 ? false : true;}
	inline bool isDummy()           const {return (flags & flag_Dummy)           == 0 ? false : true;}
	inline bool isStatic()          const {return (flags & flag_Static)          == 0 ? false : true;}
	inline bool isRebuildRequired() const {return (flags & flag_RebuildRequired) == 0 ? false : true;}
	inline bool isKnownStructure()  const {return (flags & flag_KnownStructure)  == 0 ? false : true;}
		
	// Flag gets: Chemical/Structural Nature
	inline bool isHydrogen()        const {return (flags & flag_Hydrogen)        == 0 ? false : true;}
	inline bool isBackbone()        const {return (flags & flag_Backbone)        == 0 ? false : true;}
	inline bool isCoreBackbone()    const {return (flags & flag_CoreBackbone)    == 0 ? false : true;}
	inline bool isCAlpha()          const {return (flags & flag_C_alpha)         == 0 ? false : true;}
	inline bool isSideChain()       const {return !isBackbone();}

	// Flag gets: Status/State
	inline bool isMoved()           const {return (flags & flag_Moved)           == 0 ? false : true;}


public: // TODO: may soon be private ?

	// -------------------------------------------
	// Functions that can change the state of the
	// Particle. These should only be available to
	// setup functions.
	// -------------------------------------------

	// PRIVATE flag gets: <- should be set via FFPS and then *left*

	// Flag sets: General Properties
	inline void setValid( bool _set )           { _set ? flags |= flag_Valid           : flags &= ~flag_Valid;  }
	inline void setUsed( bool _set )            { _set ? flags |= flag_Used            : flags &= ~flag_Used;  }
	inline void setDummy( bool _set )           { _set ? flags |= flag_Dummy           : flags &= ~flag_Dummy;  }
	inline void setKnownStructure( bool _set )  { _set ? flags |= flag_KnownStructure  : flags &= ~flag_KnownStructure;  }

	// Flag sets: Chemical/Structural Nature
	inline void setHydrogen( bool _set )     { _set ? flags |= flag_Hydrogen     : flags &= ~flag_Hydrogen;  }
	inline void setBackbone( bool _set )     { _set ? flags |= flag_Backbone     : flags &= ~flag_Backbone;  }
	inline void setCoreBackbone( bool _set ) { _set ? flags |= flag_CoreBackbone : flags &= ~flag_CoreBackbone;  }
	inline void setCAlpha( bool _set )       { _set ? flags |= flag_C_alpha      : flags &= ~flag_C_alpha;  }

public: 	
	// PUBLIC flag gets: Status/State -- these must be public so that forcefields etc can modify them
	inline void setStatic( bool _set )          { _set ? flags |= flag_Static          : flags &= ~flag_Static;  }
	inline void setMoved( bool _set )           { _set ? flags |= flag_Moved           : flags &= ~flag_Moved;  }
	inline void setRebuildRequired( bool _set ) { _set ? flags |= flag_RebuildRequired : flags &= ~flag_RebuildRequired;  }

private:

	// -----------------------
	// Core private functions!
	// -----------------------

	void init();
	void CopyParamsFrom(const Particle& _CopyParticle);

	/// Adds an offset to the indeces stored in the covalent atom arrays
	/// This function is used by the setup routines when adding multiple molecules
	/// into one linear array in order to maintain the correct connectivity
	void offsetInternalIndices( int offset );

	/// Adds to the covelent 1:2 list only - calc1314BondOrders() 
	/// **MUST** now be called on the parent container! Hence this function is private.
	void add12OnlyCovbond(int iatom2);

	// ---------------------
	// Core private members!
	// ---------------------

	/// indices to custom properties stored in the ValueStore (see valuestore.h)
	std::vector<size_t> m_CustomProp;

	// The internal workings of a flagset should really be private. Code can modify the flagset by
	// utilising the isXXX() and setXXX() function calls below.
	/// Flags - these set various boolean properties of the particle
	unsigned int flags;

	// General Properties
	static const unsigned int flag_Valid;           ///< All good or something wrong ?
	static const unsigned int flag_Used;            ///< Used or not ?
	static const unsigned int flag_Dummy;           ///< Real atom or dummy particle ?
	static const unsigned int flag_Static;          ///< should it be moved ?
	static const unsigned int flag_KnownStructure;  ///< Is this atom position of known structural data. 
	static const unsigned int flag_RebuildRequired; ///< Does this atom need rebuilding ?

	// Chemical/Structural Nature
	static const unsigned int flag_Hydrogen;        ///< Heavy atom or not ?
	static const unsigned int flag_Backbone;        ///< backbone or sidechain ?
	static const unsigned int flag_C_alpha;         ///< Alpha Carbon ?
	static const unsigned int flag_CoreBackbone;    ///< One of the main heacy backbone atoms?

	// State flags
	static const unsigned int flag_Moved;           ///< Has atom moved ?
};







//-------------------------------------------------
//
/// \brief  Holds whole system, a molecule, a residue or whatever (depending on usage)
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
class PD_API ParticleGroup
{
	// Nameing
	char letter; // Letter identification (if present)
	std::string name; // Full Name
	std::string l3_name; // Three Letter Name (for PDB)
	std::vector<Particle> atom; ///< Component atoms
};






//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
/// The contents of this class are public, however it is only changed through MoleculeBase
/// and associated control classes. Access elsewhere is performed through a const& and therefore
/// publicness here is fine.
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API Residue
{
public:
	size_t memuse(int level){ level--;  return sizeof(*this); }

	void info( Verbosity::Type vType = Verbosity::Silent ) const;

	char letter;

	/// Residue index
	int ir;

	// Atom indices
	int ifirst, ilast; // first & last particle belonging to this residue

	// Protein Specific indices
	int iCA, iHA, iN, iC, iO, iH, iCB; // indices of these key atoms

	// Nucleic Acid Specific indices
	// int iC1, iC2, iC3, iC4, iC5; // indices of these key atoms

	const MoleculeDefinition *param;	// Link to parent properties

	size_t size() const;

private:
	void addOffset(int offset);
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
class PD_API MoleculeRange
{
public:
	MoleculeRange() : ifirst(-1), ilast(-1), irfirst(-1), irlast(-1), chainID(' ')
	{
		name = "";
	}

	inline size_t natoms(){ return (ilast-ifirst + 1); }
	inline size_t nResidues(){ return (irlast-irfirst + 1); }
	size_t memuse(int level) { level--; return sizeof(*this); }
	std::string name; ///< The name of this sub-molecule
	int ifirst; ///< first particle belonging to this sub-molecule
	int ilast; ///< last particle belonging to this sub-molecule
	int irfirst; ///< first residue belonging to this sub-molecule
	int irlast; ///< last residue belonging to this sub-molecule
	char chainID; ///< The chain ID of the sub-molecule
};

#endif
