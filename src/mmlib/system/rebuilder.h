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

#ifndef __REBUILDER_H
#define __REBUILDER_H

#include "pickers/basicpickers.h"

class MoleculeBase;
typedef MoleculeBase_Array < Particle >  ParticleStore;

// -------------------------------------------------------------------------------------
// Note the rebuilders in this file deal with MoleculeBase. They are normally called by 
// ImportFileBase derived classes and therefore have no knowledge of a WorkSpace and are
// NOT allowed to use protocols. Special protocols that can set isRebuildRequired() are
// therefore needed if forcefields etc. are to be used in a rebuild process. Those special
// classes should NOT be placed in this file.

// ---------------------
//  Rebuilder Functions
// ---------------------

int findLevelOne( int i, const ParticleStore & atom, std::vector < int >& AlignAtom_L1, int iTemplateRangeStart = -1, int iTemplateRangeEnd = -1 );
int findLevelTwo( const ParticleStore & atom, std::vector < int >& AlignAtom_L1, std::vector < int >& AlignAtom_L2, int iTemplateRangeStart = -1, int iTemplateRangeEnd = -1 );
bool findLevelThree( const ParticleStore & atom, std::vector < int >& AlignAtom_L1, std::vector < int >& AlignAtom_L2, int& AlignAtom_L3, int iTemplateRangeStart = -1, int iTemplateRangeEnd = -1 );

/// rebuild all atoms for which the isRebuildRequired flag is set
bool rebuildMissingAtoms( MoleculeBase& mol, Verbosity::Type verbose );
bool rebuildMissingAtoms( MoleculeBase& mol, int ir, Verbosity::Type verbose );

bool rebuildAndPolymerise(MoleculeBase& mol, int start, int end);

// -------------------
//  Rebuilder Classes
// -------------------

//-------------------------------------------------
//
/// \brief Rebuilder Abstract Base class 
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class RebuilderBase
{
public:
	RebuilderBase()
	{ 
	};
	virtual bool invokeBuild(MoleculeBase& mol, Verbosity::Type verbose) = 0;
};


//-------------------------------------------------
//
/// \brief   Does nothing
///
/// \details Does exactly what is says on the tin. Nothing. 
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class NullRebuilder : public RebuilderBase
{
public:
	NullRebuilder(){};
	virtual bool invokeBuild(MoleculeBase& mol, Verbosity::Type verbose){ return true; };
};


//-------------------------------------------------
//
/// \brief   Places all unknown atoms on the origin
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class ToOriginRebuilder : public RebuilderBase
{
public:
	ToOriginRebuilder() : m_FlagBuilt(false) {};
	ToOriginRebuilder(bool _FlagBuilt) : m_FlagBuilt(_FlagBuilt) {};

	virtual bool invokeBuild(MoleculeBase& mol, Verbosity::Type verbose); ///< Invoke the builder. If setAsBuilt is flagged, atoms moved to 0,0,0 will be marked as built.

	void setFlagAsBuilt( bool _Flag ) { m_FlagBuilt = _Flag; }
	bool getFlagAsBuilt() const { return m_FlagBuilt; }

private:
	bool m_FlagBuilt; ///< Should atoms set to (0,0,0) be counted as 'built' (i.e. isRebuildRequired() set to true). Default to false.
};


//-------------------------------------------------
//
/// \brief   Builds missing atoms using the definitions in the forcefield
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class MissingAtomRebuilder : public RebuilderBase
{
public:
	MissingAtomRebuilder(){};
	virtual bool invokeBuild(MoleculeBase& mol, Verbosity::Type verbose);
};


//-------------------------------------------------
//
/// \brief   Builds missing atoms using the definitions in the forcefield
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class MainchainHydrogenRebuilder : public RebuilderBase
{
public:
	MainchainHydrogenRebuilder(){};
	virtual bool invokeBuild(MoleculeBase& mol, Verbosity::Type verbose);
	bool invokeBuild(MoleculeBase& mol, Verbosity::Type verbose, const PickResidueBase& _picker, int iTemplateRangeStart = -1, int iTemplateRangeEnd = -1);
	bool invokeBuild(MoleculeBase& mol, Verbosity::Type verbose, const PickResidueBase& _picker, const PickAtomRange& _TemplateRange );
};


//-------------------------------------------------
//
/// \brief   Represents a base for classes that own their own overridable rebuilder.
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class RebuilderOwner
{
public:
	RebuilderOwner();

	void setRebuilder( RebuilderBase& _rebuilder );
	void setRebuilderDefault();
	void disableRebuild();
	RebuilderBase& getRebuilder();

private:
	RebuilderBase* m_Rebuilder;
	MissingAtomRebuilder m_DefaultBuilder;
	ToOriginRebuilder m_DisabledBuilder;
};

#endif

