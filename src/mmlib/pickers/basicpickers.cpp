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
#include "tools/vector.h"
#include "system/fundamentals.h"
#include "system/molecule.h"
#include "basicpickers.h"

PickNothing* PickNothing::clone() const
{
	return new PickNothing();
}

bool PickNothing::matches( const Particle& particle ) const
{
	return false;
}

bool PickNothing::matches( const Residue& _res ) const
{
	return false;
}

PickEverything* PickEverything::clone() const
{
	return new PickEverything();
}

bool PickEverything::matches( const Particle& particle ) const
{
	return true;
}

bool PickEverything::matches( const Residue& _res ) const
{
	return true;
}

PickAllParticles* PickAllParticles::clone() const
{
	return new PickAllParticles();
}

bool PickAllParticles::matches( const Particle& particle ) const
{
	return true;
}



PickAllAtoms* PickAllAtoms::clone() const
{
	return new PickAllAtoms();
}

bool PickAllAtoms::matches( const Particle& particle ) const
{
	return !particle.isDummy();
}



PickHeavyAtoms* PickHeavyAtoms::clone() const
{
	return new PickHeavyAtoms(*this);
}

bool PickHeavyAtoms::matches( const Particle& particle ) const
{
	return !particle.isHydrogen() && !particle.isDummy();
}


PickResidue* PickResidue::clone() const
{
	return new PickResidue(*this);
}

bool PickResidue::matches( const Particle& particle ) const
{
	return (particle.ir == m_iRes);
}

bool PickResidue::matches( const Residue& _res ) const
{
	return _res.ir == m_iRes;
}


PickMolecule* PickMolecule::clone() const
{
	return new PickMolecule(*this);
}

bool PickMolecule::matches( const Particle& particle ) const
{
	return (particle.imol == m_iMol);
}


PickMoleculeRange* PickMoleculeRange::clone() const
{
	return new PickMoleculeRange(*this);
}

bool PickMoleculeRange::matches( const Particle& particle ) const
{
	return ((particle.imol >= m_iMol_start)&&(particle.imol <= m_iMol_end));
}



PickCAAtoms * PickCAAtoms::clone() const
{
	return new PickCAAtoms(*this);
}

bool PickCAAtoms::matches( const Particle& particle ) const
{
	return particle.isCAlpha();
}



PickBackbone * PickBackbone::clone() const
{
	return new PickBackbone(*this);
}

bool PickBackbone::matches( const Particle& particle ) const
{
	return particle.isBackbone();
}



PickCoreBackbone* PickCoreBackbone::clone() const
{
	return new PickCoreBackbone(*this);
}

bool PickCoreBackbone::matches( const Particle& particle ) const
{
	return particle.isBackbone() && !particle.isHydrogen();
}


PickHydrogens* PickHydrogens::clone() const
{
	return new PickHydrogens(*this);
}

bool PickHydrogens::matches( const Particle& particle ) const
{
	return particle.isHydrogen();
}


PickRebuildRequired* PickRebuildRequired::clone() const
{
	return new PickRebuildRequired(*this);
}

bool PickRebuildRequired::matches( const Particle& particle ) const
{
	return particle.isRebuildRequired();
}


PickSidechains* PickSidechains::clone() const
{
	return new PickSidechains(*this);
}

bool PickSidechains::matches( const Particle& particle ) const
{
	return !particle.isBackbone();
}


PickAtomIndex* PickAtomIndex::clone() const
{
	return new PickAtomIndex(*this);
}

bool PickAtomIndex::matches( const Particle& particle ) const
{
	return particle.i == iat;
}


PickAtomRange::PickAtomRange() : reversed(false)
{
	m_StartAtomIndex=0;
	m_NAtoms=0;
}

PickAtomRange::PickAtomRange(
	size_t _StartAtomIndex,
	size_t _EndAtomIndex
) : reversed(false)
{
	m_StartAtomIndex=_StartAtomIndex;
	m_NAtoms=_EndAtomIndex - m_StartAtomIndex + 1;
}

PickAtomRange::PickAtomRange(const MoleculeBase &molecule) : reversed(false)
{
	m_StartAtomIndex = 0;
	m_NAtoms=molecule.atom.size();
}

PickAtomRange::PickAtomRange(
	const MoleculeBase &molecule,
	size_t _StartAtomIndex,
	size_t _EndAtomIndex
) : reversed(false)
{
	setRange(molecule,_StartAtomIndex,_EndAtomIndex);
}

void PickAtomRange::setRange(
	const MoleculeBase &molecule,
	size_t _StartAtomIndex,
	size_t _EndAtomIndex
	)
{
	m_StartAtomIndex=_StartAtomIndex;
	m_NAtoms=_EndAtomIndex - m_StartAtomIndex + 1;
	assertRange(molecule);
}

PickAtomRange* PickAtomRange::clone() const
{
	return new PickAtomRange(m_StartAtomIndex,m_NAtoms);
}

bool PickAtomRange::operator==(const PickAtomRange& _Compare) const
{
	return (m_StartAtomIndex == _Compare.m_StartAtomIndex) && (m_NAtoms == _Compare.m_NAtoms);
}

bool PickAtomRange::operator!=(const PickAtomRange& _Compare) const
{
	return (m_StartAtomIndex != _Compare.m_StartAtomIndex) || (m_NAtoms != _Compare.m_NAtoms);
}

void PickAtomRange::assertRange(const MoleculeBase &molecule) const
{
	if(getEndAtomIndex()>=molecule.atom.size()) throw(ArgumentException("End Atom Index must smaller/equal to number of atoms in this molecule\n"));
	if(getStartAtomIndex()>(getEndAtomIndex()+1)) throw(ArgumentException("End Atom Index must be greater than Start Atom Index\n"));
}

bool PickAtomRange::matches( const Particle& particle ) const
{
	return ( (particle.i >= m_StartAtomIndex) && (particle.i < (m_StartAtomIndex+m_NAtoms)));
}


PickAtomRanges::PickAtomRanges()
	: m_PreviousMatchRangeIndex(SIZE_T_FAIL)
{
}

PickAtomRanges* PickAtomRanges::clone() const
{
	return new PickAtomRanges(*this);
}

bool PickAtomRanges::operator==(const PickAtomRanges& _Compare) const
{
	if( _Compare.size() != size() ) return false;

	for( size_t i = 0; i < m_NAtoms.size(); i++ )
		if( (m_StartAtomIndex != _Compare.m_StartAtomIndex) || (m_NAtoms != _Compare.m_NAtoms) )
			return false;

	return true;
}

bool PickAtomRanges::operator!=(const PickAtomRanges& _Compare) const
{
	return !operator==(_Compare);
}

void PickAtomRanges::assertRanges(const MoleculeBase &molecule) const
{
	ASSERT( m_StartAtomIndex.size() == m_NAtoms.size(), CodeException, "Internal array management failure!");
	for( size_t i = 0; i < m_NAtoms.size(); i++ )
	{
		if(getEndAtomIndex(i)>=molecule.atom.size()) throw(ArgumentException("End Atom Index must smaller/equal to number of atoms in this molecule\n"));
		if(getStartAtomIndex(i)>(getEndAtomIndex(i)+1)) throw(ArgumentException("End Atom Index must be greater than Start Atom Index\n"));
	}
}

bool PickAtomRanges::matches( const Particle& particle ) const
{
	for( size_t i = 0; i < m_NAtoms.size(); i++ )
	{
		if( (particle.i >= m_StartAtomIndex[i]) && (particle.i < (m_StartAtomIndex[i]+m_NAtoms[i])) )
		{
			m_PreviousMatchRangeIndex = i;
			return true;
		}
	}
	m_PreviousMatchRangeIndex = SIZE_T_FAIL;
	return false;
}

void PickAtomRanges::addRange( size_t start, size_t length, bool reversed )
{
	m_StartAtomIndex.push_back(start);
	m_NAtoms.push_back(length);
	m_Reversed.push_back( reversed );
}

void PickAtomRanges::addRange( const PickAtomRange& _range, bool reversed )
{
	m_StartAtomIndex.push_back(_range.getStartAtomIndex());
	m_NAtoms.push_back(_range.getNAtoms());
	m_Reversed.push_back( reversed );
}


PickResidueRange::PickResidueRange():
	PickAtomRange()
{
	m_StartResidueIndex=0;
	m_NResidues=0;
}

PickResidueRange::PickResidueRange(
	const MoleculeBase &molecule,
	size_t _StartRes,
	size_t _EndRes
):PickAtomRange()
{
	setRange(molecule,_StartRes,_EndRes);
}

PickResidueRange::PickResidueRange(size_t _StartRes, size_t _EndRes) :
	PickAtomRange(),
	m_StartResidueIndex(_StartRes),
	m_NResidues(_EndRes)
{
}

PickResidueRange* PickResidueRange::clone() const
{
	return new PickResidueRange(*this);
}

void PickResidueRange::setRange(const MoleculeBase &molecule, size_t _StartRes, size_t _EndRes)
{
	if( _StartRes > _EndRes ) throw OutOfRangeException("PickResidueRange requires that the end residue index is after the start index!");
	m_StartResidueIndex = _StartRes;
	m_NResidues = _EndRes - m_StartResidueIndex + 1;
	assertRangeCore(molecule);
	m_StartAtomIndex=molecule.res[_StartRes].ifirst;
	m_NAtoms=molecule.res[_EndRes].ilast - m_StartAtomIndex + 1;
	PickAtomRange::assertRange(molecule);
}

void PickResidueRange::assertRange(const MoleculeBase &molecule) const
{
	assertRangeCore(molecule);
	PickAtomRange::assertRange(molecule);
}

void PickResidueRange::assertRangeCore(const MoleculeBase &molecule) const
{
	if(getEndResIndex()>=molecule.res.size()) throw(ArgumentException("End Residue Index must smaller/equal to number of atoms in this molecule\n"));
	if(getStartResIndex()>(getStartResIndex()+1)) throw(ArgumentException("End Residue Index must be greater than Start Residue Index\n"));
}


PickResidueList::PickResidueList()
{
}

PickResidueList::PickResidueList(const MoleculeBase &_molecule)
{
	resetMolecule( _molecule );
}

PickResidueList::PickResidueList(const PickResidueRange& _MakeFromRange)
{
	m_ResList.clear();
	m_ResList.reserve(_MakeFromRange.getNRes());
	size_t end = _MakeFromRange.getEndResIndex();
	for( size_t i = _MakeFromRange.getStartResIndex(); i <= end; i++ )
	{
		m_ResList.push_back( i );
	}
}

void PickResidueList::assertInternal() const
{
	ASSERT( m_Mol != NULL, UninitialisedException, "PickResidueList: Internal molecule is NULL");
}

void PickResidueList::resetMolecule(const MoleculeBase &_molecule)
{
	m_Mol = &_molecule;
	assertInternal();
	clearList();
}

void PickResidueList::add( size_t _ResidueIndex )
{
	assertInternal();
	ASSERT( _ResidueIndex < m_Mol->nResidues(), OutOfRangeException, "PickResidueList: The supplied residue index does not lie within the scope of the current molecule");
	m_ResList.push_back(_ResidueIndex);
}

void PickResidueList::add( const Residue& _Residue )
{
	assertInternal();
	ASSERT( _Residue.ir < m_Mol->nResidues(), OutOfRangeException, "PickResidueList: The supplied residues <internal> index does not lie within the scope of the current molecule");
	m_ResList.push_back(_Residue.ir);
}

void PickResidueList::add( const PickResidueRange& resRange )
{
	assertInternal();
	for( size_t i = 0; i < resRange.getNRes(); i++ )
	{
		size_t resIndex = resRange.getStartResIndex()+i;
		ASSERT( resIndex < m_Mol->nResidues(), OutOfRangeException, "PickResidueList: The supplied residue index does not lie within the scope of the current molecule");
		if( !VectorContains(m_ResList,resIndex) )
			m_ResList.push_back(resIndex);
	}
}

void PickResidueList::clearList()
{
	m_ResList.clear();
}

bool PickResidueList::matches( const Particle& _particle ) const
{
	for( size_t i = 0; i < m_ResList.size(); i++ )
	{
		if( m_ResList[i] == _particle.ir )
			return true;
	}
	return false;
}

bool PickResidueList::matches( const Residue& _res ) const
{
	for( size_t i = 0; i < m_ResList.size(); i++ )
	{
		if( m_ResList[i] == _res.ir )
			return true;
	}
	return false;
}

void PickResidueList::invertSelection()
{
	ASSERT( m_Mol != NULL, NullInternalException, "PickResidueList: Internal Molecule is NULL" );
	std::vector<size_t> temp_ResList(m_ResList);
	m_ResList.clear();
	for( size_t ir = 0; ir < m_Mol->nResidues(); ir++ )
	{
		if( !VectorContains(temp_ResList, ir ) )
		{
			m_ResList.push_back( ir );
		}
	}
}



PickElement* PickElement::clone() const
{
	return new PickElement(*this);
}

bool PickElement::matches( const Particle& particle ) const
{
	return (particle.Z == m_Z);
}




PickAtomPDBName* PickAtomPDBName::clone() const
{
	return new PickAtomPDBName(*this);
}

bool PickAtomPDBName::matches( const Particle& particle ) const
{
	return (particle.pdbname == m_Name);
}




PickMoleculeName* PickMoleculeName::clone() const
{
	return new PickMoleculeName(*this);
}

bool PickMoleculeName::matches( const Particle& particle ) const
{
	return (particle.parentname == m_Name);
}




PickMolecule3LetterName* PickMolecule3LetterName::clone() const
{
	return new PickMolecule3LetterName(*this);
}

bool PickMolecule3LetterName::matches( const Particle& particle ) const
{
	return (particle.parentl3name == m_Name);
}































