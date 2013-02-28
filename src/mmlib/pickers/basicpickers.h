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

#ifndef __BASIC_PICKERS_H
#define __BASIC_PICKERS_H

#include "pickbase.h"

class PickAtomRange;
class MoleculeBase;

/// Picks Nothing - all functions return false
class PickNothing: public PickResidueBase
{
public:
	virtual PickNothing* clone() const;
	virtual bool matches( const Particle& particle ) const;
	virtual bool matches( const Residue& _res ) const;
};

/// Picks Everything - all functions return true
class PickEverything: public PickResidueBase
{
public:
	virtual PickEverything* clone() const;
	virtual bool matches( const Particle& particle ) const;
	virtual bool matches( const Residue& _res ) const;
};

/// Picks every particle (returns true for any particle entry)
class PickAllParticles: public PickBase
{
public:
	virtual PickAllParticles* clone() const;
	virtual bool matches( const Particle& particle ) const;
};

/// Picks every real atom (not dummy atoms)
class PickAllAtoms: public PickBase
{
public:
	virtual PickAllAtoms* clone() const;
	virtual bool matches( const Particle& particle ) const;
};

/// Pick all atoms deemed 'heavy'
class PickHeavyAtoms : public PickBase
{
public:
	virtual PickHeavyAtoms* clone() const;
	virtual bool matches( const Particle& particle ) const;
};

/// Pick only a single given residue
class PickResidue : public PickResidueBase
{
public:
	PickResidue(size_t ires){ m_iRes = ires; }
	virtual PickResidue* clone() const;
	virtual bool matches( const Particle& particle ) const;
	virtual bool matches( const Residue& _res ) const;
private:
	size_t m_iRes;
};

/// Pick only a single given molecule
class PickMolecule : public PickBase
{
public:
	PickMolecule(size_t imol){ m_iMol = imol; }
	virtual PickMolecule* clone() const;
	virtual bool matches( const Particle& particle ) const;
private:
	size_t m_iMol;
};

/// Pick molecule range
class PickMoleculeRange : public PickBase
{
public:
	PickMoleculeRange(size_t imol_start, size_t imol_end)
	{ 
		m_iMol_start = imol_start; 
		m_iMol_end   = imol_end; 
	}

	virtual PickMoleculeRange* clone() const;
	virtual bool matches( const Particle& particle ) const;
private:
	size_t m_iMol_start;
	size_t m_iMol_end;
};

/// Pick CA only
class PickCAAtoms : public PickBase
{
public:
	PickCAAtoms(){}
	virtual PickCAAtoms* clone() const;
	virtual bool matches( const Particle& particle ) const;
};

/// Pick N CA C O H and HA only
class PickBackbone : public PickBase
{
public:
	PickBackbone(){}
	virtual PickBackbone* clone() const;
	virtual bool matches( const Particle& particle ) const;
};

/// Pick N CA C and O only
class PickCoreBackbone : public PickBase
{
public:
	PickCoreBackbone(){}
	virtual PickCoreBackbone* clone() const;
	virtual bool matches( const Particle& particle ) const;
};

/// Pick Hydrogens
class PickHydrogens : public PickBase
{
public:
	PickHydrogens(){}
	virtual PickHydrogens* clone() const;
	virtual bool matches( const Particle& particle ) const;
};

/// Pick Hydrogens
class PickRebuildRequired : public PickBase
{
public:
	PickRebuildRequired(){}
	virtual PickRebuildRequired* clone() const;
	virtual bool matches( const Particle& particle ) const;
};

/// Pick only sidechains
class PickSidechains : public PickBase
{
public:
	PickSidechains(){}
	virtual PickSidechains* clone() const;
	virtual bool matches( const Particle& particle ) const;
};

class PickAtomIndex : public PickBase
{
public:
	PickAtomIndex( int _iat ) : iat(_iat) {}
	virtual PickAtomIndex* clone() const;
	virtual bool matches( const Particle& particle ) const;
protected:
	int iat;
};

/// This class definition defines a molecule and a residue range within that molecule
class PickAtomRange: public PickBase
{
public:
	PickAtomRange();
	PickAtomRange(
		size_t _StartAtomIndex,
		size_t _EndAtomIndex
	);
	PickAtomRange(const MoleculeBase &molecule);
	PickAtomRange(
		const MoleculeBase &molecule,
		size_t _StartAtomIndex,
		size_t _EndAtomIndex
	);

	virtual PickAtomRange* clone() const;

	bool operator==(const PickAtomRange& _Compare) const;
	bool operator!=(const PickAtomRange& _Compare) const;

	inline size_t getNAtoms() const { return m_NAtoms; }; ///< Returns the number of atoms in this range
	inline size_t getStartAtomIndex() const { return m_StartAtomIndex; };  ///< Returns the index of the first atom
	inline size_t getEndAtomIndex() const { return m_StartAtomIndex + m_NAtoms - 1; };  ///< Returns the index of the final atom

	/// checks that the range passed is "compatible" with this molecule
	/// and throws an exception of type ArgumentException otherwise
	void assertRange(const MoleculeBase &molecule) const;

	virtual bool matches( const Particle& particle ) const;

	bool reversed; ///< A value indicating if this should be conceptually treated as a reverse range. Even if this is true, m_StartAtomIndex will still me the N-terminal end.

protected:
	void setRange(const MoleculeBase &molecule, size_t _StartAtomIndex,	size_t _EndAtomIndex); ///< It is not desirable for setRange() to be public. The base class for range should be immutable - derived classes can make this public if they wish.
	size_t m_StartAtomIndex;
	size_t m_NAtoms;
};

class PickAtomRanges: public PickBase
{
public:
	PickAtomRanges();

	virtual PickAtomRanges* clone() const;
	bool operator==(const PickAtomRanges& _Compare) const;
	bool operator!=(const PickAtomRanges& _Compare) const;

	size_t size() const { return m_NAtoms.size(); }
	void clear() { m_NAtoms.clear(); m_StartAtomIndex.clear(); }

	inline size_t getPreviousMatchRangeIndex() const { return m_PreviousMatchRangeIndex; }
	inline size_t getNAtoms(size_t i) const { return m_NAtoms[i]; }; ///< Returns the number of atoms in this range
	inline size_t getStartAtomIndex(size_t i) const { return m_StartAtomIndex[i]; };  ///< Returns the index of the first atom
	inline size_t getEndAtomIndex(size_t i) const { return m_StartAtomIndex[i] + m_NAtoms[i] - 1; };  ///< Returns the index of the final atom
	inline bool getReversed(size_t i) const { return m_Reversed[i]; }

	/// checks that the range passed is "compatible" with this molecule
	/// and throws an exception of type ArgumentException otherwise
	void assertRanges(const MoleculeBase &molecule) const;

	virtual bool matches( const Particle& particle ) const;

	void addRange( size_t start, size_t length, bool reversed = false );
	void addRange( const PickAtomRange& _range, bool reversed = false );

private:
	mutable size_t m_PreviousMatchRangeIndex;
	std::vector<size_t> m_StartAtomIndex;
	std::vector<size_t> m_NAtoms;
	std::vector<bool> m_Reversed;
};

/// This class definition defines a molecule and a residue range within that molecule
/// NOTE! : This class derives from 'PickAtomRange' and **not** 'PickResidueBase' for efficiency reasons. 
/// MD for example only takes PickResidueRange beacuse it can treat it as an atom range ...
class PickResidueRange: public PickAtomRange
{
public:
	PickResidueRange();

	/// Initialises this class and calls assertRange() to validate the range arguments
	PickResidueRange(
		const MoleculeBase &molecule,
		size_t _StartRes,
		size_t _EndRes
	);

	virtual PickResidueRange* clone() const;

	/// Returns the number of residues in this range
	inline size_t getNRes() const { return m_NResidues; }
	/// Returns the index of the first residue in this range
	inline size_t getStartResIndex() const { return m_StartResidueIndex; }
	/// Returns the index of the last residue in this range
	inline size_t getEndResIndex() const { return m_StartResidueIndex + m_NResidues - 1; }
	/// Returns the index of the central residue in this range. If the numbber of residues is odd, the extra residue is added to the start.
	inline size_t getCentreResIndex() const { return m_StartResidueIndex + ((m_NResidues+1) / 2) - 1; }

	/// checks that the range passed is "compatible" with this molecule
	/// and throws an exception of type ArgumentException otherwise
	void assertRange(const MoleculeBase &molecule) const;

protected:
	PickResidueRange(size_t _StartRes, size_t _EndRes); // Protected constructor
	void assertRangeCore(const MoleculeBase &molecule) const;
	void setRange(const MoleculeBase &molecule,	size_t _StartRes, size_t _EndRes); ///< It is not desirable for this to be public. The base class for range should be immutable - derived classes can make this public if they wish.
	size_t m_StartResidueIndex;
	size_t m_NResidues;
};

class PickResidueList: public PickResidueBase
{
public:
	PickResidueList();
	PickResidueList(const MoleculeBase &_molecule);
	PickResidueList(const PickResidueRange& _MakeFromRange);

	virtual ~PickResidueList(){}
	virtual PickResidueList* clone() const { return new PickResidueList(*this); }

	void invertSelection(); ///< Pick whatever we dont pick now and visa versa

	void resetMolecule(const MoleculeBase &_molecule); ///< Reassign the internal molecule. This triggers clearList() internally.
	void add( size_t _ResidueIndex ); ///< Add to the list from a simple index
	void add( const Residue& _ResidueIndex ); ///< Add to the list from a Residue reference
	void add( const PickResidueRange& resRange );
	void clearList(); ///< Clear all residue definitions
	virtual bool matches( const Particle& _particle ) const;
	virtual bool matches( const Residue& _res ) const;
protected:
	std::vector<size_t> m_ResList;
	const MoleculeBase* m_Mol;
private:
	void assertInternal() const; ///< Assert that we have no internal NULL pointers.
};



class PickElement: public PickBase
{
public:
	PickElement(int _Z){ m_Z = _Z; }
	virtual PickElement* clone() const;
	virtual bool matches( const Particle& particle ) const;
private:
	size_t m_Z;
};


class PickAtomPDBName : public PickBase
{
public:
	PickAtomPDBName(const std::string &_name){ m_Name = _name; }
	virtual PickAtomPDBName* clone() const;
	virtual bool matches( const Particle& particle ) const;
private:
	std::string m_Name;
};


class PickMoleculeName: public PickBase
{
public:
	PickMoleculeName(const std::string &_name){ m_Name = _name; }
	virtual PickMoleculeName* clone() const;
	virtual bool matches( const Particle& particle ) const;
private:
	std::string m_Name;
};


class PickMolecule3LetterName: public PickBase
{
public:
	PickMolecule3LetterName(const std::string &_name){ m_Name = _name; }
	virtual PickMolecule3LetterName* clone() const;
	virtual bool matches( const Particle& particle ) const;
private:
	std::string m_Name;
};

#endif

