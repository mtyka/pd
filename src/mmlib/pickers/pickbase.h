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

#ifndef __PICKER_H
#define __PICKER_H

class Particle;
class Residue;
class MoleculeBase;

/// Arbitrary atom selection mask
class PickBase
{
public:
	PickBase(){};
	virtual ~PickBase(){};
	virtual PickBase* clone() const = 0;
	/// Returns true if the atom is valid when compared to the internal definition
	virtual bool matches( const Particle& particle ) const = 0;

	// Helpers
	void test( const MoleculeBase& _mol ) const;
	size_t countPickedAtoms( const MoleculeBase& _mol ) const;
};

/// A base class for all pickers which choose residues (and therefore atoms within those residues)
class PickResidueBase : public PickBase
{
public:
	PickResidueBase(){}
	virtual bool matches( const Residue& _res ) const = 0;

	// Helpers
	size_t countPickedResidues( const MoleculeBase& _mol ) const;
};

// --- Logical combinations of pickers --------------
//
// These allow the user to do things like this:
// SelectedIndices myselection( Pick_AND(PickHeavyAtom(), PickAtomRange(10,40) ) );
// SelectedIndices myselection( Pick_NOT( PickAtomRange(10,40) ) );
// SelectedIndices myselection( Pick_OR( PickHeavyAtom(), Pick_NOT( PickAtomRange(10,40) ) );
//
// Note that these dont own their arguments so it's the users responsibility
// to keep the original objects alive as long as needed!
//

class PickLogical_UNARY: public PickBase
{
public:
	PickLogical_UNARY( const PickBase &_Picker );
	PickLogical_UNARY( const PickLogical_UNARY &_Copy );
#ifndef SWIG
	PickLogical_UNARY& operator=( const PickLogical_UNARY &_Copy );
#endif

	virtual ~PickLogical_UNARY();
	virtual PickLogical_UNARY* clone() const = 0;

	/// Returns true if the atom is valid when compared to the internal definition
	virtual bool matches( const Particle& particle ) const = 0;

protected:
	PickBase& getPicker() { return *m_Picker; }
	const PickBase& getPicker() const { return *m_Picker; }

private:
	PickBase *m_Picker;
};


class Pick_NOT : public PickLogical_UNARY
{
public:
	Pick_NOT( const PickBase &_Picker );
	virtual Pick_NOT* clone() const;
	virtual bool matches( const Particle& particle ) const;
};


class PickLogical_BINARY: public PickBase
{
public:
	PickLogical_BINARY( const PickBase &_PickerA, const PickBase &_PickerB );
	PickLogical_BINARY( const PickLogical_BINARY &_Copy );
#ifndef SWIG
	PickLogical_BINARY& operator=( const PickLogical_BINARY &_Copy );
#endif 
	virtual ~PickLogical_BINARY();
	virtual PickLogical_BINARY* clone() const = 0;

	/// Returns true if the atom is valid when compared to the internal definition
	virtual bool matches( const Particle& particle ) const = 0;

protected:
	PickBase& getPickerA() { return *m_PickerA; }
	PickBase& getPickerB() { return *m_PickerB; }
	const PickBase& getPickerA() const { return *m_PickerA; }
	const PickBase& getPickerB() const { return *m_PickerB; }

protected:
	PickBase *m_PickerA;
	PickBase *m_PickerB;
};

class Pick_AND: public PickLogical_BINARY
{
public:
	Pick_AND( const PickBase &_PickerA,	const PickBase &_PickerB);
	virtual Pick_AND* clone() const;
	virtual bool matches( const Particle& particle ) const;
};

class Pick_OR: public PickLogical_BINARY
{
public:
	Pick_OR( const PickBase &_PickerA, const PickBase &_PickerB);	
	virtual Pick_OR* clone() const;
	virtual bool matches( const Particle& particle ) const;
};

#ifdef SWIG
%template(ObjectContainer_PickBase) ObjectContainer<PickBase>;
#endif

/// Returns true if any of the internal pickers return true.
class Pick_AND_Group : public PickBase, public ObjectContainer<PickBase>
{
public:
	Pick_AND_Group(){};
	virtual Pick_AND_Group* clone() const { return new Pick_AND_Group(*this); }
	virtual bool matches( const Particle& particle ) const;
};


/// Returns true if all of the internal pickers return true.
class Pick_OR_Group : public PickBase, public ObjectContainer<PickBase>
{
public:
	Pick_OR_Group(){};
	virtual Pick_OR_Group* clone() const { return new Pick_OR_Group(*this); }
	virtual bool matches( const Particle& particle ) const;
};

#endif

