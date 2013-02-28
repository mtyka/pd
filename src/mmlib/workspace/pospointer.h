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

// Implements a number of classes for internal use that supply pre-picked
// index/atom pointer arrays which is considerably faster than 
// asking each atom all the time if it passes any given filter (Picker)

#ifndef __POS_POINTER_H
#define __POS_POINTER_H

#include <string>
#include <vector>

#include "workspace/workspace.h" // Required for inline implementation
#include "workspace/workspace.fwd.h"
#include "system/molecule.fwd.h"
#include "system/fundamentals.fwd.h"
#include "pickers/basicpickers.h"
#include "tools/cloneholder.h"

class PickBase;
class SegmentDefBase;

//-------------------------------------------------
/// \brief The simplest possible selection - it just supplies indices with no garantuee
/// of their up-to-date-ness thus this class should be short lived
/// \details Pretty much what it says on the tin
/// \author Jon Rea 
class PD_API SelectedIndices
{
public:
	SelectedIndices(){};
	SelectedIndices(const MoleculeBase& _mol, const PickBase& _picker);

	void setPicking(const MoleculeBase & _mol, const PickBase & _picker);

	size_t size()const{ return m_Atoms.size(); } ;

#ifndef SWIG
	int operator[]( size_t index ) const;
#endif
	int getAtomIndex( size_t index ) const;

private:
	std::vector<int> m_Atoms;
};


//-------------------------------------------------
/// \brief A base class for more advanced types of prepicker compared to SelectedIndices
/// - they hold a pointer to a workspace and can provide atom pointers directly. 
/// \details NOTE - For performance reasons it is underisable for this class to contain *ANY* virtual functions
/// \author Jon Rea 
class PD_API PosPointerBase
{
protected:
	PosPointerBase();

public:
	size_t size() const; ///< The current number of matching atoms
	
#ifndef SWIG
	const Particle& operator[](size_t index ) const; ///< Public const accessor for the particle at the index in the list
	Particle& operator[](size_t index ); ///< Public accessor for the particle at the index in the list
#endif

	const Particle& at(size_t index ) const; ///< Public const accessor for the particle at the index in the list
	Particle& at(size_t index ); ///< Public accessor for the particle at the index in the list
	
	const Maths::dvector& p(size_t index ) const; ///< Public const accessor for the atom position at the index in the list
	Maths::dvector& p(size_t index ); ///< Public accessor for the atom position at the index in the list

	const Maths::dvector& f(size_t index ) const; ///< Public const accessor for the atom force at the index in the list
	Maths::dvector& f(size_t index ); ///< Public accessor for the atom force at the index in the list

	const Maths::dvector& v(size_t index ) const; ///< Public const accessor for the atom velocity at the index in the list
	Maths::dvector& v(size_t index ); ///< Public accessor for the atom velocity at the index in the list

	const Maths::dvector& op(size_t index ) const; ///< Public const accessor for the old atom position at the index in the list
	Maths::dvector& op(size_t index ); ///< Public accessor for the old atom position at the index in the list

	const Maths::dvector& of(size_t index ) const; ///< Public const accessor for the old atom force at the index in the list
	Maths::dvector& of(size_t index ); ///< Public accessor for the old atom force at the index in the list

	const Maths::dvector& ov(size_t index ) const; ///< Public const accessor for the old atom velocity at the index in the list
	Maths::dvector& ov(size_t index ); ///< Public accessor for the old atom velocity at the index in the list

	int iat( size_t index ) const; ///< The index of the refered atom

	double calcCRMS( const PickBase& _picker = PickEverything() ) const;

protected:	
	void lookupAtoms(const PickBase & _picker);

	// Range
	size_t m_Start;
	size_t m_Length;

	// Parent WSpace Pointer
	WorkSpace* m_WSpace;

	// Atom and index cache
	std::vector<int> m_Index;
	std::vector<Particle*> m_Atom;

	// Current Pos
	std::vector<Maths::dvector*> m_Pos;
	std::vector<Maths::dvector*> m_Force;
	std::vector<Maths::dvector*> m_Velocity;

	// Old Pos
	std::vector<Maths::dvector*> m_OldPos;
	std::vector<Maths::dvector*> m_OldForce;
	std::vector<Maths::dvector*> m_OldVelocity;
};


class PD_API PosPointer : public PosPointerBase
{
public:
	PosPointer(); ///< Default constructor. The atom list will be 0 length until setPicking() is called.
	PosPointer( WorkSpace& _wspace, const PickBase& _picker ); ///< Intenrnally calls setPicking() to generate our atom list.
	//PosPointer( SegmentDefBase& _seg, const PickBase& _picker ); ///< Intenrnally calls setPicking() to generate our atom list.

	/// set the current molecule and search pattern. 
	void setPicking( WorkSpace& _wspace, const PickBase & _picker );
	//void setPicking( SegmentDefBase& _seg, const PickBase & _picker );		
};


//-------------------------------------------------
/// \brief  Takes a simple position-only snapshot of a given PosPointer array.
/// \details Use store and revert to cacke positions
/// \author Jon Rea 
class PD_API PosStore : public PosPointerBase
{
public:
	PosStore();
	PosStore( WorkSpace& _mol );
	PosStore( WorkSpace& _mol, const PickBase& _picker );

	/// set a new position pointer
	void setPicking( const PickBase& _picker );
	void setPicking( WorkSpace& _wspace, const PickBase& _picker );

	const Maths::dvector& sp(size_t index ) const; ///< Public const accessor for the STORED atom position at the index in the list
	Maths::dvector& sp(size_t index ); ///< Public accessor for the STORED atom position at the index in the list

	double calcCRMSOfStoreTo( const PosStore& comp, const PickBase& _picker = PickEverything() ) const;

	/// save the structure from the workspace
	void store();
	/// load the structure into the workspace
	void revert();

	// report original picker
	const PickBase& getPicker() const { 
		return m_Picker.data(); 
	}
protected:
	bool m_DataContained;
	std::vector<Maths::dvector> m_PosCache;

	CloneHolder<PickBase> m_Picker;
};


//--------------------------------------------------------------
//  Begin pseudo-CPP file and define all inline function bodies
// -------------------------------------------------------------

inline size_t PosPointerBase::size() const
{
	return m_Atom.size(); 
}

inline const Particle& PosPointerBase::operator[](size_t index ) const
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_Atom.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_Atom[index];
}

inline Particle& PosPointerBase::operator[](size_t index )
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_Atom.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_Atom[index];
}

inline const Particle& PosPointerBase::at(size_t index ) const
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_Atom.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_Atom[index];
}

inline Particle& PosPointerBase::at(size_t index )
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_Atom.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_Atom[index];
}

inline int PosPointerBase::iat(size_t index ) const
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_Index.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return m_Index[index];
}

inline const Maths::dvector& PosPointerBase::p(size_t index ) const
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_Pos.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_Pos[index];
}

inline Maths::dvector& PosPointerBase::p(size_t index )
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_Pos.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_Pos[index];
}

inline const Maths::dvector& PosPointerBase::f(size_t index ) const
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_Force.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_Force[index];
}

inline Maths::dvector& PosPointerBase::f(size_t index )
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_Force.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_Force[index];
}

inline const Maths::dvector& PosPointerBase::v(size_t index ) const
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_Velocity.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_Velocity[index];
}

inline Maths::dvector& PosPointerBase::v(size_t index )
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_Velocity.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_Velocity[index];
}

inline const Maths::dvector& PosPointerBase::op(size_t index ) const
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_OldPos.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_OldPos[index];
}

inline Maths::dvector& PosPointerBase::op(size_t index )
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_OldPos.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_OldPos[index];
}

inline const Maths::dvector& PosPointerBase::of(size_t index ) const
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_OldForce.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_OldForce[index];
}

inline Maths::dvector& PosPointerBase::of(size_t index )
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_OldForce.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_OldForce[index];
}

inline const Maths::dvector& PosPointerBase::ov(size_t index ) const
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_OldVelocity.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_OldVelocity[index];
}

inline Maths::dvector& PosPointerBase::ov(size_t index )
{
	D_ASSERT( (m_WSpace != NULL ), ProcedureException,"Internal molecule is NULL in call (setPattern() has probably not been called!)");
	D_ASSERT( (index < m_OldVelocity.size()), OutOfRangeException, "Index given refers to a position outside of the internal array." );
	return *m_OldVelocity[index];
}

inline const Maths::dvector& PosStore::sp(size_t index ) const
{
	D_ASSERT( m_WSpace != NULL, CodeException, "PosStore wspace has not been assigned");
	D_ASSERT( m_DataContained, CodeException, "no data!");
	return m_PosCache[index];
}

inline Maths::dvector& PosStore::sp(size_t index )
{
	D_ASSERT( m_WSpace != NULL, CodeException, "PosStore wspace has not been assigned");
	D_ASSERT( m_DataContained, CodeException, "no data!");
	return m_PosCache[index];
}

#endif


