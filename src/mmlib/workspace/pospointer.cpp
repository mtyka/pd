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

#include "system/fundamentals.h"
#include "system/molecule.h"
#include "pickers/atomregex.h"
#include "pickers/basicpickers.h"

#include "pospointer.h"

// ---- SelectedIndices Implementation ---------------------

SelectedIndices::SelectedIndices(const MoleculeBase & _mol, const PickBase & _picker)
{
	setPicking( _mol, _picker );
}

void SelectedIndices::setPicking(const MoleculeBase & _mol, const PickBase & _picker)
{
	ASSERT( (&_mol != NULL ), CodeException,"Supplied Molecule/Workspace is NULL in call to SelectedIndices constructor");
	for( size_t i = 0; i < _mol.nAtoms(); i++ )
	{
		if( _picker.matches(_mol.getAtom(i)) )
		{
			m_Atoms.push_back(i);
		}
	}
}

int SelectedIndices::operator[](size_t index) const
{
	ASSERT( index < m_Atoms.size(), OutOfRangeException, "SelectedIndices call out of range");
	return m_Atoms[index];
}

int SelectedIndices::getAtomIndex(size_t index) const
{
	ASSERT( index < m_Atoms.size(), OutOfRangeException, "SelectedIndices call out of range");
	return m_Atoms[index];
}

// -------------------------------------
// --- PosPointerBase Implementation ---
// -------------------------------------

PosPointerBase::PosPointerBase() 
	: m_WSpace(NULL),
	m_Start(0),
	m_Length(0)
{
}

double PosPointerBase::calcCRMS( const PickBase& _picker ) const
{
	ASSERT( size() > 0, ArgumentException, "No atoms are present, cannot calc CRMS");
	double sum = 0.0;
	for( size_t i = 0; i < size(); i++ )
	{
		if( _picker.matches( at(i) ) )
		{
			double sqrDist = p(i).sqrdist( at(i).posRef() );
			sum += sqrDist;
		}
	}
	return sum / (double) size();
}

void PosPointerBase::lookupAtoms(const PickBase & _picker)
{
	ASSERT( (&m_WSpace != NULL ), CodeException,"Internal WorkSpace is NULL in call to PosPointer::lookupAtoms() (setPicking() has probably not been called!)");

	// Clear the current arrays
	m_Index.clear();
	m_Atom.clear();

	m_Pos.clear();
	m_Force.clear();
	m_Velocity.clear();
	
	m_OldPos.clear();
	m_OldForce.clear();
	m_OldVelocity.clear();

	// Proxies
	SnapShotAtom* cur = m_WSpace->cur.atom;
	SnapShotAtom* old = m_WSpace->old.atom;

	size_t toLength = m_Start + m_Length;
	for( size_t i = m_Start; i < toLength; i++ )
	{
		if( _picker.matches( m_WSpace->getAtom(i)) )
		{
			m_Index.push_back(i);
			m_Atom.push_back( &m_WSpace->getAtom(i) );

			m_Pos.push_back( &cur[i].p );
			m_Force.push_back( &cur[i].f );
			m_Velocity.push_back( &cur[i].v );

			m_OldPos.push_back( &old[i].p );
			m_OldForce.push_back( &old[i].f );
			m_OldVelocity.push_back( &old[i].v );
		}
	}
}

// ---------------------------------
// --- PosPointer Implementation ---
// ---------------------------------

PosPointer::PosPointer()
{
}

PosPointer::PosPointer( WorkSpace& _mol, const PickBase& _picker )
{
	setPicking(_mol,_picker);
}

//PosPointer::PosPointer( SegmentDefBase& _seg, const PickBase& _picker )
//{
//	setPicking(_seg,_picker);
//}

void PosPointer::setPicking( WorkSpace& _wspace, const PickBase& _picker )
{
	m_WSpace = &_wspace;
	m_Start = 0;
	m_Length = _wspace.nAtoms();
	lookupAtoms(_picker);
}

//void PosPointer::setPicking( SegmentDefBase& _seg, const PickBase& _picker )
//{	
//	m_WSpace = &_seg.getWorkSpace();
//	m_Start = _seg.getStartAtomIndex();
//	m_Length = _seg.getNAtoms();
//	lookupAtoms(_picker);
//}

// ---------------------------------
// --- PosStore Implementation ---
// ---------------------------------

PosStore::PosStore() 
	: m_DataContained(false)
{
	  m_Picker = new PickNothing();
}

PosStore::PosStore( WorkSpace& _mol ) 
	: m_DataContained(false)
{
	 m_Picker = new PickNothing();
	m_WSpace = &_mol;
}

PosStore::PosStore( WorkSpace& _mol, const PickBase& _picker ) 
	: m_DataContained(false)
{
	setPicking( _mol, _picker );
}

void PosStore::setPicking( const PickBase& _picker )
{
	ASSERT( m_WSpace != NULL, CodeException, "PosStore has not been assigned");
	m_DataContained = false;
	m_Picker = _picker.clone();
	lookupAtoms(_picker);
}

void PosStore::setPicking( WorkSpace& _wspace, const PickBase& _picker )
{
	m_DataContained = false;
	m_WSpace = &_wspace;
	m_Start = 0;
	m_Length = _wspace.nAtoms();
	m_Picker = _picker.clone();
	lookupAtoms(_picker);
}

void PosStore::store()
{	
	ASSERT( m_WSpace != NULL, CodeException, "PosStore has not been assigned");
	m_DataContained = true;
	m_PosCache.resize(size());
	for( size_t i = 0; i < size(); i++ )
	{		
		m_PosCache[i].setTo(p(i));
	}
}

void PosStore::revert()
{
	ASSERT( m_WSpace != NULL, CodeException, "PosStore has not been assigned");
	ASSERT( m_DataContained, CodeException, "store() must be called before revert(). store() must also be recalled following a call to setPicking()");
	for( size_t i = 0; i < size(); i++ )
	{
		p(i).setTo(m_PosCache[i]);
	}
}

double PosStore::calcCRMSOfStoreTo( const PosStore& comp, const PickBase& _picker ) const
{
	ASSERT( size() > 0, ArgumentException, "No atoms are present, cannot calc CRMS");
	ASSERT( size() == comp.size(), ArgumentException, "PosPointerBase correspondence mis-match!");
	ASSERT( m_DataContained && comp.m_DataContained, ProcedureException, "calcCRMSOfStoreTo() required both to have stored data!");

	double sum = 0.0;
	for( size_t i = 0; i < size(); i++ )
	{
		if( _picker.matches( at(i) ) )
		{
			D_ASSERT( at(i).i == comp.at(i).i, ArgumentException, "PosPointerBase correspondence mis-match!" );
			D_ASSERT( _picker.matches( comp.at(i) ), ArgumentException, "PosPointerBase correspondence mis-match!" );
			double sqrDist = m_PosCache[i].sqrdist( comp.m_PosCache[i] );
			sum += sqrDist;
		}
	}
	return sum / (double) size();
}

