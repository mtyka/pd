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
#include "system/molecule.h"
#include "pickbase.h"

void PickBase::test( const MoleculeBase& _mol ) const
{
	const ParticleStore& atom = _mol.atom;
	printf("Matches: ");
	for( size_t i = 0; i < atom.size(); i++ )
	{
		if( matches( atom[i] ) )
		{
			printf("%d, ",i);
		}
	}
	printf("\nNot Matches: ");
	for( size_t i = 0; i < atom.size(); i++ )
	{
		if( !matches( atom[i] ) )
		{
			printf("%d, ",i);
		}
	}
	printf("\n");
}

size_t PickBase::countPickedAtoms( const MoleculeBase& _mol ) const
{
	const ParticleStore& atom = _mol.atom;
	size_t count = 0;
	for( size_t i = 0; i < atom.size(); i++ )
	{
		if( matches( atom[i] ) )
		{
			count++;
		}
	}
	return count;
}

size_t PickResidueBase::countPickedResidues( const MoleculeBase& _mol ) const
{
	const ResidueStore& res = _mol.res;
	size_t count = 0;
	for( size_t i = 0; i < res.size(); i++ )
	{
		if( matches( res[i] ) )
		{
			count++;
		}
	}
	return count;
}

PickLogical_UNARY::PickLogical_UNARY(
	const PickBase &_Picker )
	: PickBase(),
	m_Picker(_Picker.clone())	
{
}

PickLogical_UNARY::PickLogical_UNARY( const PickLogical_UNARY &_Copy ) 
	: PickBase(),
	m_Picker(_Copy.getPicker().clone())	
{
}

PickLogical_UNARY::~PickLogical_UNARY()
{
	delete m_Picker;
}

PickLogical_UNARY& PickLogical_UNARY::operator=( const PickLogical_UNARY &_Copy )
{
	delete m_Picker;
	m_Picker = _Copy.getPicker().clone();
	return *this;
}



Pick_NOT::Pick_NOT(
	const PickBase &_Picker
): 
	PickLogical_UNARY(_Picker) 
{
}

Pick_NOT* Pick_NOT::clone() const
{
	return new Pick_NOT(getPicker());
}

bool Pick_NOT::matches( const Particle& particle ) const  
{
	return !getPicker().matches(particle);
}



PickLogical_BINARY::PickLogical_BINARY(
	const PickBase &_PickerA,
	const PickBase &_PickerB )
	: PickBase(),
	m_PickerA(_PickerA.clone()),
	m_PickerB(_PickerB.clone())
{
}

PickLogical_BINARY::PickLogical_BINARY( const PickLogical_BINARY &_Copy ) 
	: PickBase(),
	m_PickerA(_Copy.getPickerA().clone()),
	m_PickerB(_Copy.getPickerB().clone())
{
}

PickLogical_BINARY::~PickLogical_BINARY()
{
	delete m_PickerA;
	delete m_PickerB;
}

PickLogical_BINARY& PickLogical_BINARY::operator=( const PickLogical_BINARY &_Copy )
{
	delete m_PickerA;
	delete m_PickerB;
	m_PickerA = _Copy.getPickerA().clone();
	m_PickerB = _Copy.getPickerB().clone();
	return *this;
}


Pick_AND::Pick_AND(
	const PickBase &_PickerA,
	const PickBase &_PickerB
): 
	PickLogical_BINARY(_PickerA,_PickerB) 
{
}

Pick_AND* Pick_AND::clone() const
{
	return new Pick_AND(getPickerA(),getPickerB());
}

bool Pick_AND::matches( const Particle& particle ) const  
{
	return getPickerA().matches(particle) && getPickerB().matches(particle);
}


Pick_OR::Pick_OR(
	const PickBase &_PickerA,
	const PickBase &_PickerB
): PickLogical_BINARY(_PickerA,_PickerB) {};


Pick_OR* Pick_OR::clone() const
{
	return new Pick_OR(getPickerA(),getPickerB());
}

bool Pick_OR::matches( const Particle& particle ) const  
{
	return getPickerA().matches(particle) || getPickerB().matches(particle);
}


bool Pick_AND_Group::matches( const Particle& particle ) const
{
	for( size_t i = 0; i < size(); i++ )
	{
		if( !element(i).matches(particle ) )
		{
			return false;
		}
	}
	return true;
}


bool Pick_OR_Group::matches( const Particle& particle ) const
{
	for( size_t i = 0; i < size(); i++ )
	{
		if( element(i).matches(particle ) )
		{
			return true;
		}
	}
	return false;
}


