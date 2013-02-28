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

#ifndef _ATOM_RANGE_STORE
#define _ATOM_RANGE_STORE

#include "pickers/basicpickers.h"

class MoleculeBase;

//-------------------------------------------------
/// \brief  A specialisation of a PosStore for an atomic range as opposed to an
/// atomic picker. There may be efficiency to be had by doing this...
/// The class was written before PosStore , which is capable of the same thing with similar syntax,
/// but requires pointer lookups to store, it may therefore be slower. (unbenchmarked)
/// We should discuss wether to keep or remove this class.
/// \author Jon Rea 
class AtomRangeStore
{
public:
	AtomRangeStore();
	AtomRangeStore( MoleculeBase &_Molecule, const PickAtomRange& _Range ) ;
	void setTarget( MoleculeBase &_Molecule, const PickAtomRange& _Range );
	void store();
	void revert();

protected:
	std::vector<Maths::dvector> m_PosCache;
	PickAtomRange m_Range;
	MoleculeBase* m_Molecule;
};

#endif



