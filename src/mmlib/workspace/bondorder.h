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

#ifndef __BOND_ORDER_H
#define __BOND_ORDER_H

#include <vector>
#include "componentbase.h"
#include "system/fundamentals.h"

extern const int BondOrder_1_1_Pair;
extern const int BondOrder_1_2_Pair;
extern const int BondOrder_1_3_Pair;
extern const int BondOrder_1_4_Pair;
extern const int BondOrder_1_5_Pair;
extern const int BondOrder_1_6_Pair;

/// Holds the atom index j and its bond order
struct AtomOrder
{
	int j;
	int order;
};

/// Define bond orders for pairs of atoms whose index is far apart.
struct LongBondOrder
{
	int i;
	std::vector<AtomOrder> partner;
	size_t memuse( int level );
};


//-------------------------------------------------
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
class PD_API BondOrder : public WorkSpaceComponentBase
{
public:

	BondOrder(int _MaxBondOrder=7);
	size_t memuse(int level);

	// BondOrder queries

	/// Generic, always works, but slow
	int getBondOrder( int i, int j ) const; 

	/// Works if i and j are not sparse. Applicable to 99.9% of bonds.
	int getDiagBondOrder( int i, int j ) const;    

	/// Works if i and j are sparse. These bonds are between 
	/// highly separated regions of a polymer, e.g. long range disulphide bonds
	int getOffDiagBondOrder( int i, int j ) const; 

	/// Tests if a pair of atoms are treated by the band matrix or the off-diagonal
	/// array
	bool isOffDiag( int i, int j ) const { return abs(i-j) > m_MaxIndexDelta; } 

	// Accessors
	int getMaxIndexDelta() const   { return m_MaxIndexDelta; }
	int getMaxBondOrder() const   { return m_MaxBondOrder; }

	// get raw data pointers
	const std::vector<LongBondOrder>& getOffDiagData() const { return m_LongBondOrder; } ///< Obtain the internal data of the sparse matrix - advanced lookup.
	const std::vector<char>&          getDiagData()     const { return m_Data; } 
protected:

	/// Reinitialise to use any new workspace
	virtual void reinit( WorkSpace* _wspace ); 

	/// Calculate the max difference between the indexes of 
	/// bonded atoms within local residue space
	void calcMaxIndexDelta(); 

	/// Internally fills the member data arrays. Called by reinit()
	void calcOrder( int i, std::vector<CovalentAtom>& indexList, int order ); 

	/// Ensure we have only copy of each 'i' in the m_LongBondOrder array
	LongBondOrder& getLongBondOrder( int i ); 

private:
	void clear();

	/// bondorder above this are truncated to this value. This means that
	/// that all pairs which have a bondorder greater than m_MaxBondOrder
	/// will appear to just have a bondorder of m_MaxBondOrder
	int m_MaxBondOrder;       

	/// Essentially the width of the band matrix
	int m_MaxIndexDelta;

	/// This is a band matrix to store the bondorder between atoms 
	/// not too seperated in sequence
	std::vector<char> m_Data; 

	/// This is essentially a list of low bondorder pairs which
	/// fall outside the band matrix. such pairs are rare and thus
	/// are more efficiently stored in a list
	std::vector<LongBondOrder> m_LongBondOrder; /// OffDiag additional bonding data, e.g. disulphides are stored in a separate array.
};

#endif



