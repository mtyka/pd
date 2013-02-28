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
#include "workspace/workspace.h"
#include "bondorder.h"


const int BondOrder_1_1_Pair = 0;
const int BondOrder_1_2_Pair = 1;
const int BondOrder_1_3_Pair = 2;
const int BondOrder_1_4_Pair = 3;
const int BondOrder_1_5_Pair = 4;
const int BondOrder_1_6_Pair = 5;

size_t LongBondOrder::memuse( int level )
{
	return sizeof(LongBondOrder) + (sizeof(IndexPair) * partner.size());
}

BondOrder::BondOrder(int _MaxBondOrder) : 
	WorkSpaceComponentBase(NULL), 
	m_MaxBondOrder(_MaxBondOrder),
	m_MaxIndexDelta(-1)
{
}

size_t BondOrder::memuse(int level)
{
	size_t useage = sizeof(*this) + (sizeof(char) * m_Data.size());
	for( size_t i = 0; i < m_LongBondOrder.size(); i++ )
	{
		useage += m_LongBondOrder[i].memuse( level + 1 );
	}
	return useage;
}

int BondOrder::getBondOrder( int i, int j ) const
{
	return isOffDiag( i,j ) ? getOffDiagBondOrder(i,j) : getDiagBondOrder(i,j);
}

int BondOrder::getDiagBondOrder( int i, int j ) const 
{
	if( i == j ) return 0;
	if( i > j ) std::swap( i, j );
	int delta = j - i; // We have only made 1/2 a symmetrical matrix, i.e. it is mirrored - 'j' must therefore be the greater value
	D_ASSERT( ( delta <= m_MaxIndexDelta ), CodeException, "Incorrect use of getBondOrder(). Atom index delta is too large - isOffDiag() and the OffDiagBondOrder should have been used.");
	return m_Data[sqrmat(i, delta-1, m_MaxIndexDelta)];
}

int BondOrder::getOffDiagBondOrder( int i, int j ) const
{
	if( i > j ) 
		std::swap(i,j);
	for( size_t q = 0; q < m_LongBondOrder.size(); q++ )
	{
		if( i == m_LongBondOrder[q].i )
		{
			const std::vector<AtomOrder>& partner = m_LongBondOrder[q].partner;
			for( size_t r = 0; r < partner.size(); r++ )
			{
				if( partner[r].j == j )
				{
					return partner[r].order;
				}
			}
		}
	}
	return m_MaxBondOrder;
}

void BondOrder::clear()
{
	m_MaxIndexDelta = -1;
	m_Data.clear();
	m_LongBondOrder.clear();
}

void BondOrder::calcMaxIndexDelta()
{
	// Proxies
	ParticleStore &atom = wspace->atom;
	int natom = (int) wspace->atom.size();  

	// Init
	m_MaxIndexDelta = 0;

	// Look at all the 1-4 partners
	for(int i = 0; i < natom; i++) 
	{
		{
			std::vector<CovalentAtom>& indexList = atom[i].cov12atom;
			for(size_t nb = 0; nb < indexList.size(); nb++)
			{
				int lookup = indexList[nb].i;
				int atomIndexDelta = abs( i - lookup );
				int resDelta = abs( atom[i].ir - atom[lookup].ir );
				if( resDelta > 1 ) 
					continue;
				m_MaxIndexDelta = std::max( atomIndexDelta, m_MaxIndexDelta );
			}
		}
		{
			std::vector<CovalentAtom>& indexList = atom[i].cov13atom;
			for(size_t nb = 0; nb < indexList.size(); nb++)
			{
				int lookup = indexList[nb].i;
				int atomIndexDelta = abs( i - lookup );
				int resDelta = abs( atom[i].ir - atom[lookup].ir );
				if( resDelta > 1 ) 
					continue;
				m_MaxIndexDelta = std::max( atomIndexDelta, m_MaxIndexDelta );
			}
		}
		{
			std::vector<CovalentAtom>& indexList = atom[i].cov14atom;
			for(size_t nb = 0; nb < indexList.size(); nb++)
			{
				int lookup = indexList[nb].i;
				int atomIndexDelta = abs( i - lookup );
				int resDelta = abs( atom[i].ir - atom[lookup].ir );
				if( resDelta > 1 ) 
					continue;
				m_MaxIndexDelta = std::max( atomIndexDelta, m_MaxIndexDelta );
			}
		}
	}
	m_MaxIndexDelta+=1;
}

LongBondOrder& BondOrder::getLongBondOrder( int i )
{
	size_t q; // Defined outside for loop scope on purpose
	for( q = 0; q < m_LongBondOrder.size(); q++ )
	{
		if( m_LongBondOrder[q].i == i )
		{
			return m_LongBondOrder[q];
		}
	}
	LongBondOrder lr;
	lr.i = i;
	m_LongBondOrder.push_back( lr );
	return m_LongBondOrder[q]; // q will be == m_LongBondOrder.size()
}

void BondOrder::calcOrder( int i, std::vector<CovalentAtom>& indexList, int order )
{
	for(size_t nb = 0; nb < indexList.size(); nb++)
	{
		int j = indexList[nb].i;
		if( j <= i ) continue; // we are calculating half of a symmetrical diagonal matrix
		int delta = j - i; // Index delta
		if( delta > m_MaxIndexDelta ) // Is it bigger than the current storage width?
		{
			// This is therefore a 'sparse' bond order, positioned away from the main array diagonal
			// This happens when there are covalent bonds at large index deltas along a polymer
			AtomOrder pd;
			pd.j = j;
			pd.order = Maths::min(order,m_MaxBondOrder);
			getLongBondOrder(i).partner.push_back( pd );
		}
		else
		{
			// Main Bonded Region in the stripe along the diagonal of the theoretical full matrix
			int sq = sqrmat(i, delta-1, m_MaxIndexDelta);
			m_Data[sq] = Maths::min(order,m_MaxBondOrder);
		}
	}
}

void BondOrder::reinit( WorkSpace* _wspace )
{
	printf("  Initialising bondorder array...\n");

	clear(); // reinitialise

	wspace = _wspace;
	ParticleStore &atom = wspace->atom;
	int natom = (int) wspace->atom.size();  

	calcMaxIndexDelta(); // What width do we want?

	// allocate bondOrderData array
	// '-1' as we have no need to record the order when i == j (its always 0), so we need less space by one.
	m_Data.clear();
	m_Data.resize( m_MaxIndexDelta * (natom - 1), m_MaxBondOrder );

	// now go through atoms and set the appropriate bondOrderData entries to the
	// bondOrderData between each respective atom
	for(int i = 0; i < natom; i++) 
	{		
		calcOrder( i, atom[i].cov14atom, 3 ); // bondOrderData 3 - i.e. 1-4 pairs
		calcOrder( i, atom[i].cov13atom, 2 ); // bondOrderData 2 - i.e. 1-3 pairs
		calcOrder( i, atom[i].cov12atom, 1 ); // bondOrderData 1 - i.e. 1-2 pairs
	}
}

