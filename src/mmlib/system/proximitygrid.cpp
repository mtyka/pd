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
#include "tools/draw.h"
#include "system/molecule.h"
#include "proximitygrid.h"

const double moveTollerance = 0.5; // Each atom can have moved move by approximately 0.5A (arbritrary) before the grid is invalidated
const double defaultMaxSpacing = 1.0; // We want a grid point roughly at least every 1.0A on the longest axis.
const double defaultSpacingTollerance = defaultMaxSpacing * 1.2; // The generated spacing can be upto 20% larger than the default, before we trigger an error status

ProximityGrid::ProximityGrid( const MoleculeBase& _Mol, const PickBase& _Picker, double _InclusionCutoff )
	: m_Grid(recommendGranularity(_Mol,_Picker,_InclusionCutoff)),
	m_Mol(&_Mol),
	m_Picker(_Picker),
	m_SqrInclusionCutoff(_InclusionCutoff*_InclusionCutoff)
{
	//refresh(); -- this needs to be called prior to use anyway, so probably best not to here, it can be expensive
}

ProximityGrid::ProximityGrid( const MoleculeBase& _Mol, const PickBase& _Picker, double _InclusionCutoff, size_t _Granularity )
	: m_Grid(_Granularity),
	m_Mol(&_Mol),
	m_Picker( _Picker),
	m_SqrInclusionCutoff(_InclusionCutoff*_InclusionCutoff)
{
	//refresh(); -- this needs to be called prior to use anyway, so probably best not to here, it can be expensive
}

void ProximityGrid::setCutoff( double _InclusionCutoff )
{
	ASSERT( _InclusionCutoff > 0.0, ArgumentException, "InclusionCutoff must be >0.0");
	_InclusionCutoff *= _InclusionCutoff; // square it
	if( m_SqrInclusionCutoff == _InclusionCutoff )
	{
		refresh(); // refresh if it needs it
		return; // no point if its the same!	
	}
	m_SqrInclusionCutoff = _InclusionCutoff;
	refreshForce(); // FORCE the refresh, we have just changed the cutoff!
}

void ProximityGrid::refresh()
{
	if( !requiresRefresh() )
		return;
	refreshForce();
}

void ProximityGrid::refreshForce()
{
	storePickedPositions();
	generateGrid();
	D_ASSERT( !requiresRefresh(), CodeException, "Something smells fishy!" );
}

void Bound( const ParticleStore& atoms, const PickBase& _Picker, const Maths::dvector& worldCentre, double _Cutoff, double& m_XBound, double& m_YBound, double& m_ZBound )
{
	// Measure the X, Y and Z deviatons from the world centre
	// we kinda assume that most systems are *roughly* spherical, as each
	// grid direction has the same number of points, but will have different
	// cartesian bounds dependent on the current system.
	m_XBound = 0.0;
	m_YBound = 0.0;
	m_ZBound = 0.0;
	for( size_t iat = 0; iat < atoms.size(); iat++ )
	{
		if( _Picker.matches( atoms[iat] ) )
		{
			double x = fabs(worldCentre.x - atoms[iat].pos().x);
			if( x > m_XBound )
				m_XBound = x;

			double y = fabs(worldCentre.y - atoms[iat].pos().y);
			if( y > m_YBound )
				m_YBound = y;

			double z = fabs(worldCentre.z - atoms[iat].pos().z);
			if( z > m_ZBound )
				m_ZBound = z;
		}
	}

	// Our grid has to extend further than the scope of the atoms themselves, i.e. also by the search cutoff.
	m_XBound += _Cutoff;
	m_YBound += _Cutoff;
	m_ZBound += _Cutoff;

	return;
}

int ProximityGrid::recommendGranularity(const MoleculeBase& _Mol, const PickBase& _Picker, double _Cutoff)
{
	double m_XBound, m_YBound, m_ZBound;
	Maths::dvector worldCentre = _Mol.getCentreOfGeometry(_Picker);
	Bound( _Mol.atom, _Picker, worldCentre, _Cutoff, m_XBound, m_YBound, m_ZBound );
	return (int) 
		std::max(
		(m_XBound / defaultMaxSpacing),
		std::max( 
		(m_YBound / defaultMaxSpacing),
		std::max( 
		(m_ZBound / defaultMaxSpacing),
		2.0 ))); // any less than 2 and we will get a code violation due to divide by zero later
}

void ProximityGrid::drawGrid( IDrawProvider& _Vect ) const
{
	for( int i = 0; i < m_Grid.dim1(); i++ )
	{				
		for( int j = 0; j < m_Grid.dim2(); j++ )
		{					
			for( int k = 0; k < m_Grid.dim3(); k++ )
			{
				drawPoint( _Vect, m_Grid[i][j][k].gridPos, Colour::Orange );
			}
		}
	}
}

bool ProximityGrid::drawMatches( IDrawProvider& _Vect, const Maths::dvector& _Pos ) const
{
	const ParticleStore& atoms = m_Mol->atom;
	const GridPoint* p_list = NULL;
	if( !atomList( _Pos, p_list ) ) 
	{
		drawPoint( _Vect, _Pos, Colour::Magenta );
		return false;
	}

	const GridPoint& list = *p_list;

	drawPoint( _Vect, _Pos, Colour::Cyan );
	drawPoint( _Vect, list.gridPos, Colour::Green );

	for( size_t iat = 0; iat < atoms.size(); iat++ )
	{	
		if( m_Picker->matches( atoms[iat] ) )
		{
			bool found = false;
			for( size_t i = 0; i < list.size(); i++ )
			{
				if( list[i] == iat )
				{
					found = true;
					double dist = list.gridPos.sqrdist( atoms[iat].pos() );
					ASSERT( dist < m_SqrInclusionCutoff,
						CodeException, "Ugh!" );

					DrawingVector* dv = _Vect.request();
					dv->v1.setTo(list.gridPos);
					dv->v2.setTo( atoms[iat].pos() );
					dv->colourCode = Colour::Green;

					break;
				}
			}
			if( !found )
			{
				double dist = list.gridPos.sqrdist( atoms[iat].pos() );
				ASSERT( 
					dist >= m_SqrInclusionCutoff,
					CodeException, "Ugh!" );

				DrawingVector* dv = _Vect.request();
				dv->v1.setTo(list.gridPos);
				dv->v2.setTo( atoms[iat].pos() );
				dv->colourCode = Colour::Red;
			}
		}
	}
	return true;
}

void ProximityGrid::generateGrid()
{
	const ParticleStore& atoms = m_Mol->atom;
	m_WorldCentre = m_Mol->getCentreOfGeometry(m_Picker.data());	
	Bound( atoms, m_Picker.data(), m_WorldCentre, 
		m_SqrInclusionCutoff <= 0.0 ? 0.0 : m_SqrInclusionCutoff / m_SqrInclusionCutoff, 
		m_XBound, m_YBound, m_ZBound );

	ASSERT( 
		m_XBound / (double)(m_Grid.dim1()-1) < defaultSpacingTollerance &&
		m_YBound / (double)(m_Grid.dim2()-1) < defaultSpacingTollerance &&
		m_ZBound / (double)(m_Grid.dim3()-1) < defaultSpacingTollerance, 
		ProcedureException,
		"ProximityGrid over-stretched! Increase granularity in initial setup!");

	const double xTermMul = (1.0 / (double)(m_Grid.dim1()-1)) * 2.0 * m_XBound;
	const double xTermPlus = m_WorldCentre.x - m_XBound;
	const double yTermMul = (1.0 / (double)(m_Grid.dim2()-1)) * 2.0 * m_YBound;
	const double yTermPlus = m_WorldCentre.y - m_YBound;
	const double zTermMul = (1.0 / (double)(m_Grid.dim3()-1)) * 2.0 * m_ZBound;
	const double zTermPlus = m_WorldCentre.z - m_ZBound;

	m_GridMulX = (1.0 / (2.0 * m_XBound)) * (double)(m_Grid.dim1()-1);
	m_GridMulY = (1.0 / (2.0 * m_YBound)) * (double)(m_Grid.dim2()-1);
	m_GridMulZ = (1.0 / (2.0 * m_ZBound)) * (double)(m_Grid.dim3()-1);
	m_GridAddX = m_XBound - m_WorldCentre.x + (0.5/m_GridMulX);
	m_GridAddY = m_YBound - m_WorldCentre.y + (0.5/m_GridMulY);
	m_GridAddZ = m_ZBound - m_WorldCentre.z + (0.5/m_GridMulZ);

	// Re-Init, reseting the list, and assigning the cartsian coordinate of this grid point
	Maths::dvector pos;
	for( int i = 0; i < m_Grid.dim1(); i++ )
	{
		pos.x = ((double)i * xTermMul) + xTermPlus;
		for( int j = 0; j < m_Grid.dim2(); j++ )
		{
			pos.y = ((double)j * yTermMul) + yTermPlus;
			for( int k = 0; k < m_Grid.dim3(); k++ )
			{
				pos.z = ((double)k * zTermMul) + zTermPlus;
				m_Grid[i][j][k].clear();
				m_Grid[i][j][k].gridPos.setTo(pos);
			}
		}
	}
	
	for( size_t iat = 0; iat < atoms.size(); iat++ )
	{		
		if( m_Picker->matches( atoms[iat] ) )
		{
			for( int i = 0; i < m_Grid.dim1(); i++ )
			{				
				for( int j = 0; j < m_Grid.dim2(); j++ )
				{					
					for( int k = 0; k < m_Grid.dim3(); k++ )
					{			
						GridPoint& g = m_Grid[i][j][k];
						if( g.gridPos.sqrdist(atoms[iat].pos()) < m_SqrInclusionCutoff )
						{
							g.push_back(iat);							
						}
					}
				}
			}
		}
	}
}

bool ProximityGrid::atomList( const Maths::dvector& pos, const GridPoint*& _List ) const
{
	int gridX = (int)((pos.x + m_GridAddX) * m_GridMulX);
	int gridY = (int)((pos.y + m_GridAddY) * m_GridMulY);
	int gridZ = (int)((pos.z + m_GridAddZ) * m_GridMulZ);

	if( 
		gridX < 0 || gridX >= m_Grid.dim1() ||
		gridY < 0 || gridY >= m_Grid.dim2() ||
		gridZ < 0 || gridZ >= m_Grid.dim3() )
	{
		_List = &nullGrid;
		return false;
	}
	else
	{
		_List = &m_Grid[gridX][gridY][gridZ];
		return true;
	}
}

void ProximityGrid::storePickedPositions()
{
	const ParticleStore& atoms = m_Mol->atom;
	m_PreviousLocations.clear();
	for( size_t i = 0; i < atoms.size(); i++ )
	{
		if( m_Picker->matches( atoms[i] ) )
		{
			m_PreviousLocations.push_back( atoms[i].pos() );
		}
	}
}

bool ProximityGrid::requiresRefresh()
{
	const ParticleStore& atoms = m_Mol->atom;
	size_t i = 0;
	size_t j = 0;

	while( i < m_PreviousLocations.size() && j++ < atoms.size() )
	{
		if( m_Picker->matches( atoms[j] ) )
		{			
			if( atoms[j].pos().sqrdist( m_PreviousLocations[i] ) > Maths::sqr(moveTollerance) )
			{
				return false;
			}
			i++;
		}
	}

	bool toppedI = i == m_PreviousLocations.size();
	bool toppedJ = j == atoms.size();
	return !(toppedI && toppedJ);
}

