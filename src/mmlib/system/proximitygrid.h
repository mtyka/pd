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

#ifndef __PROXIMITY_GRID_H
#define __PROXIMITY_GRID_H

#include "maths/maths.fwd.h"
#include "system/fundamentals.fwd.h"
#include "maths/tntjama/tnt_array3d.h"
#include "tools/cloneholder.h"
#include "pickers/pickbase.h"

struct GridPoint
{
	GridPoint()
	{
		gridPos.setTo(DBL_MAX,DBL_MAX,DBL_MAX);
	}
	std::vector<int> atomIndexes;
	Maths::dvector gridPos;
	inline void clear() { atomIndexes.clear(); }
	inline void push_back( int i ) { atomIndexes.push_back(i); }
	inline size_t size() const { return atomIndexes.size(); }
	inline int operator[]( int i ) const { return atomIndexes[i];}
};

//-------------------------------------------------
/// \brief Defines a grid that can give a client a list of atoms close to a given cartesian location.
/// \details Defines a cartesian grid, centred on the molecules centre of geometry. For each grid point, 
/// a list of atom indexes is maintained. A client can thenr request atoms which are proximal to any 
/// coordinate. Internally the closest grid point is found and the corresponding list returned.
/// \author Jon Rea 
class PD_API ProximityGrid : public ICloneable
{
public:
	ProximityGrid( const MoleculeBase& _Mol, const PickBase& _Picker, double _InclusionCutoff ); ///< Initialise the grid with the recommended granularity
	ProximityGrid( const MoleculeBase& _Mol, const PickBase& _Picker, double _InclusionCutoff, size_t _Granularity ); ///< Initialise the grid with a custom granularity. Probably pointless, and maybe dangerous.
	
	virtual ProximityGrid* clone() const { return new ProximityGrid(*this); }

	void setCutoff( double _InclusionCutoff ); ///< The cutoff distance before a given atom is assigned as close to a given grid-point

	void refresh(); ///< Recalculate the grid indexes. Each grid point is assigned atom indexes which are within the _InclusionCutoff and pass the picker.
	void refreshForce(); ///< Performs the actual refresh of the Grid()
	bool requiresRefresh(); ///< Used by refresh() to decide whether to call refreshCore() which actually does the refresh. If any of the picked atoms have moved significantly, then refersh the list.

	bool atomList( const Maths::dvector& position, const GridPoint*& _List ) const; ///< Returns true if a proximal grid point is found. _List will then point to the GridPoint for the indexes of atoms that are proximal to this cartesian coordinate; the atoms in the list corresponding to the closest grid point.
	static int recommendGranularity(const MoleculeBase& _Mol, const PickBase& _Picker, double _Cutoff); ///< How many grid points do we need to cover the picked space with a sufficient grid density?

	inline const Maths::dvector& getWorldCentre() const { return m_WorldCentre; } /// the cartesian coordinate of the centre of the grid
	inline double getXBound() const { return m_XBound; } ///< +- distance of the grid on the X-axis
	inline double getYBound() const { return m_YBound; } ///< +- distance of the grid on the Y-axis
	inline double getZBound() const { return m_ZBound; } ///< +- distance of the grid on the Z-axis
	inline const TNT::Array3D<GridPoint>& getGrid() const { return m_Grid; } /// obtain a const reference to the internal grid
	const PickBase& getPicker() const { return m_Picker.data(); } ///< obtain a const reference to the internal atom picker

	// Debug / drawing provision
	void drawGrid(IDrawProvider& _Vect) const; ///< Draw the grid to a drawing provider to see where it is :-D
	bool drawMatches( IDrawProvider& _Vect, const Maths::dvector& pos ) const; ///< Draw vectors between the cartesian positions of atoms relating to the closest grid point to this cartesian location

private:
	void getGridLoc( const Maths::dvector& pos, int& gridX, int& gridY, int& gridZ ) const; ///< Get the indexers for a given cartesian location
	void generateGrid(); ///< Make our grid :-D
	void storePickedPositions(); ///< We cache the positions of the picked positions so that we can later decide if the list needs to be updated.
	
	TNT::Array3D<GridPoint> m_Grid; ///< The grid itself.
	std::vector<Maths::dvector> m_PreviousLocations; ///< The stored position list.
	GridPoint nullGrid; ///< A zero length list, returned when an off-grid request is made.
	Maths::dvector m_WorldCentre; ///< The centre of the grids world.
	CloneHolder<PickBase> m_Picker; ///< The picker which determines what atoms can be added to the grid.
	const MoleculeBase* m_Mol; ///< The molecule to which the grid refers.

	double m_SqrInclusionCutoff; ///< The square of the value passed to setCutoff()
	double m_XBound, m_YBound, m_ZBound; ///< The cartesian bounds of the grid
	double m_GridMulX, m_GridMulY, m_GridMulZ; ///< Store this frequently used value for efficiency
	double m_GridAddX, m_GridAddY, m_GridAddZ; ///< Store this frequently used value for efficiency
};

#endif

