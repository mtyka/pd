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

#ifndef __ROTAMER_FILTER
#define __ROTAMER_FILTER

#include "rotamer_applicatorbase.h" // base class

namespace Manipulator
{
	/// \brief A filter that ensures that only the most likely states are kept within the library.
	/// \details SCWRL 3.0 (original paper, page 2005, column 1, step 2, states that "rotamers are read in from
	/// highest to lowest probability until the cumulative density reaches at least 90%." I take that to mean that
	/// the <10% of low probability rotamers are removed. This class acomplishes that in any given library by
	/// filtering these improbable states.
	class PD_API RotamerCumulativeProbabilityDensityFilter : public RotamerFilterBase
	{
	public:
		RotamerCumulativeProbabilityDensityFilter( const RotamerApplicatorBase& _app, double _MinDensity );
		virtual ~RotamerCumulativeProbabilityDensityFilter(){}
		virtual RotamerCumulativeProbabilityDensityFilter* clone() const;
		virtual bool passes( const RotamerLink& _ir, size_t _rotID );
		double getMinDensity() const { return m_MinDensity; }
	private:
		double m_MinDensity;
	};

	/// Only allow CYS rotamers which are compatible with forming a disulphide bond
	class PD_API RotamerDisulphideFilter : public RotamerFilterBase
	{
	public:
		RotamerDisulphideFilter( const RotamerApplicatorBase& _app );
		virtual ~RotamerDisulphideFilter(){}
		virtual RotamerDisulphideFilter* clone() const;
		virtual bool passes( const RotamerLink& _ir, size_t _rotID );
	};

	/// Only allow rotamers that are compatible with the coordination of a given location in space.
	class PD_API RotamerCoordinationFilter : public RotamerFilterBase
	{
	public:
		RotamerCoordinationFilter( const RotamerApplicatorBase& _app );
		virtual ~RotamerCoordinationFilter(){}
		virtual RotamerCoordinationFilter* clone() const;
		virtual bool passes( const RotamerLink& _ir, size_t _rotID );
	};
}

#endif

