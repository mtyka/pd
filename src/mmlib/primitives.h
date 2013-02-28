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

#ifndef __PRIMITIVES_H
#define __PRIMITIVES_H

#include "maths/maths_vector.h"

// ----------------------------------------------
//  Header containing the most primitive of data
//  types. These *shouldn't* relate to anything
//  biological they are pure abstract concepts.
// ----------------------------------------------

//-------------------------------------------------
//
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
class PD_API IndexPair
{
 public:
	IndexPair();
	IndexPair(int _i, int _j);
	int i, j;
};


//-------------------------------------------------
//
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
class PD_API IndexTriplet
{
 public:
	IndexTriplet();
	IndexTriplet(int _i, int _a, int _j);
	int i, a, j;
};


//-------------------------------------------------
//
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
class PD_API IndexQuartet
{
 public:
	IndexQuartet();
	IndexQuartet(int _i, int _a, int _b, int _j);
	int i, a, b, j;
};


//-------------------------------------------------
/// \brief  A simple set of 8 size_t indeces
/// \details Used by RotLib
/// \author Jon Rea 
class PD_API IndexHexet
{
 public:
	IndexHexet();
	IndexHexet(size_t _a, size_t _b, size_t _c, size_t _d, size_t _i, size_t _j);
	size_t a, b, c, d, i, j;
};


//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
struct DrawingVector
{
	// the two points between which to draw a line
	Maths::fvector v1; ///< Coordinate 1
	Maths::fvector v2; ///< Coordinate 2
	int colourCode; ///< the colour of the line. See constants.h : Colour. see DAVE documentation for what the full range of integer codes are ...

	DrawingVector();
	void clear();
};


//------------------------------------------------
/// \brief  Defines a set of two related strings.
/// \author Mike Tyka & Jon Rea 
class PD_API StringPair
{
public:
	StringPair();
	StringPair(
		const std::string& _p, 
		const std::string& _q
		);
#ifndef SWIG
	bool operator==( const StringPair& sp ) const;
#endif
	std::string p;
	std::string q;
};


//-------------------------------------------------
/// \brief  Defines a set of four related strings.
/// \author Mike Tyka & Jon Rea 
class PD_API StringQuartet
{
public:
	StringQuartet();
	StringQuartet(
		const std::string& _p, 
		const std::string& _q, 
		const std::string& _r, 
		const std::string& _s );
	std::string p;
	std::string q;
	std::string r;
	std::string s;
};

/// Place the value between + and - PI
double torsionalRangeRadians( double value );

/// Place the value between + and - PI
double torsionalRangeDegrees( double value );

//-------------------------------------------------
/// \brief A class defining a torsional tollerance i.e. a torsional value +- other values
/// \details Internally deals with torsional periodicity. Works in RADIANS!
/// \author Jon Rea 
/// \todo Complete
/// \bug No Known
class TorsionalTolleranceRange
{
public:
	TorsionalTolleranceRange();
	TorsionalTolleranceRange( double _value );
	TorsionalTolleranceRange( double _value, double _plus_or_minus );
	TorsionalTolleranceRange( double _value, double _plus, double _minus );
	
	bool isInRange( double anyRadianValue ) const; ///< Says wether or not a given value is in the range, taking torsional periodicity into account.

	void setValue( double _value ) { m_Value = torsionalRangeRadians(_value); }
	void setPlus( double _plus ) { m_Plus = _plus; }
	void setMinus( double _minus ) { m_Minus = _minus; }

	double getValue() const { return m_Value; }
	double getPlus() const { return m_Plus; }
	double getMinus() const { return m_Minus; }

private:
	double m_Value;
	double m_Plus;
	double m_Minus;
};

#endif
