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
#include "primitives.h"

IndexPair::IndexPair()
{
}

IndexPair::IndexPair(int _i, int _j)
	: i(_i), j(_j) 
{
}

IndexTriplet::IndexTriplet()
{
}

IndexTriplet::IndexTriplet(int _i, int _a, int _j)
	: i(_i),a(_a),j(_j) 
{
}

IndexQuartet::IndexQuartet()
{
}

IndexQuartet::IndexQuartet(int _i, int _a, int _b, int _j)
	: i(_i), a(_a), b(_b), j(_j)
{
}

IndexHexet::IndexHexet()
{
}

IndexHexet::IndexHexet(size_t _a, size_t _b, size_t _c, size_t _d, size_t _i, size_t _j)
	: a(_a), b(_b), c(_c), d(_d), i(_i), j(_j)
{
}

DrawingVector::DrawingVector()
{
	clear();
}

void DrawingVector::clear()
{
	v1.setTo(0.0f,0.0f,0.0f);
	v2.setTo(0.0f,0.0f,0.0f);
	colourCode = 0;
}

StringPair::StringPair()
{
}

StringPair::StringPair( const std::string& _p, const std::string& _q )
	: p(_p), q(_q)
{
}

bool StringPair::operator==( const StringPair& sp ) const
{
	return (0 == sp.p.compare(p)) && (0 == sp.q.compare(q));
}

StringQuartet::StringQuartet()
{
}

StringQuartet::StringQuartet( const std::string& _p, const std::string& _q, const std::string& _r, const std::string& _s )
	: p(_p), q(_q), r(_r), s(_s)
{
}

double torsionalRangeRadians( double value )
{
	while( value > Maths::MathConst::PI )
	{
		value -= Maths::MathConst::TwoPI;
	}
	while( value < -Maths::MathConst::PI )
	{
		value += Maths::MathConst::TwoPI;
	}
	return value;
}

double torsionalRangeDegrees( double value )
{
	while( value > 180.0 )
	{
		value -= 360.0;
	}
	while( value < -180.0 )
	{
		value += 360.0;
	}
	return value;
}

TorsionalTolleranceRange::TorsionalTolleranceRange() 
	: m_Value(0.0), m_Plus(0.0), m_Minus(0.0)
{
}

TorsionalTolleranceRange::TorsionalTolleranceRange( double _value ) 
	: m_Plus(0.0), m_Minus(0.0)
{
	setValue( _value );
}

TorsionalTolleranceRange::TorsionalTolleranceRange( double _value, double _plus_or_minus ) 
{
	ASSERT( _plus_or_minus >= 0.0, ArgumentException, "Required: _plus_or_minus >= 0.0");
	setValue( _value );
	setPlus(_plus_or_minus);
	setMinus(-_plus_or_minus);
}

TorsionalTolleranceRange::TorsionalTolleranceRange( double _value, double _plus, double _minus ) 
{
	ASSERT( _minus <= _plus, ArgumentException, "Required: minus <= plus");
	ASSERT( _minus <= 0.0, ArgumentException, "Required: minus <= 0.0");
	ASSERT( _plus >= 0.0, ArgumentException, "Required: plus >= 0.0");
	setValue( _value );
	setPlus(_plus);
	setMinus(_minus);
}

bool TorsionalTolleranceRange::isInRange( double testMe ) const
{
	// m_Value is garanteed to be in the radian range +- PI

	testMe = torsionalRangeRadians( testMe );
	double max = torsionalRangeRadians(m_Value+m_Plus);
	double min = torsionalRangeRadians(m_Value+m_Minus); // + cos its a -ve number!

	if( max > min )
	{
		return min <= testMe && testMe <= max;
	}
	else
	{
		return max <= testMe && testMe <= min;
	}
}

