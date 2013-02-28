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
#include "draw.h"

void drawOrigin( IDrawProvider& _Vect, int colourCode )
{
	Maths::dvector origin(0.0,0.0,0.0);
	drawPoint( _Vect, origin, colourCode );
}

void drawPoint( IDrawProvider& _Vect, const Maths::dvector& pos, int colourCode )
{
	const float ptSize = 0.2f;
	DrawingVector* dv = _Vect.request();
	if( dv != NULL )
	{
		dv->v1.setTo( pos );
		dv->v1.x -= ptSize;
		dv->v2.setTo(pos);
		dv->v2.x += ptSize;
		dv->colourCode = colourCode;
	}
	dv = _Vect.request();
	if( dv != NULL )
	{
		dv->v1.setTo( pos );
		dv->v1.y -= ptSize;
		dv->v2.setTo(pos);
		dv->v2.y += ptSize;
		dv->colourCode = colourCode;
	}
	dv = _Vect.request();
	if( dv != NULL )
	{
		dv->v1.setTo( pos );
		dv->v1.z -= ptSize;
		dv->v2.setTo(pos);
		dv->v2.z += ptSize;
		dv->colourCode = colourCode;
	}
}

