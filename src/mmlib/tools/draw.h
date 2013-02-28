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

#ifndef __DRAW_TOOLS_H
#define __DRAW_TOOLS_H

/// \brief Draw the origin (0,0,0) as a Point via the IDrawProvider
void drawOrigin( IDrawProvider& _Vect, int colourCode );

/// \brief Draw the given cartesian coordinate as a Point via the IDrawProvider
void drawPoint( IDrawProvider& _Vect, const Maths::dvector& pos, int colourCode );

#endif

