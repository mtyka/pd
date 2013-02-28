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

#ifndef __CONSTANTS_H
#define __CONSTANTS_H

// these can be directly casted to integers for interpretation... i.e. red is 0 in DAVE
class Colour
{
	private:
		Colour();
	public:
		static const int Red;
		static const int Blue;
		static const int Yellow;
		static const int Green;
		static const int Cyan;
		static const int Magenta;
		static const int Orange;
		static const int Turquoise;
};

#endif

