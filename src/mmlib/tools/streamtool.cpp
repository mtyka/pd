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

#include "streamtool.h"

#include <iomanip> // Needed for some stream manipulation functionality e.g. setw()

using namespace std;

void setFloatFmt( std::ostream &_file, int _width, int _precision )
{
	_file.setf(ios::showpoint); // ensure the decimal point is set
	_file.setf(ios::right,ios::adjustfield); // ensure all is right aligned to the available width
	_file.setf(ios::fixed,ios::floatfield); // ensure that the floating point representation
	_file.width(_width); // set the print width
	_file.precision(_precision); // set the number of decimal palaces
}

