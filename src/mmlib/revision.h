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

#ifndef __REVISION
#define __REVISION

#include <string>

//-------------------------------------------------
/// \brief Holds PDs Revision information
/// \author Mike Tyka & Jon Rea
class Revision
{
private:
	Revision(); // You cant create an instance of this class
public:
	static std::string LongRevString();
	static std::string ShortRevString();
	static std::string LibVersion();
	static std::string LibName();
	static std::string FullStamp();
	static std::string CopyRight();
private:
	static void doInit();
	static bool init;
	static std::string longRev;
	static std::string shortRev;
	static std::string fullStamp;
};

#endif

