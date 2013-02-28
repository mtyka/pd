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
#include "tools/stringbuilder.h"
#include "revision.h"

bool Revision::init = false;
std::string Revision::longRev = "";
std::string Revision::shortRev = "";
std::string Revision::fullStamp = "";

void Revision::doInit()
{
	// Auto-make the subversion version string.
	StringBuilder build( "$Revision: 1184 $" ); // *NOTE* this string is automatically edited by subversion when revision.cpp is altered
	build.Trim("$"); // remove $ signs
	build.erase(0,9); // remove 'Revision:'
	build.Trim(); // remove remaining whitespace

	shortRev = build.toString();
	longRev = "SVN Revision: " + shortRev;

	build.clear();

	build.append( LibName() );
	build.append( ' ' );
	build.append( LibVersion() );
	build.append( ' ' );
	build.append( longRev );

	fullStamp = build.toString();

	init = true;
}

std::string Revision::FullStamp()
{
	if( !init ) doInit();
	return fullStamp;
}

std::string Revision::CopyRight()
{
	return "(c) J.Rea & M.Tyka 2003-2007";
}

std::string Revision::LibName()
{
	return "MMLib";
}

std::string Revision::LibVersion()
{
	return "v0.8b";
}

std::string Revision::LongRevString()
{
	if( !init ) doInit();
	return longRev;
}

std::string Revision::ShortRevString()
{
	if( !init ) doInit();
	return shortRev;
}

