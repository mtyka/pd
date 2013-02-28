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

#ifndef _VERBOSITY_H
#define _VERBOSITY_H

/// \brief An enumeration representing a PD-wide verbosity standard
/// \details Most objects have a member called 'Verbosity::Type OutputLevel;'
/// This defines the verbosity of that object. One then makes comparisons e.g.:
/// if( OutputLevel >= Verbosity::Normal ) Print("My wiiife!");
/// \author Jon Rea
class PD_API Verbosity
{
private:
	Verbosity();
public:
	enum Type
	{
		Silent = 0, ///< Just keep completely schtum!!
		Quiet = 2,  ///< Less output that normal. Just what I REALLY want to know!
		Normal = 4, ///< The kind of output you want to see every day
		Loud = 6,   ///< A bit more than normal. Mostly minor extra info.
		Scream = 8, ///< Debuggers paradise
		Eleven = 10 ///< Ever seen Spinal Tap? :-D God I'm whitty...
	};
};

#endif

