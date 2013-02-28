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

#ifndef __LIBPATH_H
#define __LIBPATH_H

#include <string>
#include <vector>

//-------------------------------------------------
//
/// \brief Singleton class which holds the PD library paths and provides filenames and/or file handles
/// Singleton pattern - there's only one library path store managing all the
/// library paths. Libraries can then access this centrally to find their files.
///
/// \details  You just get the instance point from the singleton, pass it the naked file name and
///           the class with try to find the file in the current directory (forst) and then in any direcotires
///           specified in PD_LIB_PATH in the order they're defined. It returns a string with the fully qualified
///           filename (with the absolute path):
///
///           Example:
///
///             std::string fullfilename = LibraryPathStore::getSingleton()->findFullFilename(filename);
/// \author Mike Tyka
///
class PD_API LibraryPathStore
{
public:
	static LibraryPathStore* getSingleton(); /// 'Singleton Pattern' Implementation
	LibraryPathStore();

	std::string findFullFilename(const std::string &filename);

private:
	void parseEnvironmentVariable();
	std::vector<std::string> m_LibraryPath;
};

#endif


