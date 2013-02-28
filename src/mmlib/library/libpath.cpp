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
#include "tools/io.h"
#include "library/libpath.h"

LibraryPathStore* LibraryPathStore::getSingleton()
{
	static LibraryPathStore inst;
	return &inst;
}

LibraryPathStore::LibraryPathStore()
{
	parseEnvironmentVariable();
}

// Parses the environment variable PD_PARAM_PATH into search paths to look for
// library files. The order of searching is: local  directory, lib path,
void LibraryPathStore::parseEnvironmentVariable()
{
	m_LibraryPath.clear();
	m_LibraryPath.push_back(".");

	char *c_env = getenv("PD_PARAM_PATH");

	// if the return value of getenv is NULL the variable is undefined.
	if(c_env == NULL)
	{
		printf("Library Paths:  PD_PARAM_PATH undefined\n");
		return;
	}

	std::string env_PD_PARAM_PATH = "";
	env_PD_PARAM_PATH = c_env;
	std::vector <std::string> temp_m_LibraryPath = chopstr(env_PD_PARAM_PATH, ";:, \t\12\15" );

	for(size_t i=0;i<temp_m_LibraryPath.size();i++)
	{
		m_LibraryPath.push_back( temp_m_LibraryPath[i] );
	}

	printf("Library Paths: \n");
	for(size_t i=0;i<m_LibraryPath.size();i++)
	{
		printf("  %s/\n",m_LibraryPath[i].c_str());
	}
}

std::string LibraryPathStore::findFullFilename(const std::string &filename)
{
	for(int i=0;i<m_LibraryPath.size();i++)
	{
		std::string composite;
		composite = m_LibraryPath[i] + "/" + filename;
		if(IO::fileExists(composite))
		{
			return composite;
		};
	}
	return filename; // return filename evwen if it fails such that error messages are fine
}

