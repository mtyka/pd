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

#ifndef __QUOTE_H
#define __QUOTE_H

#include <string>
#include <vector>






//-------------------------------------------------
//
/// \brief A Quote Engine :-) 
///
/// \details No MM library would be complete without a quote engine ;-)
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API Quote
{
public:
	static Quote* getSingleton(); /// 'Singleton Pattern' Implementation
	Quote();
	void IncludeDB( std::string _filename );
	void printQuote();
private:
	std::vector<std::string> m_Quotes;
};

#endif
