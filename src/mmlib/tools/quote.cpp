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
#include <fstream>
#include "tools/stringbuilder.h"
#include "maths/fastrandom.h"
#include "quote.h"

Quote* Quote::getSingleton()
{
	static Quote inst;
	return &inst;
}

Quote::Quote()
{
	IncludeDB("./lib/quotes.dat");
	IncludeDB("./quotes.dat");
}

void Quote::IncludeDB( std::string _filename )
{
	std::ifstream file(_filename.c_str(), std::ifstream::in);
	if( !file.is_open() ) 
	{
		// Opening the file failed. It either is locked or doesn't exist. But hey, we don't really care that much...
		return;
	}

	std::string line;
	while(std::getline(file,line))
	{
		line = trim(line," \t\n\r");
		if( line.length() > 0 ) 
			m_Quotes.push_back(line);
	}

	file.close(); // close the filestream to release the file_handle
}

void Quote::printQuote()
{
	if( m_Quotes.size() == 0 ) return;

	const int delimitWidth = 70;
	const int delimitBorder = 3; // must be >= 2
	const char borderDelimiter = ' ';
	const int width = delimitWidth - 2 * delimitBorder;

	static Maths::FastRandom f; // make ourselves a new system time based random number generator
	int index = (int)f.nextUInt() % m_Quotes.size();
	StringBuilder quote(m_Quotes[index]);
	size_t cutAuthor = quote.LastOf('-');
	std::string authorString = "";
	if( cutAuthor != SIZE_T_FAIL )
	{
		authorString = quote.toString(cutAuthor,quote.size()-cutAuthor);
		quote.erase(cutAuthor,quote.size()-cutAuthor);
	}

	// We don't allow tabs, they will bolox the alignment system.
	for( size_t i = 0; i < quote.size(); i++ )
	{
		if( quote[i] == '\t' )
		{
			quote.replace(i,' ');
		}
	}
	quote.TrimRight();
	quote.insert(0,'"');
	quote.append('"');

	// Header line
	std::cout << borderDelimiter << borderDelimiter << std::string(delimitWidth-4,'-');

	int printedWidth = width; // '= width' ensure our border is printed on the 1st line ...
	while( true )
	{
		// make sure we have space (below logic can fill upto the end of the line)
		if( printedWidth == width )
		{
			std::cout << std::endl << std::string(delimitBorder,borderDelimiter);			
			printedWidth = 0;
		}

		// Our break condition
		if( quote.size() == 0 )
		{
			std::cout << std::endl;
			break;
		}

		size_t indexer = quote.FirstOf(' ');
		if( indexer == SIZE_T_FAIL )
		{
			indexer = quote.size(); // print the remainder, either single word or a load of crap remaining
		}
		else
		{
			indexer++; // include the ' ' itself
		}

		if( indexer > width )
		{
			// This is one FAT word! - print what we can on the current line.
			while(true)
			{
				size_t spaceRemain = width-printedWidth-1; // -1 for the lines terminal '-'
				if( spaceRemain > indexer )
				{
					std::cout << quote.toString(0,indexer);
					quote.erase(0,indexer);
					printedWidth += indexer;
					break;
				}
				else
				{				
					std::cout << quote.toString(0,spaceRemain); // fill in available space
					std::cout << '-' << std::endl << std::string(delimitBorder,borderDelimiter); // hyphenate to show broken word
					quote.erase(0,spaceRemain);
					printedWidth = 0;
					indexer -= spaceRemain;
				}
			}
		}
		else if( indexer > width - printedWidth )
		{
			std::cout << std::endl << std::string(delimitBorder,borderDelimiter);
			printedWidth = 0;
			std::cout << quote.toString(0,indexer);
			quote.erase(0,indexer);
			printedWidth = indexer;
		}
		else
		{
			std::cout << quote.toString(0,indexer);
			quote.erase(0,indexer);
			printedWidth += indexer;
		}
	}

	if( authorString.length() != 0 )
	{
		std::cout << std::endl << std::string(delimitBorder+1,borderDelimiter) << authorString << std::endl;
	}

	// Footer line
	std::cout << borderDelimiter << borderDelimiter << std::string(delimitWidth-4,'-');
	std::cout << std::endl << std::endl;
}



