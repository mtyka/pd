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
#include "statclock.h"

const double StatClock::TO_SECONDS = 1.0 / ((double)CLOCKS_PER_SEC);
const double StatClock::TO_MILLI_SECONDS = 1.0 / ((double)(CLOCKS_PER_SEC) / 1000.0);

StatClock::StatClock():
	m_Name("NoName"),
	m_Running( false )
{
}

StatClock::StatClock( const std::string& name ):
	m_Name(name),
	m_Running( false )
{
}

void StatClock::Reset()
{
	m_Running = false;
	m_Times.clear();
}

void StatClock::Begin()
{
	AssertNotBegun();
	m_Running = true;
	m_Last = std::clock();
}

void StatClock::Stamp()
{
	AssertBegun();
	clock_t prev = m_Last;
	m_Last = std::clock();
	m_Times.push_back( m_Last - prev );
}

void StatClock::End()
{
	AssertBegun();
	clock_t prev = m_Last;
	m_Last = std::clock();
	m_Times.push_back( m_Last - prev );
	m_Running = false;
}

clock_t StatClock::getAverage()const
{
	AssertHasData();
	if( m_Times.size() == 1 ) return m_Times[0];
	double average = 0;
	for( size_t i = 0; i < m_Times.size(); i++ )
	{
		average += m_Times[i];
	}
	average /= (double) m_Times.size();
	return (clock_t) average;
}

clock_t StatClock::getStdDev()const
{
	AssertHasData();
	if( m_Times.size() == 1 ) return -1;
	double stddev = 0;
	for( size_t i = 0; i < m_Times.size(); i++ )
	{
		stddev += m_Times[i];
	}
	stddev /= (double) (m_Times.size() - 1);
	stddev = sqrt( stddev );
	return (clock_t) stddev;
}

clock_t StatClock::getStdErr()const
{
	AssertHasData();
	if( m_Times.size() == 1 ) return -1;
	double stderrmean = 0;
	for( size_t i = 0; i < m_Times.size(); i++ )
	{
		stderrmean += m_Times[i];
	}
	double size = (double) m_Times.size();
	stderrmean /= (size - 1.0);
	stderrmean = sqrt( stderrmean ) / sqrt( size );
	return (clock_t) stderrmean;
}

clock_t StatClock::getMax()const
{
	AssertHasData();
	if( m_Times.size() == 1 ) return m_Times[0];

	// 'std::numeric_limits' would be better, but unfortunatly doesn't come as standard...
	// clock_t Max = std::numeric_limits<std::clock_t>::Min; // Part of BOOST library
	// We will therefore have to assume that (clock_t == long):
	clock_t Max = (clock_t)LONG_MIN;

	for( size_t i = 0; i < m_Times.size(); i++ )
	{
		if( m_Times[i] > Max ) Max = m_Times[i];
	}

	return Max;
}

clock_t StatClock::getMin()const
{
	AssertHasData();
	if( m_Times.size() == 1 ) return m_Times[0];

	clock_t Min = (clock_t)LONG_MAX;

	for( size_t i = 0; i < m_Times.size(); i++ )
	{
		if( m_Times[i] < Min ) Min = m_Times[i];
	}

	return Min;
}

bool StatClock::HasData()const
{
	return (m_Times.size() != 0);
}

std::string StatClock::getName()const
{
	return m_Name;
}

void StatClock::ReportSeconds( int _verbosity )const
{
	Report( _verbosity, StatClock::TO_SECONDS, "seconds" );
}

void StatClock::ReportMilliSeconds( int _verbosity )const
{
	Report( _verbosity, StatClock::TO_MILLI_SECONDS, "milli seconds" );
}

void StatClock::Report( int _Verbosity, double _TimeFactor, std::string _TimeName )const
{
	_Verbosity = 0;// we haven't coded any other types yet.
	// 1 = Show some data
	// 2 = Show a histogram ??
	switch( _Verbosity )
	{
	case 1:
		{
			break;
		}
	case 0:
	default:
		{
			std::string name = m_Name;
			if( name.length() == 0 )
			{
				name = "NoName";
			}
			if( m_Times.size() == 1 )
			{
				std::cout << std::string("'") << name << std::string("' Time: ");
				std::cout << ((double)m_Times[0] * _TimeFactor) << _TimeName << std::endl;
			}
			else
			{
				std::cout << std::string("'") << name << std::string("' Average Time: ");
				std::cout << ((double)getAverage() * _TimeFactor) << ' ' << (char)241;
				std::cout << (double)getStdDev() * _TimeFactor << ' ' << _TimeName << std::endl;
			}
			break;
		}
	}
}

void StatClock::AssertBegun()const
{
	if( !m_Running )
	{
		throw new ProcedureException("The clock is already ticking!");
	}
}

void StatClock::AssertNotBegun()const
{
	if( m_Running )
	{
		throw new ProcedureException("The clock isn't ticking yet!");
	}
}

void StatClock::AssertHasData()const
{
	if( !HasData() )
	{
		throw new ProcedureException("The StatWatch has no data yet!");
	}
}

TextProgressBar::TextProgressBar(int _Total, int _PrintWidth):
m_Ticker( 0 ),
m_Current(0),
m_Total(_Total)
{
	m_PrintWidth = _PrintWidth;
	if( m_PrintWidth < 10 )
	{
		m_PrintWidth = 10;
	}
}

void TextProgressBar::Begin()
{
	// Fill in the relevent characters to the output stream...
	DrawProgress(); // our first draw...
}

void TextProgressBar::End()
{
	clearProgress();
}

void TextProgressBar::Reset()
{
	m_Current = 0;
}

void TextProgressBar::Reset(int _Total)
{
	m_Current = 0;
	m_Total = _Total;
}

void TextProgressBar::Reset(int _Total, int _PrintWidth)
{
	m_Current = 0;
	m_Total = _Total;
	m_PrintWidth = _PrintWidth;
}

void TextProgressBar::next()
{
	m_Current ++;
	if( m_Current > m_Total ) m_Current = m_Total;
	clearProgress();
	DrawProgress();
}

void TextProgressBar::clearProgress()
{
	for( int i = 0; i < m_PrintWidth; i++ )
	{
		std::cout << '\b';
	}
	for( int i = 0; i < m_PrintWidth; i++ )
	{
		std::cout << ' ';
	}
	for( int i = 0; i < m_PrintWidth; i++ )
	{
		std::cout << '\b';
	}
}

void TextProgressBar::DrawProgress()
{
	const int barOverhead = 11; // '|' + '|' + ' (|) 100%'

	// print Bar
	std::cout << '|';

	int drawCount = (int)( ((double)m_Current / (double)m_Total ) * (double)(m_PrintWidth-barOverhead));
	for( int i = 0; i < drawCount; i++ )
	{
		std::cout << (char)178;
	}
	for( int i = drawCount; i < m_PrintWidth-barOverhead; i++ )
	{
		std::cout << (char)176;
	}
	std::cout << '|';

	// print Ticker
	std::cout << " (";
	switch( m_Ticker )
	{
	case 0:
		std::cout << '-';
		break;
	case 1:
		std::cout << '\\';
		break;
	case 2:
		std::cout << '|';
		break;
	default:
		std::cout << '/';
		break;
	}
	m_Ticker++;
	if( m_Ticker >= 4 ) m_Ticker = 0;
	std::cout << ") ";

	// print Percentage
	int percentage = (int)( ((double)m_Current / (double)m_Total) * 100.0 );
	if( percentage < 10 ) std::cout << ' ';
	if( percentage < 100 ) std::cout << ' ';
	std::cout << percentage;
	std::cout << '%';
}

