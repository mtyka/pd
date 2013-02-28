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
#include "streamwriter.h"

StreamWriter::StreamWriter() :
  	m_Mode(Screen),
	m_Streaming(false),
	m_Appending(false),
	m_FileName(""),
	m_Stream(NULL),
	m_FileStream(NULL)
{
	streamToScreen();
}

StreamWriter::StreamWriter( const std::string& _filename, bool _Append ) :
  	m_Mode(Screen),
	m_Streaming(false),
	m_Appending(_Append),
	m_FileName(""),
	m_Stream(NULL),
	m_FileStream(NULL)
{
	streamTo(_filename,_Append);
}

StreamWriter::StreamWriter( std::ostream& _stream, bool _Append ) :
  	m_Mode(Screen),
	m_Streaming(false),
	m_Appending(_Append),
	m_FileName(""),
	m_Stream(NULL),
	m_FileStream(NULL)
{
	streamTo(_stream);
}

void StreamWriter::assertStream() const
{
	ASSERT(isStreaming(),CodeException,"Stream open!");
}

void StreamWriter::streamToScreen()
{
	ASSERT( !m_Streaming, CodeException, "Stream currently open!!");
	m_Mode = Screen;
	m_Stream = &std::cout;
}

void StreamWriter::streamTo( const std::string& _filename, bool _Append, bool _ClearFile )
{
	ASSERT( !m_Streaming, CodeException, "Stream currently open!!");

	m_Mode = File;
	m_Stream = NULL;
	m_FileStream = NULL;
	m_FileName = _filename;
	m_Appending = _Append;

	if( _ClearFile )
	{	
		std::ofstream file;
		try
		{
			// Remove the current contents of the file.
			file.open(_filename.c_str(), std::ios::out);
			file.close();
		}
		catch( std::exception ex )
		{
			file.close();
			throw ex;
		}
	}
}

void StreamWriter::streamTo( std::ostream& _stream )
{
	ASSERT( !m_Streaming, CodeException, "Stream currently open!!");
	m_Mode = Stream;
	m_Stream = &_stream;
}

void StreamWriter::beginStream()
{
	ASSERT( !m_Streaming, CodeException, "Already Streaming!!");
	m_Streaming = true;
	switch( m_Mode )
	{
	case File:
		{
			ASSERT( m_FileStream == NULL, CodeException, "Assumption failed - m_FileStream is not NULL!!");
			m_FileStream = new std::ofstream();
			m_Stream = m_FileStream;
			if( m_Appending )
			{
				m_FileStream->open(m_FileName.c_str(), std::ios::app);
			}
			else
			{
				m_FileStream->open(m_FileName.c_str(), std::ios::out);
			}
			if( !m_FileStream->is_open() )
			{
				THROW(IOException,"Could not open the filestream for writing!");
			}
			return;
		}
	case Stream:
		{
			// Nothing required - stream is under external control.
			return;
		}
	case Screen:
		{
			// Nothing required - stream is under external control.
			return;
		}
	default:
		{
			THROW(CodeException,"Unknown 'InternalStream' type!");
		}
	}
}

void StreamWriter::endStream()
{
	ASSERT( m_Streaming, CodeException, "Stream is not open!!");
	m_Streaming = false;
	switch( m_Mode )
	{
	case File:
		{
			ASSERT( m_FileStream != NULL, CodeException, "Internal member m_FileStream is not meant to be NULL!");
			if( m_FileStream->is_open() )
			{
				m_FileStream->close(); // release the file handle!
			}
			delete m_FileStream;
			m_FileStream = NULL;
			m_Stream = NULL;
			return;
		}
	case Stream:
		{
			// Nothing required - stream is under external control.
			return;
		}
	case Screen:
		{
			// Nothing required - stream is under external control.
			return;
		}
	default:
		{
			THROW(CodeException,"Unknown 'InternalStream' type!");
		}
	}
}

