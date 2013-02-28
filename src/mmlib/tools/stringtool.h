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

#ifndef __STRING_TOOL_H
#define __STRING_TOOL_H

#include <string>

// ----------------------------------------------------------------------------------
// --- String Handing Functions & Classes -------------------------------------------
// ----------------------------------------------------------------------------------

char PD_API getScreenChar( int singleDigitNumber );
char PD_API getScreenChar( unsigned int singleDigitNumber );
char PD_API getScreenChar( short singleDigitNumber );
char PD_API getScreenChar( unsigned short singleDigitNumber );
char PD_API getScreenChar( char singleDigitNumber );
char PD_API getScreenChar( unsigned char singleDigitNumber );

// Function for c-style strings
bool PD_API isSameCharIgnoringCase( char a, char b );
bool PD_API isWhitespace(char ch);	// tests for ' ' or tab
bool PD_API isZeroToNine(char ch); // is a single decimal numeric character
bool PD_API isASCIINumeric(char ch); // is an ASCII decimal numeric character
void PD_API strlc(char *text);	// converts a string to all lower case
void PD_API struc(char *text);	// converts a string to all UPPER case
void PD_API truncLF(char *buffer);
int PD_API chopString(char *data, int *beginning, int *end, int maxelements);	// chops a string according to whitespaces

// Functions for STL-style (C++) strings
// DEPRECTATED: FILE* is deprecated. Use std::ifsream and std::getline()
int PD_API getSingleLine(FILE *file, std::string& retline);
//returns true is the two strings are equal, false otherwise
bool PD_API cmpstring(const std::string& str1, const std::string& str2);
//returns true is the two strings are equal, false otherwise
bool PD_API cmpstring(const std::string& str1, const char* str2);

// string manipulation functions
void PD_API makeupper(std::string &_input);
void PD_API makelower(std::string &_input);
std::string PD_API ltrim( const std::string &str, const std::string &delimiters = " ");
std::string PD_API rtrim( const std::string &str, const std::string &delimiters = " ");
std::string PD_API trim( const std::string &str, const std::string &delimiters = " ");

std::string PD_API makeStringOfLength( const std::string& donorString, int desiredLength, char Padding, int offset = 0 );
std::vector<std::string> PD_API chopstr(const std::string& inputstring, const char* delim);
std::vector<std::string> PD_API chopstr(const std::string& inputstring, const char* delim, const char* control, const std::string& leftbracket = "", const std::string& rightbracket = "");
std::string PD_API firsttoken(const std::string& inputstring, const char* delim);

void PD_API removecomments(std::string& inputstring, const char *comment);
void PD_API removecomments(std::string& inputstring, const char *comment, std::string& commentline);

int PD_API str2float(const std::string& inputstring, float& value);
int PD_API str2double(const std::string& inputstring, double& value);
int PD_API str2int(const std::string& inputstring, int& value);

std::string PD_API float2str(float value,char *formatstring="%f");
std::string PD_API double2str(double value,char *formatstring="%lf");
std::string PD_API int2str(int value,char *formatstring="%d");
std::string PD_API long2str(long value,char *formatstring="%ld");

#endif

