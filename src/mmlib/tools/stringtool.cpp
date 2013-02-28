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
#include <vector>
#include <algorithm> // so that we can use transform() to make toupper() work
#include "tools/stringtool.h"

// add this to an integer in the range 0-9 to get a screen printable ASCII character ...
char getScreenChar( int singleDigitNumber )
{
	if( singleDigitNumber < 0 || singleDigitNumber > 9 ) return '?';
	return (char) singleDigitNumber + 48;
}

char getScreenChar( unsigned int singleDigitNumber )
{
	if( singleDigitNumber > 9 ) return '?';
	return (char) singleDigitNumber + 48;
}

char getScreenChar( short singleDigitNumber )
{
	if( singleDigitNumber < 0 || singleDigitNumber > 9 ) return '?';
	return (char) singleDigitNumber + 48;
}

char getScreenChar( unsigned short singleDigitNumber )
{
	if( singleDigitNumber > 9 ) return '?';
	return (char) singleDigitNumber + 48;
}

char getScreenChar( char singleDigitNumber )
{
	if( singleDigitNumber < 0 || singleDigitNumber > 9 ) return '?';
	return (char) singleDigitNumber + 48;
}

char getScreenChar( unsigned char singleDigitNumber )
{
	if( singleDigitNumber > 9 ) return '?';
	return (char) singleDigitNumber + 48;
}

bool isSameCharIgnoringCase( char a, char b )
{
	// A to Z == 65 to 90
	// a to z == 97 to 122
	if( a == b ) return true;
	if( abs(a-b) != (97-65) ) return false;
	return (((a >= 65 && a <= 90 ) && (b >= 97 && b <= 122 )) ||
		((b >= 65 && b <= 90 ) && (a >= 97 && a <= 122 )) );
}

bool isWhitespace(char ch)
{
	if(ch == ' ')
		return true;
	if(ch == '\t')
		return true;
	return false;
}

bool isZeroToNine(char ch)
{
	return ((ch >= 0) && (ch <= 9));
}

bool isASCIINumeric(char ch)
{
	// 0 to 9 == 48 to 57
	return ((ch >= '0') && (ch <= '9'));
}

void strlc(char *text)
{
	int i = 0;
	while(text[i] != 0) {
		if((text[i] >= 'A') && (text[i] <= 'Z'))
			text[i] -= 'A' - 'a';
		i++;
	}
}

void struc(char *text)
{
	int i = 0;
	while(text[i] != 0) {
		if((text[i] >= 'a') && (text[i] <= 'z'))
			text[i] += 'A' - 'a';
		i++;
	}
}

void truncLF(char *buffer)
{
	if(strlen(buffer) < 1)
		return;
	while((buffer[strlen(buffer) - 1] == 10) || (buffer[strlen(buffer) - 1] == 13)) {
		buffer[strlen(buffer) - 1] = 0;
	}
}

int chopString(char *data, int *beginning, int *end, int maxelements)
{
	int total = 0;
	int i = 0, j = 0;

	if(!isWhitespace(data[0])) {
		beginning[total] = i;
		goto chopString_jumpin; // goto's are nasty, but kind of neat in some low level functions.
	}
	for(i = 0; i < (int) (strlen(data) - 1); i++) {
		if(isWhitespace(data[i]) && (!isWhitespace(data[i + 1]))) {
			beginning[total] = i;

chopString_jumpin:
			// now find the corresponding end
			for(j = i; j < (int) (strlen(data) - 1); j++) {
				if(isWhitespace(data[j + 1]) && (!isWhitespace(data[j]))) {

					end[total] = j;
					break;
				}
			}
			if(!(j < (int) (strlen(data) - 1))) {
				end[total] = (int) strlen(data) - 1;
				total++;
				return total;
			}
			i = j;
			total++;
			if(total >= maxelements)
				return total;
		}
	}
	return total;
}

// Functions for STL-style (C++) strings

// Renamed from getSingleLine() as there is also a Templated std::getline() causing confustion
int getSingleLine(FILE *file, std::string& retline)
{
	char buffer[65535];
	if(fgets(&buffer[0],65535,file)==NULL) return -1;
	retline = &buffer[0];
	return 0;
}

void makeupper(std::string &_input)
{
	transform(_input.begin(),_input.end(),_input.begin(),toupper);
}

void makelower(std::string &_input)
{
	transform(_input.begin(),_input.end(),_input.begin(),tolower);
}

std::string ltrim( const std::string &str, const std::string &delimiters)
{
	int idx = str.find_first_not_of(delimiters);
	if( idx != std::string::npos )
		return str.substr(idx);
	return "";
}

std::string rtrim( const std::string &str, const std::string &delimiters)
{
	int idx = str.find_last_not_of(delimiters);
	if( idx != std::string::npos )
		return str.substr(0,idx+1);
	return str;
}

std::string trim( const std::string &str, const std::string &delimiters)
{
	return rtrim(ltrim(str,delimiters),delimiters);
}

std::string makeStringOfLength( const std::string& donorString, int desiredLength, char Padding, int offset )
{
	if( donorString.length() > desiredLength )
	{
		// truncate to 'desiredLength'
		return donorString.substr(offset,desiredLength);
	}
	std::string returnString = donorString;
	for( int i = returnString.length(); i < desiredLength; i++ )
	{
		// extend the current string using the delimiter
		returnString.push_back(Padding);
	}
	return returnString;
}

//returns true is the two strings are equal, false otherwise
bool cmpstring(const std::string& str1,const std::string& str2)
{
	return str1.compare(str2) == 0;
}
//returns true is the two strings are equal, false otherwise
bool cmpstring(const std::string& str1,const char* str2)
{
	return str1.compare(str2) == 0;
}

std::vector<std::string> chopstr(const std::string& inputstring, const char* delim)
{
	std::vector<std::string> token;
	if(inputstring.size() == 0) return token; // return empty container if string's empty

	size_t length = inputstring.size();
	size_t pos = 0;
	size_t newpos = 0;
	
	while(pos<length)
	{
		pos = inputstring.find_first_not_of(delim,pos);
		if(pos==std::string::npos) break;
		newpos = inputstring.find_first_of(delim,pos);
		if(newpos==std::string::npos) newpos = length;
		token.push_back( inputstring.substr(pos,newpos-pos) );
		pos = newpos;
	}

	return token;
}

std::vector<std::string> chopstr(const std::string& inputstring,
								 const char* delim,
								 const char* control,
								 const std::string& leftbracket,
								 const std::string& rightbracket)
{
	ASSERT(leftbracket.size() == rightbracket.size(), CodeException, "Bracket strings are of unequal length!" );

	std::string controlorbracket;
	controlorbracket = "";
	controlorbracket += control;
	controlorbracket += leftbracket;

	std::vector<std::string> token;

	size_t length = inputstring.size();
	size_t pos = 0;
	size_t posdelim,poscontrol;
	size_t newpos = 0;
	if(inputstring=="") return token; // return empty container if string's empty
	while(pos<length){
		posdelim = inputstring.find_first_not_of(delim,pos); // stop at next nonwhitespace
		pos = posdelim;
		if(pos==std::string::npos) break;

		posdelim = inputstring.find_first_of(delim,pos); // stop at next non whitespace
		poscontrol = inputstring.find_first_of(controlorbracket,pos); // or control character
		if((posdelim!=std::string::npos)&&
			(poscontrol==std::string::npos)){
				newpos = posdelim; // if we can't fins another control character just proceed to next whitespace
		}else if((posdelim==std::string::npos)&&
			(poscontrol!=std::string::npos)){
				newpos = poscontrol; // equally if we can't find another whitespace character just proceed to next deliminator
		}else if((posdelim!=std::string::npos)&&(poscontrol!=std::string::npos)){
			newpos = posdelim < poscontrol ? posdelim : poscontrol; // otherwise go to what ever comes first // otherwise go to next deliminator;
		}else{
			newpos = length; // if there are not whitespaces or deliminators left proceed to end
		}

		if(poscontrol == pos){// if we are currently on a control paraemter then just move on one
			// is the control character a bracket ??
			unsigned ichar;
			for(ichar=0;ichar<leftbracket.size();ichar++){
				if(leftbracket[ichar] == inputstring[pos]) break; // break if we have a match
			}
			// ichar now contains the index of the Type of bracket or
			// ichar = leftbracket.size() if no bracket char was found

			if(ichar >= leftbracket.size()){ // not a bracket, so just move on 1 character
				newpos = pos + 1;
			}else{// we hit on a bracket
				// now find the equivalent rightbracket of the leftbracket
				std::string rightbracketchar = rightbracket.substr(ichar,1);
				size_t bracketpos = inputstring.find_first_of(rightbracketchar,pos+1);
				if(bracketpos==std::string::npos) newpos = length; // go to end - throw an error cos we expected a bracket
				else newpos = bracketpos+1;
			}

		}
		token.push_back( inputstring.substr(pos,newpos-pos) );
		//printf("tk:'%s'\n",token[token.size()-1].c_str());
		pos = newpos;
	}
	return token;
}

std::string firsttoken(const std::string& inputstring, const char* delim)
{
	std::string token = "";

	size_t length = inputstring.size();
	size_t pos = 0;
	size_t newpos = 0;
	if(inputstring=="") return token; // return empty container if string's empty

	pos = inputstring.find_first_not_of(delim,pos);
	if(pos==std::string::npos) return token;
	newpos = inputstring.find_first_of(delim,pos);
	if(newpos==std::string::npos) newpos = length;
	token = inputstring.substr(pos,newpos-pos);
	pos = newpos;

	return token;
}

void removecomments(std::string& inputstring, const char *comment)
{
	size_t firstcomment;

	firstcomment = inputstring.find_first_of(comment,0);
	if(firstcomment==std::string::npos) return ; // no comment found

	inputstring = inputstring.substr(0,firstcomment); // truncate string to before comment
}

void removecomments(std::string& inputstring, const char *comment, std::string& commentline)
{
	size_t firstcomment;

	firstcomment = inputstring.find_first_of(comment,0);
	if(firstcomment==std::string::npos) return ; // no comment found

	commentline = inputstring.substr(firstcomment); // truncate string to after comment
	inputstring = inputstring.substr(0,firstcomment); // truncate string to before comment
}


int str2float(const std::string& inputstring, float& value)
{
	float lvalue;

	if(sscanf(inputstring.c_str(),"%f",&lvalue)!=1) return -1;
	value = lvalue;
	return 0;
}

int str2double(const std::string& inputstring, double& value)
{
	double lvalue;

	if(sscanf(inputstring.c_str(),"%lf",&lvalue)!=1) return -1;
	value = lvalue;
	return 0;
}

int str2int(const std::string& inputstring, int& value)
{
	int lvalue;

	if(sscanf(inputstring.c_str(),"%d",&lvalue)!=1) return -1;
	value = lvalue;
	return 0;
}

std::string float2str(float value,char *formatstring)
{
	char buffer[256];
	sprintf(&buffer[0],formatstring,value);
	return std::string(&buffer[0]);
}

std::string double2str(double value,char *formatstring)
{
	char buffer[256];
	sprintf(&buffer[0],formatstring,value);
	return std::string(&buffer[0]);
}

std::string int2str(int value,char *formatstring)
{
	char buffer[256];
	sprintf(&buffer[0],formatstring,value);
	return std::string(&buffer[0]);
}

std::string long2str(long value,char *formatstring)
{
	char buffer[256];
	sprintf(&buffer[0],formatstring,value);
	return std::string(&buffer[0]);
}

