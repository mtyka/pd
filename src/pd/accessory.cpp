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

#include <signal.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>

// Disable API __delcspec definitions
// Required define to allow include of any MMLib headers here ...
#define PD_API 

#include "tools/quote.h"
#include "revision.h"
#include "accessory.h"

void PrintFullpdHeader()
{
	Printpd();		
	PrintQuote();
	info();
}

void Printpd()
{
	std::cout << "PD" << std::endl;
		<< std::endl;
}

void PrintQuote()
{
	Quote* q = Quote::getSingleton();
	q->printQuote();
}

void cseed(unsigned int newseed)
{
	printf("setting random seed to %d \n",newseed);
	srand(newseed);
}

void ctimeseed()
{
	unsigned int newseed = (unsigned int)time(NULL) ;
	printf("setting random seed to %d using current system time\n",newseed);
	srand(newseed);
}

void cprint(std::string text)
{
	printf("%s\n",text.c_str());
	fflush(stdout);
}

void cflush()
{
	fflush(stdout);
}

void detail()
{
	Printpd();		
	PrintQuote();
	info();
}

void info()
{
	printf("PD %s %s %s\n",
		Revision::LibVersion().c_str(),
		Revision::CopyRight().c_str(),
		Revision::LongRevString().c_str());
	printTimeStamp();
	printHostname();
	printWorkingDirectory();
	printf("\n");
	fflush(stdout);
}

void timer(bool reset){
	static long starttime = 0;

	if((starttime == 0)){
		starttime = (long)time(NULL);
		printf("Timer initialized.\n");
		return;
	}
	if(reset){
		starttime = (long)time(NULL);
		printf("Timer reset.\n");
		return;
	}
	long timediff = (long)time(NULL) - starttime;

	printf("Timer: %ld sec (%d hrs %d mins %d secs) \n",
		timediff,
		(timediff)/3600,
		((timediff)%3600)/60,
		((timediff)%3600)%60 );
}

#ifdef WIN32
#include "winsock.h"

void printHostname(){
	printf("System: Windows ");
	WSADATA p;
	WSAStartup ((2<<8) | 2, &p);
	char hostname[256];
	if(gethostname(&hostname[0],256)!=0){
		sprintf(&hostname[0],"unknown");
	}
	printf("Hostname: %s \n",&hostname[0]);
}

#elif __linux__
#include <unistd.h>
void printHostname(){
	printf("System: Linux ");
	char hostname[256];
	if(gethostname(&hostname[0],256)!=0){
		sprintf(&hostname[0],"unknown");
	}
	printf("Hostname: %s \n",&hostname[0]);
}
#else

void printHostname(){
	printf("System: Unknown Hostname: Unknown \n");
}

#endif

// _MAX_PATH is a platform dependent value for the maximum length of a system file path
#ifdef WIN32
	#include <direct.h>
#else
	#ifndef _MAX_PATH
		#define _MAX_PATH 260
	#endif
#endif


void printWorkingDirectory()
{
	char appPathBuffer[_MAX_PATH+1];
	getcwd(appPathBuffer,_MAX_PATH); // get the current working directory of the executable
	strcat(appPathBuffer,"/");
	for(int i=0;i<strlen(appPathBuffer);i++){
		if(appPathBuffer[i] == '\\') appPathBuffer[i] = '/';
	}
	printf( "Current working directory: %s\n", appPathBuffer );
}

void printTimeStamp(){
	tm *ltime;
	time_t stime;

	// put a timestamp on this file
	stime = (int) time(NULL);
	ltime = localtime(&stime);
	printf("Compile Date and Time: %s %s\n", __DATE__, __TIME__);
	printf("Execution Date and Time: %s", asctime(ltime));
}

void printTitle( const std::string &_Title )
{
	// Character Mapping:
#ifdef WIN32
	char a = (char)205;
	char s = (char)186;
	char tr = (char)187;
	char tl = (char)201;
	char bl = (char)200;
	char br = (char)188;
#elif __linux__
	char a = (char)150;
	char s = (char)124;
	char tr = (char)92;
	char tl = (char)47;
	char bl = (char)92;
	char br = (char)47;
#else // Who Knows :-D
	char a = (char)150;
	char s = (char)124;
	char tr = (char)92;
	char tl = (char)47;
	char bl = (char)92;
	char br = (char)47;
#endif

	// Top Line
	std::cout << std::endl << ' ' << tl;
	for( size_t i = 0; i < _Title.length() + 2; i++ )
	{
		std::cout << a;
	}
	std::cout << tr << ' ' << std::endl;

	// Middle Line
	std::cout << ' ' << s << ' ' << _Title << ' ' << s << ' ' << std::endl;

	// Bottom Line
	std::cout << ' ' << bl;
	for( size_t i = 0; i < _Title.length() + 2; i++ )
	{
		std::cout << a;
	}
	std::cout << br << ' ' << std::endl << std::endl;
}







