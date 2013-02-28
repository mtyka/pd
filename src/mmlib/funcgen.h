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

#ifndef __TGEN_H
#define __TGEN_H




// Template Version Generator for 3 template integer parameters
//
// Imagine the user has a multiply templated, generative function like this:
//  
//  template <int a, int b, int c>
//  void func(int x, int y, int z){
//    .... 
//  }
//
// The purpose of the ints a,b,c is to modify the behaviour of the function
// For example they may turn on/off blocks of code
//
// Unless the template paramters are known at compile time, calling this function
// would require multiple, nested if/switch statements which are hard to maintain.
// instead the user may use an automatic generation of all the possible template parameter
// combinations. To so he must define the following class in **exactly** this way.
// 
//
//  template <int a, int b, int c>         // <-- same template as the function
//  struct func_wrap{                      // <-- choose a name, but postfixing with wrap is good
//  	func_wrap():ptr(&func<a, b, c>){}    // <-- this sets the function pointer ptr to point to the 
//                                         //     function with the template params a,b and c
//  	typedef void (*Tptr)(int,int,int);   // <-- This is a typedef which defines the function pointer type
//                                         //     note that the param list must match that of func(..)
//  	Tptr ptr;                            // <-- no need to change this line
//  	const static int limita=2;           // <-- set the limit of the template paramter a
//  	const static int limitb=2;           // <-- ditto
//  	const static int limitc=2;           // <-- ditto
//  };
//
//
//  Having done this the user may then call his function in the following way 
//    int mya, myb, myc, myx, myy, myz;
//    tcall3<func_wrap>(mya, myb, myc)( myx, myy, myz );
//  Note that mya, myb and myc are NOT known at compile time !!!
//
//  Alternate forms tcall1,tcall2, tcall4 and tcall5 are available for different lengths of template
//  lists. Longer ones are easily generatable, but remember that the number of generated functions
//  is limita*limitb*...*limitx , which grows VERY quickly and may bloat your code ! (of course this
//  depends on how big func(..) is !
//


// Main template definition
template <
	// this template parameter defines what the function is to be called and what the 
	// the parameter limits are. 
	template<int,int,int> class G,

	// these are the limits of the parameters, i.e. ta may run from 0 to limita, tb may run from 0 to limitb, ...
	int limita=-1,
	int limitb=-1,
	int limitc=-1,
	
	// these are the current template parameters
	int ta=0, 
	int tb=0,
	int tc=0
>
struct tcall3
{
	// this just makes a short cut to the function pointer type that we are calling.
	// it is defined in the user-supplied template class G
	typedef typename G<ta,tb,tc>::Tptr T;

	// This is the constructor. It compares the template parameters ta,tb,.. with the 
	// normal parameters a,b,c, .. 
	// if they match it sets the internal fptr to the appropriate function pointer, 
	// namely G<ta,tb,tc>().ptr
	// if not it counts up on of the template parameters and **calls itself** to try again
	// in this way ta counts up until it hits ta==a, then tb gets counted up until tb==b
	// finally then fc==c the first three if statements become false and the function pointer
	// is set.
	// this pointer is then passed back through the chain of callers and callees until the last
	// class ends up with having determined the right pointer.
	// because we have overloaded the type conversion operator T
	// one can call this class like this:
	// fptr = tcall3<G,...>(a,b,c)
	// or, calling the resultant function directly
	// tcall3<G,...>(a,b,c)(x,y,z)  where (x,y,z) are just the normal function parameters
	// passed to this function.
	tcall3(int a, int b, int c):fptr(G<0,0,0>().ptr){
		//printf("%d %d   %d %d\n",ta,tb,tc,a,b,c);
				 if(ta<a){ fptr=tcall3<G,limita,limitb,limitc,ta+1,tb,tc>(a,b,c); }
		else if(tb<b){ fptr=tcall3<G,limita,limitb,limitc,ta,tb+1,tc>(a,b,c); }
		else if(tc<c){ fptr=tcall3<G,limita,limitb,limitc,ta,tb,tc+1>(a,b,c); }
		else fptr=G<ta,tb,tc>().ptr;  // we've chosen our template parameters - set a pointer to the right function
	}
	// this is just for convenience and short hand writing.
	operator T(){ return fptr; }
private:
	// this stores the current pointer to the function.
	T fptr;
};

// These function are template specialisations of the above which dont call them selves any more.
// In fact they simply return a dummy pointer G<0,0,0>().ptr .
// The purpose of this is to make the compiler STOP recursing templates indefinitely.
// when it reaches the limits of each of the parameters, it instantiates one of these special 
// "stopper" specializations and thus this terminates the compiler's recursion.

// there is One stoper per parameter. note that the instantiation simply replaces
// the paramter in question with the limit (which is just an integer rather than a variable)

template <template<int,int,int> class G,					
					int limita,
					int limitb,
					int limitc,
					int tb,
					int tc>
struct tcall3<G,limita,limitb,limitc,limita,tb,tc>
{
	typedef typename G<limita,tb,tc>::Tptr T;
	tcall3(int a, int b, int c){} operator T(){ G<0,0,0>::overflow();  return G<0,0,0>().ptr; } 
};
template <template<int,int,int> class G,					
					int limita,
					int limitb,
					int limitc,
					int ta,
					int tc>
struct tcall3<G,limita,limitb,limitc,ta,limitb,tc>
{
	typedef typename G<ta,limitb,tc>::Tptr T;
	tcall3(int a, int b, int c){} operator T(){ G<0,0,0>::overflow();  return G<0,0,0>().ptr; } 
};
template <template<int,int,int> class G,					
					int limita,
					int limitb,
					int limitc,
					int ta,
					int tb>
struct tcall3<G,limita,limitb,limitc,ta,tb,limitc>
{
	typedef typename G<ta,tb,limitc>::Tptr T;
	tcall3(int a, int b, int c){} operator T(){ G<0,0,0>::overflow();  return G<0,0,0>().ptr; } 
};

// This is a special specialization: it corresponds to all the default parameters
// most importantly all limits set to -1. Instead it reads the limits out of the
// user-supplied class G itself!
// this allows the user to call the whole thing like this:
//  tcall3<G>(a,b,c)(x,y,z) - without having to supply the limits.
// of course one can skip this simplification step an call the
// main class directly by supplying limits,
//  tcall3<G,2,2,2>(a,b,c)(x,y,z)
// - although this is dangerous (wrong limits could be supplied) and thus is not recommended
// (unfortunately i have not thouhght of a way to prevent the user from doing that)

template <template<int,int,int> class G>
struct tcall3<G,-1,-1,-1,0,0,0>
{
	typedef typename G<0,0,0>::Tptr T;
	tcall3(int a, int b, int c){
		fptr=tcall3<G, G<0,0,0>::limita+1, G<0,0,0>::limitb+1, G<0,0,0>::limitc+1>(a,b,c); 
	} 
	operator T(){ return fptr; }
 private:
	T fptr;

};





/////////////////////////////////////////////////////////////
//////////////////       tcall4         /////////////////////
/////////////////////////////////////////////////////////////




// Main template definition
template <
	template<int,int,int,int> class G,
	int limita=-1,
	int limitb=-1,
	int limitc=-1,
	int limitd=-1,
	int ta=0, 
	int tb=0,
	int tc=0,
	int td=0
>
struct tcall4
{
	typedef typename G<ta,tb,tc,td>::Tptr T;

	tcall4(int a, int b, int c, int d):fptr(G<0,0,0,0>().ptr){
		//printf("%d %d   %d %d\n",ta,tb,tc,a,b,c);
				 if(ta<a){ fptr=tcall4<G,limita,limitb,limitc,limitd,ta+1,tb,tc,td>(a,b,c,d); }
		else if(tb<b){ fptr=tcall4<G,limita,limitb,limitc,limitd,ta,tb+1,tc,td>(a,b,c,d); }
		else if(tc<c){ fptr=tcall4<G,limita,limitb,limitc,limitd,ta,tb,tc+1,td>(a,b,c,d); }
		else if(td<d){ fptr=tcall4<G,limita,limitb,limitc,limitd,ta,tb,tc,td+1>(a,b,c,d); }
		else fptr=G<ta,tb,tc,td>().ptr;  
	}
	operator T(){ return fptr; }
private:
	T fptr;
};

template <template<int,int,int,int> class G,					
					int limita,
					int limitb,
					int limitc,
					int limitd,
					int tb,
					int tc,
					int td>
struct tcall4<G,limita,limitb,limitc,limitd,  limita,tb,tc,td>
{
	typedef typename G<0,0,0,0>::Tptr T;
	tcall4(int a, int b, int c,int d){} operator T(){ G<0,0,0,0>::overflow();  return G<0,0,0,0>().ptr; } 
};


template <template<int,int,int,int> class G,					
					int limita,
					int limitb,
					int limitc,
					int limitd,
					int ta,
					int tc,
					int td>
struct tcall4<G,limita,limitb,limitc,limitd,  ta,limitb,tc,td>
{
	typedef typename G<0,0,0,0>::Tptr T;
	tcall4(int a, int b, int c, int d){} operator T(){ G<0,0,0,0>::overflow();  return G<0,0,0,0>().ptr; } 
};


template <template<int,int,int,int> class G,					
					int limita,
					int limitb,
					int limitc,
					int limitd,
					int ta,
					int tb,
					int td>
struct tcall4<G,limita,limitb,limitc,limitd,ta,tb,limitc,td>
{
	typedef typename G<0,0,0,0>::Tptr T;
	tcall4(int a, int b, int c, int d){} operator T(){ G<0,0,0,0>::overflow();  return G<0,0,0,0>().ptr; } 
};

template <template<int,int,int,int> class G,					
					int limita,
					int limitb,
					int limitc,
					int limitd,
					int ta,
					int tb,
					int tc>
struct tcall4<G,limita,limitb,limitc,limitd,ta,tb,tc,limitd>
{
	typedef typename G<0,0,0,0>::Tptr T;
	tcall4(int a, int b, int c, int d){} operator T(){ G<0,0,0,0>::overflow();  return G<0,0,0,0>().ptr; } 
};


template <template<int,int,int,int> class G>
struct tcall4<G,-1,-1,-1,-1,0,0,0,0>
{
	typedef typename G<0,0,0,0>::Tptr T;
	tcall4(int a, int b, int c, int d){
		fptr=tcall4<G, G<0,0,0,0>::limita+1, 
		               G<0,0,0,0>::limitb+1, 
									 G<0,0,0,0>::limitc+1,
									 G<0,0,0,0>::limitd+1
		>(a,b,c,d); 
	} 
	operator T(){ return fptr; }
 private:
	T fptr;

};


















#endif
