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
#include "fastrandom.h"

namespace Maths
{
	// -------------------------------------------
	// Begin 'FastRandom' Singleton Implementation
	// -------------------------------------------
	FastRandom* FastRandom::GlobalFastRand = NULL; // initialise to 'NULL'
	FastRandom* FastRandom::getInstance()
	{
		if( GlobalFastRand == NULL )
		{
			GlobalFastRand = new FastRandom();
		}
		return GlobalFastRand;
	}
	// -------------------------------------------
	// End
	// -------------------------------------------

	FastRandom::FastRandom() : bitBufferIdx(32),
		REAL_UNIT_INT(1.0/((double)INT_MAX+1.0)), // The +1 ensures nextDouble doesn't generate 1.0
		REAL_UNIT_UINT(1.0/((double)UINT_MAX+1.0))
	{
		reinitialise((int)time(0));
	}

	FastRandom::FastRandom(int seed) : bitBufferIdx(32),
		REAL_UNIT_INT(1.0/((double)INT_MAX+1.0)), // The +1 ensures nextDouble doesn't generate 1.0
		REAL_UNIT_UINT(1.0/((double)UINT_MAX+1.0))
	{
		reinitialise(seed);
	}

	FastRandom::FastRandom(const FastRandom &_Clone) : 
		REAL_UNIT_INT(1.0/((double)INT_MAX+1.0)),
		REAL_UNIT_UINT(1.0/((double)UINT_MAX+1.0))
	{
		(*this) = _Clone;
	}

	FastRandom& FastRandom::operator=(const FastRandom &_Clone)
	{
		x = _Clone.x;
		y = _Clone.y; 
		z = _Clone.z; 
		w = _Clone.w;
		bitBuffer = _Clone.bitBuffer;
		bitBufferIdx = _Clone.bitBufferIdx; 
		return (*this);
	}

	FastRandom::~FastRandom()
	{
	}

	void FastRandom::reinitialise(int seed)
	{
		// The only stipulation stated for the xorshift RNG is that at least one of
		// the seeds x,y,z,w is non-zero. We fulfill that requirement by only allowing
		// resetting of the x seed
		// Jon Comment: These three numbers doublely can be anything. Ive tested, the random
		// generator still works fine.
		x = (unsigned int)seed;
		y = (unsigned)842502087;
		z = (unsigned)357987591; // these three declatations were class PD_API constants:
		w = (unsigned)273326509; // Y, Z, W respectively in the original implementation
	}

	// Generates a unsigned int. Values returned are over the full range of a unsigned int,
	// unsigned int.MinValue to unsigned int.MaxValue, including the min and max values.
	unsigned int FastRandom::nextUInt()
	{
		unsigned int t=(x^(x<<11));
		x = y; y = z; z = w;
		return (w=(w^(w>>19))^(t^(t>>8)));
	}

	// Generates a random int. Values returned are over the range 0 to int.MaxValue-1.
	// MaxValue is not generated to remain functionally equivalent to System.Random.next().
	// If you require an int from the full range, including negative values then call
	// nextUint() and cast the value to an int.
	int FastRandom::next()
	{
		unsigned int t = (x^(x<<11));
		x = y; y = z; z = w;
		w=(w^(w>>19))^(t^(t>>8));
		// Handle the special case where the value INT_MAX is generated. This is outside of 
		// the range of permitted values, so we therefore call next() to try again.
		unsigned int rtn = w & 0x7FFFFFFF;
		if(rtn==0x7FFFFFFF)
			return next();
		return (int) rtn;
	}

	// Generates a random int over the range 0 to upperBound-1, and not including upperBound.
	int FastRandom::next(int upperBound)
	{
		D_ASSERT((upperBound > 0),CodeException,"UpperBound must be >=0");
		unsigned int t=(x^(x<<11));
		x = y; y = z; z = w;

		// The explicit int cast before the first multiplication gives better performance.
		// See comments in nextDouble.
		return (int)((REAL_UNIT_INT*(int)(0x7FFFFFFF&(w=(w^(w>>19))^(t^(t>>8)))))*upperBound);
	}

	/// Generates a random int over the range lowerBound to upperBound-1, and not including upperBound.
	/// upperBound must be >= lowerBound. lowerBound may be negative.
	int FastRandom::next(int lowerBound, int upperBound)
	{
		unsigned int t=(x^(x<<11));
		x = y; y = z; z = w;

		// The explicit int cast before the first multiplication gives better performance.
		// See comments in nextDouble.
		int range = upperBound-lowerBound;
		if(range<0)
		{// If range is <0 then an overflow has occured and must resort to using long integer arithmetic instead (slower).
			// We also must use all 32 bits of precision, instead of the normal 31, which again is slower.
			return lowerBound+(int)((REAL_UNIT_UINT*(double)(w=(w^(w>>19))^(t^(t>>8))))*(double)((long)upperBound-(long)lowerBound));
		}

		// 31 bits of precision will suffice if range<=int.MaxValue. This allows us to cast to an int and gain
		// a little more performance.
		return lowerBound+(int)((REAL_UNIT_INT*(double)(int)(0x7FFFFFFF&(w=(w^(w>>19))^(t^(t>>8)))))*(double)range);
	}

	/// Generates a random double. Values returned are from 0.0 up to but not including 1.0.
	double FastRandom::nextDouble()
	{
		unsigned int t=(x^(x<<11));
		x = y; y = z; z = w;

		// Here we can gain a 2x speed improvement by generating a value that can be cast to
		// an int instead of the more easily available unsigned int. If we then explicitly cast to an
		// int the compiler will then cast the int to a double to perform the multiplication,
		// this final cast is a lot faster than casting from a unsigned int to a double. The extra cast
		// to an int is very fast (the allocated bits remain the same) and so the overall effect
		// of the extra cast is a significant performance improvement.
		return (REAL_UNIT_INT*(int)(0x7FFFFFFF&(w=(w^(w>>19))^(t^(t>>8)))));
	}

	// Generates random bool.
	// Increased performance is achieved by buffering 32 random bits for
	// future calls. Thus the random number generator is only invoked once
	// in every 32 calls.
	bool FastRandom::nextBool()
	{
		if(bitBufferIdx==32)
		{
			// Generate 32 more bits.
			unsigned int t=(x^(x<<11));
			x = y; y = z; z = w;
			bitBuffer = w=(w^(w>>19))^(t^(t>>8));

			// Reset the idx that tells us which bit to read next.
			bitBufferIdx = 1;
			return (bitBuffer & 0x1)==1;
		}

		bitBufferIdx++;
		return ((bitBuffer>>=1) & 0x1)==1;
	}
} // namespace Maths


