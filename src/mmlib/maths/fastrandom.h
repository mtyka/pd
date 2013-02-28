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

#ifndef __FAST_RANDOM
#define __FAST_RANDOM

namespace Maths
{






//-------------------------------------------------
//
/// \brief Provides Random Numbers 
///
/// \details 
/// Adapted from an article on www.codeproject.com
/// http://www.codeproject.com/csharp/FastRandom.asp
/// Jon Rea, December 2005
///
/// The original owner has stated that the software is distributed under 'something like' the MIT licence:
/// http://www.opensource.org/licenses/mit-license.php
/// I quote a response in the original article:
/// 'I'm happy for the code to be reused in whatever context so if you want I can send a modified version 
/// of the source with extra license blurb and get them to update the zip file on here.'
///
/// Key points:
/// 1) Based on a simple and fast xor-shift pseudo random number generator (RNG) specified in:
/// Marsaglia, George. (2003). Xorshift RNGs.
/// http://www.jstatsoft.org/v08/i14/xorshift.pdf
///
/// This particular implementation of xorshift has a period of 2^128-1. See the above paper to see
/// how this can be easily extened if you need a longer period.
///
/// 2) Allows fast re-initialisation with a seed, unlike System.Random which accepts a seed at construction
/// time which then executes a relatively expensive initialisation routine. This provides a vast speed improvement
/// if you need to reset the pseudo-random number sequence many times, e.g. if you want to re-generate the same
/// sequence many times. An alternative might be to cache random numbers in an array, but that approach is limited
/// by memory capacity and the fact that you may also want a large number of different sequences cached. Each sequence
/// can each be represented by a single seed value (int) when using FastRandom.
///
/// \note 
/// A further performance improvement can be obtained by declaring local variables as static, thus avoiding
/// re-allocation of variables on each call. However care should be taken if multiple instances of
/// FastRandom are in use or if being used in a multi-threaded environment.
///
/// \author  Jon Rea 
///
///
	class PD_API FastRandom
	{
		// -------------------------------------------
		// Begin 'FastRandom' Singleton Implementation
		// -------------------------------------------
	public:
		static FastRandom* getInstance();
	private:
		static FastRandom* GlobalFastRand;
		// -------------------------------------------
		// End
		// -------------------------------------------

	private:
		const double REAL_UNIT_INT, REAL_UNIT_UINT;
		unsigned int x, y, z, w;
		unsigned int bitBuffer;     ///< Buffer 32 bits in bitBuffer, return 1 at a time
		int bitBufferIdx;           ///< keep track of how many have been returned with bitBufferIdx.
	public:

		FastRandom();
		FastRandom(int seed);
		~FastRandom();
		FastRandom(const FastRandom& _Clone);

#ifndef SWIG
		FastRandom& operator=(const FastRandom &_Clone);
#endif 

	  /// Allows fast re-initialisation with a seed. Unlike srand() which accepts a seed at construction time which then executes a relatively expensive initialisation routine. This provides a vast speed improvement if you need to reset the pseudo-random number sequence many times, e.g. if you want to re-generate the same sequence many times. An alternative might be to cache random numbers in an array, but that approach is limited by memory capacity and the fact that you may also want a large number of different sequences cached. Each sequence can each be represented by a single seed value (int) when using FastRandom.
		void reinitialise(int seed);
		unsigned int nextUInt();
		int next();
		int next(int upperBound);
		int next(int lowerBound, int upperBound);
		double nextDouble();
		bool nextBool();

		/// returns two normally distributed numbers with a given sigma and mean 0
		template < class T >
		void nextNormal(T & nr1, T & nr2, T sigma)
		{
			T zeta1, zeta2;
			T x1, x2;

			zeta1 = (T) nextDouble(); // uniform random number 0..1
			zeta2 = (T) nextDouble();

			double l = sqrt(-2.0 * log(zeta1));
			x1 =  l * cos( (T)(-MathConst::TwoPI) * zeta2 ); // normal number -$..$
			x2 =  l * sin( (T)(-MathConst::TwoPI) * zeta2 ); // normal number -$..$

			nr1 = x1 * sigma;
			nr2 = x2 * sigma;
		}

		/// returns a normally distributed number with a given sigma and mean 0
		template < class T >
		T nextNormal( T sigma )
		{
			return ( sqrt(-2.0 * log( (T) nextDouble() )) * cos(-2.0 * (T) MathConst::PI * (T) nextDouble() ) ) * sigma;
		}

	}; // class PD_API FastRandom
}; // namepace Maths

#endif

