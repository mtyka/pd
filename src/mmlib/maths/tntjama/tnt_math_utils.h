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

#ifndef MATH_UTILS_H
#define MATH_UTILS_H

/* needed for fabs, sqrt() below */
#include <cmath>



namespace TNT
{
	/**
	@returns hypotenuse of real (non-complex) scalars a and b by
	avoiding underflow/overflow
	using (a * sqrt( 1 + (b/a) * (b/a))), rather than
	sqrt(a*a + b*b).
	*/
	template <class Real>
	Real hypot(const Real &a, const Real &b)
	{

		if (a== 0)
			return abs(b);
		else
		{
			Real c = b/a;
			return fabs(a) * sqrt(1 + c*c);
		}
	}
} /* TNT namespace */



#endif
/* MATH_UTILS_H */
