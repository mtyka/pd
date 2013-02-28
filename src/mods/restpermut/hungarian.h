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

#ifndef __HUNGARIAN_H
#define __HUNGARIAN_H

#include <valarray>

// Description of this solution to the linear ssignment problem
//
// http://www.public.iastate.edu/~ddoty/HungarianAlgorithm.html
//


// Hungarian Algorithm for Linear assignement - the original 
// algorithm due to Kuhn and Munkres. Its pretty slow ( O(n^4) )
// 
// Kuhn, H. W.: The Hungarian method for the assignment problem. 
// Naval Research Logistics Quarterly 2, 83-97 (1955).
int LinearAssignmentHungarian(
	std::valarray<int> &costm, 
	std::vector<int> &row2col, 
	std::vector<int> &col2row
);


/// A much faster algorithm for linear assignement: 
///  Reference: Jonker, R. and Volgenant, A., 1987. A shortest augmenting
///  path algorithm for dense and sparse linear assignment problems. 
///  Computing 38, pp. 325–340
int LinearAssignmentJVC(
	std::valarray<int> &costm,
	std::vector<int> &row2col,
	std::vector<int> &col2row
);

#endif



