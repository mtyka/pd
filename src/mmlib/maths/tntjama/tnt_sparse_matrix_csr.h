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

/*
*
* Template Numerical Toolkit (TNT)
*
* Mathematical and Computational Sciences Division
* National Institute of Technology,
* Gaithersburg, MD USA
*
*
* This software was developed at the National Institute of Standards and
* Technology (NIST) by employees of the Federal Government in the course
* of their official duties. Pursuant to title 17 Section 105 of the
* United States Code, this software is not subject to copyright protection
* and is in the public domain. NIST assumes no responsibility whatsoever for
* its use by other parties, and makes no guarantees, expressed or implied,
* about its quality, reliability, or any other characteristic.
*
*/


#ifndef TNT_SPARSE_MATRIX_CSR_H
#define TNT_SPARSE_MATRIX_CSR_H

#include "tnt_array1d.h"

namespace TNT
{


	/**
	Read-only view of a sparse matrix in compressed-row storage
	format. Neither array elements (nonzeros) nor sparsity
	structure can be modified. If modifications are required,
	create a new view.

	<p>
	Index values begin at 0.

	<p>
	<b>Storage requirements:</b> An (m x n) matrix with
	nz nonzeros requires no more than ((T+I)*nz + M*I)
	bytes, where T is the size of data elements and
	I is the size of integers.


	*/
	template <class T>
	class Sparse_Matrix_CompRow{

	private:
		Array1D<T> val_; // data values (nz_ elements)
		Array1D<int> rowptr_; // row_ptr (dim_[0]+1 elements)
		Array1D<int> colind_; // col_ind (nz_ elements)

		int dim1_; // number of rows
		int dim2_; // number of cols

	public:

		Sparse_Matrix_CompRow(const Sparse_Matrix_CompRow &S);
		Sparse_Matrix_CompRow(int M, int N, int nz, const T *val,
			const int *r, const int *c);



		inline const T& val(int i) const { return val_[i]; }
		inline const int& row_ptr(int i) const { return rowptr_[i]; }
		inline const int& col_ind(int i)const { return colind_[i];}

		inline int dim1() const {return dim1_;}
		inline int dim2() const {return dim2_;}
		int NumNonzeros()const {return val_.dim1();}


		Sparse_Matrix_CompRow& operator=(
			const Sparse_Matrix_CompRow &R);



	};

	/**
	Construct a read-only view of existing sparse matrix in
	compressed-row storage format.

	@param M the number of rows of sparse matrix
	@param N the number of columns of sparse matrix
	@param nz the number of nonzeros
	@param val a contiguous list of nonzero values
	@param r row-pointers: r[i] denotes the begining position of row i
	(i.e. the ith row begins at val[row[i]]).
	@param c column-indices: c[i] denotes the column location of val[i]
	*/
	template <class T>
	Sparse_Matrix_CompRow<T>::Sparse_Matrix_CompRow(int M, int N, int nz,
		const T *val, const int *r, const int *c): val_(nz,val),
		rowptr_(M, r), colind_(nz, c), dim1_(M), dim2_(N) {}


}
// namespace TNT

#endif
