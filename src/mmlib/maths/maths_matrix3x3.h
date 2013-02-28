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

#ifndef __MATHS_MATRIX3X3_H
#define __MATHS_MATRIX3X3_H

namespace Maths
{
	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	///
	/// \author Mike Tyka  
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class PD_API matrix3x3
	{
	public:
		double r[3][3]; // nine values

		matrix3x3()
		{
			setToNull();
		}

		void setTo(const matrix3x3 & mat); ///< copy function
		void setToTranspose(const  matrix3x3 & mat); ///< set to transpose of another matrix
		void setToijk(
			const Maths::dvector & i,
			const Maths::dvector & j,
			const Maths::dvector & k); ///< create a rotation matrix from 3 orthogonal ijk vectors
		void setToStarVector(const Maths::dvector & vec);///< create a star matrix (matrix crossproduct)
		void setTo(
			double r00, double r01, double r02, ///< set values individually
			double r10, double r11, double r12,
			double r20, double r21, double r22);
		void setToIdentity(); ///< set to identity matrix
		void setToNull(); ///< set to Null matrix

		/// create simple rotation matrix around X-axis
		void setToXrot(double angle);

		/// create simple rotation matrix around Y-axis
		void setToYrot(double angle);

		/// create simple rotation matrix around Z-axis
		void setToZrot(double angle);

		/// create rotation matrix around x, y and z axes 
		void setToXYZrot(double xangle, double yangle, double zangle);

		/// create simple rotation matrix around arbitrary
		void setToAxisRot(const Maths::dvector & axis, double angle);

		/// \brief set to a random uniformly distributed rotation matrix
		/// \details
		/// Algorithm used:James Arvo, "Fast Random Rotation Matrices", in
		/// Graphics Gems III, pages 117-120, edited by David Kirk, Academic
		/// Press, New York, 1992.
		void setToRandomRot();

		void add(const matrix3x3 & mat);
		void sub(const matrix3x3 & mat);

		void postmul(const matrix3x3 & mmat); ///< matrix POST multiplication i.e. C = AB : C.setTo(A); C.mul(&B);
		void premul(const matrix3x3 & mmat);  ///< matrix PRE multiplication i.e. C = BA : C.setTo(A); C.mul(&B);

		void mul(double k); ///< scalar multiplication
		void div(double k); ///< scalar division
		void orthonormalize(); ///< othonormalises a rotation matrix
		double normalise(); ///< normalises the matrix so that the internal sum = 1
		void transpose(); ///< sets it's own transpose
		void invert(); ///< inverts matrix (does nothing if det A = 0 )
		double determinant() const; ///< returns determinant
		int diagonaliseSymetric(
			double &lambda1, double &lambda2, double &lambda3,
			Maths::dvector &v1, Maths::dvector &v2, Maths::dvector &v3);
		void info() const; ///< prints values
		void printMathematicaFormat() const;
	};
}

#endif
