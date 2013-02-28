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

#ifndef __MATHS_VECTOR_H
#define __MATHS_VECTOR_H

#include "maths.h"

namespace Maths
{





//-------------------------------------------------
//
/// \brief Tempalte 3D vector class (template so it can be double or float) 
///
/// \details 
/// Function: provides all function needed to manipulate vectors
/// Templatisation allows double and single precision doubles to be used
/// from the same code and interconverted
///
/// \author Mike Tyka  
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	template < class T >
	class Tvector
	{
	public:
		T x, y, z; ///< Tvector coordinates (x,y,z) of the Tvector
		Tvector(void);
		Tvector(T, T, T); ///< basic constructor
		Tvector(const Tvector < T > &v1); ///< universal copy constructor
		template < class U >
		Tvector(const Tvector < U > &v1); ///< universal copy constructor

		void info() const; ///< prints the Tvector

		void setTo(const Tvector < T > &v2)
		{
			x = v2.x;
			y = v2.y;
			z = v2.z;
		}
#ifndef SWIG
		Tvector < T > &operator=(const Tvector < T > &v2)
		{
			x = v2.x;
			y = v2.y;
			z = v2.z;
			return *this;
		}
#endif

		void setTo(T _x, T _y, T _z)
		{
			x = _x;
			y = _y;
			z = _z;
		}

		void setToCentrePointBetween( const Maths::dvector& v1, const Maths::dvector& v2 )
		{
			x = (v1.x + v2.x) / 2.0;
			y =	(v1.y + v2.y) / 2.0;
			z = (v1.z + v2.z) / 2.0;
		}

		void setToRandomUnit()
		{
			double t, r;

			z = frand() * 2.0 - 1.0;
			t = frand() * MathConst::TwoPI;
			r = sqrt(1 - sqr(z));
			x = r * cos(t);
			y = r * sin(t);
		}

		bool isReal() const 
		{
			if(!isNumber(x)) return false;
			if(!isNumber(y)) return false;
			if(!isNumber(z)) return false;
			return true;
		}

		inline void zero()
		{
			x = 0;
			y = 0;
			z = 0;
		}

		inline void inv()
		{
			x = -x;
			y = -y;
			z = -z;
		}

		inline T innerdot() const 
		{
			return (T) (x * x + y * y + z * z);
		}

		inline T mag() const 
		{
			return (T) sqrt(x * x + y * y + z * z);
		}

		inline void add(T _x, T _y, T _z)
		{
			x += _x;
			y += _y;
			z += _z;
		}

		inline void sub(T _x, T _y, T _z)
		{
			x -= _x;
			y -= _y;
			z -= _z;
		}

		inline void add(const Tvector < T > &v2)
		{
			x += v2.x;
			y += v2.y;
			z += v2.z;
		}

		inline void sub(const Tvector < T > &v2)
		{
			x -= v2.x;
			y -= v2.y;
			z -= v2.z;
		}

		inline void mul(T f)
		{
			x *= f;
			y *= f;
			z *= f;
		}

		inline void div(T f)
		{
			x /= f;
			y /= f;
			z /= f;
		}

		inline T scalarProduct(const Tvector < T > &v2) const; // calculates scalar product

		inline T angleWith(const Tvector < T > &v2) const; // outputs an acute angle (rad)

		inline void vectorProduct(const Tvector < T > *v2, Tvector < T > *product); 

		inline void diff(const Tvector < T > &v1, const Tvector < T > &v2)
		{
			x = v1.x - v2.x;
			y = v1.y - v2.y;
			z = v1.z - v2.z;
		};

		inline void sum(const Tvector < T > &v1, const Tvector < T > &v2)
		{
			x = v1.x + v2.x;
			y = v1.y + v2.y;
			z = v1.z + v2.z;
		};

		inline void crossProduct(const Tvector < T > &v1, const Tvector < T > &v2);

		inline void crossProduct(const Tvector < T > &v1, const Tvector < T > &v2, T scalar);

		inline void mulmat(const matrix3x3 & mmat);

		inline T dist(const Tvector < T > &v2) const;

		inline T sqrdist(const Tvector < T > &v2) const;

		inline void unify()
		{
			T magnitude = mag();
			if(magnitude != 0) 
			{
				x /= magnitude;
				y /= magnitude;
				z /= magnitude;
			}
		}

		inline void rotateX(T angle);
		inline void rotateY(T angle);
		inline void rotateZ(T angle);

		inline void rotateAxis(
			T angle, 
			const Tvector < T > &offset, 
			const Tvector < T > &axis
			);

		inline void rotateYat(
			T angle, 
			const Tvector < T > &rotCentre
			);

		inline void interpolate_2(
			const Tvector < T > &v1, T p1, 
			const Tvector < T > &v2, T p2
			);

		inline void interpolate_3(
			const Tvector < T > &v1, T p1,
			const Tvector < T > &v2, T p2, 
			const Tvector < T > &v3, T p3
			);

		inline void interpolate_4(
			const Tvector < T > &v1, T p1,
			const Tvector < T > &v2, T p2,
			const Tvector < T > &v3, T p3, 
			const Tvector < T > &v4, T p4
			);
	};

	// Tvector non-member classes

	template < class T >
	inline T dotProduct(const Tvector < T > &v1, const Tvector < T > &v2){
		return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
	}

	template < class T >
	inline T sqrdist(const Tvector < T > &v1,const Tvector < T > &v2){
		T g = (T) (v1.x - v2.x);
		T h = (T) (v1.y - v2.y);
		T j = (T) (v1.z - v2.z);
		return (g * g + h * h + j * j);
	}

	template < class T >
	inline T dist(const Tvector < T > &v1,const Tvector < T > &v2){
		return sqrt(sqrdist(v1,v2));
	}
#ifndef SWIG
	template < class T >
	inline bool operator== (const Tvector < T > &v1,const Tvector < T > &v2){
		if((v1.x==v2.x)&&(v1.y==v2.y)&&(v1.z==v2.z)) return true;
		return false;
	}
#endif

	template < class T >
	inline void rotateVectors(const Tvector < T > *zero, matrix3x3 * rmat, const Tvector < T > *inarray,
		Tvector < T > *outarray, int n){
			for(int i = 0; i < n; i++) {
				outarray[i].x = inarray[i].x - zero->x;
				outarray[i].y = inarray[i].y - zero->y;
				outarray[i].z = inarray[i].z - zero->z;
				outarray[i].mulmat(rmat);
				outarray[i].add(zero);
			}
	}

	// ----------------------
	// Tvector class member functions
	// ----------------------

	template < class T >
	Tvector < T >::Tvector(void)
	{
	}

	template < class T >
	Tvector < T >::Tvector(T px, T py, T pz)
	{
		x = px;
		y = py;
		z = pz;
	}

	template < class T >
	template < class U >
	inline Tvector < T >::Tvector(const Tvector < U > &v1)
	{
		x = (T) v1.x;
		y = (T) v1.y;
		z = (T) v1.z;
	}

	template < class T >
	inline Tvector < T >::Tvector(const Tvector < T > &v1)
	{
		x = v1.x;
		y = v1.y;
		z = v1.z;
	}

	template < class T >
	inline void Tvector < T >::info() const
	{
		printf("x:% 7.3lf y:% 7.3lf z:% 7.3lf mag:%7.3lf\n", (double) x, (double) y, (double) z, mag());
	};

	template < class T >
	inline T Tvector < T >::scalarProduct(const Tvector < T > &v2) const 
	{
		return (x * (T) v2.x + y * (T) v2.y + z * (T) v2.z);
	}

	template < class T >
	inline void Tvector < T >::crossProduct(const Tvector < T > &v1, const Tvector < T > &v2)
	{
		x = (v1.y * v2.z - v1.z * v2.y);
		y = (v1.z * v2.x - v1.x * v2.z);
		z = (v1.x * v2.y - v1.y * v2.x);
	}

	template < class T >
	inline void Tvector < T >::crossProduct(const Tvector < T > &v1, const Tvector < T > &v2, T scalar)
	{
		x = (v1.y * v2.z - v1.z * v2.y) * scalar;
		y = (v1.z * v2.x - v1.x * v2.z) * scalar;
		z = (v1.x * v2.y - v1.y * v2.x) * scalar;
	}

	template < class T >
	inline void Tvector < T >::vectorProduct(const Tvector < T > *v2, Tvector < T > *product)
	{
		product->x = (y * v2->z - z * v2->y); // finding the Tvector<U> product
		product->y = (z * v2->x - x * v2->z); // the quotient makes it a unit Tvector<U>
		product->z = (x * v2->y - y * v2->x); // so might be taken away for economy reasons
	}

	template < class T >
	inline T Tvector < T >::angleWith(const Tvector < T > &v2) const 
	{
		T scalar = scalarProduct(v2) / (mag() * v2.mag()); //simple scalar product calculation
		if (scalar > 1)
			return 0;
		if (scalar < -1)
			return (T) MathConst::PI;
		return acos(scalar);
	}

	template < class T >
	inline T Tvector < T >::dist(const Tvector < T > &v2)const 
	{

		T g = (T) (v2.x - x);
		T h = (T) (v2.y - y);
		T j = (T) (v2.z - z);
		return sqrt(g * g + h * h + j * j);
	}

	template < class T >
	inline T Tvector < T >::sqrdist(const Tvector < T > &v2)const 
	{
		T g = (T) (v2.x - x);
		T h = (T) (v2.y - y);
		T j = (T) (v2.z - z);
		return (g * g + h * h + j * j);
	}

	template < class T >
	inline void Tvector < T >::mulmat(const matrix3x3 & mmat)
	{
		T nx, ny, nz;
		nx = (T) mmat.r[0][0] * x + (T) mmat.r[0][1] * y + (T) mmat.r[0][2] * z;
		ny = (T) mmat.r[1][0] * x + (T) mmat.r[1][1] * y + (T) mmat.r[1][2] * z;
		nz = (T) mmat.r[2][0] * x + (T) mmat.r[2][1] * y + (T) mmat.r[2][2] * z;
		x = nx;
		y = ny;
		z = nz;
	}

	template < class T >
	inline void Tvector < T >::rotateY(T angle)
	{
		T cangle = cos(angle);
		T sangle = sin(angle);
		T tx = x * cangle + z * sangle;
		T tz = -x * sangle + z * cangle;
		x = tx;
		z = tz;
	}

	template < class T >
	inline void Tvector < T >::rotateX(T angle)
	{
		T cangle = cos(angle);
		T sangle = sin(angle);
		T ty = y * cangle + z * sangle;
		T tz = -y * sangle + z * cangle;
		y = ty;
		z = tz;
	}

	template < class T >
	inline void Tvector < T >::rotateZ(T angle)
	{
		T cangle = cos(angle);
		T sangle = sin(angle);
		T tx = x * cangle + y * sangle;
		T ty = -x * sangle + y * cangle;
		x = tx;
		y = ty;
	}


	template < class T >
	inline void Tvector < T >::rotateAxis(
		T angle, 
		const Tvector < T > &offset, 
		const Tvector < T > &aaxis
		)
	{
		sub(offset); //translate so offset becomes 0

		Tvector axis(aaxis); //make a copy of the axis

		Tvector xaxis(1, 0, 0); //create an xaxis

		Tvector xzplane(axis);
		xzplane.y = 0; //obtain angle in the xzplane
		T xzangle = xaxis.angleWith(xzplane);
		if(xzplane.z < 0)
			xzangle *= -1;

		rotateY(xzangle); //rotate the point and the axis
		axis.rotateY(xzangle);

		Tvector xyplane(axis);
		xyplane.z = 0; //repeat with xyplane
		T xyangle = xaxis.angleWith(xyplane);
		if(xyplane.y < 0)
			xyangle *= -1;

		rotateZ(xyangle);

		rotateX(angle); //now rotate round axis (X)
		rotateZ(-xyangle); //rotate back in backwards direction
		rotateY(-xzangle); //

		add(offset); //translate back
	}

	template < class T >
	inline void Tvector < T >::rotateYat(
		T angle, 
		const Tvector < T > &rotCentre)
	{
		sub(rotCentre);

		T tx = x * cos(angle) + z * sin(angle);
		T tz = -x * sin(angle) + z * cos(angle);

		x = tx;
		z = tz;

		add(rotCentre);
	}


	template < class T >
	inline void Tvector < T >::interpolate_2(
		const Tvector < T > &v1, T p1, 
		const Tvector < T > &v2, T p2)
	{
		x = v1.x * p1 + v2.x * p2;
		y = v1.y * p1 + v2.y * p2;
		z = v1.z * p1 + v2.z * p2;
	}

	template < class T >
	inline void Tvector < T >::interpolate_3(
		const Tvector < T > &v1, T p1,
		const Tvector < T > &v2, T p2, 
		const Tvector < T > &v3, T p3)
	{
		x = v1.x * p1 + v2.x * p2 + v3.x * p3;
		y = v1.y * p1 + v2.y * p2 + v3.y * p3;
		z = v1.z * p1 + v2.z * p2 + v3.z * p3;
	}

	template < class T >
	inline void Tvector < T >::interpolate_4(
		const Tvector < T > &v1, T p1,
		const Tvector < T > &v2, T p2,
		const Tvector < T > &v3, T p3, 
		const Tvector < T > &v4, T p4 )
	{
			x = v1.x * p1 + v2.x * p2 + v3.x * p3 + v4.x * p4;
			y = v1.y * p1 + v2.y * p2 + v3.y * p3 + v4.y * p4;
			z = v1.z * p1 + v2.z * p2 + v3.z * p3 + v4.z * p4;
	}



}

#endif

