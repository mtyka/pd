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

#ifndef __MATHS_H
#define __MATHS_H

#include <math.h> // Reqired as template bodies must be in this file, and functions like pow() are used
#include "maths/maths.fwd.h"


namespace Maths
{
	/// Double Precision Mathematical Constants
	class PD_API MathConst
	{
	private:
		MathConst(){}
	public:
		static const double PI;
		static const double TwoPI;
		static const double FourPI;
		static const double PIOver180;
		static const double OneEightyOverPI;
	};

	// returns true if 'a' and 'b' are equal to '_DecPlaces' decimal places...
	template <class T>
	bool SigFigEquality( T _a, T _b, int _DecPlaces )
	{
		return fabs(_a-_b) < (pow(10.0,(double)(-_DecPlaces)));
	}

	template < class T >
	inline T min(T a, T b)
	{
		return ((a) < (b) ? (a) : (b));
	}

	template < class T >
	inline T min(T a, T b, T c)
	{
		return ((a) < (b) ? min( a, c) : min( b, c));
	}

	template < class T >
	inline T minmax(T a, T b, T c)
	{
		if(b>c)return c;
		if(b<a)return a;
		return b;
	}

	template < class T >
	inline T max(T a, T b)
	{
		return ((a) > (b) ? (a) : (b));
	}

	template <class T>
	inline T sqr(T x)
	{
		return x*x;
	}

	template <class T>
	inline T cube(T x)
	{
		return x*x*x;
	}

	inline double RadToDeg(double rad)
	{
		return rad * MathConst::OneEightyOverPI;
	}

	inline double DegToRad(double deg)
	{
		return deg * MathConst::PIOver180;
	}

	template <class T>
	inline T SqrDistanceBetween( T x1, T y1, T z1, T x2, T y2, T z2 )
	{
		return sqr((x1)-(x2)) + sqr((y1)-(y2)) + sqr((z1)-(z2));
	}

	template <class T>
	inline T DistanceBetween( T x1, T y1, T z1, T x2, T y2, T z2 )
	{
		return sqrt(sqr((x1)-(x2)) + sqr((y1)-(y2)) + sqr((z1)-(z2)));
	}

	// used to 'simulate' 2D arrays using linear arrays
#define sqrmat(y,x,a) ((y)*(a) + (x))

#define periodic(x,limit) (0<(x)?((x)%(limit)): ((((limit)*10000)+x)%(limit)) ) // keeps argument in periodic limit

	/// returns true is its a valid double, false if its a NaN, +-INF,QNAN etc..
	bool isNumber(float number); 

	/// returns true is its a valid double, false if its a NaN, +-INF,QNAN etc..
	bool isNumber(double number);

	/// returns a float NAN
	float fNAN();

	/// returns a double NAN
	double dNAN();

	/// returns the exponential for of a number as in 
	/// prefix E exponent
	template < class T >
	void getExpForm(T number, T & prefix, T & exponent){
		if(number == 0) {
			prefix = 0;
			exponent = 0;
		} else {
			exponent = floor(log10(fabs(number)));
			prefix = number / pow((T) 10, exponent);
		}
	}

	/// returns a binary random decision
	PD_API bool brand();

	/// returns a random number evenly distributed between 0.0 .. 1.0 (not including!)
	PD_API double frand();

	/// returns a random number evenly distributed between 0.0 .. 1.0 (including!)
	PD_API double frand01(); 

	/// returns one normally distributed number with a given sigma and mean 0
	template < class T >
	inline T nrand( T sigma )
	{
		return ( sqrt(-2.0 * log( (T) frand() )) * cos((T)(-MathConst::TwoPI) * (T) frand() ) ) * sigma;
	}

	/// returns two normally distributed numbers with a given sigma and mean 0
	template < class T >
	inline void nrand(T & nr1, T & nr2, T sigma = 1.0){

		T zeta1, zeta2;
		T x1, x2;

		zeta1 = (T) frand(); // uniform random number 0..1
		zeta2 = (T) frand();

		double l = sqrt(-2.0 * log(zeta1));
		x1 =  l * cos( (T)(-MathConst::TwoPI) * zeta2 ); // normal number -$..$
		x2 =  l * sin( (T)(-MathConst::TwoPI) * zeta2 ); // normal number -$..$

		nr1 = x1 * sigma;
		nr2 = x2 * sigma;
	}

	/// chooses one of N according to probability provided in a flaot array the probabilities need not to be normalised
	template < class T >
	int chooseByProbability(T * problist, int N){
		T sum = 0.0;
		T randomnr;
		int i;

		for(i = 0; i < N; i++)
			sum += problist[i];
		randomnr = (T) frand01() * sum;
		for(i = 0; i < N; i++) {
			randomnr -= problist[i];
			if(randomnr < 0.0)
				return i;
		}
		return N - 1;
	}

	/// returns 'i!'
	PD_API unsigned factorial(unsigned i);

	/// returns 'ln(i!)'
	PD_API double log_factorial(unsigned i);

	/// returns 'ln(i!)'
	PD_API double log_factorial_sterling(unsigned i);

	/// returns the surface area of a sphere with radius R
	PD_API double sphereSurfaceArea( double r ); 

	///calculates the overlap of two spheres of sizes ri, rj at a distance dij
	PD_API double sphereSurfaceOverlap(double ri, double rj, double dij); 

	///calculates the derivative overlap of two spheres of sizes ri, rj at a distance dij
	PD_API double sphereSurfaceOverlapDeriv(double ri, double rj, double dij); 




	/// returns the value of x^2 + a1*x + a2 = 0
	double calcQuadraticPolynomial(
		double a1, 
		double a2, 
		double x
		);

	/// solves the quadratic equation 	x^2 + a1*x + a2 = 0 returning the real
	/// roots in z1 and z2. The number of real roots is returned directly.
	int solveQuadraticPolynomial(
		double a1, 
		double a2, 
		double &z1, 
		double &z2
		);

	/// returns the value of x^3 + a1*x^2 + a2*x + a3 = 0;
	double calcCubicPolynomial(
		double a1, 
		double a2, 
		double a3, 
		double x
		);

	/// solves the quadratic equation x^3 + a1*x^2 + a2*x + a3 = 0 returning any real
	/// roots in z1,z2 and z3. The number of real roots is returned directly.
	int solveCubicPolynomial(
		double a1, 
		double a2, 
		double a3, 
		double &z1, 
		double &z2, 
		double &z3
		);

}

#include "maths/maths_matrix3x3.h"
#include "maths/maths_vector.h"

#ifdef SWIG
%template(dvector)	Maths::Tvector<double>;
%template(fvector)	Maths::Tvector<float>;
#endif

namespace Maths
{
	template < class T >
	T calcTorsionAngle(const Tvector < T > &i, const Tvector < T > &a, 
		const Tvector < T > &b, const Tvector < T > &j)
	{
		Maths::dvector r12, r23, r34;
		Maths::dvector A, B, C;
		T rA, rB, rC;
		T sin_phi, cos_phi, phi;
		r12.diff(i, a);
		r23.diff(a, b);
		r34.diff(b, j);
		A.crossProduct(r12, r23);
		B.crossProduct(r23, r34);
		C.crossProduct(r23, A);
		rA = A.mag();
		rB = B.mag();
		rC = C.mag();
		cos_phi = dotProduct(A, B) / (rA * rB);
		sin_phi = dotProduct(C, B) / (rC * rB);
		phi = -atan2(sin_phi, cos_phi);
		return phi;
	}

	// this function tries to superimpose three vectors onto three target vectors the best it can.
	// giving perfect alignment of the o->a Maths::dvector and as good as good an alignment on the a->b Maths::dvector
	// as possible.
	template < class T, class U >
	void superimpose(
		const Tvector < T >& targeto, 
		const Tvector < T >& targeta, 
		const Tvector < T >& targetb,
		const Tvector < U >& veco, 
		const Tvector < U >& veca, 
		const Tvector < U >& vecb, 
		matrix3x3 & rmat)
	{
		Tvector < T > toa, oa, tab, ab; // create all the nessessary vectors
		toa.diff(targeta, targeto);
		tab.diff(targetb, targeta);
		oa.diff(veca, veco);
		ab.diff(vecb, veca);

		if( (toa==oa) && (tab == ab) )
		{
			rmat.setToIdentity();
			return;
		}

		Tvector < T > toa_x_oa;
		T angle = toa.angleWith(oa);
		if(angle != 0.0)
		{
			toa_x_oa.crossProduct(toa, oa);
			rmat.setToAxisRot(toa_x_oa, -angle); // construct rotation to superimpose oa and toa
		}
		else
		{
			rmat.setToIdentity();
		}

		ab.mulmat(rmat);

		Tvector < T > toa_x_tab, toa_x_ab;
		Tvector < T > doublecross;

		toa_x_tab.crossProduct(toa, tab);
		toa_x_ab.crossProduct(toa, ab);

		doublecross.crossProduct(toa_x_tab, toa_x_ab);
		angle = toa_x_tab.angleWith(toa_x_ab);

		if(angle != 0.0)
		{
			if(dotProduct(doublecross, toa) > 0) angle *= -1.0;
			matrix3x3 rmat2;
			rmat2.setToAxisRot(toa, angle); // construct rotation to superimpose toa_x_tab and toa_x_ab
			rmat.premul(rmat2); // add the second r
		}
	}

	// this is like the previous function but only aligns on two atoms.
	// this of course leaves an ambiguity about the rotation around those two
	// atoms - this is left random and will depend on the previous orientations
	template < class T, class U >
	void superimpose(
		Tvector < T > const &targeto, 
		Tvector < T > const &targeta,
		Tvector < U > const &veco,    
		Tvector < U > const &veca, 
		matrix3x3 & rmat)
	{
		Tvector < T > toa, oa, tab, ab; // create all the nessessary vectors
		toa.diff(targeta, targeto);
		oa.diff(veca, veco);

		Tvector < T > toa_x_oa;
		toa_x_oa.crossProduct(toa, oa);
		T angle = toa.angleWith(oa);

		rmat.setToAxisRot(toa_x_oa, -angle); // construct rotation to superimpose oa and toa
	}


	template < class T >
	int leastSquaresFit(T * x, T * y, int N, T & a, T & b, T & R){
		if(N <= 1)
			return -1;

		T sumx = 0;
		T sumy = 0;
		T meanx = 0;
		T meany = 0;
		T sumxx = 0, sumyy = 0, sumxy = 0;

		int i;
		for(i = 0; i < N; i++) {
			sumx += x[i];
			sumy += y[i];
		}
		meanx = sumx / (T) N;
		meany = sumy / (T) N;
		for(i = 0; i < N; i++) {
			sumxx += sqr(x[i] - meanx);
			sumyy += sqr(y[i] - meany);
			sumxy += (x[i] - meanx) * (y[i] - meany);
		}

		b = sumxy / sumxx;
		a = meany - b * meanx;

		R = sqr(sumxy) / (sumxx * sumyy);

		return 0;
	}


	void leastSquaresFit(const std::vector<double> &x, const std::vector<double> &y, double & a, double & b, double & R);


	/// A simple graph of lnear segments which allows simple interpolation. Some simple forcefields
	/// can use this.
	class PD_API PolygonGraph
	{
	public:
		PolygonGraph(int n, double *x, double *y);
		~PolygonGraph();
		double InterpolateFromX(double x);
		double InterpolateFromX(double x, double *dydx);
	protected:
		Maths::dvector * point;
		int npoints;
	};
}

#endif
