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
#include "exception.h"

namespace Maths
{
	const double MathConst::PI = 3.14159265358979323846;
	const double MathConst::TwoPI = 2.0 * PI;
	const double MathConst::FourPI = 4.0 * PI;
	const double MathConst::PIOver180 = PI / 180.0;
	const double MathConst::OneEightyOverPI = 180.0 / PI;

	void leastSquaresFit(const std::vector<double> &x, const std::vector<double> &y, double & a, double & b, double & R)
	{
		if(x.size() != y.size()) throw(ArgumentException("The two arrays given to leastSquaresFit must be of equal length"));
		if(x.size() <= 1)          throw(ArgumentException("At least 2 data points must be given for leastSquaresFit to work"));

		size_t N = x.size();

		double sumx = 0.0;
		double sumy = 0.0;
		double sumxx = 0.0; 
		double sumyy = 0.0; 
		double sumxy = 0.0;

		for(size_t i = 0; i < N; i++) 
		{
			sumx += x[i];
			sumy += y[i];
		}

		double meanx = sumx / (double) N;
		double meany = sumy / (double) N;

		for(size_t i = 0; i < N; i++) 
		{
			sumxx += sqr(x[i] - meanx);
			sumyy += sqr(y[i] - meany);
			sumxy += (x[i] - meanx) * (y[i] - meany);
		}

		b = sumxy / sumxx;
		a = meany - b * meanx;

		R = sqr(sumxy) / (sumxx * sumyy);
	}

	PolygonGraph::PolygonGraph(int n, double *x, double *y)
	{
		npoints = n;
		point = new dvector[npoints];
		for(int i = 0; i < npoints; i++) 
		{
			point[i].setTo(x[i], y[i], 0);
		}
	}

	PolygonGraph::~PolygonGraph()
	{
		if(point != NULL)
			delete[]point;
	}

	double PolygonGraph::InterpolateFromX(double x)
	{
		double dydx;
		return InterpolateFromX(x, &dydx);
	}

	double PolygonGraph::InterpolateFromX(double x, double *dydx)
	{
		int i;
		for(i = 0; i < npoints; i++) {
			if(x < point[i].x)
				break;
		}
		i--;
		if(i >= (npoints - 1)) {
			*dydx = 0;
			return 0;
		}

		double dx, dy;
		dx = point[i + 1].x - point[i].x;
		dy = point[i + 1].y - point[i].y;

		*dydx = dy / dx;

		return point[i].y + dy * (x - point[i].x) / dx;
	}


	// -------------------------------------------------------------------
	// Standalone functions
	// -------------------------------------------------------------------

	// The surface area of a sphere
	double sphereSurfaceArea( double r )// returns the surface area of a sphere with radius R
	{
		return Maths::MathConst::FourPI * r * r;
	}

	// Sphere surface overlap
	double sphereSurfaceOverlap(double ri, double rj, double dij)//calculates the overlap of two spheres of sizes ri, rj at a distance dij
	{
		return (Maths::MathConst::TwoPI*(ri) * ((ri) - (dij)/2.0 - ((ri)*(ri) - (rj)*(rj))/(2.0*dij)));
	}

	//calculates the derivative overlap of two spheres of sizes ri, rj at a distance dij
	double sphereSurfaceOverlapDeriv(
		double ri, 
		double rj,
		double dij
		)
	{
		return ( Maths::MathConst::PI*(ri)*((((ri)*(ri) - (rj)*(rj))/((dij)*(dij))) - 1));
	}

	double calcQuadraticPolynomial(
		double a1, 
		double a2, 
		double x
		)
	{
		return (sqr(x) + a1 * x + a2);
	}

	int solveQuadraticPolynomial(
		double a1, 
		double a2, 
		double &z1, 
		double &z2
		)
	{
		double D = sqr(a1) - 4 * a2;
		if(D < 0)
			return 0;
		D = sqrt(D);
		z1 = (-a1 - D) / 2.0;
		z2 = (-a1 + D) / 2.0;

		return 2;
	}

	double calcCubicPolynomial(
		double a1, 
		double a2, 
		double a3, 
		double x
		)
	{
		return (cube(x) + a1 * sqr(x) + a2 * x + a3);
	}

	int solveCubicPolynomial(
		double a1, 
		double a2, 
		double a3, 
		double &z1, 
		double &z2, 
		double &z3
		)
	{
		double Q = (sqr(a1) - 3.0 * a2) / 9.0;
		double R = (2 * cube(a1) - 9.0 * a1 * a2 + 27.0 * a3) / 54.0;
		double Qcubed = Q * Q * Q;
		double D = Qcubed - R * R;

		// Three double roots
		if(D >= 0) 
		{
			double theta = acos(R / (sqrt(Q) * Q));
			double sqrtQ = sqrt(Q);
			z1 = -2.0 * sqrtQ * cos(theta / 3.0) - a1 / 3.0;
			z2 = -2.0 * sqrtQ * cos((theta + MathConst::TwoPI) / 3.0) - a1 / 3.0;
			z3 = -2.0 * sqrtQ * cos((theta + MathConst::FourPI) / 3.0) - a1 / 3.0;
			return 3;
		} 
		else // one root
		{			
			double e = pow((double) (sqrt(-D) + fabs(R)), (double) (1.0 / 3.0));
			if(R > 0)
				e *= -1.0;
			z1 = (e + Q / e) - a1 / 3.0;
			return 1;
		}
	}

	unsigned factorial(unsigned i)
	{
		unsigned result = 1;
		unsigned a;
		for(a = 2; a <= i; a++) 
		{
			result *= a;
		}
		return result;
	}

	double log_factorial(unsigned i)
	{
		double result = 0;
		unsigned a;
		for(a = 2; a <= i; a++) 
		{
			result += log((double)a);
		}
		return result;
	}

	double log_factorial_sterling(unsigned i)
	{
		return double(i)*log(double(i)) - double(i);
	}

	bool is_f_finite(double number)
	{
		if((1.0 == number) && (1.0 != number))
			return false;
		if((!(number >= 0.0)) && (!(number <= 0.0)))
			return false;

		return true;
	}

	bool isNumber(float number)
	{
		return (number==number);
	}

	bool isNumber(double number)
	{
		return (number==number);
	}

	float fNAN()
	{
		int iNAN = 0x7800000; // 0 1111 111 0000 0000 0000 0000 0000 0000
		float *pfNAN = (float *) (&iNAN);
		return *pfNAN;
	}

	double dNAN()
	{
		int iNAN = 0x7800000; // 0 1111 111 0000 0000 0000 0000 0000 0000
		double *pfNAN = (double *) (&iNAN);
		return *pfNAN;
	}

	bool brand()
	{
		return ( 0 ==
			rand() % 2 // either 0 or 1 with equal probability ...
			);
	}

	double frand() // returns a random number evenly distributed between 0.0 .. 1.0 (not including!)
	{
		double zeta1;
		zeta1 = (double) (1.0 + (double) rand()) / ((double) (RAND_MAX) + 1.0); //uniform random number 0..1
		return zeta1;
	}

	double frand01() // returns a random number evenly distributed between 0.0 .. 1.0 (including!)
	{ 
		double zeta1;
		zeta1 = (double) ((double) rand()) / ((double) (RAND_MAX)); //uniform random number 0..1
		return zeta1;
	}

	// very plain straight forward distance matrix calculation
	int calcDistanceMatrix(double *dmatrix, dvector * atom, int atoms)
	{
		for(int i = 0; i < atoms; i++) 
		{
			for(int j = i + 1; j < atoms; j++) 
			{
				double Dist_ij = sqrt(sqr((atom[i].x - atom[j].x)) + sqr((atom[i].y - atom[j].y)) + sqr((atom[i].z - atom[j].z)));
				dmatrix[sqrmat(i, j, atoms)] = Dist_ij;
				dmatrix[sqrmat(j, i, atoms)] = Dist_ij;
			}
			dmatrix[sqrmat(i, i, atoms)] = 0;
		}
		return 0;
	}

	// very plain straight forward distance matrix calculation
	int calcDistanceMatrix(double *dmatrix, dvector * atom, int atoms, double Cutoff)
	{
		const double sqrcutoff = sqr(Cutoff);
		for(int i = 0; i < atoms; i++) 
		{
			for(int j = i + 1; j < atoms; j++) 
			{
				double sqrdistij = sqr((atom[i].x - atom[j].x)) + sqr((atom[i].y - atom[j].y)) + sqr((atom[i].z - atom[j].z));
				if(sqrdistij < sqrcutoff) 
				{
					double Dist_ij = sqrt(sqrdistij);
					dmatrix[sqrmat(i, j, atoms)] = Dist_ij;
					dmatrix[sqrmat(j, i, atoms)] = Dist_ij;
				} 
				else 
				{
					dmatrix[sqrmat(i, j, atoms)] = Cutoff; // these are invalid
					dmatrix[sqrmat(j, i, atoms)] = Cutoff; // these are invalid
				}
			}
			dmatrix[sqrmat(i, i, atoms)] = 0;
		}
		return 0;
	}

	void createVertexSphere(dvector * vertex, double radius, int number)
	{
		int i;
		double x, y, z, t, r;

		for(i = 0; i < number; i++) 
		{
			z = ((double) (rand() % 100000) - 50000) / 50000;
			t = ((double) (rand() % 100000) / 50000) * Maths::MathConst::TwoPI;

			r = sqrt(1 - sqr(z));
			x = r * cos(t);
			y = r * sin(t);

			vertex[i].setTo(x * radius, y * radius, z * radius);
		}

	}

	void createEvenVertexSphere(dvector * vertex, double radius, int number)
	{
		double x, y, z, t, r;

		int nvertex = number * 10;

		dvector *nvert = new dvector[nvertex]; // reserve space for 10 times more
		bool *validate = new bool[nvertex];

		for(int i = 0; i < nvertex; i++) 
		{
			z = ((double) (rand() % 100000) - 50000) / 50000;
			t = ((double) (rand() % 100000) / 50000) * Maths::MathConst::TwoPI;

			r = sqrt(1 - sqr(z));
			x = r * cos(t);
			y = r * sin(t);

			nvert[i].setTo(x * radius, y * radius, z * radius);
			validate[i] = true;
		}

		double ilimit = 2 * r / sqrt((double) number);
		double limit = ilimit * 0.9;

		int done = 1, j, k;

		printf("ilimit is %lf\n", ilimit);

		for(int i = 0; i < 100; i++) {

			for(j = 0; j < (nvertex); j++) 
			{
				if(validate[j] == false )
					continue;
				for(k = 0; k < (nvertex); k++) 
				{
					if(validate[k] == false )
						continue;
					if(k == j)
						continue;
					if(nvert[j].dist(nvert[k]) < limit) 
					{
						validate[j] = false; // mark vertex invalid
						nvertex--; // decrease number of vertexes left by 1
						break;
					}
					if(nvertex <= number) 
					{
						done = 0;
						break;
					}
				}
				if(nvertex <= number) 
				{
					done = 0;
					break;
				}
			}
			printf("REMARK limit %lf nvertex %d \n", limit, nvertex);
			if(done == 0)
				break;
			limit *= 1.221;
		}


		// now put points back into original array

		int pos = 0;
		for(int i = 0; i < nvertex; i++) {
			if(validate[i] == true) {
				vertex[pos].setTo(nvert[i]);

				printf("% 6.4lf % 6.4lf % 6.4lf \n", vertex[pos].x, vertex[pos].y, vertex[pos].z);
				pos++;
				if(pos >= number)
					break;
			}
		}

		delete[]nvert;
		delete[]validate;
	}


	void printSphFile(int number)
	{
		dvector *stdpoint;

		stdpoint = new dvector[number];

		createEvenVertexSphere(&stdpoint[0], 1, number);
		for(int i = 0; i < number; i++)
			printf("% 6.4lf % 6.4lf % 6.4lf \n", stdpoint[i].x, stdpoint[i].y, stdpoint[i].z);

		delete[]stdpoint;
	}

	int loadXYZcoords(char *filename, dvector * point, int number)
	{
		FILE *file;
		int i;
		char buffer[100];
		double x, y, z;

		file = fopen(filename, "r");
		if(file == NULL) 
		{
			printf("ERROR - file not found \n");
			return -1;
		}

		for(i = 0; (i < number) && (!feof(file)); i++) 
		{
			fgets(&buffer[0], 100, file);
			sscanf(&buffer[0], "%lf %lf %lf", &x, &y, &z);
			point[i].setTo(x, y, z);
		}
		// printf("read %d coordinates from %s \n",i,filename);
		fclose(file);
		return 0;
	}

	double calcNumericalSASA(dvector * atom, double *radius, int atoms, double probeRadius, char *SPfile, int SurfacePoints, // unit sphere surface point file
		double *SASA)
	{
		dvector *stdpoint;
		dvector point;
		double *dmatrix;
		int validpoints, valid;
		double distjp;
		int i, j, p;
		double SAi;
		double SASAtot = 0;

		stdpoint = new dvector[SurfacePoints];
		dmatrix = new double[sqr(atoms)];

		if(SPfile == NULL)
			createVertexSphere(&stdpoint[0], 1, SurfacePoints);
		else
			loadXYZcoords(SPfile, &stdpoint[0], SurfacePoints); // load coords

		calcDistanceMatrix(&dmatrix[0], atom, atoms);

		for(i = 0; i < atoms; i++) 
		{
			if(radius[i] <= 0.0001) 
			{
				validpoints = 0;
			} 
			else 
			{
				validpoints = SurfacePoints;
				for(p = 0; p < SurfacePoints; p++) 
				{
					point = stdpoint[p]; // get point
					point.mul((radius[i] + probeRadius)); // multiply by relevant radius + proberadius
					point.add(atom[i]); // move to coordinate of respective atom

					valid = 1;
					for(j = 0; j < atoms; j++) 
					{
						// check if atom is out of reach anyway
						if(j == i)
							continue;
						if(radius[j] <= 0.0001)
							continue;
						if(dmatrix[sqrmat(i, j, atoms)] > (radius[i] + radius[j] + 2 * probeRadius))
							continue;
						// if not, see if it occludes the point of interest
						// distjp = atom[j].dist(point);

						distjp = sqr(atom[j].x - point.x) + sqr(atom[j].y - point.y) + sqr(atom[j].z - point.z);

						if(distjp < sqr(radius[j] + probeRadius)) 
						{
							valid = 0;
							break;
						}
					}
					if(valid != 1)
						validpoints--;
				}
			}

			// printf("SASA %d %lf\n",i,((double)validpoints/(double)SurfacePoints) );

			SAi = ((double) validpoints / (double) SurfacePoints)
				* sphereSurfaceArea(radius[i] + probeRadius);
			SASA[i] = SAi;
			SASAtot += SAi;
			// printf("Atom %d: Surface Area: %lf %c\n",i,SAi,'%');
		}

		delete[]stdpoint;
		delete[]dmatrix;
		return SASAtot;
	}
}

