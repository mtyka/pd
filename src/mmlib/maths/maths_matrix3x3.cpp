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
#include "maths.h"

namespace Maths{
	void matrix3x3::setTo(const matrix3x3 & mat){
		r[0][0] = mat.r[0][0];
		r[0][1] = mat.r[0][1];
		r[0][2] = mat.r[0][2];
		r[1][0] = mat.r[1][0];
		r[1][1] = mat.r[1][1];
		r[1][2] = mat.r[1][2];
		r[2][0] = mat.r[2][0];
		r[2][1] = mat.r[2][1];
		r[2][2] = mat.r[2][2];
	}

	void matrix3x3::setToTranspose(const matrix3x3 & mat){
		r[0][0] = mat.r[0][0];
		r[0][1] = mat.r[1][0];
		r[0][2] = mat.r[2][0];
		r[1][0] = mat.r[0][1];
		r[1][1] = mat.r[1][1];
		r[1][2] = mat.r[2][1];
		r[2][0] = mat.r[0][2];
		r[2][1] = mat.r[1][2];
		r[2][2] = mat.r[2][2];
	}

	void matrix3x3::setToijk(const dvector & i,
													 const dvector & j,
													 const dvector & k){ // create a rotation matrix from 3 orthogonal ijk vectors
		r[0][0] = i.x;
		r[1][0] = i.y;
		r[2][0] = i.z;
		r[0][1] = j.x;
		r[1][1] = j.y;
		r[2][1] = j.z;
		r[0][2] = k.x;
		r[1][2] = k.y;
		r[2][2] = k.z;
	}

	// makes a start notation matrix out of a dvector
	void matrix3x3::setToStarVector(const dvector & vec){
		r[0][0] = 0;
		r[0][1] = -vec.z;
		r[0][2] = vec.y;
		r[1][0] = vec.z;
		r[1][1] = 0;
		r[1][2] = -vec.x;
		r[2][0] = -vec.y;
		r[2][1] = vec.x;
		r[2][2] = 0;
	}



	void matrix3x3::setToIdentity(){
		r[0][0] = 1;
		r[0][1] = 0;
		r[0][2] = 0;
		r[1][0] = 0;
		r[1][1] = 1;
		r[1][2] = 0;
		r[2][0] = 0;
		r[2][1] = 0;
		r[2][2] = 1;
	}

	void matrix3x3::setToNull(){
		r[0][0] = 0;
		r[0][1] = 0;
		r[0][2] = 0;
		r[1][0] = 0;
		r[1][1] = 0;
		r[1][2] = 0;
		r[2][0] = 0;
		r[2][1] = 0;
		r[2][2] = 0;
	}

	void matrix3x3::setTo(double r00, double r01, double r02,
												double r10, double r11, double r12,
												double r20, double r21, double r22){
			r[0][0] = r00;
			r[0][1] = r01;
			r[0][2] = r02;

			r[1][0] = r10;
			r[1][1] = r11;
			r[1][2] = r12;

			r[2][0] = r20;
			r[2][1] = r21;
			r[2][2] = r22;
	}

	void matrix3x3::setToXrot(double angle){
		double c = cos(angle);
		double s = sin(angle);

		r[0][0] = 1;
		r[0][1] = 0;
		r[0][2] = 0;

		r[1][0] = 0;
		r[1][1] = c;
		r[1][2] = -s;

		r[2][0] = 0;
		r[2][1] = s;
		r[2][2] = c;
	}

	void matrix3x3::setToYrot(double angle){
		double c = cos(angle);
		double s = sin(angle);

		r[0][0] = c;
		r[0][1] = 0;
		r[0][2] = s;

		r[1][0] = 0;
		r[1][1] = 1;
		r[1][2] = 0;

		r[2][0] = -s;
		r[2][1] = 0;
		r[2][2] = c;
	}

	void matrix3x3::setToZrot(double angle){
		double c = cos(angle);
		double s = sin(angle);

		r[0][0] = c;
		r[0][1] = -s;
		r[0][2] = 0;

		r[1][0] = s;
		r[1][1] = c;
		r[1][2] = 0;

		r[2][0] = 0;
		r[2][1] = 0;
		r[2][2] = 1;
	}


	/// create rotation matrix around x, y and z axes 
	void matrix3x3::setToXYZrot(double xangle, double yangle, double zangle){
		setToXrot(xangle);
		matrix3x3 rot; 
		rot.setToYrot(yangle);
		postmul(rot);
		rot.setToZrot(zangle);
		postmul(rot);
	}


	void matrix3x3::setToAxisRot(const dvector & axis, double angle)
	{
		double invmag = 1.0 / axis.mag();
		double c = cos(angle);
		double s = sin(angle);
		double t = 1.0 - c;
		double tx = t * axis.x * invmag;

		double sx = s * axis.x * invmag;
		double sy = s * axis.y * invmag;
		double sz = s * axis.z * invmag;
		double tyz = t * axis.y * invmag * axis.z * invmag;

		r[0][0] = tx * axis.x * invmag + c;
		r[0][1] = tx * axis.y * invmag - sz;
		r[0][2] = tx * axis.z * invmag + sy;

		r[1][0] = tx * axis.y * invmag + sz;
		r[1][1] = t * sqr(axis.y * invmag) + c;
		r[1][2] = tyz - sx;

		r[2][0] = tx * axis.z * invmag - sy;
		r[2][1] = tyz + sx;
		r[2][2] = t * sqr(axis.z * invmag) + c;
	}

	void matrix3x3::setToRandomRot(){
		double theta = Maths::MathConst::TwoPI * frand();
		double psi = Maths::MathConst::TwoPI * frand();
		double z = frand();

		dvector V(  cos(psi) * sqrt(z),
								sin(psi) * sqrt(z),
								sqrt( 1.0 - z));

		setTo( 2.0*V.x*V.x - 1.0, 2.0*V.x*V.y      ,2.0*V.x*V.z,
					 2.0*V.y*V.x      , 2.0*V.y*V.y - 1.0,2.0*V.y*V.z,
					 2.0*V.z*V.x      , 2.0*V.z*V.y      ,2.0*V.z*V.z - 1.0);
		matrix3x3 mat;
		mat.setTo( cos(theta), sin(theta), 0.0,
							-sin(theta), cos(theta), 0.0,
							 0.0       , 0.0       , 1.0 );
		postmul(mat);
	}


	void matrix3x3::add(const matrix3x3 & mat){
		r[0][0] += mat.r[0][0];
		r[0][1] += mat.r[0][1];
		r[0][2] += mat.r[0][2];
		r[1][0] += mat.r[1][0];
		r[1][1] += mat.r[1][1];
		r[1][2] += mat.r[1][2];
		r[2][0] += mat.r[2][0];
		r[2][1] += mat.r[2][1];
		r[2][2] += mat.r[2][2];
	}


	void matrix3x3::sub(const matrix3x3 & mat){
		r[0][0] -= mat.r[0][0];
		r[0][1] -= mat.r[0][1];
		r[0][2] -= mat.r[0][2];
		r[1][0] -= mat.r[1][0];
		r[1][1] -= mat.r[1][1];
		r[1][2] -= mat.r[1][2];
		r[2][0] -= mat.r[2][0];
		r[2][1] -= mat.r[2][1];
		r[2][2] -= mat.r[2][2];
	}


	void matrix3x3::postmul(const matrix3x3 & mmat){

		matrix3x3 re;

		re.r[0][0] = r[0][0] * mmat.r[0][0] + r[0][1] * mmat.r[1][0] + r[0][2] * mmat.r[2][0];
		re.r[1][0] = r[1][0] * mmat.r[0][0] + r[1][1] * mmat.r[1][0] + r[1][2] * mmat.r[2][0];
		re.r[2][0] = r[2][0] * mmat.r[0][0] + r[2][1] * mmat.r[1][0] + r[2][2] * mmat.r[2][0];

		re.r[0][1] = r[0][0] * mmat.r[0][1] + r[0][1] * mmat.r[1][1] + r[0][2] * mmat.r[2][1];
		re.r[1][1] = r[1][0] * mmat.r[0][1] + r[1][1] * mmat.r[1][1] + r[1][2] * mmat.r[2][1];
		re.r[2][1] = r[2][0] * mmat.r[0][1] + r[2][1] * mmat.r[1][1] + r[2][2] * mmat.r[2][1];

		re.r[0][2] = r[0][0] * mmat.r[0][2] + r[0][1] * mmat.r[1][2] + r[0][2] * mmat.r[2][2];
		re.r[1][2] = r[1][0] * mmat.r[0][2] + r[1][1] * mmat.r[1][2] + r[1][2] * mmat.r[2][2];
		re.r[2][2] = r[2][0] * mmat.r[0][2] + r[2][1] * mmat.r[1][2] + r[2][2] * mmat.r[2][2];
		this->setTo(re);
	}




	void matrix3x3::premul(const matrix3x3 & mmat){

		matrix3x3 re;

		re.r[0][0] = mmat.r[0][0] * r[0][0] + mmat.r[0][1] * r[1][0] + mmat.r[0][2] * r[2][0];
		re.r[1][0] = mmat.r[1][0] * r[0][0] + mmat.r[1][1] * r[1][0] + mmat.r[1][2] * r[2][0];
		re.r[2][0] = mmat.r[2][0] * r[0][0] + mmat.r[2][1] * r[1][0] + mmat.r[2][2] * r[2][0];

		re.r[0][1] = mmat.r[0][0] * r[0][1] + mmat.r[0][1] * r[1][1] + mmat.r[0][2] * r[2][1];
		re.r[1][1] = mmat.r[1][0] * r[0][1] + mmat.r[1][1] * r[1][1] + mmat.r[1][2] * r[2][1];
		re.r[2][1] = mmat.r[2][0] * r[0][1] + mmat.r[2][1] * r[1][1] + mmat.r[2][2] * r[2][1];

		re.r[0][2] = mmat.r[0][0] * r[0][2] + mmat.r[0][1] * r[1][2] + mmat.r[0][2] * r[2][2];
		re.r[1][2] = mmat.r[1][0] * r[0][2] + mmat.r[1][1] * r[1][2] + mmat.r[1][2] * r[2][2];
		re.r[2][2] = mmat.r[2][0] * r[0][2] + mmat.r[2][1] * r[1][2] + mmat.r[2][2] * r[2][2];
		this->setTo(re);
	}





	void matrix3x3::mul(double k){
		r[0][0] *= k;
		r[0][1] *= k;
		r[0][2] *= k;

		r[1][0] *= k;
		r[1][1] *= k;
		r[1][2] *= k;

		r[2][0] *= k;
		r[2][1] *= k;
		r[2][2] *= k;
	}

	void matrix3x3::div(double k){
		mul(1 / k);
	}


	void matrix3x3::orthonormalize(){
		dvector X(r[0][0], r[1][0], r[2][0]);
		dvector Y(r[0][1], r[1][1], r[2][1]);
		dvector Z;

		X.unify();
		Z.crossProduct(X, Y);
		Z.unify();
		Y.crossProduct(Z, X);
		Y.unify();

		r[0][0] = X.x;
		r[1][0] = X.y;
		r[2][0] = X.z;

		r[0][1] = Y.x;
		r[1][1] = Y.y;
		r[2][1] = Y.z;

		r[0][2] = Z.x;
		r[1][2] = Z.y;
		r[2][2] = Z.z;
	}


	void matrix3x3::transpose(){
		matrix3x3 T;
		T.setTo(*this);
		setToTranspose(T);
	}

	void matrix3x3::invert(){
		matrix3x3 in;
		double det = determinant();
		if(det == 0)
			return; // matrix is singular

		in.r[0][0] = r[1][1] * r[2][2] - r[2][1] * r[1][2];
		in.r[0][1] = -r[0][1] * r[2][2] + r[2][1] * r[0][2];
		in.r[0][2] = r[0][1] * r[1][2] - r[1][1] * r[0][2];

		in.r[1][0] = -r[1][0] * r[2][2] + r[2][0] * r[1][2];
		in.r[1][1] = r[0][0] * r[2][2] - r[2][0] * r[0][2];
		in.r[1][2] = -r[0][0] * r[1][2] + r[1][0] * r[0][2];

		in.r[2][0] = r[1][0] * r[2][1] - r[2][0] * r[1][1];
		in.r[2][1] = -r[0][0] * r[2][1] + r[2][0] * r[0][1];
		in.r[2][2] = r[0][0] * r[1][1] - r[1][0] * r[0][1];

		in.div(det);

		setTo(in);
	}

	double matrix3x3::determinant() const{ // returns the determiannt of the 3x3 matrix
		return r[0][0] * (r[1][1] * r[2][2] - r[2][1] * r[1][2])
				 - r[0][1] * (r[1][0] * r[2][2] - r[2][0] * r[1][2])
				 + r[0][2] * (r[1][0] * r[2][1] - r[2][0] * r[1][1]);
	}

	double matrix3x3::normalise(){ // normalises so that sum of the matrix is 9
		double factor =
			fabs(r[0][0]) + fabs(r[1][0]) + fabs(r[2][0]) +
			+fabs(r[0][1]) + fabs(r[1][1]) + fabs(r[2][1]) + +fabs(r[0][2]) + fabs(r[1][2]) + fabs(r[2][2]);
		factor /= 9.0;
		if(factor == 0)
			return 1.0;

		double invfactor = 1.0 / factor;

		this->mul(invfactor);
		return factor;
	}

	// diagonalises a symetric 3x3 matrix by finding the three
	// eigen values and associated eigenvectors
	// this function *enforces* symetricality by crosscopying the correspondent entries
	int matrix3x3::diagonaliseSymetric(
		double &lambda1,
		double &lambda2,
		double &lambda3,
		dvector &v1,
		dvector &v2,
		dvector &v3
	)
	{
			int lam; // counts through the lambdas
			double lambda[3]; // eigenvalues

			double det1; // determinant
			double det2; // determinant
			double det3; // determinant
			int roots; // nr of roots of characteristic equation
			double a1, a2, a3; // polynomial coefficients
			matrix3x3 A, thisSave;
			dvector v, va, vb, vc, w;
			double factor = 1.0;

			thisSave.setTo(*this); // save this
			factor = this->normalise();

			// enforce symetricality - WARNING: this obliterates any non symetric matrix !
			r[1][0] = r[0][1];
			r[2][0] = r[0][2];
			r[2][1] = r[1][2];

			// to calculate the eigenvalues, solve det(A-lambdaI)=0

			// construct characteristic polynomial by determining the coefficients from
			// the determinant formula
			a1 = -(r[0][0] + r[1][1] + r[2][2]);
			a2 =
				-(r[2][1] * r[1][2] + r[0][1] * r[1][0] + r[0][2] * r[2][0] - r[0][0] * r[2][2] - r[0][0] * r[1][1] -
				r[1][1] * r[2][2]);
			a3 = -(r[0][0] * r[1][1] * r[2][2]
			- r[0][0] * r[2][1] * r[1][2]
			- r[2][2] * r[0][1] * r[1][0]
			- r[0][2] * r[2][0] * r[1][1]
			+ r[0][1] * r[2][0] * r[1][2] + r[2][0] * r[1][0] * r[2][1]);

			// find the roots of the polynomial, which will equal the three lambda's (not nessessarily different !)
			roots = solveCubicPolynomial(a1, a2, a3, lambda[0], lambda[1], lambda[2]);
			/*
			printf("coeffs : %e %e %e \n",a1,a2,a3);
			printf("lambdas: %e %e %e \n",lambda[0]*factor,lambda[1]*factor,lambda[2]*factor);
			printf("Polynomial solutions: %e %e %e \n",
			calcCubicPolynomial(a1,a2,a3,lambda[0]),
			calcCubicPolynomial(a1,a2,a3,lambda[1]),
			calcCubicPolynomial(a1,a2,a3,lambda[2]));
			*/

			// a symetric matrix *should* have three double eigenvalues - if not, something's seriously wrong
			if(roots == 1) {
				//printf("MATH ERROR: maths.cpp: diagonaliseSymetric(..): \"Characteristic Polynomial with only 1 double root\" \n");
				return -1;
			}
			// give lambdas to the caller
			lambda1 = lambda[0];
			lambda2 = lambda[1];
			lambda3 = lambda[2];

			// now find eigenvectors by solving (A-lambdaI)V=0 for each eigenvalue lambda
			for(lam = 0; lam <= 2; lam++) {
				A.setTo(*this);
				A.r[0][0] -= lambda[lam];
				A.r[1][1] -= lambda[lam];
				A.r[2][2] -= lambda[lam];

				// printf("A-lambdaI= \n");
				// A.info();
				det1 = A.r[0][0] * A.r[1][1] - A.r[0][1] * A.r[1][0];
				det2 = A.r[1][0] * A.r[2][1] - A.r[1][1] * A.r[2][0];
				det3 = A.r[2][0] * A.r[0][1] - A.r[2][1] * A.r[0][0];

				//printf("determinants = \n%.20lf\n%.20lf\n%.20lf\n%.20lf\n",
				// det,det ,det2,det3);

				// check if the eigenvector is along a particular axis
				if((fabs(A.r[0][0]) < 1e-12) && (fabs(A.r[1][0]) < 1e-12) && (fabs(A.r[2][0]) < 1e-12))
					v.setTo(1, 0, 0);
				else if((fabs(A.r[0][1]) < 1e-12) && (fabs(A.r[1][1]) < 1e-12) && (fabs(A.r[2][1]) < 1e-12))
					v.setTo(0, 1, 0);
				else if((fabs(A.r[0][2]) < 1e-12) && (fabs(A.r[1][2]) < 1e-12) && (fabs(A.r[2][2]) < 1e-12))
					v.setTo(0, 0, 1);
				else {/*
					  // double sum1,sum2,sum3;
					  // sum1 = + A.r[1][0] + A.r[2][0];
					  // sum2 = A.r[0][1] + A.r[1][1] + A.r[2][1];
					  // sum3 = A.r[0][2] + A.r[1][2] + A.r[2][2];
					  // printf("A - lambdaI:\n");

					  //printf("%lf %lf %lf\n",sum1,sum2,sum3);
					  */
					/* if((1+fabs(det)) < 1e-100){ // if subdeterminant is 0 then we dont need to solve
					//printf("determinant == 0\n");

					v.setTo((-A.r[0][1]-A.r[0][2])/A.r[0][0],1,1);
					w.setTo(v);
					w.mulmat(A);
					if(w.mag()>0.0000000001){
					v.setTo(1,1,(-A.r[2][0]-A.r[2][1])/A.r[2][2]);
					w.setTo(v);
					w.mulmat(A);
					if(w.mag()>0.000000001){
					v.setTo(1,(-A.r[1][0]-A.r[1][2])/A.r[1][1],1);
					w.setTo(v);
					w.mulmat(A);
					}
					}

					}else
					*/// if subdeterminant is not 0, then we can solve the top two
					// linear equations by 2x2 matrix inversion

					v.setTo((A.r[1][1] * (-A.r[0][2]) - A.r[0][1] * (-A.r[1][2])),
						(-A.r[1][0] * (-A.r[0][2]) + A.r[0][0] * (-A.r[1][2])), det1);
					if((fabs(v.x) < 1e-12) && (fabs(v.y) < 1e-12)) {
						v.setTo((A.r[2][1] * (-A.r[1][2]) - A.r[1][1] * (-A.r[2][2])),
							(-A.r[2][0] * (-A.r[1][2]) + A.r[1][0] * (-A.r[2][2])), det2);
						if((fabs(v.x) < 1e-12) && (fabs(v.y) < 1e-12)) {
							v.setTo((-A.r[0][0] * (-A.r[2][2]) + A.r[2][0] * (-A.r[0][2])),
								(A.r[0][1] * (-A.r[2][2]) - A.r[2][1] * (-A.r[0][2])), det3);
							if((fabs(v.x) < 1e-12) && (fabs(v.y) < 1e-12)) {
								v.setTo((-A.r[0][1] - A.r[0][2]) / A.r[0][0], 1, 1);
								// v.setTo(1,1,(-A.r[2][0]-A.r[2][1])/A.r[2][2]);
								/*
								w.setTo(v);
								w.mulmat(A);
								if(w.mag()>0.0000000001){
								v.setTo(1,1,(-A.r[2][0]-A.r[2][1])/A.r[2][2]);
								w.setTo(v);
								w.mulmat(A);
								if(w.mag()>0.000000001){
								v.setTo(1,(-A.r[1][0]-A.r[1][2])/A.r[1][1],1);
								w.setTo(v);
								w.mulmat(A);
								}}
								*/
							}
						}
					}
				}

				//v.info();
				if(v.mag() < 1e-100) {
					v.info();
					printf("MATHS ERROR: trivial solution (0,0,0) found in Diag..Symetric for lambda %.20lf(...)\n", lambda[lam]);
				}

				v.unify();
				// printf("RESULT:-------------------------------\n");
				//v.info();
				if(lam == 0)
					v1.setTo(v);
				else if(lam == 1)
					v2.setTo(v);
				else if(lam == 2)
					v3.setTo(v);



				//v1.info();
				w.setTo(v);
				w.mulmat(A);
				if(w.mag() > 0.001) {
					printf("MATH ERROR: Invalid solution found in DiagonaliseSymetric(..)\n");
					printf("Eigenvector: for lambda=%.20lf : ",lambda[lam]);
					v.info();
					printf("Result (should be 0,0,0):");
					w.info();
					printf("Original matrix (rescaled by factor= %.20lf): \n", factor);
					mul(factor);
					info();
					div(factor);
					printf("Original matrix (rescaled) minus lambda (also rescaled) \n");
					A.info();
					printf("det: %.20lf \n", det1);
					printf("Vs..: \n");
					va.info();
					vb.info();
					vc.info();
				}
			}

			// now rescale eigenvalues back to what they should have been
			// without the normalisation
			// the vectors are left as they are, since they are redundant
			// in terms of their length
			lambda1 *= factor;
			lambda2 *= factor;
			lambda3 *= factor;

			this->setTo(thisSave);

			return 0;
	}


	void matrix3x3::info() const{
		printf(" / %13.10e %13.10e %13.10e \\\n", r[0][0], r[0][1], r[0][2]);
		printf(" | %13.10e %13.10e %13.10e |\n", r[1][0], r[1][1], r[1][2]);
		printf(" \\ %13.10e %13.10e %13.10e /\n\n", r[2][0], r[2][1], r[2][2]);
	}

	void matrix3x3::printMathematicaFormat() const{
		printf("{");
		for(int i = 0; i < 3; i++) {
			printf("{");
			for(int j = 0; j < 3; j++) {
				double prefix, exponent;
				getExpForm(r[i][j], prefix, exponent);
				printf("%lf*10^%d ", prefix, (int) exponent);
				if(j < 2)
					printf(",");
			}
			if(i < 2)
				printf("},");
			else
				printf("}");
		}
		printf("}\n");
	}



}
