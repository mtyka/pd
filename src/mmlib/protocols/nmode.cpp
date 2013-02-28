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

// Nmode.cpp/h  - provides functionality for Hessian & Covariance matrices & associated
// normal mode calculations
//
//               Matrix3Nx3N                    
//          .---------^---------.               
//       Hessian            CovarianceMatrix    
//          |                                   
//     Hessian_Numerical
//
//
//               EigenSystem
//          .----------^---------.
//          |                    |
//   EigenSystem_Hessian  EigenSystem_Covariance
//
//
//
//

#include "global.h"

// Matrix Library providing the code for large matrix diagonalisation
//#include "maths/tntjama/tnt.h" - included in header
#include "maths/tntjama/jama_eig.h"

#include "tools/rdstdout.h"
#include "workspace/workspace.h"
#include "workspace/neighbourlist.h"
#include "forcefields/ffparam.h"
#include "forcefields/forcefield.h"

// Own header include
#include "protocols/nmode.h"

// Namespace includes
using namespace Physics;
using namespace TNT;
using namespace JAMA;
using namespace Maths;

namespace Protocol {


  //-------------------------------------------------------------------
	//
	// Matrix3Nx3N 

	void Matrix3Nx3N::printMathematica(){
		int i,j;
		for(i=0;i<( matrix->m_ );i++){
			printf("{");
			for(j=0;j<(  matrix->n_ );j++){
				double prefix,exponent;
				getExpForm( (double)(matrix->data_[sqrmat(i,j, matrix->m_)]), prefix, exponent);
				printf("%lf*10^%lf ",prefix,exponent);
				if(j<(matrix->n_-1)) printf(",");
			}
			if(i<(matrix->m_-1)) printf("},");
			else           printf("}");
		}
		printf("}\n");
	}





  //-------------------------------------------------------------------
	//
	// Hessian 

	void Hessian::massWeight(WorkSpace &wspace){
		int i,j;

		if( matrix == NULL ){
			throw( ProcedureException("Cannot reweight empty Hessian"));
		}
		if( matrix->m_ != (wspace.atom.size() * 3)||
				matrix->n_ != (wspace.atom.size() * 3) ){
			throw( ProcedureException("Number fo particles in workspace * 3 does not equal widht of hessian"));
		}

		// mass weight the matrix
		for(i = 0; i < (wspace.atom.size() * 3); i++) {
			//printf("%e \n", wspace.atom[int (i / 3)].mass);
			for(j = 0; j < (wspace.atom.size() * 3); j++) {
				matrix->data_[sqrmat(i, j, (int)wspace.atom.size() * 3)] /=
					(float)(sqrt(wspace.atom[int (i / 3)].mass * wspace.atom[int (j / 3)].mass));

			}
		}
	}



	
	void Hessian::writeRaw(const std::string &filename){  // writes raw data to a file
		FILE *file=NULL;

		file = fopen(filename.c_str(), "wb" );
		if(file==NULL){
			throw( IOException(std::string("Error writing to file '") +
												 filename + 
												 std::string("' in Hessian::writeRaw(const std::string &filename)")) );
		}
		int dim;

		dim = matrix->m_;
		fwrite((void*)&dim, sizeof(int) , 1, file);
		dim = matrix->n_;
		fwrite((void*)&dim, sizeof(int) , 1, file);

		int i;
		double data;
		for(i = 0; i<(matrix->m_*matrix->n_) ; i++ ){
			data = matrix->data_[i];
			fwrite((void*)&data, sizeof(LINALG_real), 1  , file );
		}
		for(i = 0; i<x0.nAtoms() ; i++ ){
			data = x0.atom[i].p.x;
			fwrite((void*)&data, sizeof(LINALG_real), 1  , file );
			data = x0.atom[i].p.y;
			fwrite((void*)&data, sizeof(LINALG_real), 1  , file );
			data = x0.atom[i].p.z;
			fwrite((void*)&data, sizeof(LINALG_real), 1  , file );
		}


		fclose(file);		
	}

	void Hessian::readRaw(const std::string &filename){   // reads raw data back froma a file
		FILE *file=NULL;

		file = fopen(filename.c_str(), "rb" );
		if(file==NULL){
			throw( IOException(std::string("Error opening file '") +
												 filename + 
												 std::string("' in Hessian::readRaw(const std::string &filename)")) );
		}
		
		int m,n;

		fread((void*)&m, sizeof(int) , 1, file);
		fread((void*)&n, sizeof(int) , 1, file);
		if( m != n ){
			throw( ParseException("Dimension of hessian do not appear to be square, i.e. m != n\n - cannot continue reading in data "));
		}
		Matrix3Nx3N::createEmpty(m,n);

		printf("reading in Covariance matrix from %s (%d x %d)\n",filename.c_str(), m, n );
		int i;
		double data;
		for(i = 0; i<(matrix->m_*matrix->n_) ; i++ ){
			fread((void*)&data, sizeof(LINALG_real), 1  , file );
			matrix->data_[i] = (float)(data);
		}

		printf("reading in Average positions %s (%d)\n",filename.c_str(), m );
		x0 = SnapShot( matrix->m_ );
		for(i = 0; i< matrix->m_ ; i++ ){
			fread((void*)&data, sizeof(LINALG_real), 1  , file );
			x0.atom[i].p.x = data;
			fread((void*)&data, sizeof(LINALG_real), 1  , file );
			x0.atom[i].p.y = data;
			fread((void*)&data, sizeof(LINALG_real), 1  , file );
			x0.atom[i].p.z = data;
		}

		fclose(file);
	}





	




  //-------------------------------------------------------------------
	//
	// Hessian_Numerical 

	#define MAX_NMEASURES 10
	struct dFstruct{
		double    dFx[MAX_NMEASURES * 2 + 1];
		double    dFy[MAX_NMEASURES * 2 + 1];
		double    dFz[MAX_NMEASURES * 2 + 1];
	};

	int Hessian_Numerical::runcore(){
		//proxies
		WorkSpace &wspace = getWSpace();
		SnapShot &cur = getWSpace().cur;
		SnapShot &old = getWSpace().old;
		size_t nAtoms = getWSpace().atom.size();

		// other variables
		int i, j;

		double disp;
		dvector orig;
		int ci;
		double dxarray[MAX_NMEASURES * 2 + 1];
		dFstruct *dF = new dFstruct[nAtoms];
		double dFxdx;
		double dFydx;
		double dFzdx;
		double a, R;
	
		if( (FiniteDiffSteps/2) >= MAX_NMEASURES ){
			throw(ArgumentException("The maximum number of force measurements during the finite difference calculation of the numerical hessian is 20"));
		}

		int fd_points_before = FiniteDiffSteps / 2;
		int fd_points_after = FiniteDiffSteps / 2 + FiniteDiffSteps % 2;


		int pos = 0;
		int nmes;
		long starttime = (long)time(NULL);

		if(OutputLevel){
			printf("--- Numerical determination of Hessian -----------\n");
			printf("Degrees of Freedom: %d \n",nAtoms * 3);
			printf("Force evaluations:  %d \n",nAtoms * FiniteDiffSteps * 3+1);
			printf("Finite difference steps:  %d ", FiniteDiffSteps); 
			printf("Finite difference ladder: "); 
			for(int p = -fd_points_before; p <= fd_points_after; p++) {
				if( p==0 ){ printf("0, "); continue; }
				printf("%d*StepSize, ", p);
			}
			printf("\n");
		}

		createEmpty(wspace);

		wspace.nlist().calcNewList();
		ff->calcForces();

		// print residual force magnitude
		a = 0;
		for(i = 0; i < nAtoms; i++)
			a += cur.atom[i].f.mag();
		a /= nAtoms;

		if(OutputLevel){
			printf("Average residual force magnitude:  %6.4e (N)   %lf  (kcal/A)\n", (double) a, (double) (a / 6.9525E-11));
		}
		// save original forces and original positions
		for(i = 0; i < nAtoms; i++) {
			old.atom[i].f.setTo(cur.atom[i].f);
			old.atom[i].p.setTo(cur.atom[i].p);
		}

		x0 = wspace.save();
		e0 = wspace.ene.epot;

		for(i = 0; i < nAtoms; i++) {
			if((OutputLevel)&&(UpdateScr > 0)){
				if((i%UpdateScr)==0){
					printf("Progress: %d/%d \n",i,nAtoms);
				}
			}
			for(ci = 0; ci < 3; ci++) {
				
				// first of all restore all the atom positions
				for(j = 0; j < nAtoms; j++) {	
					cur.atom[j].p.setTo(old.atom[j].p);
				}

				dxarray[0] = 0;
				for(j = 0; j < nAtoms; j++) {
					dF[j].dFx[0] = old.atom[j].f.x;
					dF[j].dFy[0] = old.atom[j].f.y;
					dF[j].dFz[0] = old.atom[j].f.z;
				}
				nmes = 0;
				for(disp = (-StepSize * (double) fd_points_before); disp <= (StepSize * (double) fd_points_after); disp += StepSize) {
					if( fabs(disp) < fabs(StepSize/2.0) ){
						dxarray[nmes] = disp * Physics::PhysicsConst::Angstrom;
						for(j = 0; j < nAtoms; j++) {
							dF[j].dFx[nmes] = 0; 
							dF[j].dFy[nmes] = 0; 
							dF[j].dFz[nmes] = 0; 
						}
					}else{
						cur.atom[i].p.setTo(old.atom[i].p);
						switch (ci) {
						case 0:
							cur.atom[i].p.x += (double)disp;
							break;
						case 1:
							cur.atom[i].p.y += (double)disp;
							break;
						case 2:
							cur.atom[i].p.z += (double)disp;
							break;
						}
						ff->calcForces();

						dxarray[nmes] = disp * Physics::PhysicsConst::Angstrom;
						for(j = 0; j < nAtoms; j++) {
							dF[j].dFx[nmes] = (cur.atom[j].f.x - old.atom[j].f.x);
							dF[j].dFy[nmes] = (cur.atom[j].f.y - old.atom[j].f.y);
							dF[j].dFz[nmes] = (cur.atom[j].f.z - old.atom[j].f.z);
						}
					}

					nmes++;
					
				}
				for(j = 0; j < nAtoms; j++) {

					leastSquaresFit(&dxarray[0], &dF[j].dFx[0], nmes, a, dFxdx, R);
					matrix->data_[pos] = (LINALG_real) -dFxdx;

					pos++;
					leastSquaresFit(&dxarray[0], &dF[j].dFy[0], nmes, a, dFydx, R);
					matrix->data_[pos] = (LINALG_real)-dFydx;

					pos++;
					leastSquaresFit(&dxarray[0], &dF[j].dFz[0], nmes, a, dFzdx, R);
					matrix->data_[pos] = (LINALG_real)-dFzdx;

					pos++;
				}
			}
		}

		// finally restore all the atom positions and recalculate the zero point energy
		for(j = 0; j < nAtoms; j++) {	
			cur.atom[j].p.setTo(old.atom[j].p);
		}
		wspace.nlist().calcNewList();
		ff->calcForces();

		if(ForceSymmetry){

			if(OutputLevel){
				printf("Forcing hessian symmetry (taking average of H[i,j] and H[j,i])\n");
			}
			// force hessian to be symetric
			// (because of numerical inaccuracies H[i,j] will differ from H[j,i]
			// so set both to average
			LINALG_real hij;
			for(i = 0; i < (nAtoms * 3); i++) {
				for(j = i; j < (nAtoms * 3); j++) {
					hij = matrix->data_[sqrmat(i, j, (int)nAtoms * 3)] + 
								matrix->data_[sqrmat(j, i, (int)nAtoms * 3)];
					hij /= 2;
					matrix->data_[sqrmat(i, j, (int)nAtoms * 3)] = hij;
					matrix->data_[sqrmat(j, i, (int)nAtoms * 3)] = hij;

				}
			}
		}

		long timediff = (long)time(NULL) - starttime;

		printf("Time taken: %ld sec  (%d hrs %d mins %d secs) \n",
				  timediff,
			   (timediff)/3600,
				((timediff)%3600)/60,
				((timediff)%3600)%60		);

	
		delete[]dF;
		return 0;
	}


	







  //-------------------------------------------------------------------
	//
	//  CovarianceMatrix


	void	CovarianceMatrix::reset(){
		int i,j;
		ndata = 0;
		int matcnt=0;
		for(i=0;i<getWSpace().atom.size();i++){
			x0.atom[i].p.zero();
			for(j = 0; j < getWSpace().atom.size(); j++) {
				matrix->data_[matcnt] = 0.0; matcnt++;
				matrix->data_[matcnt] = 0.0; matcnt++;
				matrix->data_[matcnt] = 0.0; matcnt++;
			}																		 
			for(j = 0; j < getWSpace().atom.size(); j++) {
				matrix->data_[matcnt] = 0.0; matcnt++;
				matrix->data_[matcnt] = 0.0; matcnt++;
				matrix->data_[matcnt] = 0.0; matcnt++;
			}																			 
			for(j = 0; j < getWSpace().atom.size(); j++) {
				matrix->data_[matcnt] = 0.0; matcnt++;
				matrix->data_[matcnt] = 0.0; matcnt++;
				matrix->data_[matcnt] = 0.0; matcnt++;
			}
		}
		
	}


	void	CovarianceMatrix::setcurdata(){
		int i,j;
		SnapShotAtom *atom = getWSpace().cur.atom;
		if( matrix == NULL ){
			throw( ProcedureException("Hessian not initialised"));
		}
		if( matrix->m_ != (getWSpace().atom.size() * 3)||
				matrix->n_ != (getWSpace().atom.size() * 3) ){
			throw( ProcedureException("Number fo particles in workspace * 3 does not equal width of hessian"));
		}
		int matcnt=0;
		for(i = 0; i < getWSpace().atom.size(); i++) {
			x0.atom[i].p.add(atom[i].p);
			for(j = 0; j < getWSpace().atom.size(); j++) {
				matrix->data_[matcnt] += (float)(atom[i].p.x * atom[j].p.x);matcnt++;
				matrix->data_[matcnt] += (float)(atom[i].p.x * atom[j].p.y);matcnt++;
				matrix->data_[matcnt] += (float)(atom[i].p.x * atom[j].p.z);matcnt++;
			}

			for(j = 0; j < getWSpace().atom.size(); j++) {
				matrix->data_[matcnt] += (float)(atom[i].p.y * atom[j].p.x);matcnt++;
				matrix->data_[matcnt] += (float)(atom[i].p.y * atom[j].p.y);matcnt++;
				matrix->data_[matcnt] += (float)(atom[i].p.y * atom[j].p.z);matcnt++;
			}

			for(j = 0; j < getWSpace().atom.size(); j++) {
				matrix->data_[matcnt] += (float)(atom[i].p.z * atom[j].p.x);matcnt++;
				matrix->data_[matcnt] += (float)(atom[i].p.z * atom[j].p.y);matcnt++;
				matrix->data_[matcnt] += (float)(atom[i].p.z * atom[j].p.z);matcnt++;
			}
		}
		
		ndata++;
	}


	void	CovarianceMatrix::finish(bool mass_weight,
																 bool remove_correlations ){
		int i,j;

		printf("Finish CovarianceMatrix calculation: Measurements: %d \n", ndata);
		
		for(i = 0; i < (matrix->m_*matrix->n_); i++){
			matrix->data_[i] /= (float) ndata;
		}
		
		for(i = 0; i <  getWSpace().atom.size(); i++){
			x0.atom[i].p.mul(1.0/(float) ndata);
		}

		// cov(xi,xj) = <xixj> - <xi><xj>
		for(i = 0; i < getWSpace().atom.size(); i++) 
		{
				for(j = 0; j < getWSpace().atom.size(); j++) 
				{
					matrix->data_[sqrmat(i * 3 + 0, j * 3 + 0, (int)getWSpace().atom.size() * 3)] -= (float)(x0.atom[i].p.x * x0.atom[j].p.x);
					matrix->data_[sqrmat(i * 3 + 0, j * 3 + 1, (int)getWSpace().atom.size() * 3)] -= (float)(x0.atom[i].p.x * x0.atom[j].p.y);
					matrix->data_[sqrmat(i * 3 + 0, j * 3 + 2, (int)getWSpace().atom.size() * 3)] -= (float)(x0.atom[i].p.x * x0.atom[j].p.z);
																																															
					matrix->data_[sqrmat(i * 3 + 1, j * 3 + 0, (int)getWSpace().atom.size() * 3)] -= (float)(x0.atom[i].p.y * x0.atom[j].p.x);
					matrix->data_[sqrmat(i * 3 + 1, j * 3 + 1, (int)getWSpace().atom.size() * 3)] -= (float)(x0.atom[i].p.y * x0.atom[j].p.y);
					matrix->data_[sqrmat(i * 3 + 1, j * 3 + 2, (int)getWSpace().atom.size() * 3)] -= (float)(x0.atom[i].p.y * x0.atom[j].p.z);
																																															
					matrix->data_[sqrmat(i * 3 + 2, j * 3 + 0, (int)getWSpace().atom.size() * 3)] -= (float)(x0.atom[i].p.z * x0.atom[j].p.x);
					matrix->data_[sqrmat(i * 3 + 2, j * 3 + 1, (int)getWSpace().atom.size() * 3)] -= (float)(x0.atom[i].p.z * x0.atom[j].p.y);
					matrix->data_[sqrmat(i * 3 + 2, j * 3 + 2, (int)getWSpace().atom.size() * 3)] -= (float)(x0.atom[i].p.z * x0.atom[j].p.z);
			}
		}

		// adjust units to SI
		for(i = 0; i < (matrix->m_*matrix->n_); i++){
			matrix->data_[i] *= (float)(sqr(PhysicsConst::Angstrom));
		}
	
		// mass weight the matrix if normal frequencies are required
		if(mass_weight){
			massWeight();
		}

		// remove correlations - i.e. remove off diagonal covariance 
		// elements by setting them to 0
		if(remove_correlations){
			for(i = 0; i < (getWSpace().atom.size() * 3); i++) {
				for(j = 0; j < (getWSpace().atom.size() * 3); j++) {
					if(i!=j){
						matrix->data_[sqrmat(i, j, (int)getWSpace().atom.size() * 3)] = 0;
					}
				}
			}
		}

	}


	void CovarianceMatrix::massWeight(){
		int i,j;
		for(i = 0; i < (getWSpace().atom.size() * 3); i++) {
			for(j = 0; j < (getWSpace().atom.size() * 3); j++) {
				matrix->data_[sqrmat(i, j, (int)getWSpace().atom.size() * 3)] *=
					(float)(sqrt(getWSpace().atom[int (i / 3)].mass * getWSpace().atom[int (j / 3)].mass));
			}
		}
	}

  //-------------------------------------------------------------------
	//
	//  EigenSystem


	void EigenSystem::diagonalise(TNT::Array2D < LINALG_real > &matrix,
																Verbosity::Type OutputLevel){
	
		long sttime;
		if(OutputLevel){
			printf("Diagonalising hessian ... ");
			sttime = (long) time(NULL);
		}
		double average=0;
		double count=0;
		for ( int i = 0; i < matrix.m_ ; i ++ ){
			for ( int j = 0; j < matrix.m_ ; j ++ ){
					average+= matrix.data_[sqrmat(i, j, (int)matrix.m_)];
					count++;
			}
		}
		average/=count;
		average = sqrt(average);
		average = 1.0/average;
		for ( int i = 0; i < matrix.m_ ; i ++ )
		{
			for ( int j = 0; j < matrix.m_ ; j ++ )
			{
					matrix.data_[sqrmat(i, j, (int)matrix.m_)] *= (float)(average); 
					matrix.data_[sqrmat(j, i, (int)matrix.m_)] *= (float)(average);
			}
		}

		Eigenvalue < LINALG_real >eigenmymat(matrix);
		long entime = (long) time(NULL);

		if(OutputLevel){
			printf("Time taken: %d secs\n", entime - sttime);
		}

		eigenmymat.getRealEigenvalues(eigenvalue_real);
		eigenmymat.getImagEigenvalues(eigenvalue_imag);

		eigenmymat.getV(eigenvector);

		for(int i = 0; i < eigenvalue_real.dim1() ; i ++ ) eigenvalue_real[i] /= (LINALG_real)(average);
		for(int i = 0; i < eigenvalue_imag.dim1() ; i ++ ) eigenvalue_real[i] /= (LINALG_real)(average);

		for ( int i = 0; i < matrix.m_ ; i ++ ){
			for ( int j = 0; j < matrix.m_ ; j ++ ){
					matrix.data_[sqrmat(i, j, (int)matrix.m_)] /= (float)(average); 
					matrix.data_[sqrmat(j, i, (int)matrix.m_)] /= (float)(average);
			}
		}

	}

	void EigenSystem::printEigenValues(){
		printf("N      Real      Imaginary  \n");
		for(int e=0;e<getWSpace().atom.size()*3;e++){
			printf("%4d\t%10.8e\t%10.8e\n", e, eigenvalue_real[e],eigenvalue_imag[e]);
	  }
	}
	

	/// This calculates the vibrational entropy of each eigen frequency using the relation ship:
	///
	///            Vib. Energy Level  E(i) = v*h*i
	///            Partition function Z = 1/(1-exp(-v*h/kT))
	///            Entropy            S = k ln (Z)
	///

	double EigenSystem::calcVibEntropy(double temp){
		int imode;
		double _vibentropy = 0;
		double _modeentropy;


		printf("--- Normal Modes ---------------------------- \n");
		for(imode = 0; imode < nmode_freq.size() ; imode++) {
			_modeentropy = log(1.0 / (1.0 - exp((double)(-nmode_freq[imode]) * PhysicsConst::planck / (PhysicsConst::kB * temp))));
			_vibentropy += _modeentropy;
			printf("Mode %4d  % 10.8e (Hz)   %10.3lf cm-1    %lf kcal/mol\n",
				(int) imode, (double) (nmode_freq[imode]),
				(double) (nmode_freq[imode] / (PhysicsConst::clight * 100)), (double) (-_modeentropy * (PhysicsConst::J2kcal * PhysicsConst::kB * PhysicsConst::Na * temp)));
		}
		_vibentropy *= PhysicsConst::kB;		// S = kB*ln(Z)(S = kB ln (Qvib)
		printf("--- Quantum Entropy ( S = kB.ln(Qvib) )---- \n");
		printf("Temperature:   %8.3lf  K\n",temp);
		printf("ZeroEnergy:    % 10.5lf kcal/mol\n",
			 getWSpace().ene.epot * PhysicsConst::J2kcal * PhysicsConst::Na);
		printf("Entropy:       % 10.5lf kcal/mol/K \n",
			 _vibentropy * PhysicsConst::J2kcal * PhysicsConst::Na);
		printf("Entropy(-TS):  % 10.5lf kcal/mol \n",
			 -_vibentropy * temp * PhysicsConst::J2kcal * PhysicsConst::Na);
		
		return _vibentropy;
	}

	double EigenSystem::calcVibEntropyFull(double temp){
		int imode;
		double _vibentropy = 0;
		double _modeentropy;
		double homegaoverkt;

		printf("--- Normal Modes ---------------------------- \n");
		for(imode = 0; imode < nmode_freq.size(); imode++) {
			homegaoverkt = (nmode_freq[imode]) * PhysicsConst::planck / (PhysicsConst::kB * temp);
			_modeentropy =
				homegaoverkt/(exp(homegaoverkt) - 1.0) - log(1.0 - exp(-homegaoverkt));
			_vibentropy += _modeentropy;
			printf("Mode %4d  % 10.8e (Hz)   %10.3lf cm-1    %lf kcal/mol\n",
				(int) imode, (double) (nmode_freq[imode]),
				(double) (nmode_freq[imode] / (PhysicsConst::clight * 100)), (double) (-_modeentropy * (PhysicsConst::J2kcal * PhysicsConst::kB * PhysicsConst::Na * temp)));
		}
		_vibentropy *= PhysicsConst::kB;
		printf("--- Quantum Entropy ----------------------- \n");
		printf("Temperature:   %8.3lf  K\n",temp);
		printf("ZeroEnergy:    % 10.5lf kcal/mol\n",
			 getWSpace().ene.epot * PhysicsConst::J2kcal * PhysicsConst::Na);
		printf("Entropy:       % 10.5lf kcal/mol/K \n",
			 _vibentropy * PhysicsConst::J2kcal * PhysicsConst::Na);
		printf("Entropy(-TS):  % 10.5lf kcal/mol \n",
			 -_vibentropy * temp * PhysicsConst::J2kcal * PhysicsConst::Na);

		return -_vibentropy * temp * PhysicsConst::J2kcal * PhysicsConst::Na;
	}

	double EigenSystem::calcVibEntropyClassically(
		double temp, 
		size_t skipFirst
	){
		int imode;
		double _vibentropy = 0;
		double _modeentropy;


		printf("--- Normal Modes ---------------------------- \n");
		for(imode = skipFirst; imode < nmode_freq.size(); imode++) {
			_modeentropy = log(1.0 / ((double)(nmode_freq[imode]) * PhysicsConst::planck / (PhysicsConst::kB * temp)));	// classical partition function is just kT/hv
			_vibentropy += _modeentropy;
			printf("Mode %4d  % 10.8e (Hz)   %10.3lf cm-1    %lf kcal/mol\n",
				(int) imode, (double) (nmode_freq[imode]),
				(double) (nmode_freq[imode] / (PhysicsConst::clight * 100)), (double) (-_modeentropy * (PhysicsConst::J2kcal * PhysicsConst::kB * PhysicsConst::Na * temp)));
		}
		_vibentropy *= PhysicsConst::kB;		// S = kB*ln(Z)
		printf("--- Classical Entropy ----------------------- \n");
		printf("Temperature:   %8.3lf  K\n",temp);
		printf("ZeroEnergy:    % 10.5lf kcal/mol\n",
			 getWSpace().ene.epot * PhysicsConst::J2kcal * PhysicsConst::Na);
		printf("ln(qvib):      % 10.7e \n",_vibentropy /PhysicsConst::kB);
	
		printf("Entropy:       % 10.5lf kcal/mol/K \n",
			 _vibentropy * PhysicsConst::J2kcal * PhysicsConst::Na);
		printf("Entropy(-TS):  % 10.5lf kcal/mol \n",
			 -_vibentropy * temp * PhysicsConst::J2kcal * PhysicsConst::Na);

		//printf("Epo and Entropy   % 10.5lf  % 10.5lf \n",
		//	 getWSpace().ene.epot * PhysicsConst::J2kcal * PhysicsConst::Na,
		//	 -_vibentropy * temp * PhysicsConst::J2kcal * PhysicsConst::Na);

		return -_vibentropy * temp * PhysicsConst::J2kcal * PhysicsConst::Na;
	}

	void EigenSystem::applyAllEigenModes(int firstn){
		if((firstn < 0)||(firstn>(getWSpace().atom.size()*3))) firstn = (getWSpace().atom.size()*3);
		for(int i=0;i < firstn;i++){
			applyEigenMode(i);
		}		
	}

	void EigenSystem::applyEigenMode(int e){
		// create a trajectory displaying the normal mode motions :D
		argcheck_range("eigenvector index",e,0,((int)getWSpace().atom.size() * 3)-1);

		int i, v;
		double r;
		for(i = 0; i < getWSpace().atom.size(); i++)
			getWSpace().old.atom[i].p = getWSpace().cur.atom[i].p;

		for(r = 0; r < (Maths::MathConst::PI * 2.0); r += 0.80) {
			for(v = 0; v < getWSpace().atom.size() * 3; v++) {
				getWSpace().cur.atom[v / 3].p.x = (double)(
					getWSpace().old.atom[v / 3].p.x +
					eigenvector.data_[v * (int)getWSpace().atom.size() * 3 + e] * getWSpace().atom.size() * .01 * sin(r));
				v++;
				getWSpace().cur.atom[v / 3].p.y = (double)(
					getWSpace().old.atom[v / 3].p.y +
					eigenvector.data_[v * (int)getWSpace().atom.size() * 3 + e] * getWSpace().atom.size() * .01 * sin(r));
				v++;
				getWSpace().cur.atom[v / 3].p.z = (double)(
					getWSpace().old.atom[v / 3].p.z +
					eigenvector.data_[v * (int)getWSpace().atom.size() * 3 + e] * getWSpace().atom.size() * .01 * sin(r));
			}
			getWSpace().outtra.append();
		}

		for(i = 0; i < getWSpace().atom.size(); i++)
			getWSpace().cur.atom[i].p = getWSpace().old.atom[i].p;

	}




  //-------------------------------------------------------------------
	//
	//  EigenSystem_Hessian


	void EigenSystem_Hessian::calcEigenVectors(
										Hessian &hessian,
										Verbosity::Type OutputLevel)
	{
		//proxies
		WorkSpace &wspace = getWSpace();

		int i, v, e;
		int iv, iw;


		diagonalise(*(hessian.matrix), OutputLevel);

		if(OutputLevel) printf("Normalising vectors ... (%d,%d)\n",
			eigenvector.dim1(),
			eigenvector.dim2());
		// normalise vectors
		for(v = 0; v < eigenvector.dim1(); v++) {
			double sumsq = 0;
			for(e = 0; e < eigenvector.dim2(); e++) {
				sumsq += sqr(eigenvector.data_[e * eigenvector.dim1() + v]);
			}
			sumsq = sqrt(sumsq);
			for(e = 0; e < eigenvector.dim2(); e++)
				eigenvector.data_[e * eigenvector.dim1() + v] /= (float)(sumsq);
		}
	

		// Now filter out those modes which correspond purely to translational
		// & rotational modes - do this by analysising the pattern the normal mode
		// vectors make, i.e. do they all lie in a plane (after removel of total translation)
		// and does the normal distance from the rotation normal correlate with
		// the atoms's distance from the centre of geometry
		// score them according to these criteria and skimm off the top scoring
		// 6 eigenmodes (by marking them)

		int      *rottransmode = new int[ eigenvector.dim2()];

		if(RemoveTraRot) {
			if(OutputLevel) printf("Removing translational & rotational modes ... \n");

			for(i = 0; i < getWSpace().atom.size() * 3; i++)
				rottransmode[i] = 0;

			dvector *eigenvec = new dvector [getWSpace().atom.size()];
			double   *nmode_score = new double[30];
			int      *nmode_score_index = new int[30];
			int nmode_score_cnt = 0;
			double   *nvdist = new double[getWSpace().atom.size()];
			double   *nvmag = new double[getWSpace().atom.size()];


			// check first 30 modes (which are ordered by size of their eigenvalues,
			// i.e. only check the 30 with the lowest eigenfrequencies)
			for(v = 0; v < (Maths::min((int)getWSpace().atom.size(), 10) * 3); v++) {

				// save in a dvector  array (one 'dvector ' per atom now);
				for(e = 0; e < getWSpace().atom.size() * 3; e++) {
					switch (e % 3) {
						case 0:
							eigenvec[e / 3].x = eigenvector.data_[e * (int)getWSpace().atom.size() * 3 + v];
						case 1:
							eigenvec[e / 3].y = eigenvector.data_[e * (int)getWSpace().atom.size() * 3 + v];
						case 2:
							eigenvec[e / 3].z = eigenvector.data_[e * (int)getWSpace().atom.size() * 3 + v];
					}
				}

				// now for each triplet, first calculate the angle to each other triplet
				// if the average angle is close to 0, it's an indication that it is a purely translational dvector 


				dvector netmove(0, 0, 0);
				dvector cog(0, 0, 0);
				dvector naverage(0, 0, 0);
				double dotsum = 0;
				double nva, nvb, nvR;
				dvector ev;

				dvector cog_single;
				cog_single = getWSpace().getCentreOfGeometry();
				cog.setTo(cog_single);

				// remove netmovement
				for(iv = 0; iv < getWSpace().atom.size(); iv++)
					netmove.add(eigenvec[iv]);
				netmove.div((double)getWSpace().atom.size());

				netmove.info();

				for(iv = 0; iv < getWSpace().atom.size(); iv++)
					eigenvec[iv].sub(netmove);

				for(iv = 0; iv < getWSpace().atom.size(); iv++) {
					for(iw = iv + 1; iw < getWSpace().atom.size(); iw++) {
						ev.crossProduct(eigenvec[iv], eigenvec[iw]);
						
						if((iv == 0) && (iw == 0))
							naverage.setTo(ev);
						else {
							double dot = dotProduct(ev, naverage);
							if(dot < 0)
								ev.mul(-1);
							naverage.add(ev);
						}
					}
				}
				naverage.info();
				naverage.unify();

				for(iv = 0; iv < getWSpace().atom.size(); iv++) {
					dotsum += sqr(dotProduct(eigenvec[iv], naverage));	// calculate how well the average normal fits the data

					ev.setTo(getWSpace().cur.atom[iv].p);
					ev.sub(cog);		// calc coordinate relative to CM

					// nvdist[iv] is the distance perpendicular to average normal
					nvdist[iv] = sqrt(ev.innerdot() - sqr(naverage.x * ev.x + naverage.y * ev.y + naverage.z * ev.z));
					nvmag[iv] = eigenvec[iv].mag();
					//printf("-->%lf %lf          \t      %lf\n",xj.mag(),nvdist[iv],nvmag[iv]);
				}

				leastSquaresFit(nvdist, nvmag, getWSpace().atom.size(), nva, nvb, nvR);

				printf("Anglesum: %lf   %e   %lf    %lf    %lf \n", 
					1 / (dotsum / sqr(nvR)), dotsum, nva, nvb, nvR);

				nmode_score[nmode_score_cnt] =  dotsum *fabs((double) eigenvalue_real[v]) / (sqr(nvR)) ;
				nmode_score_index[nmode_score_cnt] = v;
				nmode_score_cnt++;
			}

			printf("Sorting modes \n");
			// sort them by score and mark the ones analysed to be rot. modes
			qcksort(nmode_score, nmode_score_index, nmode_score_cnt);
			printf("deleteing memory \n");
			for(v = 0; v < 6; v++) {
				printf("%d \n", nmode_score_index[v]);
				rottransmode[nmode_score_index[v]] = 1;
			}

			delete[]nmode_score;
			delete[]nmode_score_index;
			delete[]nvdist;
			delete[]nvmag;
			delete[]eigenvec;

		}
		// print final eigenvector matrix:
		/*
		for(int e=0;e<getWSpace().atom.size()*3;e++){
		for(int v=0;v<getWSpace().atom.size()*3;v++){
		printf("% 10.8e\t", eigenvector.data_[e*getWSpace().atom.size()*3 + v] );
		}
		printf("\n");
		}
		*/

		// save normal mode frequencies  (these will only be correct if the matrix
		// provided was mass weighted !!

		if(OutputLevel) printf("calculating normal mode frequencies ... \n");

		nmode_freq.clear();
		for(e = 0; e < eigenvalue_real.dim1(); e++) {
			if((RemoveTraRot) && (rottransmode[e] != 0))
				continue;		// ignore purely rotational/translational modes
			nmode_freq.push_back(sqrt(fabs((double) eigenvalue_real[e] / 
													(4.0 * sqr(Maths::MathConst::PI)))));
		}

		// save the 0-point energy
		E0 = getWSpace().ene.epot;
		delete[]rottransmode;	
	}


	void EigenSystem_Hessian::calcHessian(Hessian &hessian){
		hessian = Hessian( getWSpace () );

		int m = eigenvector.dim1();
		int n = eigenvector.dim2();
		TNT::Array2D< LINALG_real > temp(m,n);
		int matcnt=0,i,j;
		for(i=0;i<m;i++){
			for(j=0;j<n;j++){
				temp.data_[j*m + i] = eigenvector.data_[matcnt] * eigenvalue_real[j];
				matcnt++;
			}
		}

		*hessian.matrix = TNT::matmult(eigenvector,temp);
	}

	double EigenSystem_Hessian::getE0()
	{
		return E0*Physics::PhysicsConst::Na * Physics::PhysicsConst::J2kcal; 
	}











  //-------------------------------------------------------------------
	//
	//  EigenSystem_Covariance


	void EigenSystem_Covariance::calcEigenVectors(
		CovarianceMatrix &covmat,
		Verbosity::Type OutputLevel
	)
	{
		int e;

		WorkSpace &wspace = getWSpace();

		diagonalise(*(covmat.matrix), OutputLevel);

		nmode_freq.clear();
		for(e = 0; e < getWSpace().atom.size() * 3; e++) {
			//if((RemoveTraRot) && (rottransmode[e] != 0))
			//	continue;		// ignore purely rotational/translational modes
			printf("%f \n", eigenvalue_real[e]);
			nmode_freq.push_back(sqrt(PhysicsConst::kB * Temperature / (fabs(eigenvalue_real[e]) * 4.0 * sqr(Maths::MathConst::PI))));
		}

	}




	void EigenSystem_Covariance::calcHessian(
		Hessian &hessian,
		double Temperature
	)
	{
		hessian = Hessian( getWSpace() );

		printf("calculating Hessian from Eigenvectors and Eigenvalues \n");
		int m = eigenvector.dim1();
		int n = eigenvector.dim2();
		TNT::Array2D< LINALG_real > temp(m,n);
		int matcnt=0,i,j;
		for(i=0;i<m;i++){
			for(j=0;j<n;j++){
				// note that the eigen values from a covariance matrix are not the same as the
				// ones from a hessian matrix - they are related by kBT/lambda
				temp.data_[j*m + i] = (float)(eigenvector.data_[matcnt] * 
															(PhysicsConst::kB * Temperature / eigenvalue_real[j]));
				matcnt++;
			}
		}

		*(hessian.matrix) = TNT::matmult(eigenvector,temp);
	}


}

