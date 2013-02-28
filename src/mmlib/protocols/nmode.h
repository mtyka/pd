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

#ifndef __NMODE_H
#define __NMODE_H

// Essential Headers
#include "maths/tntjama/tnt.h" // Provides a class member
#include "protocols/protocolbase.h" // Provides a base class
#include "workspace/componentbase.h" // Provides a base class
#include "workspace/workspace.fwd.h"


// float or double settign for linear algebra/matrix operations. Double is default and works for
// most purposes, but the define allows for an easy switch to float should it be required.

#ifndef LINALG_real
  #define LINALG_real float  
#endif

namespace Protocol
{
	//-------------------------------------------------
	//
	/// \brief Abstract base class for various matrices of degrees of freedom of a
	///        workspace, such as hessians and covariance matrices
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	///
	/// \author Mike Tyka 
	///

	class PD_API Matrix3Nx3N {
	public:
			Matrix3Nx3N()
		  {
			  matrix = NULL;
		  };

		  virtual ~Matrix3Nx3N(){
			  delete matrix;
		  };


			void createEmpty(size_t m, size_t n){
			  delete matrix;
			  matrix = new TNT::Array2D < LINALG_real > (m, n);
		  }


		  /// print the matrix in Mathematica Format
		  void printMathematica();					

		  /// \brief 
		  /// This, although data, is public so that other classes can efficiently
		  /// access it (and modify it if they need to) 
		  TNT::Array2D < LINALG_real  > *matrix;

	};


	//-------------------------------------------------
	//
	/// \brief  Class For Second derivative Matrices of a forcefield
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
	class PD_API Hessian: public Matrix3Nx3N{
	public:
			Hessian(): Matrix3Nx3N() {}
			Hessian( WorkSpace &wspace ): Matrix3Nx3N() {
				 createEmpty(wspace);
			}

		  virtual ~Hessian(){
		  }


		  /// mass weight the matrix (for normal mode calcs)
		  void massWeight( WorkSpace &wspace );							 

		  /// two more properties that are not technically part of the hessian itself
		  /// but nevertheless often useful:
			
			/// the potential energy offset at x==x0
		  double             e0;       

		  /// x0, the phase space position at which the matrix was 'taken'
		  SnapShot    x0;        

			void createEmpty(WorkSpace &wspace){
			  delete matrix;
			  matrix = new TNT::Array2D < LINALG_real > (wspace.atom.size() * 3, wspace.atom.size() * 3);
			  x0 = wspace.save();
			  e0 = wspace.ene.epot;
		  }


		  /// writes raw data to a file
		  void writeRaw(const std::string &filename);  

		  /// reads raw data back froma a file
		  void readRaw(const std::string &filename);   
	};


	//-------------------------------------------------
	//
	/// \brief  calculates a Hessian Matrix using numerical method. 
	///
	/// \details 
	/// Hessian_Numerical will take each atom in turn move it a small amount (StepSize)
	/// in each coordinate direction and calculate the change in the force.
	/// and thus calculate the cartesian second derivate
	/// matrix. This is very costly (taking at least ~3*natom energy evaluations)
	/// but unversal - i.e. the forcefields themselves dont need to supply a second
	/// derivative calculation.
	/// The calculation involves moving each degree of freedom in turn a small amount (StepSize)
	/// or multiples of that amount and obtain the gradient in the forces with respect to that degree of freedom.
	/// Depending on the forcefield used, a single step can be sufficient (FiniteDiff=1) or sometimes multiple steps backwards and forwards
	/// (controlled by the parameter FiniteDiffSteps >= 2) may be necessary (due to numerical inaccuracies).
	/// in the latter case a linear fit is used to the multiple data points.
	/// ForceSymmetry averages the two half matrices obtained to ensure that the Hessian is positive definite
	// /(as is should be in theory, but may not be exactly due to numerical inaccuracies).
	///
	/// \author Mike Tyka 
	///
	class PD_API Hessian_Numerical: public Hessian, public Protocol::ProtocolBase{
	public:
		Hessian_Numerical(WorkSpace &_wspace, Physics::Forcefield & _ff):
		  ProtocolBase(_ff),
			  Hessian()
		  {	
			  Steps = 1;
			  ForceSymmetry=true;   
			  StepSize = 0.0002;    // this is usually accurate enough
				FiniteDiffSteps = 2;
		  };

		  virtual ~Hessian_Numerical(){};
		  virtual Hessian_Numerical* clone() const 
		  { 
			  return new Hessian_Numerical(*this); 
		  }

			/// Run the numerical calculation
		  virtual int runcore();   

		  // Parameters

			/// force symmetry? (by averaging Hij and Hji elements) - true by default
			bool   ForceSymmetry;				

		  /// the Step size for the numerical calculation; 0.0002 Angstrom by default
		  double StepSize;	        

			/// How many finite difference measurements? 1 = 0, StepSize; 2= -StepSize, 0, StepSize; 3 = -StepSize, 0, d, 2*StepSize; etc..
			int    FiniteDiffSteps;   
	protected:
		// discarded for now.
		virtual void info() const{};						
		virtual void infoLine() const{};				
		virtual void infoLineHeader() const{};	

	};


	//-------------------------------------------------
	//
	/// \brief  is a marriage between a Hessian and a Monitor and calculates the Covariance matrix of a run 
	///
	/// \details CovarianceMatrix derives from:
	///            Hessian - Its not technically a Hessian but it has the same machinery so here this is a useful derivation.
	///            MonitorBase - It accumulates data through a run or a trajectory replay
	///            WorkSpaceOperatorBase - CovarianceMatrix is permanaetly attached to a workspace, so 
	///                                    WorkSpaceOperatorBase provides the machinery for that.
	///
	/// \author Mike Tyka  
	///
	///
	class PD_API CovarianceMatrix: public Hessian, 
																 public Monitors::MonitorBase,
																 public WorkSpaceOperatorBase {
	public:
		CovarianceMatrix(WorkSpace & _wspace):
				WorkSpaceOperatorBase( _wspace ),
			  Hessian(),
			  MonitorBase()
		  {	
			  reset();
			  state = Gathering;
		  };

		  virtual ~CovarianceMatrix(){};

		  /// resets the state to 0 (i.e. clears all the sums & sets state to Gathering
		  virtual void		reset();  

		  /// adds another data point
		  virtual void		setcurdata();

		  /// finishes the calculation by taking the appropriate averages etc.
		  void    finish(bool mass_weight=true,
			bool remove_correlations=false);
		  void massWeight();
	private:
		enum _state { Gathering, Finished } state;
		int							ndata;
	};


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
	class PD_API EigenSystem: public WorkSpaceOperatorBase{
	public:
		EigenSystem( WorkSpace &wspace ): WorkSpaceOperatorBase( wspace )
		{
		};

		double  calcVibEntropy(double temp);
		double  calcVibEntropyFull(double temp);
		double  calcVibEntropyClassically(double temp, size_t skipFirst=0);

		void applyAllEigenModes(int firstn=-1);
		void applyEigenMode(int e);
		void printEigenValues();

	protected: 

		void diagonalise(TNT::Array2D < LINALG_real > &matrix,
			Verbosity::Type OutputLevel = Verbosity::Silent);


		/// This holds the eigenvalue's real parts
		TNT::Array1D < LINALG_real >eigenvalue_real;     

		/// And this the imaginary parts
		TNT::Array1D < LINALG_real >eigenvalue_imag;     

		/// Eigenvector matrix (vectors in columns!)
		TNT::Array2D < LINALG_real >eigenvector;         


		/// Eigen frequencies (only meaningful when originating matrix was massweighted) 
		std::vector <double> nmode_freq;          
	};


	//-------------------------------------------------
	//
	/// \brief  diagonalize a Hessian matrix
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
	class PD_API EigenSystem_Hessian: public EigenSystem {
	public:
		EigenSystem_Hessian( WorkSpace &wspace): EigenSystem( wspace ) 
		{
			RemoveTraRot = false;
		}
		
	  // Functions
		virtual void calcEigenVectors(Hessian &hessian,
				                    Verbosity::Type OutputLevel=Verbosity::Silent);
		virtual void calcEigenValues(Hessian &hessian, 
	    			                Verbosity::Type OutputLevel=Verbosity::Silent)
		{
			calcEigenVectors(hessian,OutputLevel);
		}
					
		/// remove the 6 Degrees of rotational/translational freedom ?
		bool  RemoveTraRot;  
		
		/// back calculate what the hessian had been from the normal vectors and eigenvalues
		void calcHessian(Hessian &hessian );         
		double getE0();
	private:
		double E0;
		double vib;
	};


	//-------------------------------------------------
	//
	/// \brief  diagonalize a Covariance matrix
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
	class PD_API EigenSystem_Covariance: public EigenSystem {
	public:
		EigenSystem_Covariance( WorkSpace &wspace): EigenSystem( wspace ) {};

		// Functions
		virtual void calcEigenVectors(CovarianceMatrix &cov,
				                    Verbosity::Type OutputLevel=Verbosity::Silent);
		virtual void calcEigenValues(CovarianceMatrix &cov,
				    			          Verbosity::Type OutputLevel=Verbosity::Silent)
		{
			calcEigenVectors(cov,OutputLevel);
		}

		void calcHessian(Hessian &hessian,
			double Temperature );         

		double Temperature;
	private:
	};

} // namespace Protocol


#endif
