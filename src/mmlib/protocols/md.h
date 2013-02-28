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

#ifndef __MD_H
#define __MD_H

// Essential Headers
#include "protocols/protocolbase.h" // Provides a base class
#include "protocols/temperature.h"  // Provides a class member
#include "workspace/workspace.fwd.h"

namespace Protocol{

	class PD_API LocalREMD;
	class PD_API MolecularDynamics;

	/// \brief Molecular Dynamics Algorithms
	/// \details 
	/// This class implements a number of basic Molecular Dynamics algorithms
	/// including thermostating, barstating etc.
	/// The algorithms include:
	///
	///  Newtonian Molecular Dynamics, using the integrators Verlet, Velocity Verlet and Beeman
	///  Thermostats implemented: Berendsen and Andersen 
	///  Barostats implemented: Berendsen 
	///  Langevin Dynamics
	///
	/// 
	///
	/// Molecular Dynamics
	/// M. Levitt and H. Meirovitch, Integrating the Equations of Motion,
	/// J. Mol. Biol., 168, 617-620 (1983)
	///
	/// MolecularDynamics::Beeman Integrator
	/// D. MolecularDynamics::Beeman, Some Multistep Methods for Use in Molecular Dynamics
	/// Calculations, J. Comput. Phys.,20,130-139 (1976)
	///
	/// Velocity Verlet Integration
	/// W.C. Swope, H.C. Andersen, P.H. Berens, K.R. Wilson. A computer
	/// simulation method for the calculation of equilibrium constants
	/// for the formation of physical clusters of molecules: application
	/// to small water clusters. J. Chem. Phys. 76, 637-649 (1982)
	///
	/// Langevin Dynamics:
	/// Gunsteren, W.F. van, MolecularDynamics::Berendsen H.J.C. Algorithms for brownian
	/// dynamics Mol. Phys., 45 (3) 637-647 (1982)
	/// Paterlini M.G., Ferguson D.M. Constant Temperature simulations
	/// using the Langevin equation with velocity Verlet integration
	/// Chem. Phys. 236, 243-252 (1998)
	///
	/// Thermostats: MolecularDynamics::Berendsen
	/// H. J. C. MolecularDynamics::Berendsen, J. P. M. Postma, W. F. van Gunsteren, A. DiNola and
	/// J. R. Haak, Molecular Dynamics with Coupling to an External Bath,
	/// J. Chem. Phys., 81, 3684-3690 (1984)
	///
	/// Thermostats: Andersen
	/// Andersen H.C. Molecular dynmics simulation at constant pressure and/or
	/// Temperature. J. Chem. Phys. 72(4) 2384 (1980)
	///
	/// \author Michael Tyka

	class PD_API MolecularDynamics : public RangedProtocolBase
	{
	public: 

		friend class PD_API LocalREMD;
		MolecularDynamics( Physics::Forcefield& _ff );
		MolecularDynamics( Physics::Forcefield& _ff, const PickAtomRange& _def );
		virtual ~MolecularDynamics();

		virtual MolecularDynamics* clone() const 
		{ 
			return new MolecularDynamics(*this); 
		}

		/// prints a little block of parameter information
		virtual void info() const; 
		
		/// usual run mode
		virtual int  runcore(); 

		enum ThermostatType { NoThermostat, Andersen, Berendsen};
		enum BarostatType   { NoBarostat, BerendsenBaro };
		enum IntegratorType { Verlet, VelocityVerlet, Beeman, Langevin };

		/// type of MD Integrator algorithm {Verlet|VelocityVerlet|Beeman|Langevin}
		IntegratorType Integrator;     
		
		/// integration Timestep in (seconds)
		double Timestep;							 
		
		/// randomize velocities at the start of the simulation ?
		bool RandVel;									 
		
public:

		/// How often to remove total translational and rotational momentum
		int UpdateRemoveTotalMomentum;		
		
		/// Type of Thermostat {NoThermostat|Andersen|Berendsen}
		ThermostatType Thermostat;		 
		
		/// Target Pressure (only if pressure thermostat is active) [Atm-1]
		double TargetPressure;				 

		/// Type of barostat {NoThermostat|Berendsen}
		BarostatType Barostat;		

		/// Bulk solvent isothermal compressibility [Atm-1], ( default is water, 0.000046)
		/// The default is appropriate for water
		double	Compress;

		/// tau parameter (in seconds) (for Berendsen Thermostat)
		double BerendsenTau;					 
		
		/// tau parameter (in seconds) (for Berendsen Barostat)
		double BerendsenPressureTau;	 
		
		/// collision rate (~0.1) for andersen Thermostat
		double AndersenRate;					 
		
		/// friction coefficient for langevin dynamics [ per second ] ( default = 10E12 ( which is 10.0 per picosecond) )
		double FricCoeff;							

		/// Apply langevin dynamics to hydrogen atoms ?
		bool LangevinOnHydrogens;							

		/// Move system after each step such that center of mass is at 0,0,0. This is 
		/// only recommended in Implicit solvent or other Infinite Space simulations,
		/// not in periodic boundary conditions
		bool CentreAfterMove;

	protected:

		/// sets default parameter values
		void settodefault();								

		/// prints a line of current energies
		virtual void infoLine() const;			

		/// prints the headers for the above function
		virtual void infoLineHeader() const;

		/// prints the Final set of statistics after run 
		virtual void printFinalStatistics(size_t steps);

	private:
		int							m_StepOffset;

		/// Current Target Kinetic Energy 
		double					m_TargetEkin;   

		/// current Temperature (calculated by calcKineticEnergy)
		double					m_CurTemp;			

		/// current Pressure (calculated by calcKineticEnergy)
		double					m_CurPress;			

		/// current Step number
		int							Step;						


		/// Potential energy statistics
		Maths::Statistics m_StatEpot;	

		/// Total energy statistics
		Maths::Statistics m_StatEtot;	

		/// temperature statistics
		Maths::Statistics m_StatTemp;	

		/// pressure statistics
		Maths::Statistics m_StatPres;	

		/// density statistics
		Maths::Statistics m_StatDens;	

		/// volume statistics
		Maths::Statistics m_StatVolu;	
		std::vector<double> m_TotalEnergyDrift_x;
		std::vector<double> m_TotalEnergyDrift_y;
																		 
		Maths::dvector linvel;
		Maths::dvector angvel;
		Maths::dvector linmom;
		Maths::dvector angmom;

		void setup(); 
		int run_core();

		//Statistical calculations
		void calcTotalLinearVelocity(int istart, int iend);
		void calcTotalAngularVelocity(int istart, int iend);
		void calcTotalLinearMomentum(int istart, int iend);
		void calcTotalAngularMomentum(int istart, int iend);
		void zeroTotalLinearMomentum(int istart, int iend);
		void zeroTotalAngularMomentum(int istart, int iend);



		/// \brief Calculate Kinetic Energy, Instantaneous Temperature and Pressure
		/// \details
		///
		///  The kinetic energy is calculated from the velocity v of each atom i
		///  using K(i) = 0.5 * m(i) * v(i) ^ 2
		///  where m(i) is atom i's mass and v(i) is it's velocity
		///  The kinetic energy is stored in the workspace's Energy structure wspace->ene.ekin
		///  
		///  Then, it calculates the instantaneuos temperature T from the kinetic energy K = 2NkBT/2
		///
		///  Finally, it calculates the instantaneous pressure from P = (2.0K + W)/(3V)
		///
		///  Where W is the internal virial which is basically the derivative of the internal 
		///  energy with respect to the volume, i.e. dH/dV
		///
		///  dH/dV = W/(3.0*V)
		///  For pairwise interactions W is: (see virial calculation in the forcefield components)
		///  W = Sum[ Dist * |Force|  ] 
		///
		///  The first term  (2.0K/3V) is basically the ideal gas contribution PV = nRT.
		
		void calcKineticEnergy();
		
		
		void setKineticEnergyTo(double Edes);
		void setInitialSpeeds();
		void setInitialSpeeds(double tgtTemp);
		void LocalREMD();
		void LocalREMD(double temperature);

		void applyThermostat();
		void applyBerendsenThermostat(double TargetTemperature);
		void applyAndersonThermostat(double TargetTemperature);

		void applyBarostat();
		void applyBerendsenBarostat();

		void equaliseVelocities();

		//Particle Protocols
		void calcOldPositions();
		void reCalculateNewPositions();

		void applyForces();
		void applyForces_BeemanIntegration();
		void applyForces_VelocityVerletIntegration();
		void applyForces_VerletIntegration();
		void applyForces_LangevinIntegration();
		void applyForces_LangevinIntegration_NonHydrogens();
	};
}

#endif



