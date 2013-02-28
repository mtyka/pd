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

#ifndef __FORCEFIELDBASE_H
#define __FORCEFIELDBASE_H

#include <typeinfo>
#include <vector>

#include "verbosity.h"
#include "tools/cloneholder.h"
#include "pickers/pickbase.h"          // provides a member variable
#include "workspace/pospointer.h"   // provides a member variable
#include "workspace/componentbase.h"   // provides a member variable
#include "workspace/workspace.fwd.h"
#include "workspace/snapshot.fwd.h"
#include "forcefields/forcefield.fwd.h"

namespace Protocol
{
	class ProtocolBase;
}

class Object;
class FFParamSet;
class PickBase;

namespace Physics
{
	//-------------------------------------------------
	//
	/// \brief  A container for all the general Physics Constants
	///
	/// \details
	/// Various Physical Constants - Description & [Units] given for each. 
	/// These are from The NIST Reference on Constants, Units and Uncertainty (2006)
	/// http://physics.nist.gov/constants
	///
	/// \author Mike Tyka & Jon Rea 
	///
	class PD_API PhysicsConst
	{
	private:
		PhysicsConst();
	public:
		static const double kcal2J;          ///< Multiplier to convert Joules to kcal [kcal/J]
		static const double J2kcal;          ///< Multiplier to convert kcals to Joule [J/kcal]
		static const double Na;              ///< Avogadros's Number - i.e. no of particles per mol [1/mol]
		static const double kcal2JDivNa;     ///< Conversion Factor from kcal/mol to Joules/molecule
		static const double kB;              ///< Boltzmann Constant [J/K]
		static const double R_kcal;          ///< The Gas Constant in kcal/mol/K
		static const double planck;          ///< Plank Constant [Js]
		static const double clight;          ///< Speed of Light [m/s]
		static const double amu;             ///< Atomic mass unit (mass of Carbon-12 divided by 12) [kg]
		static const double Angstrom;        ///< Length of an Angstrom in meters [m/A]
		static const double invAngstrom;     ///< Number of A in a meter [A/m]
		static const double e_charge;        ///< Charge of electron - Unit of Charge [C]
		static const double sqr_e_charge;    ///< Square of electronic charge [C^2]
		static const double e0;              ///< Permittivity of Space [C^2/J m]
		static const double _4pi_e0;         ///< Electrostatic  Units (premultiplier for electrostatic  interactions) [Jm/C^2]
		static const double econv;           ///< Energy of two unit charges at 1 Angstrom distance [kcal/mol]
		static const double halfeconv;       ///< half of the above
		static const double econv_joule;     ///< as above but in in joules
		static const double halfeconv_joule; ///< and half of that
		static const double Bar2Pa;          ///< Pressure: Bar to Pascal

		/// print a list of fundamental constants and their names, units and meaning
		static void printConstants();
	};


	/// \class ForcefieldBase
	/// \brief Base class to all forcefield components - defines common interface.
	/// \author M.Tyka
	///
	/// A Forcefield is a functor class that implements an energy and force 
	/// calculation, i.e. a forcefield component. It’s interface is defined 
	/// using ForcefieldBase which defines a number of pure virtual functions 
	/// which must be overloaded by derived classes. Most importantly calcEnergies() 
	/// and calcForces() calculate energies and forces respectively (the latter 
	/// does both). A container version of ForcefieldBase, Forcefield, acts as a
	/// collector of individual forcefield components to make a complete forcefield.
	/// The collector mirrors the functions described above and simply calls the 
	/// equivalents in each of the stored forcefield components. By adding components
	/// into the container, various forecfields can be created from the user interface
	/// giving the user full control and flexibility over the simulation parameters.
	/// Forcefield components can also be made ’Passive’ in which energies are 
	/// calculated but not added to the total energy to be used in Free Energy
	/// MoveBase calculations. The forcefield parameters are supplied through 
	/// a auxillary class called FFParam which in turn reads in force field defintion
	/// files containing all the parameters. A new format has been developed which 
	/// was mainly inspired by the CHARMM format, but has several improvements
	/// inspired by other formats. FFParam however can also read in CHARMM .prm 
	/// and .top files for compatibility.  

	class PD_API ForcefieldBase
		: public Object,
		public WorkSpaceOperatorBase
	{
		friend class PD_API Forcefield;
	public:

		enum AtomicVerbosity
		{ 
			Summary, 
			Atomwise, 
			Detailed 
		};

		ForcefieldBase( WorkSpace &newwspace )
			: WorkSpaceOperatorBase( newwspace )
		{
			setDefaults();
		}

		/// virtual destructor
		virtual ~ForcefieldBase(){}; 

		virtual ForcefieldBase* clone() const = 0;

private:
		void setDefaults()
		{
			name = "Forcefield";
			ShortName = "FF";
			Active = true;
			Passive = false;
			OutputLevel = Verbosity::Normal;
			needsetup = true;
			m_CheckSum = 0;
		}
public:

		// setup functions -------------------------------------------------
		
		/// \brief This virtual function will be called automatically before any
		/// simulation. Derived classes can use it for setting up things before 
		/// the actual simulation begins, for example custom arrays over all the atoms,
		/// load simulation parameters or anyhting else.
		/// The function *MUST* be safe to be run multiple times! 
		virtual void setup() = 0;

		/// checks if the forcefield still 'points' to the wspace it was 
		/// last setup to run on. If not it will call setup() to ensure the setup 
		/// is up to date!
		virtual int ensuresetup(WorkSpace &newwspace);

		/// This is a special checksum function that returns the state of essential
		/// variables who's change must trigger a re-setup (i.e. calling of setup)
		/// or in other words setting of needsetup to true;
		/// see examples in other forcefields to see how to implement this function.
		/// Basically you want to return a unique number dependant on the special
		/// variables you want "checked"
		virtual unsigned long calcCheckSum() { return getWSpace().getCheckSum(); } 

	private:
		unsigned long m_CheckSum;
	public:
		// Energy/Force functions -------------------------------------------------
		
		/// calculate forces and energies (used during simulation)
		/// derived classes must define this and implement the force calculation
		virtual void calcForces() = 0; 

		/// calculate energies only (for monte carlo for example)
		/// By default this calls calcForces() but if calculation of forces
		/// creates a significant overhead then a specialised implementation
		/// can be useful.
		virtual void calcEnergies();   

		///display energies verbosely
		virtual void calcEnergiesVerbose(AtomicVerbosity level); 
	                                       
		// inspector functions -------------------------------------------------

		/// prints a little block of parameter information               		
		virtual void info() const;           

		/// prints a line of current energies                       		
		virtual void infoLine() const;      

		/// prints the headers for the above function                       		
		virtual void infoLineHeader() const; 

		void activate();
		void deactivate();

		/// Short name appears as captions on columns etc..
		std::string ShortName;

		/// \brief Controls if Energy&Forces are added to respective sums
		/// false by default, a 'Passive' forcefield calculates the energy but does not add
		/// energy or force components to wspace->epot/wspace->atom[].f which
		/// means it has no effect on the behaviour of the system.
		/// This functionality allows forcefields to be used in Free Energy MoveBase
		/// (FEP) situations where the canonical exponential average of a potential energy
		/// difference is obtained.
		bool Passive;

		/// produce running screen output or not. Default is false.
		Verbosity::Type OutputLevel;

		/// potential energy of a forcefield (component) in SI units (J/molecule)
		double epot;

		/// returns a pointer to epot (for setting up FEP calcultions from python)
		double *epot_ptr(){ return &epot; }

		void forceSetup(){ needsetup = true; };

	protected:
		/// \brief (de)activetes forcefield component
		/// true by default, when false the forecfield (component) is deactived
		/// and will not perform any calculation
		bool Active;

		/// Encapsulated forcefields are allowed to call setup on their encapsulation
		//void setupCoreCall( ForcefieldBase& ff, WorkSpace& wspace );

		/// \brief Indicates state of forcefield: if true setup() must be run.
		/// If this is set to true it indicates that the forcefield needs to run setup()
		/// before being allowed to run. This is true when class PD_API is born, and when radical
		/// modifications are done
		bool needsetup;
	};


	////-------------------------------------------------
	///// \brief   Does nothing. Acts as an inert placeholder in situation 
	/////          where a Forcefield is required (e.g. syntactically) but
	/////          no operation is desired.
	/////
	///// \author Mike Tyka 
	/////
	//	class NullForcefield: public Physics::ForcefieldBase
	//	{
	//	 public:
	//		
	//		NullForcefield(){}
	//		virtual void setup(){ return 0; };
	//		virtual void calcForces(){}
	//		virtual void calcEnergies();   
	//		virtual void calcEnergiesVerbose(Verbosity::Type level) {}; 
	//		virtual void infoLine() const {};       
	//		virtual void infoLineHeader() const {}; 
	//	};

}

#ifdef SWIG
%template(ObjectContainer_ForcefieldBase) ObjectContainer<Physics::ForcefieldBase>;
#endif

namespace Physics
{
	//-------------------------------------------------
	//
	/// \brief  This is a special container holding forcefield components to make a forcefield
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	///
	/// \author Mike Tyka & Jon Rea 
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class PD_API Forcefield: 
		public ForcefieldBase, 
		public ObjectContainer<ForcefieldBase> 
	{
	public:
		friend class Protocol::ProtocolBase;

		Forcefield( WorkSpace &newwspace)
			: ForcefieldBase( newwspace )
		{
			name = "Forcefield";
			ShortName = "ff";
			max_force = -1;
		}

		virtual Forcefield* clone() const { return new Forcefield(*this); }

		virtual void setup();

		virtual void calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level);
		virtual void calcEnergies();
		virtual void calcForces();

		virtual void info();
		virtual void infoLine() const;
		virtual void infoLineHeader() const;

		virtual void addPreAddition( ForcefieldBase  &newelement );

		/// calculate energy and display total potential energy
		void printEnergyShort();

		/// calculate energy and display summary
		void printEnergySummary();

		/// calculate energy verbosely and display atomwise energies
		void printEnergyByAtom();

		/// calculate energy verbosely and display interactions individually
		void printEnergyDetailed();


		/// List the runtime types of the forcefields that this forcefield contains.
		void listForcefields();  

	protected:
		int ensuresetup(WorkSpace &newwspace);

	private:
		virtual void info() const {} // "dispose" of this inherited function

		/// if max_force > 0 Forcefield will truncate forces if they exceed max_force
		double max_force; 
	};


#ifndef SWIG

	template< typename T >
	int hasFFComponent( Forcefield& ff )
	{
		try
		{
			int found = 0;
			for( size_t i = 0; i < ff.size(); i++ )
			{
				T* ptr = dynamic_cast<T*>(&ff.element(i)); // Safe downcast
				if( ptr != NULL )
				{
					found++;
				}
			}
			return found;
		}
		catch( std::bad_typeid ex )
		{
			THROW(CodeException,"Invalid type int the obtainFFComponent() downcast call!");
		}
	}

	/// Find me my desired forcefield component at runtime. Returns true ONLY if one instance of the forcefield is found.
	/// If true is returned, the _ObtainedPointer will be assigned to the first instance found in the container.
	/// Returns the number of components of that type Fills the argument pointer with the first forcefield object of the templated type, fills NULL if there is no such component.
	
	template< typename T >
	bool obtainFFComponent( Forcefield& ff, T*& _ObtainedPointer ) 
	{		
		try
		{
			_ObtainedPointer = NULL; // Firstly nullify the pointer to be returned
			size_t found = 0;
			for( size_t i = 0; i < ff.size(); i++ )
			{
				T* ptr = dynamic_cast<T*>(&ff.element(i)); // Safe downcast
				if( ptr != NULL )
				{
					if( found == 0 )
					{
						// We have found the 1st instance
						_ObtainedPointer = ptr;
					}
					found++;
				}
			}
			if( found == 1 )
			{
				return true; // All is well :-D
			}
			else
			{
				_ObtainedPointer = NULL; // If there is 0 instances or more than one, this is an error condition!
				return false;
			}
		}
		catch( std::bad_typeid ex )
		{
			_ObtainedPointer = NULL; // Ensure we dont return anything dubious.
			THROW(CodeException,"Invalid type int the obtainFFComponent() downcast call!");
		}
	}

	/// Pushes all instances of a particular desired forcefield component into the supplied std::vector<T>
	/// Returns the number that have been added.
	/// Returns the number of components of that type Fills the argument pointer with the first 
	/// forcefield object of the templated type, fills NULL if there is no such component.
	template< typename T >
	size_t obtainFFComponents( Forcefield& ff, std::vector<T*>& _Vector ) 
	{		
		try
		{
			size_t found = 0;
			for( size_t i = 0; i < ff.size(); i++ )
			{
				T* ptr = dynamic_cast<T*>(&ff.element(i)); //safe downcast
				if( ptr != NULL )
				{
					_Vector.push_back( ptr );
					found++;
				}
			}
			return found;
		}
		catch( std::bad_typeid ex )
		{
			THROW(CodeException,"Invalid type int the obtainFFComponent() downcast call!");
		}
	}
#endif

	// some little extra classes that can be "crossed" with forcefields to make
	// special versions such as thermodynamic integration etc

	// Type of atom in annihilation-type calculations
	enum LambdaAtomGroup { Independent, Couple, Decouple };








	//-------------------------------------------------
	//
	/// \brief Extension for forcefields which are used to make TI calculations.
	/// \details This class provides a number of properties which can be inherited by a forcefield
	///         to implement thermodynamic integration
	///
	///
	///
	/// \author Mike Tyka 
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class PD_API FF_Extension_TI{
	 protected: // constructor is protected - the user shouldnt instantiate this by itself
		 FF_Extension_TI(const WorkSpace &wspace):
		  m_wspace(&wspace)
		 {

			lambda = 1.0;

			setupDefaultLambdaGroups();
		 }
	 public:
		void   setLambda(double _lambda){ lambda = _lambda; }
		double getLambda() const { return lambda; }
		double get_dEdlambda() const { return dEdlambda; }

		void printLambdaGroups(){
			printf("Lambda Group:\n Atom No.  Name    Group\n");
			for(size_t i = 0; i < m_wspace->nAtoms(); i++){
				printf(" %6d   %4s   %s\n",i,m_wspace->atom[i].pdbname.c_str(), getGroupString( mAtomGroup[i] ).c_str() );
			}
		}

		void setLambdaGroup(const PickBase &_Picker, LambdaAtomGroup group ){
			int count=0;
			for(size_t i = 0; i < m_wspace->nAtoms(); i++){
				if( _Picker.matches( m_wspace->atom[i] ) ){
					mAtomGroup[i] = group;
					count++;
				}
			}	
			printf("Set %d atoms to LambdaGroup = %s. Use 'printLambdaGroup()' to print the full list\n",
			        count, getGroupString(group).c_str() );
		}

	 protected:

		double lambda;
		double dEdlambda;

		void setupDefaultLambdaGroups() {
		 for(size_t i=0; i < m_wspace->nAtoms(); i++ )
		 { 
		 		mAtomGroup.push_back( Independent ); 
		 }
		}

		std::string getGroupString( LambdaAtomGroup group ){
				switch( group ){
					case Independent:   return("Independent"); break;
					case Couple:        return("Couple"); break;
					case Decouple:      return("Decouple"); break;
				}
				return "";
		}

		std::vector <LambdaAtomGroup> mAtomGroup;  // takes a value of 0, 1 or 2


		void info() const{ // prints a little block of parameter information
			printf("INFO: Lambda:                   %6.3lf\n", lambda);
		}
	 private:
	  const WorkSpace *m_wspace;
	};

	double calc_q_IdealGas(unsigned N, double mass, double T, double V, bool dist);

} // namespace Physics

#endif

