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

#include "workspace/workspace.h"
#include "workspace/neighbourlist.h"
#include "pickers/basicpickers.h"

#include "forcefield.h"

// namespace includes
using namespace Maths;

namespace Physics
{
	// Various Physical Constants - Description & [Units] given for each. (as Doxygen comments)
	// These are from The NIST Reference on Constants, Units and Uncertainty (2006)
	// http://pConstanthysics.nist.gov/constants

	const double PhysicsConst::kcal2J = 4184.0;                                         // Multiplier to convert Joules to kcal [kcal/J]
	const double PhysicsConst::J2kcal = 1.0 / kcal2J;                                   // Multiplier to convert kcals to Joule [J/kcal]
	const double PhysicsConst::Na = 6.0221418E23;                                       // Avogadros's Number - i.e. no of particles per mol [1/mol]
	const double PhysicsConst::kcal2JDivNa = kcal2J / Na;                               // Conversion Factor from kcal/mol to Joules/molecule
	const double PhysicsConst::kB = 1.38065E-23;                                        // Boltzmann Constant [J/K]
	const double PhysicsConst::R_kcal = kB * Na * J2kcal;                               // The Gas Constant in kcal/mol/K (0.5957 kcal/mol at T=300)
	const double PhysicsConst::planck = 6.626069E-34;                                   // Plank Constant [Js] 
	const double PhysicsConst::clight = 2.99792458E8;                                   // Speed of Light [m/s]
	const double PhysicsConst::amu = 1.66053878E-27;                                    // Atomic mass unit (mass of Carbon-12 divided by 12) [kg]
	const double PhysicsConst::Angstrom = 1E-10;                                        // Length of an Angstrom in meters [m/A]
	const double PhysicsConst::invAngstrom = 1.0/Angstrom;                              // Number of A in a meter [A/m]
	const double PhysicsConst::e_charge = 1.60217645E-19;                               // Charge of electron - Unit of Charge [C]
	const double PhysicsConst::sqr_e_charge = e_charge*e_charge;                        // Square of electronic charge [C^2]
	const double PhysicsConst::e0 = 8.854187817E-12;                                    // Permittivity of Space [C^2/J m]
	const double localPhysics_4pi = 3.14159265358979323846 * 4.0;                       // Somehow, Maths::FourPI contains 0.0 when used here ?!? (I've made a local version)
	const double PhysicsConst::_4pi_e0 = localPhysics_4pi * e0;                         // Electrostatic Units (premultiplier for electrostatic interactions) [Jm/C^2]
	const double PhysicsConst::econv = (J2kcal*Na*sqr_e_charge)/(_4pi_e0*Angstrom);     // Energy of two unit charges at 1 Angstrom distance [kcal/mol]
	const double PhysicsConst::halfeconv = econv*0.5;                                   // half of the above
	const double PhysicsConst::econv_joule = (sqr_e_charge)/(_4pi_e0*Angstrom);         // as above but in in joules
	const double PhysicsConst::halfeconv_joule = econv_joule*0.5;                       // and half of that
	const double PhysicsConst::Bar2Pa = 100000.0;                                       // Pressure: Bar to Pascal

	void PhysicsConst::printConstants()
	{
		 printf("PhysicsConst::kcal2J             %10.8e     Multiplier to convert Joules to kcal [kcal/J]\n",                              PhysicsConst::kcal2J         );
		 printf("PhysicsConst::J2kcal             %10.8e     Multiplier to convert kcals to Joule [J/kcal]\n",                              PhysicsConst::J2kcal         );
		 printf("PhysicsConst::Na                 %10.8e     Avogadros's Number - i.e. no of particles per mol [1/mol]\n",                  PhysicsConst::Na             );
		 printf("PhysicsConst::kcal2JDivNa        %10.8e     Conversion Factor from kcal/mol to Joules/molecule\n",                         PhysicsConst::kcal2JDivNa    );
		 printf("PhysicsConst::kB                 %10.8e     Boltzmann Constant [J/K]\n",                                                   PhysicsConst::kB             );
		 printf("PhysicsConst::R_kcal             %10.8e     The Gas Constant in kcal/mol/K (0.5957 kcal/mol at T=300)\n",                  PhysicsConst::R_kcal         );
		 printf("PhysicsConst::planck             %10.8e     Plank Constant [Js] \n",                                                       PhysicsConst::planck         );
		 printf("PhysicsConst::clight             %10.8e     Speed of Light [m/s]\n",                                                       PhysicsConst::clight         );
		 printf("PhysicsConst::amu                %10.8e     Atomic mass unit (mass of Carbon-12 divided by 12) [kg]\n",                    PhysicsConst::amu            );
		 printf("PhysicsConst::Angstrom           %10.8e     Length of an Angstrom in meters [m/A]\n",                                      PhysicsConst::Angstrom       );
		 printf("PhysicsConst::invAngstrom        %10.8e     Number of A in a meter [A/m]\n",                                               PhysicsConst::invAngstrom    );
		 printf("PhysicsConst::e_charge           %10.8e     Charge of electron - Unit of Charge [C]\n",                                    PhysicsConst::e_charge       );
		 printf("PhysicsConst::sqr_e_charge       %10.8e     Square of electronic charge [C^2]\n",                                          PhysicsConst::sqr_e_charge   );
		 printf("PhysicsConst::e0                 %10.8e     Permittivity of Space [C^2/J m]\n",                                            PhysicsConst::e0             );
		 printf("PhysicsConst::_4pi_e0            %10.8e     Electrostatic Units (premultiplier for electrostatic interactions) [Jm/C^2]\n",PhysicsConst::_4pi_e0        );
		 printf("PhysicsConst::econv              %10.8e     Energy of two unit charges at 1 Angstrom distance [kcal/mol]\n",               PhysicsConst::econv          );
		 printf("PhysicsConst::halfeconv          %10.8e      half of the above\n",                                                         PhysicsConst::halfeconv      );
		 printf("PhysicsConst::econv_joule        %10.8e      as above but in in joules\n",                                                 PhysicsConst::econv_joule    );
		 printf("PhysicsConst::halfeconv_joule    %10.8e      and half of that\n",                                                          PhysicsConst::halfeconv_joule);
		 printf("PhysicsConst::Bar2Pa             %10.8e     Pressure: Bar to Pascal\n",                                                    PhysicsConst::Bar2Pa         );
	}


	// --------------
	// ForcefieldBase
	// --------------

	int ForcefieldBase::ensuresetup(WorkSpace &newwspace)
	{
		// first check that the new pointers provided match the old ones - if not
		// we're potentially dealing with a new wspace/ffps and we MUST call setup(..)!
		if(&newwspace != &getWSpace() ) needsetup = true;		

		// Check if checksum has changed
		unsigned long _CheckSum = calcCheckSum();
		if(_CheckSum != m_CheckSum) needsetup = true; // we need to setup again if checksum has changed	

		int result = 0;
		if( needsetup )
		{
			printf("Setting up forcefield:  %s \n",name.c_str());
			setup();
		}

		// remember latest checksum 
		m_CheckSum	= calcCheckSum();
		return 0;
	}

	void ForcefieldBase::calcEnergies()
	{ 
		if(!Active) return;
		calcForces(); 
	}

	void ForcefieldBase::calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level)
	{
		if(!Active) return;
		calcForces();
		printf(" %-16s%10.3lf kcal/mol\n", name.c_str(), double (epot) * PhysicsConst::J2kcal * PhysicsConst::Na);
	}	

	void ForcefieldBase::info() const
	{				
		printf(" %s  ", name.c_str());
		if(!Active) printf("[InActive]  ");
		if(!OutputLevel) printf("[!OutputLevel]  ");
		if(Passive) printf("[Passive]  ");
		printf("\n");
	}

	void ForcefieldBase::infoLine() const 
	{			
		if(OutputLevel)
		{
			if(Passive)printf("( ");
			printf("\t%8.2lf", double (epot) * PhysicsConst::J2kcal * PhysicsConst::Na);
			if(Passive)printf(" )");
		}
	}

	void ForcefieldBase::infoLineHeader() const 
	{			
		if(OutputLevel)
		{
			if(Passive)printf("( ");
			printf("\t%8s", ShortName.c_str());
			if(Passive)printf(" )");
		}
	}

	void ForcefieldBase::activate()
	{
		Active = true;
	}

	void ForcefieldBase::deactivate()
	{
		Active = false;
		epot = 0.0; // this should really be reset, as the component stops calculating it once deactivated...
	}


	// ------------
	//  Forcefield
	// ------------
	

	void Forcefield::addPreAddition( ForcefieldBase  &newelement )
	{
		// Check if the forcefieldComponent about to be added uses the same
		// workspace as this forcefield container, by matching memory addresses.
		if( ! (&newelement.getWSpace() ==  &getWSpace() )){
			throw(ArgumentException(
				"You cannot add this force field component (" 
				+ newelement.name +
				") to this forcefield because\nthey use different workspaces.\n"));
		}
	}

	void Forcefield::setup()
	{
		for(size_t i = 0; i < size(); i++)
		{
			element(i).setup();
		}
		needsetup = false; // flag that we're set up
	}

	int Forcefield::ensuresetup(WorkSpace &newwspace)
	{
		if(&newwspace != &getWSpace()){
			throw( ProcedureException("This forcefield (" + name + 
			  ") was initialised with a different workspace than you intend to use it with." ) );
		}

		// else just check that each forcefield is setup in turn
		for(size_t i = 0; i < size(); i++)
		{
			if((element(i).ensuresetup(getWSpace())) < 0)
			{
				throw ProcedureException("Unknown error occured during forcefield setup");
			}
			element(i).needsetup = false;
		}
		
		needsetup = false; // if all setups were successful then flag that we're ok
		return 0;
	}

	void Forcefield::calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level)
	{
		printf(" -----------------------------------\n");
		printf(" --- Verbose Energy Summary ---\n");
		printf(" -----------------------------------\n");
		getWSpace().zeroForces();
		// loop through all the initiated forcefields and make them calculate their energies
		size_t nforcefields = size();
		for(size_t i = 0; i < nforcefields; i++)
			element(i).calcEnergiesVerbose(level);
		if(level != ForcefieldBase::Summary) 
		{
			printf("\n");
			for(unsigned int i = 0; i < nforcefields; i++)
				element(i).calcEnergiesVerbose(ForcefieldBase::Summary);
		}
		printf(" -----------------------------------\n");
		printf(" Total energy: %10.3lf kcal/mol\n\n", (double) getWSpace().ene.epot * PhysicsConst::J2kcal * PhysicsConst::Na);
	}
	
	void Forcefield::calcEnergies()
	{
		getWSpace().zeroForces();
		size_t nforcefields = size();
		for(size_t i = 0; i < nforcefields; i++)
			element(i).calcEnergies();
	}
	
	void Forcefield::calcForces(){
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array
		getWSpace().zeroForces();
		size_t nforcefields = size();
		for(size_t i = 0; i < nforcefields; i++)
		{
			element(i).calcForces();
		}
		if(max_force > 0) 
		{
			int forcetrunc = 0;
			for(int i = 0; i < getWSpace().atom.size(); i++) 
			{
				double magn = atom[i].f.mag();
				if(magn > max_force) 
				{ 
					atom[i].f.mul(max_force / magn); // truncate excessive forces!
					forcetrunc++;
				}
			}
			// if(forcetrunc > 0) printf("Force truncation: %d \n",forcetrunc);
		}
		// check that none of the component energies isnt a NAN
		if(!isNumber(getWSpace().ene.epot_vdw))
			getWSpace().ene.epot = dNAN();
		if(!isNumber(getWSpace().ene.epot_elec))
			getWSpace().ene.epot = dNAN();
	}

	void Forcefield::info() 
	{ 
		ensuresetup( getWSpace() );
		// prints a little block of parameter information
		size_t nforcefields = size();
		printf("    Forcefields: %d      \n", nforcefields);
		printf("                         \n");
		for(size_t i = 0; i < nforcefields; i++)
			element(i).info();
	}

	void Forcefield::infoLine() const 
	{ 
		// prints a line of current energies
		size_t nforcefields = size();
		for(size_t i = 0; i < nforcefields; i++)
			element(i).infoLine();
	}

	void Forcefield::infoLineHeader() const 
	{
		// prints the headers for the above function
		size_t nforcefields = size();
		for(size_t i = 0; i < nforcefields; i++)
			element(i).infoLineHeader();
	}

	void Forcefield::printEnergyShort()
	{
		if(ensuresetup(getWSpace())!=0) return;
		getWSpace().nlist().calcNewList();
		calcEnergies();
		printf("epot: %.7lf \n", getWSpace().ene.epot * PhysicsConst::Na * PhysicsConst::J2kcal );
	}

	void Forcefield::printEnergySummary()
	{
		if(ensuresetup(getWSpace())!=0) return;
		getWSpace().nlist().calcNewList();
		calcEnergiesVerbose(ForcefieldBase::Summary);
	}

	void Forcefield::printEnergyByAtom()
	{
		if(ensuresetup(getWSpace())!=0) return;
		getWSpace().nlist().calcNewList();
		calcEnergiesVerbose(ForcefieldBase::Atomwise);
	}

	void Forcefield::printEnergyDetailed()
	{
		if(ensuresetup(getWSpace())!=0) return;
		getWSpace().nlist().calcNewList();
		calcEnergiesVerbose(ForcefieldBase::Detailed);
	}

	void Forcefield::listForcefields()
	{
		std::cout << "Forcefield conatins:" << std::endl;
		for( size_t i = 0; i < size(); i++ )
		{
			std::cout << element(i).name << std::endl;
		}
	}









	// some classic theoretical functions
	double calc_q_IdealGas(unsigned N, double mass, double T, double V, bool dist)
	{
		double q = log( cube(sqrt( 2 * MathConst::PI * mass * PhysicsConst::amu *
			                 PhysicsConst::kB * T / sqr(PhysicsConst::planck))) * V * 1E-30);
		q *= double(N);
		if(!dist) q -= log_factorial(N);
		printf("Partition Function: IdealGas of %d %s particle(s):\n", N, dist ? "distinguishable" : "indistinguishable" );
		printf("  Mass:   %10.7lf amu \n", mass );
		printf("  Temp:   %10.7lf K \n", T);
		printf("  Volume:   %10.7lf Angstrom^3\n", V );
		printf("  ln(q) =   %15.9lf \n",q);
		return q;
	}
} // namespace Physics


