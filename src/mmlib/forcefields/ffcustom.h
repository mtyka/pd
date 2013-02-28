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

#ifndef __FFCUSTOM_H
#define __FFCUSTOM_H

#include "forcefields/forcefield.h" // Provides a base class

namespace Physics{

	class PD_API FF_Custom;






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
	class PD_API CustomForce
	{
	public:
		friend class PD_API FF_Custom;

		enum CustomForceType
		{
			INACTIVE = 0,
			HARMONIC,
			BELLSHAPED,
			VSHAPED,
			INVREP, // 1/x repulsive potential
			LINEARRESTRAINT
		}; // linear distance restraint

		CustomForce()
			: Type(INACTIVE), 
			i(-1), 
			j(-1), 
			Active(false) 
		{
			absv.setTo(0, 0, 0);
		}

		CustomForce(int _i, int _j)
			: i(_i), 
			j(_j), 
			Active(true) 
		{
			absv.setTo(0, 0, 0);
		}

		CustomForce(int _i, const Maths::dvector& _absv)
			: i(_i), 
			j(-1), 
			Active(true) 
		{
			absv.setTo(_absv);
		}

		virtual void setup()
		{
		}

		void calcForceSingle(double Dist_ij, double &potential, double &force);

		void info(const WorkSpace &wspace) const;
		bool compareAtoms(const CustomForce& cf);

		void setToHarmonic(double _tdist, double _epsilon)
		{
			Type = HARMONIC;
			tdist = _tdist;
			epsilon = _epsilon / PhysicsConst::Na * PhysicsConst::kcal2J;
		}

		void setToBellShaped(double _tdist, double _epsilon, double _gamma)
		{
			Type = BELLSHAPED;
			tdist = _tdist;
			epsilon = _epsilon / PhysicsConst::Na * PhysicsConst::kcal2J;
			gamma = _gamma;
		}

		void setToVShaped(double _tdist, double _epsilon, double _gamma, double _beta)
		{
			Type = VSHAPED;
			tdist = _tdist;
			epsilon = _epsilon / PhysicsConst::Na * PhysicsConst::kcal2J;
			gamma = _gamma;
			beta = _beta;
		}

		void setToLinearRestraint(double _tdist, double _gamma, double _tdist2, double _delta)
		{
			Type = LINEARRESTRAINT;
			tdist = _tdist;
			gamma = _gamma / PhysicsConst::Na * PhysicsConst::kcal2J;
			tdist2 = Maths::max(tdist, _tdist2);
			delta = _delta;
		}

		void activate()
		{
			if(Type != INACTIVE)
				Active = true;
		}

		void deactivate()
		{
			Active = false;
		}

		void setAbsPosTo( const Maths::dvector& pos )
		{
			absv.setTo(pos);
		}

		void printTextLine(WorkSpace &wspace,FILE * file); // prints a text description of itself

		int ID; // ID number for grouping etc..
	private:
		void calcEnergyVerbose(FF_Custom * cff, ForcefieldBase::AtomicVerbosity level);
		void calcEnergy(FF_Custom * cff);
		void calcForce(FF_Custom * cff);

		double rdist; // double distance
		int i,j; // acts between atoms i and j
		Maths::dvector absv; // if j == NULL, absv defines an absolute partner point
		CustomForceType Type; // types to be defined
		bool Active; // 0 = inactive, else Active

		double tdist; // target distance
		double tdist2; // second distance
		double epsilon; // welldepth Type parameter
		double alpha, beta, // additional parameters to be used n different ways
			gamma, delta; // by different forcefield types

	};


	//-------------------------------------------------
	//
	/// \brief Custom interaction forcefield 
	///
	/// \details 
	/// This defines a forcefield which can be used to define individual
	/// interactions between atom pairs with different potential
	/// types. It is the basis for all kinds of constraint & itemised
	/// forcefields
/*!
Custom forces can also be read in from a file:

\section custom_force_definition_file Custom Forcefield Definition file '.con'} 


\verbatim

File syntax:


General Parameters ..
CUSTOMFORCELIST
   1 'forcedefinition' per line
   ....
ENDCUSTOMFORCELIST



Force definition syntax:


ResidueNr1  AtomName1  ResidueNr2   AtomName2  Forcetype [Param1 Param2 ...]
Forcetype = {H | B | V | LR}


H   Harmonic}


Parameters:  tdist  epsilon 
      tdist:     equilibrium distance in    Angstrom
      epsilon:   forceconstant in           kcal/mol/A^2
Formula      V = epsilon (dist - tdist)^2

Shape:       o      o
      o      o
      o      o
       o    o
       o    o
	o  o
 0 ------oo------------


B   BellShaped}


Parameters:  tdist  epsilon  gamma

      tdist:     equilibrium distance in    Angstrom
      epsilon:   welldepth in               kcal/mol
      gamma      wellwidth                  -- (dimensionless)   (the higher the value the tighter)

Formula      V =  epsilon + epsilon/(1 + gamma*sqr(distij - tdist));

Shape:       epsilon   - ooooooo          ooooooooooo
		|         o      o
		|          o    o
		|          o    o
		|           o  o
                          0   `------------oo----------------------


V   Vshaped


Parameters:  tdist  epsilon  gamma  beta 

	tdist:   point at which extrapolated inear part meets x axis in Angstroms
	epsilon: gradient of linear part kcal/mol/A
	gamma:   sharpness -- (dimensionless) (the higher the value the tighter)
	beta:    reversal should take 1.0 or -1.0, -1.0 mirrows the function around the y axis

 Formula      V = epsilon * gamma * ln(1.0 + exp(beta*(distij - tdist)/gamma));

 Shape:                 -                   o    
			|                  o
			|                 o   <-- gradient = epsilon
			|                o
			|              oo   <-- curvature dependant on gamma 
		    0   `-oooooooooooo ----------------------


A symmetric V-Shaped well can be generated by placing two of these potentials back-to-back:

 Shape:                 - o                         o    
			|  o                       o
			|   o                     o   <-- gradient = epsilon
			|    o                   o
			|     oo               oo   <-- curvature dependant on gamma 
		    0   `------- ooooooooooooo ----------------------

LR  Linear Restraint}

 Parameters:  tdist  epsilon    tdist2  epsilon2
	      tdist:     interceptpoint 1           Angstrom
	      epsilon:   gradient of linear part1   kcal/mol/A
	      tdist:     interceptpoint 2           Angstrom
	      epsilon:   gradient of linear part2   kcal/mol/A

 Formula      
		     /
		    |   -gamma*(distij - tdist)     distij < tdist
	      V =  <    0                           tdist < distij < tdist2
		    |    gamma*(distij - tdist2);           tdist2 < distij
		     \

 Shape:                 -  o                o    
			|   o              o
     epsilon = gradient -->  o            o   <-- gradient = epsilon2
			|     o          o
			|      o        o   <-- sharp corner, not so good for MD
		    0   `-------oooooooo---------------------
				^       ^
			    tdist      tdist2

\endverbatim


*/
	///
	/// \author Mike Tyka  
	///
	/// \todo STATE OF DEVELOPMENT
	///
	///
	class PD_API FF_Custom : public ForcefieldBase
	{
	public:
		friend class PD_API CustomForce;

		FF_Custom( WorkSpace &newwspace, Verbosity::Type _verbose): 
		ForcefieldBase( newwspace )
		{
			epot = 0.0;
			verbose = _verbose;
			name = "customff";
			name_verbose = "Custom FF";
			name_infoheader = "Ecustom";
		}

		virtual FF_Custom* clone() const { return new FF_Custom(*this); }

		virtual void setup()
		{
			needsetup = false;
		}

		virtual void calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level);
		virtual void calcEnergies();
		virtual void calcForces();
		void calcForceSingle(double Dist_ij, double &totpot, double &totforce);

		virtual void info() const; // prints a little block of parameter information
		virtual void infoLine() const; // prints a line of current energies
		virtual void infoLineHeader() const; // prints the headers for the above function

		CustomForce& getForce( size_t i )
		{
			return force[i];
		}

		// adds/removed modifies forces
		void addForce(const CustomForce &cf)
		{
			force.push_back(cf);
		}

		void clear()
		{
			force.clear();
		}

		void empty()
		{
			force.clear();
		}

		size_t size() const { return force.size(); }

		int findForce(const CustomForce& cf);

		// these (de)activate all or subsets of the custom forces (according to ID)
		void activateforces();
		void activateforces(int ID);
		void activateforcesnot(int ID);
		void deactivateforces();
		void deactivateforces(int ID);
		void deactivateforcesnot(int ID);

		// loads&parses a file containing the restraint definitions
		int loadTextFile(const char *filename);
		int saveTextFile(const char *filename);
		int load(const char *filename); // loads the forcefield from a binary form
		int save(const char *filename); // saves the forcefield in a binary form

		std::string name_verbose;
		std::string name_infoheader;

	protected:
		std::vector < CustomForce > force;
		Verbosity::Type verbose;
	};

}

#endif

