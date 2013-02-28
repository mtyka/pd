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

// MMLib Headers
#include "workspace/workspace.h"
#include "workspace/rotbond.h"
#include "forcefields/forcefield.h"

// Self Header
#include "minimise.h"

using namespace Physics;
using namespace Maths;

namespace Protocol
{
	Minimisation::Minimisation(Physics::Forcefield & _ff):
		PickedProtocolBase(_ff)
	{
		name = "Minimiser";
		settodefault();
	}

	Minimisation::Minimisation( Physics::Forcefield& _ff, const PickBase& _Picker ) :
		PickedProtocolBase( _ff, _Picker )
	{
		settodefault();
	}

	Minimisation* Minimisation::clone() const 
	{ 
		return new Minimisation(*this); 
	}

	void Minimisation::settodefault()
	{
		Step = 0;
		Algorithm = ConjugateGradients;
		StepSize = 2E5;
		SlopeCutoff = (double) -1.0;
		Strictness = 0.0; 
	}

	void Minimisation::info() const
	{
		printf("UpdateScr:      %d\n", UpdateScr);
		printf("UpdateTra:      %d\n", UpdateTra);
		printf("Algorithm:      " );

		switch (Algorithm) 
		{
			case SteepestDescent: printf("SteepestDescent"); break;
			case ConjugateGradients: printf("ConjugateGradients"); break;
			default: printf("ERROR");
		}

		printf("\n");
		printf("Steps           %d\n", Steps );
		printf("StepSize        %e\n", StepSize );
		printf("SlopeCutoff     %e\n", SlopeCutoff);
		printf("Strictness      %e\n", Strictness );
	}

	void Minimisation::infoLine() const 
	{ 
		// prints a line of current energies
	
		if( m_OldEnergy >= getWSpace().ene.epot ){
			printf("%8d\t%3.1e\t%6.2lf\t%6.2lf\t %22.10lf \t",
				Step, 
				m_StepMultiplier, 
				getWSpace().ene.cRMS, 
				getWSpace().ene.dRMS, 
				(double)getWSpace().ene.epot * PhysicsConst::J2kcal * PhysicsConst::Na
				);
		}else{
			printf("%8d\t%3.1e\t%6.2lf\t%6.2lf\t[%22.10lf]\t",
				Step, 
				m_StepMultiplier, 
				getWSpace().ene.cRMS, 
				getWSpace().ene.dRMS, 
				(double)getWSpace().ene.epot * PhysicsConst::J2kcal * PhysicsConst::Na
				);
		}
		
		mon.printCurData();
		ff->infoLine();
		printf("\n");
	}

	void Minimisation::infoLineHeader() const 
	{ 
		// prints the headers for the above function
		printf("%8s\t%7s\t%6s\t%7s\t %22s \t", "Step", "Stepmul", "CRMS", "dRMS", "Epot");
		mon.printHeader();
		ff->infoLineHeader();
		printf("\n");
	}

	int Minimisation::runcore()
	{
		m_StepMultiplier = 1.0;
		Step = 0;
		ff->calcForces();
		m_OldEnergy = getWSpace().ene.epot;

		if(OutputLevel){
			if((UpdateScr>0)&&(UpdateScr<Steps))
			infoLineHeader();
		}
		int breakcount = 0;
		double localOldEnergy;

		// set the last lowest found to the current structure
		getWSpace().old = getWSpace().cur;

		starttime = (int) time(NULL);
		for(Step = 0; Step < Steps; Step++) 
		{ 
			// Main principal loop for simulation
			getWSpace().Step = Step; // let particle system know about the current Step nr
			refreshNeighborList();
			ff->calcForces();
			if(!isNumber(getWSpace().ene.epot)) {
				if(OutputLevel)
					printf("Simulation unstable. Terminating ... \n");
				return -1;
			}

			// Do some statistical analysis
			runmonitors();

			if(Step == 0)
				localOldEnergy = double(getWSpace().ene.epot);
			if((Step > 0) && (SlopeCutoff > 0.0)) 
			{
				// if we're meant to quit past a certain energy gradient
				if(((localOldEnergy - double(getWSpace().ene.epot)) >= 0.0) &&
					((localOldEnergy - double(getWSpace().ene.epot)) < SlopeCutoff))
					breakcount++;
				localOldEnergy = getWSpace().ene.epot;
			}
			if(breakcount > 10)
				break;

			switch (Algorithm) 
			{
				case Minimisation::SteepestDescent:
					doSteepestDescentStep();
					break;
				case Minimisation::ConjugateGradients:
					doConjugateGradientStep();
					break;
				default:
					THROW(CodeException,"Unknown Minimisation Type");	
					break;
			}

			// Save the coordinates in the trajectory as required
			if( (Step > 0) && (UpdateTra > 0) && ((Step % UpdateTra) == 0) )
			{
				getWSpace().outtra.append();
			}

			// Display any information as required
			if( (OutputLevel) && (Step > 0) && (UpdateScr > 0) && ((Step % UpdateScr) == 0) )
			{
				infoLine();
			}
		}
		endtime = (int) time(NULL);

		// set final structure to the last lowest found
		getWSpace().cur = getWSpace().old;

		// update the energies
		ff->calcForces();

		return Step;
	}

	int Minimisation::doSteepestDescentStep()
	{
		// Set up proxies to workspace to make code more readable.
		size_t length = getPosPointer().size();
		dvector f;

		if(Step != 0) 
		{
			if((double(getWSpace().ene.epot) <= m_OldEnergy)) 
			{ 
				// energy did go down on last Step
				m_StepMultiplier *= 1.2;
				for(size_t i = 0; i < length; i++) 
				{ 
					// save positions & forces					
					getPosPointer().op(i).setTo(getPosPointer().p(i)); // save current position & force
					getPosPointer().of(i).setTo(getPosPointer().f(i));
					m_OldEnergy = getWSpace().ene.epot; // accept new energy as new reference point
				}
			} 
			else 
			{ 
				// we went up the other side
				m_StepMultiplier *= 0.5 / 1.2;
				for(size_t i = 0; i < length; i++) 
				{ 
					// save positions & forces					
					getPosPointer().p(i).setTo(getPosPointer().op(i)); // restore old position
					getPosPointer().f(i).setTo(getPosPointer().of(i)); // restore old forces
				}
				m_OldEnergy += 1;
				return 0;
			}
		} 
		else 
		{
			for(size_t i = 0; i < length; i++) 
			{ 
				// save positions & forces				
				getPosPointer().op(i).setTo(getPosPointer().p(i));
				getPosPointer().of(i).setTo(getPosPointer().f(i));				
			}
			m_OldEnergy = getWSpace().ene.epot; // accept new energy as new reference point
		}

		for(size_t i = 0; i < length; i++) 
		{
			
			f.setTo(getPosPointer().f(i));
			f.mul(StepSize * m_StepMultiplier);
			getPosPointer().p(i).add(f);
		}

		return 0;
	}

	int Minimisation::doConjugateGradientStep()
	{
		PosPointer& atoms = getPosPointer();

		// Set up proxies to workspace to make code more readable.
		size_t length = atoms.size();

		dvector f;
		double fmagold, fmagnew, gamma;

		if(Step != 0) 
		{
			if(((m_StepMultiplier < 0.05) &&
				((double(getWSpace().ene.epot) - m_OldEnergy) < ( Strictness* Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na ) )) ||
				(double(getWSpace().ene.epot) < m_OldEnergy)) 
			{ 
				// energy did go down on last Step
				m_StepMultiplier *= 1.2;

				fmagold = 0;
				fmagnew = 0;
				for(size_t i = 0; i < length; i++) 
				{ 
					// calculate the inner dots of the current force;					
					dvector temp;
					temp.diff( atoms.of(i), atoms.f(i) );
					fmagold += atoms.of(i).innerdot();
					fmagnew += temp.scalarProduct(atoms.f(i)) ; // atoms.f(i).innerdot();
				}

				if(fmagold != 0) gamma = fmagnew / fmagold;
				else             gamma = 0;
				
				if( gamma < 0 )  gamma = 0;

				for(size_t i = 0; i < length; i++) 
				{ 
					// save positions & forces
					atoms.v(i).setTo(atoms.ov(i)); // set to old direction
					atoms.v(i).mul(gamma); // multiply the old direction with gamma
					atoms.v(i).add(atoms.f(i));// and add current force to make current direction
					atoms.op(i).setTo(atoms.p(i)); // save position
					atoms.of(i).setTo(atoms.f(i)); // save old forces
					atoms.ov(i).setTo(atoms.v(i)); // save old directions
				}

				m_OldEnergy = getWSpace().ene.epot;
			} 
			else 
			{ 
				// we went up the other side so reduce stepsize and reset positions. 
				m_StepMultiplier *= 0.5 / 1.2;
				for(size_t i = 0; i < length; i++) 
				{ 
					// save positions & forces					
					atoms.p(i).setTo(atoms.op(i)); // restore old position
				}
			}
		} 
		else 
		{ 
			// this block is for Step == 0 - its just a standard SD Step
			for(size_t i = 0; i < length; i++) 
			{ 
				// save positions & forces				
				atoms.op(i).setTo(atoms.p(i)); // save position (old position = current position)
				atoms.of(i).setTo(atoms.f(i)); // save old forces
				atoms.ov(i).setTo(atoms.f(i)); // save old directions, equal to old force
				atoms.v(i).setTo(atoms.f(i)); // save old directions, equal to old force
			}
			m_OldEnergy = getWSpace().ene.epot;
		}

		for(size_t i = 0; i < length; i++) 
		{			
			f.setTo(atoms.v(i));
			f.mul(StepSize * m_StepMultiplier);
			atoms.p(i).add(f);
		}

		return 0;
	}
} // namespace 'Protocol'


