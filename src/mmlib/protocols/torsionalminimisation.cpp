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
#include "workspace/pospointer.h"
#include "forcefields/forcefield.h"
#include "torsionalminimisation.h"

// namespace includes
using namespace Physics; // for forcefields and physical constants
using namespace Maths;

namespace Protocol
{
	//// -----------------------------------------------------------------------------------------------------------
	//// Begin: 'TorsionalMinimisation'
	//// Implementation of a Torsional Minimisation that works within the confines of small section 
	//// (allowing for breaks).
	//// -----------------------------------------------------------------------------------------------------------

	const double TorsionalMinimisation::RESET_PREV_SQUARE = 1E30;

	TorsionalMinimisation::TorsionalMinimisation( Physics::Forcefield &_ff ):
		RangesProtocolBase( _ff )
	{
		settodefault();
		initRotations(); // initialise the rotations that will be used during simulation
	}

	TorsionalMinimisation::TorsionalMinimisation( Physics::Forcefield &_ff, const PickAtomRange& _newRange ):
		RangesProtocolBase( _ff, _newRange )
	{
		settodefault();	
		initRotations(); // initialise the rotations that will be used during simulation
	}

	TorsionalMinimisation::TorsionalMinimisation( Physics::Forcefield &_ff, const PickAtomRanges& _newRange ):
		RangesProtocolBase( _ff, _newRange )
	{
		settodefault();	
		initRotations(); // initialise the rotations that will be used during simulation
	}

	void TorsionalMinimisation::settodefault()
	{
		UpdateNList = 1; // Anything else is really too lax

		// "public changable values"
		Steps = 60; // most energy is eeked out by 60 Steps
		RotationExcldMode = Default; // None, //(RotationExcludeMode)(ImproperSingles | Omega),
		SlopeCutoff = -1.0;
		InitialCapFactor = 1.0; // Change this to take the initial steric state into account (0.1 for very harsh starting conditions)
		
		// Leave the rest as Allans values
		CapRise = 1.2;
		CapDrop = 0.5;
		IntendedMaxAngleCap = Maths::DegToRad(1.0); //established defaults - might not be optimum for different types of structure?
		IntendedMaxSCAngleCap = Maths::DegToRad(5.0);
		SCAngleFactor = 4.33E7; // Was 0.003 = allans SCAngleFactor, but now we are not using AllanForceFactor = -14392580000.0; (mikes to allans force units conversion factor)
	};

	void TorsionalMinimisation::initRotations()
	{
		RotBond& rotbond = getWSpace().rotbond(); // local accessor for the duration of this function call
		if( rotbond.size() == 0 )
		{
			THROW(ProcedureException,"Procedural error during TorsionalMinimisation::initRotations(), the number of rotatable bonds is '0'!");
		}

		m_InitPickerSerial = getPickerSerial(); // we need to cache this - if the segment is changed, then runcore *MUST* recall initRotations()
		const PickAtomRanges& pickRanges = getRanges();
		m_BestStore.setPicking( getWSpace(), pickRanges );

		ChiRotDefs.clear();
		PhiPsiRotDefs.clear();

		const ParticleStore& atom = getWSpace().atom;

		// find out how many we have of each Type first...
		for( int i = 0; i < rotbond.size(); i++) // look in the rotatable bond array ...
		{
			bool matchFirst = pickRanges.matches( atom[rotbond[i].i] );
			if( !matchFirst ) continue;
			size_t firstMatchRange = pickRanges.getPreviousMatchRangeIndex();
			bool matchSecond = pickRanges.matches( atom[rotbond[i].j] );
			if( !matchSecond ) continue;
			size_t secondMatchRange = pickRanges.getPreviousMatchRangeIndex();
			if( firstMatchRange != secondMatchRange )
			{
				// rotbonds bridges range, we assume that there are dual-branch
				continue;
			}

			if( (rotbond[i].Type == RotatableBond::Single) )
			{
				if( (0 == (RotationExcldMode & ImproperSingles)) || // if 'Exclude ImproperSingles' isn't flagged
					!rotbond.isDelocalised(i) )
				{
					setupRotDefChi( rotbond[i] );
				}
			}
			else if( (rotbond[i].Type == RotatableBond::Phi) ||
				(rotbond[i].Type == RotatableBond::Psi) ||
				( (0 == (Omega &RotationExcldMode)) // exclude 'Omega' isn't flagged
				&& (rotbond[i].Type == RotatableBond::Omega))
				)
			{
				setupRotDefPhiPsi( rotbond[i], firstMatchRange );
			}
		}
	}

	void TorsionalMinimisation::setupRotDefPhiPsi( RotatableBond &rotBond, size_t forRange )
	{
		RotationDefinition_TM rotDef;
		rotDef.rangeLink = forRange;

		size_t rangeStart = getStartAtom(forRange);
		size_t rangeEnd = getEndAtom(forRange);
		bool reverseIt = getReversed(forRange);

		// Force the reverse flag if we are dealing with ranges which sit at the terminii
		int natom = getWSpace().nAtoms();
		if( rangeStart == 0  )
		{
			reverseIt = (rotBond.i <= (natom / 2));
		}
		else if( rangeEnd == natom-1 )
		{
			reverseIt = (rotBond.j <= (natom / 2));
		}

		bool swapIndexers = rotBond.i < rotBond.j; // check for a rear-facing rotbond
		if( reverseIt ) 
			swapIndexers = !swapIndexers; // If we are looking reerward, then it needs to flip again

		// Assumptions about the atom ordering are made in this code.
		// We also need to make sure our bond is pointing in the correct direction!
		if( swapIndexers )
		{
			rotDef.A = rotBond.i;
			rotDef.B = rotBond.j;			
			rotDef.atom1 = &getWSpace().cur.atom[rotBond.ip].p;
			rotDef.atom2 = &getWSpace().cur.atom[rotBond.i ].p;
			rotDef.atom3 = &getWSpace().cur.atom[rotBond.j ].p;
			if( rotBond.jp >= rangeStart )
			{
				rotDef.atom4 = &getWSpace().cur.atom[rotBond.jp].p;
			}
			else
			{
				// We have a problem ....
				// Not completely sure about this one....
				// Required for the N-terminal residue of a dual-branch broken conformer
				// There may be other cases, but then .j+1 may not be apropriate?!?
				rotDef.atom4 = &getWSpace().cur.atom[rotBond.j+1].p;
			}
		}
		else
		{
			rotDef.A = rotBond.j;
			rotDef.B = rotBond.i;
			rotDef.atom1 = &getWSpace().cur.atom[rotBond.jp].p;
			rotDef.atom2 = &getWSpace().cur.atom[rotBond.j ].p;
			rotDef.atom3 = &getWSpace().cur.atom[rotBond.i ].p;
			rotDef.atom4 = &getWSpace().cur.atom[rotBond.ip].p;
		}

		if( reverseIt )
		{
			rotDef.startAtomIndex = rangeStart;
			rotDef.endAtomIndex = rotDef.A - 1;
			if( rotBond.Type == rotBond.Omega )
				rotDef.arseAtom = rotDef.A + 1;
		}
		else
		{
			rotDef.startAtomIndex = rotDef.B + 1;
			rotDef.endAtomIndex = rangeEnd;
		}

		rotDef.Value = 0.0;
		rotDef.PrevValue = 0.0;
		rotDef.molBase = &getWSpace();
		rotDef.ReciprocalLen = 1.0 / rotDef.atom2->dist(*rotDef.atom3);

		PhiPsiRotDefs.push_back(rotDef);
	}

	void TorsionalMinimisation::setupRotDefChi( RotatableBond &rotBond )
	{
		RotationDefinition_TM rotDef;
		rotDef.rangeLink = SIZE_T_FAIL; // irrelent for Chi angles, but initialise just in case

		rotDef.A = rotBond.j;
		rotDef.B = rotBond.i;
		rotDef.atom1 = &getWSpace().cur.atom[rotBond.jp].p;
		rotDef.atom2 = &getWSpace().cur.atom[rotBond.j ].p;
		rotDef.atom3 = &getWSpace().cur.atom[rotBond.i ].p;
		rotDef.atom4 = &getWSpace().cur.atom[rotBond.ip].p;

		ASSERT( rotBond.segmentstart[1] == -1 && rotBond.segmentend[1] == -1, 
			CodeException, "Assumed that sidechains never have multiple rotatable segments");

		rotDef.startAtomIndex = rotBond.segmentstart[0] + 1;
		rotDef.endAtomIndex = rotBond.segmentend[0];

		if( rotBond.segmentstart[1] != -1 || rotBond.segmentend[1] != -1 )
		{
			// only counting the 1st segment is hackyness, but rather easier.
			THROW(CodeException,"Code assumption is invalid!");
		}

		rotDef.Value = 0.0;
		rotDef.PrevValue = 0.0;
		rotDef.molBase = &getWSpace();
		rotDef.ReciprocalLen = 1.0 / rotDef.atom2->dist(*rotDef.atom3);

		ChiRotDefs.push_back(rotDef);
	}

	double TMleastSquaresFit(const std::vector<double> &x, const std::vector<double> &y)
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
			sumxy += (x[i] - meanx) * (y[i] - meany);
		}

		return sumxy / sumxx;
	}

	int TorsionalMinimisation::runcore()
	{
		if( m_InitPickerSerial != getPickerSerial() )
		{
			initRotations(); // reinitialisation must occur, the atomic range definition has changed.
		}

		m_BestStore.store();

		// Important to reinitialise these to ensure that we can call Run again without affecting the minimisation outcome.
		for( size_t i = 0; i < ChiRotDefs.size(); i++ )
		{
			ChiRotDefs[i].Value = 0.0;
			ChiRotDefs[i].PrevValue = 0.0;
		}
		for( size_t i = 0; i < PhiPsiRotDefs.size(); i++ )
		{
			PhiPsiRotDefs[i].Value = 0.0;
			PhiPsiRotDefs[i].PrevValue = 0.0;
		}

		getWSpace().Step = Step = 0; // let WorkSpace know about the current Step #

		// set the initial Step size
		// should be 100% if we are not in steric clash!
		AngleCap = IntendedMaxAngleCap * InitialCapFactor; //if require a gentle start - depends a bit on steric condition
		SCAngleCap = IntendedMaxSCAngleCap * InitialCapFactor;

		previous_epot = DBL_MAX;
		PrevSquare = RESET_PREV_SQUARE;
		Uphill = 0;

		// Ene monitoring for grad cutoff
		const size_t gradWindowSize = 20;
		std::vector<double> xGradMarker(gradWindowSize);
		std::vector<double> yGradMarker(gradWindowSize);
		size_t currentGradMarker = 0;

		Hamiltonian best_epot; // record the best energy, store the conformation when we beat it ...
		best_epot.epot = DBL_MAX; // Big number, the first step of minimisation will beat this ...

		if(OutputLevel) infoLineHeader();

		// Main simulation inner-loop
		for( Step = 0; Step < Steps; Step++ )
		{
			if( AngleCap < 0.00001 )
			{
				// The minimisation has stalled, there is no point in doing more steps
				break;
			}

			getWSpace().Step = Step;// let particle system know about the current Step #
			//getWSpace().nlist().calcNewList();
			refreshNeighborList();

			ff->calcForces();
			if(!isNumber(getWSpace().ene.epot))
			{
				if(OutputLevel)
				{
					printf("Torsional CG Minimisation Unstable. Terminating... \n");
				}
				return -1;
			}
			if( getWSpace().ene.epot < best_epot.epot )
			{
				best_epot = getWSpace().ene;
				m_BestStore.store(); // record the best over the minimisation
			}

			ApplyXYZAbsoluteDotVector(); // 1. calculate the required torsional magnitudes.
			ActionConjugateGradient(); // 2. Normalise the magnitudes to give required changes in radians and apply.

			if(Step > 2) // until you have had 2 rounds, you cant assess "up-hill" or "down-hill"
			{
				if (getWSpace().ene.epot <= previous_epot)
				{
					Uphill = 0; // we are going "downhill"
					AngleCap *= CapRise; // if move was downhill then increase Step size (1.2 is an empiricly "good" weighting)
					SCAngleCap *= CapRise;
					if(AngleCap > IntendedMaxAngleCap) AngleCap = IntendedMaxAngleCap;
					if(SCAngleCap > IntendedMaxSCAngleCap) SCAngleCap = IntendedMaxSCAngleCap;
					//previous_epot = getWSpace().ene.epot; // : if Enabled here, disable after if(Step > 2)...
				}
				else
				{
					if( Uphill > 20 ) // The minimisation has stalled, there is no point in doing more steps
					{
						break;
					}

					Uphill++; // if move was uphill then halve Step size
					AngleCap *= CapDrop;
					SCAngleCap *= CapDrop;

					// reinitialise these
					for( size_t i = 0; i < ChiRotDefs.size(); i++ )
					{
						ChiRotDefs[i].PrevValue = 0.0;
					}
					for( size_t i = 0; i < PhiPsiRotDefs.size(); i++ )
					{
						PhiPsiRotDefs[i].PrevValue = 0.0;
					}

					PrevSquare = RESET_PREV_SQUARE;
				}
			}
			previous_epot = getWSpace().ene.epot; // : if Enabled here, disable in if(getWSpace().ene.epot <= previous_epot)...

			// Energy change over a number of Steps.
			// If no improvement over a gradient Cutoff then quit, our work is done...
			if( SlopeCutoff > 0.0 )
			{
				if( currentGradMarker < gradWindowSize )
				{
					xGradMarker[currentGradMarker] = (double)Step;
					yGradMarker[currentGradMarker++] = getWSpace().ene.epot;
				}
				else
				{
					double b = TMleastSquaresFit( xGradMarker, yGradMarker );
					if( SlopeCutoff > -b )
					{
						break;
					}
					xGradMarker[currentGradMarker % gradWindowSize] = (double)Step;
					yGradMarker[currentGradMarker % gradWindowSize] = getWSpace().ene.epot;
					currentGradMarker++;
				}
			}

			runmonitors();

			// Save the coordinates in the trajectory as required
			if((UpdateTra > 0) && ((Step % UpdateTra) == 0) )
			{
				getWSpace().outtra.append();
			}

			// Display any information as required
			if( OutputLevel && (UpdateScr > 0) && ((Step % UpdateScr) == 0) )
			{
				// now print info line
				infoLine();
				ff->infoLine();
				printf("\n");
			}
		}
		
		m_BestStore.revert();
		getWSpace().ene = best_epot; // its quicker to memcpy the energies than to ff->calcForces() ...

		return Step; // success
	} // end procedure

	void TorsionalMinimisation::ApplyXYZAbsoluteDotVector()
	{
		int AtomA;
		double VectorX, VectorY, VectorZ;
		double Delta2X, Delta2Y, Delta2Z, NewDeltaX, NewDeltaY, NewDeltaZ;

		for( size_t i = 0; i < PhiPsiRotDefs.size(); i++ )
		{
			PhiPsiRotDefs[i].Value = 0.0; // initialise to zero per "ApplyXYZAbsoluteDotVector()" call

			VectorX = PhiPsiRotDefs[i].atom3->x - PhiPsiRotDefs[i].atom2->x;
			VectorY = PhiPsiRotDefs[i].atom3->y - PhiPsiRotDefs[i].atom2->y;
			VectorZ = PhiPsiRotDefs[i].atom3->z - PhiPsiRotDefs[i].atom2->z;

			for (AtomA = PhiPsiRotDefs[i].startAtomIndex; AtomA <= PhiPsiRotDefs[i].endAtomIndex; AtomA++)
			{
				Delta2X = getWSpace().cur.atom[AtomA].p.x - PhiPsiRotDefs[i].atom3->x; //vector C to AtomA
				Delta2Y = getWSpace().cur.atom[AtomA].p.y - PhiPsiRotDefs[i].atom3->y;
				Delta2Z = getWSpace().cur.atom[AtomA].p.z - PhiPsiRotDefs[i].atom3->z;
				NewDeltaX = VectorY * Delta2Z - VectorZ * Delta2Y; //calculate cross vector product
				NewDeltaY = VectorZ * Delta2X - VectorX * Delta2Z;
				NewDeltaZ = VectorX * Delta2Y - VectorY * Delta2X;
				PhiPsiRotDefs[i].Value += (
					NewDeltaX * getWSpace().cur.atom[AtomA].f.x + 
					NewDeltaY * getWSpace().cur.atom[AtomA].f.y + 
					NewDeltaZ * getWSpace().cur.atom[AtomA].f.z)
					* PhiPsiRotDefs[i].ReciprocalLen; //dot product evaluates move
			} //next atomA
		} //next i

		for( size_t i = 0; i < ChiRotDefs.size(); i++ )
		{
			ChiRotDefs[i].Value = 0.0; // initialise to zero

			VectorX = ChiRotDefs[i].atom2->x - ChiRotDefs[i].atom3->x;
			VectorY = ChiRotDefs[i].atom2->y - ChiRotDefs[i].atom3->y;
			VectorZ = ChiRotDefs[i].atom2->z - ChiRotDefs[i].atom3->z;

			for(AtomA = ChiRotDefs[i].startAtomIndex; AtomA <= ChiRotDefs[i].endAtomIndex; AtomA++)
			{
				Delta2X = getWSpace().cur.atom[AtomA].p.x - ChiRotDefs[i].atom2->x; //vector C to AtomA
				Delta2Y = getWSpace().cur.atom[AtomA].p.y - ChiRotDefs[i].atom2->y;
				Delta2Z = getWSpace().cur.atom[AtomA].p.z - ChiRotDefs[i].atom2->z;
				NewDeltaX = VectorY * Delta2Z - VectorZ * Delta2Y; //calculate cross vector product
				NewDeltaY = VectorZ * Delta2X - VectorX * Delta2Z;
				NewDeltaZ = VectorX * Delta2Y - VectorY * Delta2X;
				ChiRotDefs[i].Value -= (
					NewDeltaX * getWSpace().cur.atom[AtomA].f.x + 
					NewDeltaY * getWSpace().cur.atom[AtomA].f.y + 
					NewDeltaZ * getWSpace().cur.atom[AtomA].f.z)
					* ChiRotDefs[i].ReciprocalLen; //dot product evaluates move
			} //next AtomA
		} //next i
	} //end procedure

	void TorsionalMinimisation::ActionConjugateGradient()// adjust rotation values by conjugate gradient
	{
		double Square, GradRatio, Value, MaxSoFar, Multiplier;

		Square = 0;
		MaxSoFar = 0;
		for( size_t i = 0; i < PhiPsiRotDefs.size(); i++ )
		{
			Square += PhiPsiRotDefs[i].Value * PhiPsiRotDefs[i].Value; // PHI & PSI
		}

		for( size_t i = 0; i < ChiRotDefs.size(); i++ )
		{
			Square += ChiRotDefs[i].Value * ChiRotDefs[i].Value; // PHI & PSI
		}

		GradRatio = Square / (PrevSquare * CapRise); // CapRise assumes StepSize already mul'CapRise' cos of overstep
		PrevSquare = Square;

		// PHI, PSI Defs
		for( size_t i = 0; i < PhiPsiRotDefs.size(); i++ )
		{
			PhiPsiRotDefs[i].Value += GradRatio * PhiPsiRotDefs[i].PrevValue; //actual move is this + ratio*last
			PhiPsiRotDefs[i].PrevValue = PhiPsiRotDefs[i].Value;
			Value = fabs(PhiPsiRotDefs[i].Value);
			if(Value > MaxSoFar) MaxSoFar = Value;
		}
		Multiplier = AngleCap / MaxSoFar;
		for( size_t i = 0; i < PhiPsiRotDefs.size(); i++ )
		{
			PhiPsiRotDefs[i].Value *= Multiplier;
			PhiPsiRotDefs[i].rotateByValue();
		}

		// CHI Defs
		for( size_t i = 0; i < ChiRotDefs.size(); i++ )
		{
			ChiRotDefs[i].Value += GradRatio * ChiRotDefs[i].PrevValue; //actual move is this + ratio*last
			ChiRotDefs[i].PrevValue = ChiRotDefs[i].Value;
			ChiRotDefs[i].Value *= SCAngleFactor;
			if( ChiRotDefs[i].Value > SCAngleCap ) ChiRotDefs[i].Value = SCAngleCap;
			else if( ChiRotDefs[i].Value < -SCAngleCap ) ChiRotDefs[i].Value = -SCAngleCap;
			ChiRotDefs[i].rotateByValue();
		} //end procedure
	}

	void TorsionalMinimisation::info() const
	{
		printf("tmin.UpdateScr: %d\n", UpdateScr);
		printf("tmin.UpdateTra: %d\n", UpdateTra);
		printf("tmin.Steps %d\n", Steps );
		printf("tmin.InitialCapFactor %e\n", InitialCapFactor );
		printf("tmin.SlopeCutoff %e\n", SlopeCutoff);
		printf("tmin.CapRise %e\n", CapRise );
		printf("tmin.CapDrop %e\n", CapDrop );
		printf("tmin.IntendedMaxAngleCap %e\n", IntendedMaxAngleCap );
		printf("tmin.IntendedMaxSCAngleCap %e\n", IntendedMaxSCAngleCap );
	}

	void TorsionalMinimisation::infoLine() const 
	{
		printf("%8d\t%7.6lf\t%7.6lf\t%6.2lf\t%6.2lf\t%15.10lf\t", Step, AngleCap, SCAngleCap, getWSpace().ene.cRMS, getWSpace().ene.dRMS, double(getWSpace().ene.epot) * PhysicsConst::J2kcal * PhysicsConst::Na);
		mon.printCurData();
		ff->infoLine();
		printf("\n");
	}

	void TorsionalMinimisation::infoLineHeader() const 
	{
		printf("%8s\t%6s\t%7s\t%7s\t%7s\t%10s\t", "Step#", "AngCap", "SCAngCap", "CRMS", "dRMS", "Epot");
		mon.printHeader();
		ff->infoLineHeader();
		printf("\n");
	}

	// -----------------------------------------------------------------------------------------------------------
	// End: 'TorsionalMinimisation'
	// -----------------------------------------------------------------------------------------------------------
}

