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
#include "system/molecule.h"
#include "workspace/workspace.h"
#include "breakablebonded.h"

namespace Physics
{
	FF_BreakableBonded::FF_BreakableBonded(WorkSpace &newwspace) : 
		FF_Bonded(newwspace)
	{
	}

	void FF_BreakableBonded::createBreak( int _AtomIndexA, int _AtomIndexB )
	{
		if( _AtomIndexA < 0 ) return;
		if( _AtomIndexB < 0 ) return;
		if( _AtomIndexA == _AtomIndexB ) return;			
		m_BrokenBonds.push_back(BreakDef(_AtomIndexA,_AtomIndexB));
		needsetup = true;
	}

	void FF_BreakableBonded::clearBreaks()
	{
		needsetup = true;
		m_BrokenBonds.clear();
	}

	void FF_BreakableBonded::validateBonds()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		if(OutputLevel) printf("Evaluating bond break list...\n");
		for( size_t k = 0; k < m_BrokenBonds.size(); k++ )
		{
			// Is the bond a real bond?
			bool valid = isBonded( m_BrokenBonds[k].i, m_BrokenBonds[k].j );
			m_BrokenBonds[k].valid = valid;

			std::vector<int>& status = m_BrokenBonds[k].atomList;
			status.clear();

			if( valid )
			{
				// Atom proxy
				ParticleStore& atom = wspace.atom;
				
				// Get the shortest list of atoms from one side of the bond. 			
				int top, bottom;
				if( m_BrokenBonds[k].i < m_BrokenBonds[k].j)
				{
					bottom = m_BrokenBonds[k].i;
					top = m_BrokenBonds[k].j;
				}
				else
				{
					bottom = m_BrokenBonds[k].j;
					top = m_BrokenBonds[k].i;
				}				

				bool secondRun = false;
				while(true)
				{
					int ignore;
					if( !secondRun )
					{
						// begin by looking at atom indexes larger than that of the bonds largest atom index (top)
						ignore = bottom;
						status.push_back(top);
						for( size_t q = 0; q < atom[top].cov12atom.size(); q++ )
						{
							int index = atom[top].cov12atom[q].i;
							if( ignore == index ) continue; // don't cross the bond itself!
							status.push_back(index);
						}
					}
					else
					{
						// If looking that way from the bond includes more than half
						// the total atoms, look the other way and try again!
						status.clear(); // reset the list for the next run!
						ignore = top;
						status.push_back(bottom);
						for( size_t q = 0; q < atom[bottom].cov12atom.size(); q++ )
						{
							int index = atom[bottom].cov12atom[q].i;
							if( ignore == index ) continue; // don't cross the bond itself!
							status.push_back(index);
						}
					}

					int curindex = 0;
					while( curindex++ < status.size()-1 )
					{
						int search = status[curindex];
						for( size_t q = 0; q < atom[search].cov12atom.size(); q++ )
						{
							int index = atom[search].cov12atom[q].i;
							if( index == ignore ) 
							{
								// Our bond is in a ring structure!!	
								if(OutputLevel) 
								{
									printf("WARNING: Bond-break defined within a ring structure! Between: Atoms %d and %d\n",
										m_BrokenBonds[k].i, 
										m_BrokenBonds[k].j);
									printf("Could be within chemical ring or a disulphide linked region.\n");
									printf("Side effects are possible depending upon angle and improper definitions within the forcefield!\n");
								}

								// Bugger ...
								// We now have to define some sensible atoms for the atoms on this "side" of the bond
								// - we now **have** to assume that any angles and torsions will be defined 
								// accross the bond itself and with one of the first shell of covelent partners
								
								status.clear();
								status.push_back(top);
								for( size_t q = 0; q < atom[top].cov12atom.size(); q++ )
								{
									int index = atom[top].cov12atom[q].i;
									if( ignore == index ) continue; // don't cross the bond itself!
									status.push_back(index);
								}
						
								goto NEXTBOND; // Ewww :-D - but the only decent way to break out of several loops.
							}
							if( !VectorContains(status,index) ) status.push_back(index);						
						}							
					}

					if( secondRun )
					{
						// We only have two runs, one for each side
						break;
					}
					else if(status.size() < (atom.size() / 2))
					{
						// We have found the shortest side in the first pass :-D
						break;
					}
					else
					{
						// In some rare cases, a second run will be required
						secondRun = true;
					}
				}

				if(OutputLevel) printf("Breaking bond between: Atoms %d and %d:\n",m_BrokenBonds[k].i, m_BrokenBonds[k].j);
			}
			else
			{
				if(OutputLevel) printf("Invalid bond exclustion: Atoms %d and %d are not bonded.\n",m_BrokenBonds[k].i, m_BrokenBonds[k].j);
			}

			NEXTBOND: // Need this to allow a 'goto' breaker above

			// now sort our array, because we will be using std::binary_search below
			std::sort(status.begin(),status.end());
		}
	}

	void FF_BreakableBonded::setup()
	{
		FF_Bonded::setup();

		validateBonds();

		if( m_BrokenBonds.size() > 0 )
		{
			// Now remove the relevent FF_Bonded terms from the base-class generated list!
			// It is easier and less error prone to remove them post-setup rather than try and re-write
			// the setup functions or alter the base class in any way.

			// NOTE: The '(bool (FF_BreakableBonded::*)(Bond) const)' is an explicit cast, and only required because
			// until recently (when I found it ;-)) there was a compiler bug in the Intel compiler! They have fixed it 
			// in the latest release, however this hack will be kept as the code works on all compilers, including the
			// older intel compilers.

			removematches(bond, (bool (FF_BreakableBonded::*)(Bond) const) &FF_BreakableBonded::shouldRemove, this );
			removematches(angle, (bool (FF_BreakableBonded::*)(Angle) const) &FF_BreakableBonded::shouldRemove, this );
			removematches(torsion,(bool (FF_BreakableBonded::*)(Torsion) const) &FF_BreakableBonded::shouldRemove, this );
			removematches(improper,(bool (FF_BreakableBonded::*)(Torsion) const) &FF_BreakableBonded::shouldRemove, this );
		}
	}

	bool FF_BreakableBonded::shouldRemove( const Physics::Bond bond ) const
	{
		ParticleStore& atom = getWSpace().atom;
		for( int k = 0; k < m_BrokenBonds.size(); k++ )
		{
			if( !m_BrokenBonds[k].valid ) continue;
			if( ( bond.i == m_BrokenBonds[k].i && bond.j == m_BrokenBonds[k].j ) || 
				( bond.j == m_BrokenBonds[k].i && bond.i == m_BrokenBonds[k].j ) )
			{
				if(OutputLevel) 
				{
					printf("   Del bond    (");
					atom[bond.i].info(OutputLevel,0,false);
					printf(", ");
					atom[bond.j].info(OutputLevel,0,false);
					printf(")\n");
				}
				return true;
			}
		}
		return false;
	}

	bool FF_BreakableBonded::shouldRemove( const Physics::Angle angle ) const
	{
		ParticleStore& atom = getWSpace().atom;
		// Remove all angle terms that span accross the bond in question
		for( int k = 0; k < m_BrokenBonds.size(); k++ )
		{
			if( !m_BrokenBonds[k].valid ) continue;

			// Find if all the atoms are on the same side of the bond
			// ***NOTE, you can only perform binary_search on a sorted container - this
			// is currently done when they are first generated.***
			const std::vector<int>& vect = m_BrokenBonds[k].atomList;
			
			// Identify which side of the fence each angle atom lies on
			bool iUpStream = std::binary_search( vect.begin(), vect.end(), angle.i );
			bool jUpStream = std::binary_search( vect.begin(), vect.end(), angle.j );
			bool aUpStream = std::binary_search( vect.begin(), vect.end(), angle.a );
			
			if( iUpStream != jUpStream || 
				iUpStream != aUpStream )
			{
				if(OutputLevel) 
				{
					printf("   Del angle   (");
					atom[angle.i].info(Verbosity::Normal,0,false);
					printf(", ");
					atom[angle.j].info(Verbosity::Normal,0,false);
					printf(", ");
					atom[angle.a].info(Verbosity::Normal,0,false);
					printf(")\n");
				}
				return true;
			}			
		}
		return false;
	}

	bool FF_BreakableBonded::shouldRemove( const Physics::Torsion torsion ) const
	{
		// Remove all torsion terms that span accross the bond in question
		// This should work for both torsions and imporpers, even through their 
		// 'shape definitions' are different
		ParticleStore& atom = getWSpace().atom;
		for( int k = 0; k < m_BrokenBonds.size(); k++ )
		{
			if( !m_BrokenBonds[k].valid ) continue;

			// Find if all the atoms are on the same side of the bond
			// NOTE, you can only perform binary_search on a sorted container - this
			// is currently done when they are first generated.
			const std::vector<int>& vect = m_BrokenBonds[k].atomList;

			// Identify which side of the fence each torsion atom lies on
			bool iUpStream = std::binary_search( vect.begin(), vect.end(), torsion.i );
			bool jUpStream = std::binary_search( vect.begin(), vect.end(), torsion.j );
			bool aUpStream = std::binary_search( vect.begin(), vect.end(), torsion.a );
			bool bUpStream = std::binary_search( vect.begin(), vect.end(), torsion.b );

			if( iUpStream != jUpStream || 
				iUpStream != aUpStream || 
				iUpStream != bUpStream )
			{
				if(OutputLevel) 
				{
					printf("   Del torsion (");
					atom[torsion.a].info(Verbosity::Normal,0,false);
					printf(", ");
					atom[torsion.i].info(Verbosity::Normal,0,false);
					printf(", ");
					atom[torsion.j].info(Verbosity::Normal,0,false);
					printf(", ");
					atom[torsion.b].info(Verbosity::Normal,0,false);
					printf(")\n");
				}
				return true;
			}			
		}
		return false;
	}
}

