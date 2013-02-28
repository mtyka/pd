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
#include "tools/vector.h"
#include "system/molecule.h"
#include "genpolymer.h"
#include "rebuilder.h"

using namespace Maths;

int findLevelOne( int i, const ParticleStore& atom, std::vector < int >& AlignAtom_L1, int iTemplateRangeStart, int iTemplateRangeEnd )
{
	if( iTemplateRangeStart == -1 ) iTemplateRangeStart = 0;
	if( iTemplateRangeEnd == -1 ) iTemplateRangeEnd = atom.size() - 1;

	D_ASSERT( iTemplateRangeStart >= 0 && iTemplateRangeStart < atom.size(), OutOfRangeException, "findLevelOne()" );
	D_ASSERT( iTemplateRangeEnd >= 0 && iTemplateRangeEnd < atom.size(), OutOfRangeException, "findLevelOne()" );
	D_ASSERT( iTemplateRangeStart <= iTemplateRangeEnd, ArgumentException, "findLevelOne() NOT iTemplateRangeStart <= iTemplateRangeEnd!" );

	AlignAtom_L1.clear();
	int j;
	for(int nj = 0; nj < atom[i].cov12atom.size(); nj++)
	{
		j = atom[i].cov12atom[nj].i;

		if( j < iTemplateRangeStart || j > iTemplateRangeEnd )
			continue;

		if(!atom[j].isRebuildRequired())
		{
			AlignAtom_L1.push_back(j);
			if( AlignAtom_L1.size() == 3 ) break;
		}
	}
	return AlignAtom_L1.size();
}

int findLevelTwo( const ParticleStore & atom, std::vector < int >& AlignAtom_L1, std::vector < int >& AlignAtom_L2, int iTemplateRangeStart, int iTemplateRangeEnd )
{
	if( iTemplateRangeStart == -1 ) iTemplateRangeStart = 0;
	if( iTemplateRangeEnd == -1 ) iTemplateRangeEnd = atom.size() - 1;

	D_ASSERT( iTemplateRangeStart >= 0 && iTemplateRangeStart < atom.size(), OutOfRangeException, "findLevelTwo()" );
	D_ASSERT( iTemplateRangeEnd >= 0 && iTemplateRangeEnd < atom.size(), OutOfRangeException, "findLevelTwo()" );
	D_ASSERT( iTemplateRangeStart <= iTemplateRangeEnd, ArgumentException, "findLevelTwo() NOT iTemplateRangeStart <= iTemplateRangeEnd!" );

	AlignAtom_L2.clear();
	// now go through all the neighbours of the atoms found
	for(size_t nj = 0; nj < AlignAtom_L1.size(); nj++)
	{
		int j = AlignAtom_L1[nj];

		if( j < iTemplateRangeStart || j > iTemplateRangeEnd )
			continue;

		for(size_t nk = 0; nk < atom[j].cov12atom.size(); nk++)
		{
			int k = atom[j].cov12atom[nk].i;

			if( k < iTemplateRangeStart || k > iTemplateRangeEnd )
				continue;

			if(!atom[k].isRebuildRequired() && !VectorContains(AlignAtom_L2,k))
			{
				AlignAtom_L2.push_back(k);
				if( AlignAtom_L1.size() + AlignAtom_L2.size() == 3 )
				{
					return AlignAtom_L2.size();
				}
			}
		}
	}
	return AlignAtom_L2.size();
}

bool findLevelThree( const ParticleStore & atom, std::vector < int >& AlignAtom_L1, std::vector < int >& AlignAtom_L2, int& AlignAtom_L3, int iTemplateRangeStart, int iTemplateRangeEnd )
{
	if( iTemplateRangeStart == -1 ) iTemplateRangeStart = 0;
	if( iTemplateRangeEnd == -1 ) iTemplateRangeEnd = atom.size() - 1;

	D_ASSERT( iTemplateRangeStart >= 0 && iTemplateRangeStart < atom.size(), OutOfRangeException, "findLevelTwo()" );
	D_ASSERT( iTemplateRangeEnd >= 0 && iTemplateRangeEnd < atom.size(), OutOfRangeException, "findLevelTwo()" );
	D_ASSERT( iTemplateRangeStart <= iTemplateRangeEnd, ArgumentException, "findLevelTwo() NOT iTemplateRangeStart <= iTemplateRangeEnd!" );

	// now go through all the neighbours of the atoms found
	for(size_t nj = 0; nj < AlignAtom_L2.size(); nj++)
	{
		int j = AlignAtom_L2[nj];

		if( j < iTemplateRangeStart || j > iTemplateRangeEnd )
			continue;

		for(size_t nk = 0; nk < atom[j].cov12atom.size(); nk++)
		{
			int k = atom[j].cov12atom[nk].i;

			if( k < iTemplateRangeStart || k > iTemplateRangeEnd )
				continue;

			if(!atom[k].isRebuildRequired() && !VectorContains(AlignAtom_L1,k) && !VectorContains(AlignAtom_L2,k))
			{
				AlignAtom_L3 = k;
				return true;
			}
		}
	}
	return false;
}

void rebuild( int i, ParticleStore & atom, int r0, int r1, int r2 )
{
	matrix3x3 rmat;

	superimpose(
		atom[r0].pos(),
		atom[r1].pos(),
		atom[r2].pos(),
		atom[r0].posGeom(),
		atom[r1].posGeom(),
		atom[r2].posGeom(),
		rmat);

	atom[i].pos().setTo(atom[i].posGeom());
	atom[i].pos().sub(atom[r0].posGeom());
	atom[i].pos().mulmat(rmat);
	atom[i].pos().add(atom[r0].pos());

	// make sure the new coordinates are numerically 'real':
	if(!atom[i].pos().isReal())
	{
		printf("Rebuild FAILED for: ");
		atom[i].info();

		printf("Using template atoms:\n");
		atom[r0].info();
		atom[r1].info();
		atom[r2].info();

		printf("Reference Positions:\n");
		atom[r0].posGeom().info();
		atom[r1].posGeom().info();
		atom[r2].posGeom().info();

		printf("Atom Positions:\n");
		atom[r0].pos().info();
		atom[r1].pos().info();
		atom[r2].pos().info();

		THROW(ProcedureException,"Maths Failure, non-real-number solution!");
	}

	atom[i].setRebuildRequired( false ); // It's rebuilt!
}

void rebuildReport( const Particle& atom, int& hydrogenRebuild, Verbosity::Type verbose )
{
	if( atom.isHydrogen() )
	{
		hydrogenRebuild++;
	}
	else
	{
		if( verbose )
		{
			printf("   Building: ");
			atom.info();
		}
	}
}

bool rebuildInnerLoop( MoleculeBase& mol, int errcount, size_t iBuildRangeStart, size_t iBuildRangeEnd, Verbosity::Type verbose, int iTemplateRangeStart = -1, int iTemplateRangeEnd = -1 )
{
	if( verbose ) 
		printf("Entering rebuild loop:\n");

	ParticleStore &atom = mol.atom;

	// if rebuild is requested and there were unidentified atoms, attempt rebuild
	std::vector < int > AlignAtom_L1;
	std::vector < int > AlignAtom_L2;
	int AlignAtom_L3 = -1;
	while(errcount > 0)
	{
		int hydrogenRebuild = 0;

		for(size_t i = iBuildRangeStart; i <= iBuildRangeEnd; i++)
		{
			if(!atom[i].isRebuildRequired())
				continue; // skip atoms that have been allocated by another source

			// find three atoms being as close as possible

			// start with 1-2 atoms
			int found_L1 = findLevelOne( i, atom, AlignAtom_L1, iTemplateRangeStart, iTemplateRangeEnd );
			if( found_L1 == 3 )
			{
				rebuildReport( atom[i], hydrogenRebuild, verbose );
				rebuild( i, atom, AlignAtom_L1[0], AlignAtom_L1[1], AlignAtom_L1[2] );
				continue;
			}
			else if( found_L1 == 0 )
			{
				continue;
			}

			// Now have a look at 1-3 atoms
			int found_L2 = findLevelTwo( atom, AlignAtom_L1, AlignAtom_L2, iTemplateRangeStart, iTemplateRangeEnd );
			if( found_L2 == 0 )
			{
				continue;
			}
			else if( found_L2 == 1 && found_L1 == 2 )
			{
				rebuildReport( atom[i], hydrogenRebuild, verbose );
				rebuild( i, atom, AlignAtom_L1[0], AlignAtom_L1[1], AlignAtom_L2[0] );
				continue;
			}
			else if( found_L2 == 2 )
			{
				rebuildReport( atom[i], hydrogenRebuild, verbose );
				rebuild( i, atom, AlignAtom_L1[0], AlignAtom_L2[0], AlignAtom_L2[1] );
				continue;
			}

			// You can only use level three partners if the relationship is not dependent on the rotation
			// of the torsion between BC (A -- B -- C -- D). We have to assume that if rebuild is not possible
			// in the above atom-finding calls, that it must be called here.
			if(findLevelThree( atom, AlignAtom_L1, AlignAtom_L2, AlignAtom_L3, iTemplateRangeStart, iTemplateRangeEnd ))
			{
				rebuildReport( atom[i], hydrogenRebuild, verbose );
				rebuild( i, atom, AlignAtom_L1[0], AlignAtom_L2[0], AlignAtom_L3 );
				continue;
			}
		}

		// now recount the errors
		int errcount2 = 0;
		for(size_t i = iBuildRangeStart; i <= iBuildRangeEnd; i++)
		{
			if(atom[i].isRebuildRequired())
				errcount2++;
		}
		if(errcount == errcount2)
		{
			if( verbose ) 
				printf(" No further atoms could be rebuilt - aborting - unresolved missing atoms remain!\n");
			for(size_t i = iBuildRangeStart; i <= iBuildRangeEnd; i++)
			{
				if(atom[i].isRebuildRequired())
				{
					if( verbose ) 
					{
						printf("\tCannot build: ");
						atom[i].info();
					}
				}
			}
			return false; // resign when no further improvement is reached
		}

		int tot = errcount - errcount2;
		if( verbose ) 	
			printf(" Rebuilt %d atoms (%d Heavy and %d Hydrogens)\n", tot, tot-hydrogenRebuild, hydrogenRebuild);
		errcount = errcount2;	// go round again if there was an improvement
	}

	return true;
}

bool rebuildMissingAtoms( MoleculeBase& mol, int ir, Verbosity::Type verbose )
{
	ParticleStore &atom = mol.atom;

	int iStart = mol.res[ir].ifirst;
	int iLast = mol.res[ir].ilast;

	int missing = 0;
	for(size_t j = iStart; j <= iLast; j++)
	{
		if(atom[j].isRebuildRequired())
		{
			missing++;
		}
	}

	if( missing == 0 )
	{
		return true;
	}
	else
	{
		return rebuildInnerLoop( mol, missing, iStart, iLast, verbose );
	}
}

bool rebuildMissingAtoms( MoleculeBase& mol, Verbosity::Type verbose )
{
	ParticleStore & atom = mol.atom;
	const ResidueStore& res = mol.res;
	int errcount = 0;

	// Print a report if we are being verbose ...
	// At the same time count what there is to be rebuilt
	for(size_t i = 0; i < res.size(); i++)
	{
		size_t jStart = res[i].ifirst;
		size_t jEnd =  res[i].ilast;
		size_t jLength = res[i].ilast - res[i].ifirst + 1;
        if( jLength <= 3 && !res[i].param->hasBackLink() && !res[i].param->hasFrwdLink() )
		{
			// We will get into this code path if we are something like water
			// very small and lacking forward and back links

			static Maths::dvector cog; // centre of mass
			cog.setTo(0,0,0);
			int validCount = 0;
			for(size_t j = jStart; j <= jEnd; j++)
			{
				if(!atom[j].isRebuildRequired())
				{
					validCount++;
					cog.add(mol.atomxyz(j));
				}
			}		
			if( validCount == jLength )
			{
				continue;
			}
			else if ( validCount == 0 )
			{
				printf("Cannot rebuild residue %d. No positional information is present!",i);
				return false;
			}
			else
			{
				cog.div((double)validCount);
				// Superimpose the centre of the geoPos onto the current centre of mass!
				static Maths::dvector cogGeo; // centre of mass
				cogGeo.setTo(0,0,0);
				for(size_t j = jStart; j <= jEnd; j++)
				{
					cogGeo.add(atom[j].posGeom());
				}
				cogGeo.div((double)jLength);
				for(size_t j = jStart; j <= jEnd; j++)
				{
					mol.atomxyz(j).setTo(atom[j].posGeom()); // set to the stored position of geometry
					mol.atomxyz(j).sub(cogGeo); // place small molecule centre on the origin
					mol.atomxyz(j).add(cog); // push the atom back out where its supposed to be (i.e. where the original valid positions lay.
				}
			}
		}
		else
		{
			int missing = 0;
			for(size_t j = jStart; j <= jEnd; j++)
			{
				if(atom[j].isRebuildRequired())
				{
					missing++;
				}
			}
			if( missing == 0 )
			{
				continue;
			}
			else if( missing == jLength )
			{
				if( errcount == 0 )
				{
					printf("The following things will be rebuilt:\n");
				}
				errcount += missing;
				if( verbose )
				{
					printf("   Whole Residue: ");
					res[i].info();
				}
			}
			else
			{
				for(size_t j = jStart; j <= jEnd; j++)
				{
					if(atom[j].isRebuildRequired())
					{
						if( errcount == 0 )
						{
							printf("The following things will be rebuilt:\n");
						}
						errcount++;
						if( verbose )
						{
							if( atom[j].cov12atom.size() > 1 )
							{
								printf("   Heavy-atom: ");
								atom[j].info();
							}
							else
							{
								// Again, here we only want to print if the valency > 1 as otherwise we get loads of info about
								// rebuilding missing 'H' atoms. Much more interesting to report missing heavy atoms!
								// H atom printing should maybe be included in a higher order Verbosity::Type setting - at the mo, we only have
								// true and false...
							}
						}
					}
				}
			}
		}
	}

	if( errcount == 0 )
	{
		return true;
	}
	else
	{
		return rebuildInnerLoop( mol, errcount, 0, atom.size()-1, verbose );
	}
}

bool rebuildAndPolymerise(MoleculeBase& mol, int start, int end)
{
	int sectionStart = mol.res[start].ifirst;
	int sectionEnd = mol.res[end].ilast;

	// Reassign our residues positions back to their forcefield defaults.
	// This should be quick and painless :-D
	for( int i = start; i <= end; i++ )
	{
		int atomStart = mol.res[i].ifirst;
		int atomEnd = mol.res[i].ilast;
		int atomCount =  atomEnd - atomStart + 1;
		const MoleculeDefinition* molDef = mol.res[i].param;
		ASSERT( molDef->atom.size() == atomCount, CodeException, "MoleculeBase lengths oddly dont match!");
		for( size_t j = 0; j < atomCount; j++ )
		{
			mol.atomxyz(atomStart+j).setTo( molDef->atom[j].pos() );
		}
	}

	for(int ir = end; ir > start; ir--)
	{
		const MoleculeDefinition* resdef = mol.res[ir].param;
		const MoleculeDefinition* resdefPrev = mol.res[ir-1].param;

		int iBackAtom = mol.findParticle(ir-1,resdef->backName);
		int iFrwdAtom = mol.findParticle(ir,resdefPrev->frwdName);

		if(iBackAtom < 0)
		{
			printf("\nERROR: Cannot find BackLink Atom in previous residue of residue %d\n", ir);
			return false;
		}
		if(iFrwdAtom < 0)
		{
			printf("\nERROR: Cannot find Forward Link Atom in next residue of residue %d\n", ir);
			return false;
		}

		dvector resoffset;
		matrix3x3 resrot;

		superimpose(
			mol.atomxyz(iBackAtom),
			resdefPrev->frwdpos,
			resdef->backpos,
			mol.atomxyz(iFrwdAtom),
			resrot);

		dvector subt(resdef->backpos);
		dvector addt(mol.atomxyz(iBackAtom));

		for(size_t i = sectionStart; i <= sectionEnd; i++)
		{
			if(mol.atom[i].ir < ir) continue;
			mol.atomxyz(i).sub(subt);
			mol.atomxyz(i).mulmat(resrot);
			mol.atomxyz(i).add(addt);
		}

		dvector backpos;
		backpos = resdef->backpos;
		backpos.sub(subt);
		backpos.mulmat(resrot);
		backpos.add(addt);
	}

	for(size_t i = sectionStart; i <= sectionEnd; i++)
	{
		mol.atom[i].setRebuildRequired(false); // we have rebuilt the atom
	}

	return true;
}

bool ToOriginRebuilder::invokeBuild(MoleculeBase& mol, Verbosity::Type outLevel )
{
	if( outLevel >= Verbosity::Normal ) 
		std::cout << "Zeroing the coordinates of atoms requiring rebuild." << std::endl;
	// Place all unknown atoms on the cartesian origin.
	for( size_t i = 0; i < mol.nAtoms(); i++ )
	{
		if( mol.atom[i].isRebuildRequired() )
		{
			if( m_FlagBuilt ) mol.atom[i].setRebuildRequired(false);
			mol.atom[i].pos().setTo(0,0,0);
		}
	}
	return true;
}

bool MissingAtomRebuilder::invokeBuild(MoleculeBase& mol, Verbosity::Type verbose)
{
	return rebuildMissingAtoms(mol,verbose);
}

bool MainchainHydrogenRebuilder::invokeBuild(MoleculeBase& mol, Verbosity::Type verbose)
{
	return invokeBuild( mol, verbose, PickEverything() );
}

bool MainchainHydrogenRebuilder::invokeBuild(MoleculeBase& mol, Verbosity::Type verbose, const PickResidueBase& _picker, const PickAtomRange& _TemplateRange )
{
	return invokeBuild( mol, verbose, _picker, _TemplateRange.getStartAtomIndex(), _TemplateRange.getEndAtomIndex() );
}

bool MainchainHydrogenRebuilder::invokeBuild(MoleculeBase& mol, Verbosity::Type verbose, const PickResidueBase& _picker, int iTemplateRangeStart, int iTemplateRangeEnd )
{
	for( size_t i = 0; i < mol.nResidues(); i++ )
	{
		if( _picker.matches( mol.res[i] ) )
		{
			const Residue& res = mol.res[i];
			for( int j = res.ifirst; j <= res.ilast; j++ )
			{
				Particle& part = mol.atom[j];
				if( part.isBackbone() && part.isHydrogen() )
				{
					part.setRebuildRequired(true);
					if( !rebuildInnerLoop( mol, 1, j, j, verbose, iTemplateRangeStart, iTemplateRangeEnd ) )
						return false;
				}
			}
		}
	}
	return true;
}

RebuilderOwner::RebuilderOwner() : 
	m_DefaultBuilder(),
	m_DisabledBuilder(false)
{
	setRebuilderDefault();
}

void RebuilderOwner::setRebuilder( RebuilderBase& _rebuilder )
{
	m_Rebuilder = &_rebuilder;
}

void RebuilderOwner::setRebuilderDefault()
{
	m_Rebuilder = &m_DefaultBuilder;
}

void RebuilderOwner::disableRebuild()
{
	m_Rebuilder = &m_DisabledBuilder;
}

RebuilderBase& RebuilderOwner::getRebuilder() 
{ 
	return *m_Rebuilder; 
}


