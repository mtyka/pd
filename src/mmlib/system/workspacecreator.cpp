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

#include "system/workspacecreator.h"
#include "system/system.h"
#include "workspace/workspace.h"

WorkspaceCreatorBase::WorkspaceCreatorBase(System &_system)
{
	if (_system.nMolecules() == 0)
	{
		throw(ArgumentException("You are trying to create an empty workspace!\nDid you forget to add any molecules to your system?\n"));
	}
	system = &_system;
}

Sequence::BioSequence& WorkspaceCreatorBase::WSCall_GetSeq(WorkSpace &wspace) const
{
	return wspace.m_Sequence;
}

/// Private WorkSpace Accessors
void WorkspaceCreatorBase::WSCall_Allocate(WorkSpace &wspace, int _natom) const
{
	wspace.Allocate(_natom);
}

void WorkspaceCreatorBase::WSCall_Append(WorkSpace &wspace, const System& sysspec) const
{
	wspace.append(sysspec);
}

const FFParamSet &WorkspaceCreatorBase::ffps() const 
{ 
	return system->ffps(); 
}

void WorkspaceCreatorBase::create(WorkSpace &wspace) const
{	
	// first and foremost check that all atoms have been
	// asigned parameters
	printf("Checking if molecules have forcefield parameters attributed ..\n");

	if(system->setup()<0)
	{
		THROW(ProcedureException,"System setup failed");
	}

	printf("Creating WorkSpace ..\n");
	allocateParticles(wspace);

	printf("Initialising workspace components ..\n");
	wspace.reinitAll();

	printf("Done creating system.\n");
}

void WorkspaceCreatorBase::setAtomPositionRedirection(WorkSpace &wspace) const 
{
	// this redirects things such that the wspace.atom[i].pos() actually refers to
	// wspace.cur.atom[i].p
	for(size_t i=0; i< wspace.nAtoms(); i++ ){
		wspace.atom[i].setPosPointer( wspace.cur.atom[i].p );
	}
}

void WSfull::allocateParticles(WorkSpace &wspace) const
{
	// Append copies of all particles from the system to the workspace
	WSCall_Append(wspace,*system);

	// Ensure the correct chainID is maintained - this MUST be called following WSCall_Append()
	for(size_t imol=0;imol<system->nMolecules();imol++)
	{		
		wspace.mol[imol].chainID = system->getMolecule(imol).getChainID();
	}

	// Set the originating indexes
	for(size_t imol=0;imol<system->nMolecules();imol++)
	{
		// Allocate the internal indexes
		for(size_t iatom=0;iatom<system->getMolecule(imol).nAtoms();iatom++)
		{
			wspace.isysmol.push_back(imol);
			wspace.isysatom.push_back(iatom);	
		}
	}

	// Create the 'cur' and 'old' linear arrays
	printf("  Allocating space for particles ..\n");
	WSCall_Allocate(wspace,wspace.atom.size());

	// Set the positions in the cur array
	for(size_t i=0; i<wspace.nAtoms(); i++ )
	{
		wspace.old.atom[i].p.setTo(0,0,0); // initialise
		wspace.cur.atom[i].p.setTo(wspace.atom[i].pos());
	}

	// Now redirect the wspace.atom[i].pos() to actually point to the cur array !
	// the cur array is static in memory (it can't change!) such that 
	// this should be a safe "connection"

	setAtomPositionRedirection(wspace);
}
