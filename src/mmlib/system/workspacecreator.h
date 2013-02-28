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

#ifndef __WORKSPACE_CREATOR_H
#define __WORKSPACE_CREATOR_H

#include "workspace/workspace.fwd.h"
namespace Sequence
{
	class PD_API BioSequence;
}

class PD_API System;
class PD_API FFParamSet;






//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
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
class PD_API WorkspaceCreatorBase
{
public:
	friend class PD_API WorkSpace;
	WorkspaceCreatorBase(System &_system);

	// Public accessors
	const FFParamSet &ffps() const;

protected:
	Sequence::BioSequence& WSCall_GetSeq(WorkSpace &wspace) const;
	void WSCall_Allocate(WorkSpace &wspace, int _natom) const;
	void WSCall_ReInitialiseAll(WorkSpace &wspace) const;
	void WSCall_Append(WorkSpace &wspace, const System& sysspec) const;

	void create(WorkSpace &wspace) const;
	virtual void allocateParticles(WorkSpace &wspace) const = 0;

	void setAtomPositionRedirection(WorkSpace &wspace) const;

	System *system;
};


//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
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
class PD_API WSfull: public WorkspaceCreatorBase
{
public:
	friend class PD_API WorkSpace;
	WSfull(System &system):WorkspaceCreatorBase(system) {}

protected:
	virtual void allocateParticles(WorkSpace &wspace) const;
};

#endif

