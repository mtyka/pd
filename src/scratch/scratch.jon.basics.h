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

#ifndef __SCRACTH_JON_BASICS_H
#define __SCRACTH_JON_BASICS_H

#include "global.h"
#include "forcefields/forcefield.h"
#include "protocols/montecarlo.h"
#include "workspace/workspace.fwd.h"

Physics::Forcefield createffs( WorkSpace& wspace, bool useBreakableFF = false, bool summary = true);
Physics::Forcefield createffts( WorkSpace& wspace, bool useBreakableFF = false, bool summary = true);
Physics::Forcefield createffVac( WorkSpace& wspace, bool useBreakableFF = false, bool summary = true);
Physics::Forcefield createff(WorkSpace& wspace, bool useBreakableFF = false, double dielec = 1.0, bool summary = true);

double getMeRMS( const std::vector<Maths::dvector>& native, const std::vector<Maths::dvector>& conformer );

#endif

