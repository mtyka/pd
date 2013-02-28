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
#include "forcefields/forcefield.h"
#include "workspace/snapshot.h"
#include "protocols/md.h"
#include "energy.h"

using namespace Physics;
using namespace Maths;

namespace Protocol
{
	int Energy::runcore()
	{
		if(ensureFFSetup()!=0) return -1;
		refreshNeighborList();
		ff->calcEnergies();
		return 1;
	}
/*
	int ReadPickles::runcore(){

		size_t is=0;
		
		FILE *file;
		file = fopen(filename.c_str(),"r");
		if(file==NULL){ printf("Cannot find file \n"); return -1; }

		printf("Loading all pickles and saving in trajectory \n");

		while(!feof(file)){
			if((is>0)&&((is%100)==0)) printf("%d...",is);
			SnapShot psp;
			psp.readMIME(file);
			getWSpace().load(psp);
			getWSpace().outtra.append();
			is++;
		}
		printf("Read %d snapshots \n",is);
		fclose(file);
		return 0;
	}
*/

} // namespace 'Protocol'


