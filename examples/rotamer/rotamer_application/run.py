# PD is a free, modular C++ library for biomolecular simulation with a 
# flexible and scriptable Python interface. 
# Copyright (C) 2003-2013 Mike Tyka and Jon Rea
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from pd import *
cseed(4)
info()
timer()

ffps = FFParamSet()

## Use that Charmm forcefield
ffps.readLib( after("-ff","amber03aa.ff")  )

## test simple loading of PDB file
sim = System( ffps );

## Create a Protein helix using evry amino acid
sim.add( NewProteinHelix(ffps,"*A-(CDEFGHIKLMNPQRSTVW)-Y*") );

# create workspace
wspace = WorkSpace( sim )

wspace.printPDB("rotamer_perturb_random_start.pdb");
tra = OutTra_BTF("rotamer_perturb_random_end",wspace)
wspace.addTra(tra)


rot = RotamerLibrary(ffps);

rot.convertLib("../../../param/rotlib/shetty/scl-B30-occ1.0-rmsd1.0-prop20.0",RotLibConvert_Shetty());
#2A
rot.writeLib("shetty.rotamer");

#	}
#	else
#	{
#rot.readLib("rotlib/shetty.rotamer");
#	}

timeMe = StatClock()
timeMe.Begin();

app = RandomRotamerApplicator ( wspace, rot, ApplyCartesian );
##	if( _testSterics )
##		app.addFilter_DefaultSteric(); // Add the linear-time clash filter for all rotamer states

app.test(250);  ## 250 random rotamer applications over the entire PickedResidueList
app.addFilter_DefaultSteric();  ## Add the linear-time clash filter for all rotamer states
app.test(250);  ## 250 random rotamer applications over the entire PickedResidueList

timeMe.End();
timeMe.ReportMilliSeconds();

wspace.printPDB("rotamer_perturb_random_end.pdb");
