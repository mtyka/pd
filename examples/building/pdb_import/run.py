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
import sys

ffps = FFParamSet()

ffps.readLib("amber03aa.ff");
ffps.readLib("tip3.ff");
ffps.readLib("default.alias");
ffps.readLib("default.class");

stem = "../../pdb/"

names = [
 	"bba1.pdb",       ## Correctly complains about PYA
 	"water.pdb",      ## Correctly throws an exception - MOL is not a molecule in amber03aa.ff
 	"trpzip1.pdb",    ## DIES !
	"sh3.pdb",        ## Good
 	"acyltra.pdb",    ## DIES !
	"1fsd.pdb",       ## Good
 	"villin.pdb",     ## DIES !
	"protg.pdb",      ## Good
  "proteinAZ.pdb",  ## DIES !
	"ubifrag.pdb",    ## DIES ! 
	"trpzip.pdb",     ## Good
	"trpcage.pdb",    ## Good
	"2ci2.pdb"        ## Good
]

## Mass - try loading every PDB

for pdbname in names:
	try:
		sim = PDB_In(ffps,stem + pdbname);
		sim.load();                          ## Load everything from model 1
	except:
		sys.stderr.write( "Failed on PDB " + pdbname + " !!! \n"  )
		continue
	
	## Write System to see what we've loaded from the PDB
	sim.info();
	sim.save( OutputFile_PDB( "output_" + pdbname ) )



