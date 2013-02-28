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

ffps = FFParamSet()

ffps.readLib("amber03aa.ff");
ffps.readLib("tip3.ff");
ffps.readLib("bombykol.ff");

pdbname = "../../pdb/1dqe.pdb"



# Demonstrating Content filters and specific loading of PDB contents 

sim = PDB_In(ffps,pdbname);      # specify the PDB file to be loaded
sim.setFilter( Polypeptide );       # set the filtering mode - PolyPeptide is defined in the Forcefield file
sim.load();                         # scan and load the PDB
sim.info();                         # print some summarising information info
sim.save( OutputFile_PDB( "loaded_polypeptide.pdb" ) )

sim = PDB_In(ffps,pdbname);      # specify the PDB file to be loaded
sim.setFilter( Water );             # Only load the water (we loaded tip3.ff earlier, so that will be used to model the water)
sim.load();                         # scan and load the PDB
sim.info();                         # print some summarising information info
sim.save( OutputFile_PDB( "loaded_water.pdb" ) )

sim = PDB_In(ffps,pdbname);      # specify the PDB file to be loaded
sim.setFilter( Water | Polypeptide );   # Load protein *and* water  (we loaded tip3.ff earlier, so that will be used to model the water)
sim.load();                         # scan and load the PDB
sim.info();                         # print some summarising information info
sim.save( OutputFile_PDB( "loaded_water_and_protein.pdb" ) )

## This unfortunately loads both molecules ! How do i load just the first one ??? I Wanna be able to specify the exact residue number
sim = PDB_In(ffps,pdbname);      # specify the PDB file to be loaded
sim.setFilter( SmallMolecule );   # Load protein *and* water  (we loaded tip3.ff earlier, so that will be used to model the water)
sim.load();                         # scan and load the PDB
sim.info();                         # print some summarising information info
sim.save( OutputFile_PDB( "loaded_small.pdb" ) )


sim = PDB_In(ffps,pdbname);      # specify the PDB file to be loaded
sim.loadAll();                      # scan and load the PDB
sim.info();                         # print some summarising information info
sim.save( OutputFile_PDB( "loaded_all.pdb" ) )

sim = PDB_In(ffps,pdbname);      # specify the PDB file to be loaded
sim.load('A');                      # scan and load the PDB
sim.info();                         # print some summarising information info
sim.save( OutputFile_PDB( "loaded_chain_A_only.pdb" ) )

#sim = PDB_In(ffps,pdbname);      # specify the PDB file to be loaded
#sim.loadExplicit('A',"BOM",47,' '); # Load BOM 47 with no icode from chain A of model 1
#sim.info();                         # print some summarising information info
#sim.save( OutputFile_PDB( "loaded_specific_molecule.pdb" ) )


sim = PDB_In(ffps,pdbname);      # specify the PDB file to be loaded
# create a custom new class "custom1" and add residues called ALA, GLU and ARG to it.
simclass = sim.getClass()
simclass.addClassMember("custom1","ALA"); # ALA is a member of class 'custom1'
simclass.addClassMember("custom1","GLU"); # GLU is a member of class 'custom1'
simclass.addClassMember("custom1","ARG"); # ARG is a member of class 'custom1'
sim.setFilter("custom1");                       # this loads only from the custom gilter
sim.load();                         # scan and load the PDB
sim.info();                         # print some summarising information info
sim.save( OutputFile_PDB( "loaded_custom.pdb" ) )

