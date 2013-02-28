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

pdbname = "../../pdb/protg.pdb"

sim = PDB_In(ffps, pdbname);

###  Ok, so the sequence of ProteinG is MQYKLVINGKTLKGETTTKAVDAETAEKAFKQYANDNGVDGVWTYDDATKTFTVTE
### and we're going to map over it a mildly homologous sequence
###
###                                    MTTFKLIINGKTLKGEITIEAVDAAEAEKIFKQYANDNGIDGEWTYDDATKTFTVTE
###
### We'll let PD do the alignment and byild the model. It will reuse as much of the sidechain information to be obtained
### from the original PDB although ot all clashes can be resolved (have a look at the resultant PDB file)
###

bioseq = BioSequence();
bioseq.setTo("*M-(TTFKLIINGKTLKGEITIEAVDAAEAEKIFKQYANDNGIDGEWTYDDATKTFTVT)-E*");   # Map the particles defined in the PDB over this sequence!
sim.loadExplicit('A',Polypeptide,bioseq); # Load polypeptide from chain C using a sequece override
sim.info();                               # print some summarising information info
sim.save( OutputFile_PDB( "seq.pdb" ) )   # save a PDB file with the coordinates
