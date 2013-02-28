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

## PD Example input script for a simple Monte Carlo simulation
##
##
from pd import *
ctimeseed()      ## Set a random number seed using system time
timer()          ## Start a generic timer


## Set up a simple system containing a short peptide 
## and a simple forcefield contianing electrostatics
## and vdw forces

ffps = FFParamSet("amber03aa.ff")

sim = System(ffps)


mol =  NewProteinHelix(ffps,"*A-(APGIKKQ)-A*")

sim.add(mol)
sim.info()
rotbond = RotBond()
wspace = WorkSpace( sim )
wspace.setRotatableBondList( rotbond )
wspace.info()

tra = OutTra_NAMD("output",wspace)
wspace.addTra(tra)

ff = Forcefield(wspace)

bonds = FF_Bonded(wspace)
nb    = FF_NonBonded(wspace)

nb.Cutoff = 12.0
nb.InnerCutoff =  9.0
nb.VdwCutoff =  5.0
nb.VdwInnerCutoff =  4.0
nb.DDDielectric = True

ff.add( bonds )
ff.add( nb )

ff.printEnergySummary()


## Set up the evaluator. In this case this is simply a 
## single energy evaluation
## 

ene = Energy(ff)  ## simple evaluator (just the energy)

## Create a set of moves to be executed. Each will be called
## at every move. 

moveset = MoveSet(wspace)
moveset.add( SidechainTorsionalMove( wspace, 1, 2, 0.1 ) )
moveset.add( BackboneTorsionalMove( wspace, 0.5, 0.8, 40 ) )

## Now create the Monte Carlo simulation using the moveset and the
## evaluator. MonteCarlo is simply another type of Protocol.

mc = MonteCarlo(ene,moveset)
mc.Steps = 200               ## Set the number of steps
mc.UpdateScrAcc=True;        ## These set which events result in the
mc.UpdateScrRej=True;        ## printing of a line on the screen. (acceptances/rejections)
mc.UpdateTraAcc=True;        ## ditto for saving structures to the trajectory
mc.UpdateTraRej=True;

## Launch it - above we used an Energy Calculation as the evaluator
## - thus this will be a standard metropolis monte carlo algorithm
mc.run()

timer()

