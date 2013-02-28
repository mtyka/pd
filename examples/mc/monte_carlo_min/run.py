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

## Now lets try a different evaluator - a minimisation. It will be called
## instead of the simple energy evaluation. This is called Monte Carlo
## with minimisation (Scheraga et al. 1983)
## Note that paramters to do with the evaluator (or the "subprotocol")
## are set in the appropriate object, not the monte carlo object (mc)
ene = Minimisation(ff)  ## simple evaluator (just the energy)
ene.Steps = 200         ## 200 steps of minimisation
ene.UpdateScr = 20      ## Screen update frequency
ene.UpdateTra = 20      ## Screen update frequency

## Create a set of moves to be executed. Each will be called
## at every move. 

moveset = MoveSet(wspace)
moveset.add( SidechainTorsionalMove( wspace, 0.5, 0.8, 40 ) )
moveset.add( BackboneTorsionalMove( wspace, 0.5, 0.8, 40 ) )

## Now create the Monte Carlo simulation using the moveset and the
## evaluator. MonteCarlo is simply another type of Protocol.

mc = MonteCarlo(ene,moveset)
mc.Steps = 200               ## Set the number of steps
mc.UpdateScrAcc=True;        ## These set which events result in the
mc.UpdateScrRej=True;        ## printing of a line on the screen. (acceptances/rejections)
mc.UpdateTraAcc=True;        ## ditto for saving structures to the trajectory
mc.UpdateTraRej=True;

mc.run()


timer()

