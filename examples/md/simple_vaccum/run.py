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

###
###  This example script demonstrates the usage of simple MolecularDynamics and demonstrates 
###  most of the options it provides. It demonstrates how to set the integrators and set temperature control 
###  as well as related paramters. This script concentrates on MolecularDynamics  in vaccum and
###  implicit solvent. A different MD example covers MolecularDynamics under periodic boundary conditions.
###  \section test test test
###


from pd import *
cseed(4)
info()
timer()

ffps = FFParamSet()
ffps.readLib("amber03aa.ff")

## test simple loading of PDB file
sim = PDB_In(ffps, "../../pdb/trpcage.pdb");
sim.loadAll();

# create workspace from the system sim
wspace = WorkSpace( sim )

# print loaded system - this can be compared later
wspace.printPDB("inout.pdb")

# create a few common forcefield components

ff = Forcefield(wspace)
# set up a simple vaccuum forcefield
bonds = FF_Bonded(wspace)
nb    = FF_NonBonded(wspace)
nb.Cutoff = 20.0
nb.InnerCutoff =  16.0
ff.add( bonds )
ff.add( nb )

## add all the forcefield components to our forcefield
ff.add( bonds )
ff.add( nb )

## print energies as summary
ff.printEnergySummary()

## also show parameters
ff.info()


## remember starting structurei                     ## run it !


## Here are some alternative sets of paramters to demonstrate the range of things
## This class can do.

starting = wspace.save()

## preminimise to prevent the MD to explode straight away 
min = Minimisation(ff)
min.Steps = 25
min.UpdateScr  = 1 
min.UpdateTra  = 0
min.run()

## restore and try an MD run (constant energy)
## this should conserve energy within a margin
## but definitely without drift
wspace.load(starting)
ff.printEnergySummary()


md = MolecularDynamics(ff)   ## create the MD protocol
md.Steps = 2000              ## set the number of simulation steps
md.Timestep   = 0.5E-15      ## set the timestep to 1 femto second             
md.UpdateScr = 10            ## Print a line to the screen every 10 steps
md.UpdateTra = 50            ## Save structures to each trajectory every 50 steps
md.UpdateNList = 10          ## Update the Neighbour List every 10 steps
md.setTargetTemp( 300 )      ## Set the start temperature to 300 Kelvin
md.run()                     ## run it !

ff.printEnergySummary()
wspace.load(starting)        ## restore starting structure


## Here are some alternative sets of paramters to demonstrate the range of things
## This class can do.

md.Integrator = MolecularDynamics.Beeman             ## BeeMan integrator (excellent energy conservation, default)
md.Integrator = MolecularDynamics.VelocityVerlet     ## Velocity version of Verlet algorithm
md.Integrator = MolecularDynamics.Verlet             ## The original Verlet integrator

# Temperature control
md.setTargetTemp = ConstantProfile(300)           ## Keep temperature at a constant 300
md.setTargetTemp = LinearProfile(300)     ## Linearly Ramp temperature from 300 to 10 Kelvin
md.setTargetTemp = ExponentialProfile(200,1)      ## Do an exponential Ramp of the temperature from 200 to 1K

## Regulate the temperature using the Berendsen thermostat
md.Thermostat = MolecularDynamics.Berendsen 
md.BerendsenTau =  100e-15          ## the tau parameter (in seconds) 10-1000fs is usual. The lower this parameter
                                    ## the more tightly the temperature will be forced to be close to the desired temperature
																		
## Regulate the temperature using the Andersen thermostat
md.Thermostat = MolecularDynamics.Andersen
md.AndersenRate = 0.1              ## Set the collision rate for Andersen thermostat (0.1 is reasonable)

md.RandVel = False                  ## Next time md.run() is called, do NOT reinitialise 
                                    ## the velocities to the starting temperature
## Langevin Dynamics

md = MolecularDynamics(ff)   ## create the MD protocol
md.Steps = 2000              ## set the number of simulation steps
md.Timestep   = 0.5E-15      ## set the timestep to 1 femto second             
md.UpdateScr = 10            ## Print a line to the screen every 10 steps
md.UpdateTra = 50            ## Save structures to each trajectory every 50 steps
md.UpdateNList = 10          ## Update the Neighbour List every 10 steps
md.LangevinOnHydrogens = True ## Do we want Hydrogens to be treated with stochastic 
                             ## simulation or using velcoity verlet (newtonian)
md.setTargetTemp(300)
md.run()


timer()

