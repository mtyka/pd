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
cseed(3)
info()
timer()

ffps = FFParamSet("amber03aa.ff")
ffps.readLib("tip3.ff")

## Load the protein structure and the water molecules from PDB
sim = PDB_In(ffps,"start.pdb");
sim.load(); 

## Create a periodic box with dimensions 42, 44 and 47 Angstrom respectively for x,y and z
## (since this information is not available from the PDB file)
boundary = PeriodicBox(42, 44,  47)
boundary.info()

sim.info()
sim.printPDB("inout.pdb")

wspace = WorkSpace( sim )
wspace.info()

## set the space dimensions/geometry in the workspace
wspace.setSpace(boundary)

## The cutoff should be 12 unless overwritten from command line  
UserCutoff = float(after("-cutoff",12))

ff = Forcefield(wspace)
bonds = FF_Bonded(wspace)
nb    = FF_NonBonded(wspace)
nb.ForceSwitch = True 
nb.EnergySwitch = False

## Make the inner cutoffs 2 Angstrom smaller then the absolute cutoffs
nb.Cutoff =          UserCutoff 
nb.InnerCutoff =     UserCutoff - 2.0 
nb.VdwCutoff =       UserCutoff 
nb.VdwInnerCutoff =  UserCutoff - 2.0 
# Add the bonded and nonbonded force field components
ff.add(bonds)
ff.add(nb)

## print an energy summary
ff.printEnergySummary()
ff.info()

min = Minimisation(ff)                      ## Preminimise the system such that the MD doesnt blow up
min.Steps = 50
min.UpdateScr  =    1
min.UpdateNList = 10
min.UpdateTra  = 0
min.run()

timer(0)                                    ## Reset the timer
md = MolecularDynamics(ff)
md.Steps =   250                            ## Run for 250 steps (0.25 pico seconds)
md.Timestep   = 1E-15                       ## Timestep to the usual 1 femtosecond
md.UpdateScr =  10                          ## Update screen and
md.UpdateMon = 2                            ## Update Monitors every 2 steps
md.UpdateTra =  2000                        ## Dump trajectory every 100 steps
md.UpdateNList = 20                         ## Realculate new neighbor list every 10 steps
md.Integrator = MolecularDynamics.Langevin  ## Use langevin dynamics to move system forward
md.UpdateRemoveTotalMomentum = False        ## Removing the total momentum would disturb things 
                                            ## under periodic boundary conditions
md.setTargetTemp( 300 )                     ## Set the running temperaure to 300 K
md.run()                                    ## and go!
timer()

