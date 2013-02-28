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
bonds = FF_Bonded(wspace)
gb    = FF_GeneralizedBorn_Still( wspace)
gb.Cutoff = 12.0
gb.InnerCutoff =  9.0
gb.FastMode = 1
sasa  = FF_SASA_LCPO( wspace )
sasa.GlobalASP = 0.005

## add all the forcefield components to our forcefield
ff.add( bonds )
ff.add( gb )
ff.add( sasa )

## print energies as summary
ff.printEnergySummary()

## also show parameters
ff.info()

## try some minimisation 
min = Minimisation(ff)
min.Steps = 25
min.UpdateScr  = 1 
min.UpdateTra  = 0
min.run()

## restore and try an MD run (constant TEMP)
## this should keep its temperature 
md.Steps = 2000 
md.Timestep   = 1.00E-15
md.UpdateScr = 100
md.UpdateTra = 100
md.UpdateNList = 10
md.Integrator = MolecularDynamics.Langevin
md.LangevinOnHydrogens = True
md.setTargetTemp( 300 )
md.run()

timer()

