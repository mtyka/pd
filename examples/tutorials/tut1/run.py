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
cseed(4)   # set random number seed to current system time
timer()    # start a timer

ffps = FFParamSet()
ffps.readLib("amber03aa.ff")

## test simple loading of PDB file
sim = PDB_In(ffps, "../../pdb/trpcage.pdb");
sim.loadAll();
sim.info();

# create workspace from the system sim
wspace = WorkSpace( sim )

# print loaded system - this can be compared later
wspace.printPDB("inout.pdb")

## Create an output trajectory (tra format)
tra = OutTra_BTF("output",wspace)
wspace.addTra(tra)

## create a GB/SA Forcefield using the AMBER
## 03 Parameter set
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

ff.printEnergySummary()  ## calculate energies and show them summarised
ff.info()         ## show parameters


## Run some Minimisation
min = Minimisation(ff)
min.Steps = 50
min.UpdateScr  =  10
min.UpdateNList = 10
min.run()

## Run some Molecular Dynamics
md = MolecularDynamics(ff)
md.Steps = 2000
md.UpdateScr = 100
md.UpdateTra = 10
md.UpdateNList = 10
md.Integrator = MolecularDynamics.Langevin
md.setTargetTemp(300)

# Add a atom-atom distance monitor to the simulation 
# (we will monitor atoms with indices 102 and 263)
dist = AtomDistanceMonitor(wspace, 102, 263)
md.addMonitor(dist)
md.run()

## put the collected data into a histogram with bin width 0.05
dist.getHistogram(0.05).print_count()

## print time taken
timer()

