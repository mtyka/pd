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
## we need the gromacs module for the gromacs file formats
from gromacs import *

cseed(3)
info()
timer()

## create a water box
ffps = FFParamSet("tip3.ff")
sim = System(ffps)

## create a new molecule (water, TIP3)
mol =  NewMolecule(ffps,"TIP3")
boundary = PeriodicBox(18.1500)

## and fill the box with it
sim.solvate_N(mol, boundary ,200)

## now make a workspace 
wspace = WorkSpace( sim )
wspace.info()
wspace.printPDB("inout.pdb")
wspace.setSpace(boundary)

## Now add various trajectories. Note that they're all independent
## and can be added in any combination.

## BTF trajectory format 
tra = OutTra_BTF("output",wspace)
wspace.addTra(tra)

## NAMD Format. Note that this trajectory automatically creates
## Both the PSF and DCD file needed 
tra1 = OutTra_NAMD("output",wspace);
wspace.addTra(tra1)

## Gromacs also uses a two  part trajectory but here we represented them 
## separately, so you need to add the GRO and the XTC part.
## This  was because gromacs also has a full precision format called
## TRR which will be added later.

tra2 = OutTra_GRO("output",wspace);
wspace.addTra(tra2)
tra3 = OutTra_XTC("output",wspace);
wspace.addTra(tra3)

## A PDB makeshift trajectory (creates a multi-model PDB which for example
## VMD can read.
tra4 = OutTra_PDB("output",wspace);
wspace.addTra(tra4)


## Run some simple MD to demonstrate.

ff = Forcefield(wspace)

bonds = FF_Bonded(wspace)
nb    = FF_NonBonded(wspace)

nb.ForceSwitch = True 
nb.EnergySwitch = False

nb.Cutoff =     8.00
nb.InnerCutoff =    6.00
nb.VdwCutoff =    8.00
nb.VdwInnerCutoff =  6.00
ff.add(bonds)
ff.add(nb)

ff.info()
ff.printEnergySummary()

timer()

## Minimise for a few steps first
min = Minimisation(ff)
min.Steps = 50
min.UpdateScr  =   20
min.UpdateTra  = 0
min.run()

#### MD to relax the unequilibrarted starting structure
md = MolecularDynamics(ff)
md.Steps =     2500
md.Timestep   = float(after("-timestep","1.0E-15")) 
md.UpdateScr =  10          ## Update screen 
md.UpdateTra =  50          ## Trajectory dumps occur every 50 steps
md.UpdateNList = 10         ## Calculate new neighbor list every 10 steps
md.Integrator = MolecularDynamics.Langevin
md.UpdateRemoveTotalMomentum = False
md.setTargetTemp(300)
md.Barostat = MolecularDynamics.BerendsenBaro

md.run()

