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

def createff(psystem,ffps):
  ff = Forcefield(psystem,ffps)

  bonds = BondedForcefield()
  bonds.doimproper = 0
  ff.add( bonds ) 
  
  gbsa    = GB_Forcefield()
  gbsa.fastmode = 1
  ff.add( gbsa )

  sasa  = SASA_Forcefield()
  sasa.globalASP = 0.009
  ff.add( sasa )   
  
  return ff  

ffps = ForcefieldParameterSet()
ffps.readlib("lib/amber03aa.ff")

sim = SysSpec(ffps)
sim.add( ProteinPDB(ffps,"trpcage.pdb",'A',0)  )

psystem = ParticleSystem( sim )
psystem.printPDB("read-in.pdb")
psystem.info()

# Legacy
#tra = OutTraTrajectory("outtra.tra",psystem)

# New Tra System
tra = Traj_Out("outtra.tra",psystem)

psystem.addtra(tra)

ff = createff(psystem,ffps)

epot_summary(psystem,ff)

flushed_print("------------------------------")
flushed_print("First Torsional Minimisation")
flushed_print("------------------------------\n")

min1 = TorsionalMinimisation(psystem,ff)
min1.InitialCapFactor = 0.1 # steric considerations
min1.steps = 100
min1.update_scr  = 10
min1.update_tra  = 10
min1.update_mon  = 10
min1.run()

flushed_print("\n\n------------------------------")
flushed_print("Second Cartesian Minimisation")
flushed_print("------------------------------\n")

min2 = Minimisation(psystem,ff)
min2.steps = 400
min2.update_scr  = 10
min2.update_tra  = 10
min2.update_mon  = 10
min2.run()
   
epot_summary(psystem,ff)