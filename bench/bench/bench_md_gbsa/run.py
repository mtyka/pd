from pd import *
cseed(4)
info()
timer()

ffps = FFParamSet()
ffps.readlib("../../../lib/amber03aa.ff")

## test simple loading of PDB file
sim = System(ffps)
mol =  NewProteinPDB(ffps,"trpcage.pdb")  
sim.add(mol)

# create workspace
wspace = WorkSpace( sim )

# print loaded system - this can be compared later
wspace.printPDB("inout.pdb")

# adda trajecotry
tra = Traj_Out("langevin",wspace)
wspace.addTra(tra)

# create a few common forcefields
ff = Forcefield(wspace)

bonds = BondedForcefield()
nb    = NonBondedForcefield()

nb.Cutoff = 12.0
nb.InnerCutoff =  9.0

gb    = GB_Forcefield( nb )
gb.FastMode = 1
sasa  = SASA_Forcefield()
sasa.GlobalASP = 0.005
ff.add( bonds )
ff.add( gb )
ff.add( sasa )

## print energies as summary
ff.printEnergySummary()

## also show parameters
ff.info()


## remember starting structure
starting = wspace.save()

## try some minimisation 
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
md = MolecularDynamics(ff)
md.Steps = 2000 
md.Timestep   = 0.5E-15
md.UpdateScr = 10
md.UpdateTra = 10
md.UpdateNList = 10
md.TargetTemp = Temp(300)
md.run()

## restore and try an MD run (constant TEMP)
## this should keep its temperature 
wspace.load(starting)
ff.printEnergySummary()
md.Steps = 2000 
md.Timestep   = 1.00E-15
md.UpdateScr = 100
md.UpdateTra = 100
md.UpdateNList = 10
md.Integrator = MolecularDynamics.Langevin
md.LangevinOnHydrogens = True
md.TargetTemp = Temp(300)
md.run()


timer()

