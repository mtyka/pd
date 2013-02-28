from pd import *
cseed(3)
info()
timer()

ffps = FFParamSet("../../../lib/amber03aa.ff")
ffps.readlib("tip3.ff")

sim = PDB_In(ffps,"start.pdb");
sim.load(); 

boundary = PeriodicBox(42, 44,  47)
boundary.info()

sim.info()
sim.printPDB("inout.pdb")

nlist = NeighbourList_PeriodicBox()

wspace = WorkSpace( sim )
wspace.info()
wspace.setSpace(boundary)
wspace.setNeighbourList(nlist)

cut = float(after("-cutoff",12))

ff = Forcefield(wspace)
bonds = BondedForcefield()
nb    = NonBondedForcefield_PeriodicBox()
nb.ForceSwitch = True 
nb.EnergySwitch = False
nb.Cutoff =    cut 
nb.InnerCutoff =   cut - 2.0 #10.00
nb.VdwCutoff =   cut #12.10
nb.VdwInnerCutoff =  cut - 2.0 #10.00

ff.add(bonds)
ff.add(nb)

ff.printEnergySummary()
ff.info()

min = Minimisation(ff)
min.Steps = 50
min.UpdateScr  =    5
min.UpdateNList = 10
min.UpdateTra  = 0
min.run()

timer(0)
md = MolecularDynamics(ff)
md.Steps =   250 
md.Timestep   = 1E-15 
md.UpdateScr =  10         ## Update screen and
md.UpdateMon = 2
md.UpdateTra =  2000          ## Trajectory every 100 steps
md.UpdateNList = 20         ## Calculate new neighbor list every 10 steps
md.Integrator = MolecularDynamics.Langevin
md.UpdateRemoveTotalMomentum = False
md.TargetTemp = Temp(300)
md.run()
timer()

