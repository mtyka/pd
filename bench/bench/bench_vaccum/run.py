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

# create a few common forcefields
ff = Forcefield(wspace)

bonds = BondedForcefield()
nb    = NonBondedForcefield()
nb.Cutoff = 12.0
nb.InnerCutoff =  9.0
ff.add( bonds )
ff.add( nb )

## print energies as summary
ff.printEnergySummary()
## and dy detail (useful for when things DO break)
ff.printEnergyByAtom()
## also show parameters
ff.info()


