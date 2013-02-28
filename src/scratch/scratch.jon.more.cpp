// PD is a free, modular C++ library for biomolecular simulation with a 
// flexible and scriptable Python interface. 
// Copyright (C) 2003-2013 Mike Tyka and Jon Rea
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "global.h"

#include "../pd/accessory.h"
#include "tools/statclock.h"
#include "tools/draw.h"
#include "maths/fastrandom.h"

#include "forcefields/ffbonded.h"
#include "forcefields/breakablebonded.h"
#include "forcefields/nonbonded.h"
#include "forcefields/gbff.h"
#include "forcefields/lcpo.h"
#include "forcefields/ffsoftvdw.h"
#include "forcefields/pops.h"
#include "forcefields/restraint_positional.h"

#include "monitors/basicmonitors.h"

#include "system/proximitygrid.h"
#include "workspace/workspace.h"
#include "workspace/bondorder.h"

#include "library/angleset.h"
#include "library/rotamerlib.h"
#include "library/rotamer_shetty.h"
#include "library/rotamer_dunbrack.h"

#include "manipulators/movebase.h"
#include "manipulators/basicmoves.h"

#include "manipulators/rotamer_applicatorbase.h"
#include "manipulators/rotamer_scwrl.h"

#include "protocols/minimise.h"
#include "protocols/dualffminimiser.h"
#include "protocols/md.h"
#include "protocols/scpack.h"

#include "fileio/tra.h"
#include "fileio/pdb.h"

#include "filters/basicfilters.h"

#include "arcus.h"

#include "scratch.jon.more.h"

using namespace std;
using namespace Physics;
using namespace Protocol;
using namespace IO;
using namespace Manipulator;

namespace Protocol
{
	MonteCarloLowKeep::MonteCarloLowKeep(
		ProtocolBase&  _evaluator,
		MoveSet& _moveset ): 
	MonteCarlo(
			_evaluator,
			_moveset
			), m_BestEne(DBL_MAX)
	{
		m_StartEne = getWSpace().ene.epot;
	}

	int MonteCarloLowKeep::runcore()
	{
		m_BestEne = DBL_MAX; 
		return MonteCarlo::runcore();
	}

	bool MonteCarloLowKeep::accept(double enenew, double eneold) const
	{
		const double acc = ( m_StartEne + (fabs(m_StartEne) * 0.3) );
		static const double conv = 1.0 / (Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na);
		if( enenew < m_BestEne ) m_BestEne = enenew;
		if( enenew < acc )
		{
			printf("*Valid*, Cur:%14.3lf Cutoff:%8.3lf Start:%8.3lf Best:%8.3lf\n", 
				enenew * conv, acc * conv, m_StartEne * conv, m_BestEne * conv);
			return true;
		}
		else
		{
			printf("InValid, Cur:%14.3lf Cutoff:%8.3lf Start:%8.3lf Best:%8.3lf\n", 
				enenew * conv, acc * conv, m_StartEne * conv, m_BestEne * conv);
			return false;
		}

		// Metropolis
		//return MonteCarlo::accept( enenew, eneold );
	}
}

void Test_Conformer()
{
	// set a standard random number seed for this test.
	Maths::FastRandom::getInstance()->reinitialise(100); 

	// Load our anglsest
	Library::AngleSet as;
	as.loadFromFile("lib/angleset/default.loop.angleset");
	as.info();

	// Load our forcefield
	FFParamSet ffps;
	ffps.readLib("lib/amber03aa.ff");

	// Evaluation system
	PDB_In sys( ffps, "./pd_double/test.pdb" );
	sys.load();

	// Wspace
	WorkSpace wspace( sys );
	wspace.info();
	wspace.addStdTra("impFeb");
	wspace.outtra.append();

	SegmentDef def1;
	def1.init(wspace,37,44,SegBreakStart);

	rebuildExtdFromStart(wspace,def1.getStartResIndex(),def1.getEndResIndex());
	wspace.outtra.append();
	
	rebuildExtdFromEnd(wspace,def1.getStartResIndex(),def1.getEndResIndex());
	wspace.outtra.append();
	
	rebuildExtdFromStart(wspace,def1.getStartResIndex(),def1.getCentreResIndex());
	rebuildExtdFromEnd(wspace,def1.getCentreResIndex()+1,def1.getEndResIndex());
	wspace.outtra.append();

	// A random conformer builder - lets perturb our helix!
	for( int i = 0; i <= 2; i++ )
	{
		SegmentDef def2;
		SegBreakType breakType = (SegBreakType)i;
		wspace.outtra.append();
		def2.init(wspace,37,44,breakType);
		wspace.outtra.append();

		ConfBuild_RandomSingle conf(as,def2);
		//conf.setBuildMode( ConfBuilderBase::Backbone );

		wspace.outtra.append();
		conf.MinimiseIdealised();
		conf.revertIdealised(); // we want to perform the build using the idealised structure
		wspace.outtra.append();
		conf.setRandCount(100);
		while( conf.next() )
		{
			wspace.outtra.append();
		}
	}

	return;
}

//void TestSoftSteric( WorkSpace& wspace )
//{
//	Forcefield stericFF(wspace);
//
//	BondedForcefield *bonds = new BondedForcefield(wspace);
//	stericFF.addWithOwnership(*bonds);
//
//	SoftVDWForcefield *sff = new SoftVDWForcefield(wspace);
//	stericFF.addWithOwnership(*sff);
//
//	stericFF.setup();
//
//	Minimisation min2( stericFF );
//	min2.Steps = 501;
//	min2.StepSize = 2E1; // 1.0E7
//	min2.UpdateNList = 1;
//	min2.UpdateScr = 1;
//	min2.UpdateTra = 1;
//	min2.UpdateMon = 0;
//	min2.run();
//
//	MolecularDynamics md( stericFF );
//	md.Steps = 6001;
//	md.Integrator = MolecularDynamics::Beeman;
//	md.Thermostat = MolecularDynamics::NoThermostat;
//	md.Timestep   = 1.00E-15;
//	md.UpdateScr = 1;
//	md.UpdateNList = 10;
//	md.UpdateTra = 1;
//	md.TargetTemp = new Temp(300);
//	md.run();
//
//	return;
//}
//
//void TestSoftSteric()
//{
//	FFParamSet ffps;
//	ffps.readLib("lib/amber03aa.ff");
//
//	System sys(ffps);
//	sys.add(NewProteinHelix(ffps,"*A-(AAAAPAAAA)-A*" ));
//
//	WorkSpace wspace = WorkSpace( sys );
//	wspace.addStdTra( "softvdw" );
//
//	TestSoftSteric( wspace );
//}
//
//void Test_BrokenBonded()
//{
//	StatClock clockMe;
//	clockMe.Begin();
//
//	FFParamSet ffps;
//	ffps.readLib("lib/amber03aa.ff");
//
//	System sys(ffps);
//	sys.add(NewProteinHelix(ffps,"*A-(AAAAPAAAA)-A*" ));
//
//	WorkSpace wspace = WorkSpace( sys );
//	PDB_Out out1("imported",wspace);
//	wspace.info();
//
//	wspace.addStdTra("imported");
//
//	Forcefield ff = createff(wspace,true);
//	FF_BreakableBonded& broken = (FF_BreakableBonded&)ff.element(0);
//	broken.createBreak(84,86);
//
//	Minimisation min2( ff );
//	min2.Steps = 501;
//	min2.StepSize = 2E1;
//	min2.UpdateNList = 5;
//	min2.UpdateScr = 100;
//	min2.UpdateTra = 5;
//	min2.UpdateMon = 0;
//	min2.run();
//
//	MolecularDynamics md(ff);
//	md.Steps = 16001;
//	md.Timestep   = 1.00E-15;
//	md.UpdateScr = 100;
//	md.UpdateNList = 10;
//	md.UpdateTra = 20;
//	md.TargetTemp = new Temp(3000);
//	md.run();
//
//	clockMe.End();
//
//	clockMe.ReportMilliSeconds(10);
//
//	return;
//}
//

//
//void Test_CoiledCoil()
//{
//	FFParamSet ffps;
//	ffps.readLib("lib/amber03aa.ff");
//
//	//PDB_In sim(ffps,"coiled_coil/zinc_finger.pdb");
//	PDB_In sim(ffps,"coiled_coil/coiled_coil.pdb");
//	sim.setAlignerDirect();
//
//	std::string seqString = "*Y-(TEQERQIREKKKKLSALEQTISLV)-K*";
//	Sequence::BioSequence seq( ffps );
//	seq.setTo(seqString);
//	sim.loadExplicit('A',Library::Polypeptide,seq);
//	sim.loadExplicit('B',Library::Polypeptide,seq);
//	
//	sim.printPDB("coiled_coil/pre-packed.pdb");
//
//	WorkSpace wspace(sim);
//	wspace.addStdTra("coiled_coil/built");
//
//	// Forcefield config
//	Forcefield ff = createff(wspace);
//	ff.printEnergySummary();
//	Forcefield ffs = createffs(wspace);
//	ffs.printEnergySummary();
//
//	MCPackSideChains( wspace, ffs, ff, "lib/dict.pdb" );
//
//	sim.printPDB("coiled_coil/final.pdb");
//
//	return;
//}

//void TestPops()
//{
//	// set a standard random number seed for this test.
//	Maths::FastRandom::getInstance()->reinitialise(100); 
//
//	// Load our anglsest
//	Library::AngleSet as;
//	//as.loadFromFile("lib/angleset/big.all.angleset");
//
//	as.info();
//	FFParamSet ffps;
//	//ffps.readlib("f:\\_scratch\\lib\\amber03aa.ff");
//	ffps.readLib("c:\\_source\\_scratch\\lib\\amber03aa.ff");
//
//	//PDB_In sim(ffps,"1bw8.pdb");
//	//PDB_In sim(ffps,"f:\\_scratch\\1hme.pdb");
//	//PDB_In sim(ffps,"f:\\_scratch\\1DK0A.pdb");
//	PDB_In sim(ffps,"c:\\temp.pdb");
//	sim.disableRebuild();
//	sim.setFilter(Library::Polypeptide);
//	sim.UseSEQRES = false;
//	//sim.load('A');	
//	sim.load();	
//
//	WorkSpace wspace(sim);
//	wspace.addStdTra("pops");
//	wspace.outtra.append();
//
//	Pops pop;
//	pop.readDat("c:\\_source\\_scratch\\lib\\pops\\pops-dna.dat"); 
//	//pop.setTo(wspace, Pops::AllAtom);
//	pop.setTo(wspace, Pops::Coarse);
//	pop.calc();		
//	pop.detail();
//
//	//SegmentDef def2(wspace,28,33,false);
//	SegmentDef def2(wspace,57,62,false);
//
//	pop.calc();
//	pop.detail();
//
//	std::vector<Maths::dvector> native;
//	native.resize(def2.getNRes());
//	getMePos( wspace, def2, native );
//
//	std::vector<Maths::dvector> conformer;
//	conformer.resize(def2.getNRes());
//	
//	// A random conformer builder - lets perturb our helix!
//	ConfBuild_RandomSingle conf(as,def2);
//	//conf.MinimiseIdealised();
//	conf.revertIdealised(); // we want to perform the buildup on the idealised structure
//	conf.setBuildMode( Manipulator::ConfBuilderBase::Backbone );
//	wspace.outtra.append();
//
//	// Join filter
//	Physics::BondedForcefield bonds(wspace);
//	bonds.ensuresetup( wspace );
//	OmegaGroupFilter joinFilter;
//	joinFilter.setTo( bonds, def2.getEndResIndex() );
//	joinFilter.setTollerance( 3.0 );
//
//	// 2) Surface clash filter
//	SurfacePicker surfaceBox;
//	surfaceBox.calcBoundingBox( def2 );
//	ClashFilter filter;
//	filter.setOLapFac( 0.7 ); // quite heavy clashes are allowed - chain crosses definatly are not
//	filter.setForPicker( 
//		Pick_AND( 
//			Pick_AND( 
//				Pick_NOT( PickResidue(def2.getEndResIndex()) ),
//				PickResidueRange(def2)), 
//			PickCoreBackbone()
//			)
//		);
//	filter.setAgainstPicker( Pick_AND( Pick_NOT(PickResidueRange(def2)), surfaceBox) );
//	filter.setMolecule( wspace );
//
//	conf.setRandCount(100000);
//	while( conf.next() )
//	{
//		if( joinFilter.passes() )
//		{
//			if( filter.passes() )
//			{	
//				//getMePos( wspace, def2, conformer );
//				//double rms = getMeRMS( native, conformer );
//				//if( rms < 2.5 )
//				{
//
//					pop.calc();
//					//pop.detail();
//					double fract = pop.sasaFraction( def2 );
//					double sasa = pop.SASA( def2 );
//					printf("%8.3lf %8.3lf Res SASA: ",sasa,fract);
//					for( size_t i = def2.getStartResIndex(); i <= def2.getEndResIndex(); i++ )
//					{
//						printf("%8.3lf ", pop.resFraction(i));
//					}	
//					getMePos( wspace, def2, conformer );
//					double rms = getMeRMS( native, conformer );
//					printf("%8.3lf\n", rms );
//					wspace.outtra.append();
//
//					//printf(".");
//					//wspace.outtra.append();
//				}				
//			}
//		}
//	}
//
//	//for( size_t i = 0; i < wspace.natom(); i++ )
//	//{		
//	//	double sasa =  pop.atomSASA(i);
//	//	if( sasa != 0.0 ) 
//	//	{
//	//		wspace.getAtom(i).info(0,4,false);	
//	//		printf(": Atom SASA: %8.3lf\n",sasa);
//	//	}
//	//}
//
//	//for( size_t i = 0; i < wspace.nres(); i++ )
//	//{
//	//	printf("Res SASA: %8.3lf\n", pop.resSASA(i));
//	//}		
//
//	return;
//}
//
//void Test_Mike()
//{
//	FFParamSet ffps;
//	ffps.readLib("lib/amber03aa.ff");
//
//	PDB_In sim(ffps,"pep4_A_84_237_ss.pdb");
//	sim.setFilter(Library::Polypeptide);
//	Sequence::BioSequence seq( ffps );
//	seq.setTo("*A-CYX-(T)-CYX-(SYSWPPA)-CYX-(Y)-CCYX");
//	sim.loadExplicit('I',Library::Polypeptide,seq);
//	
//	//sim.printPDB("output.pdb");
//
//	WorkSpace wspace( sim );
//	TestSoftSteric( wspace );
//
//	//wspace.addStdTra("checky_out.pdb");
//	//wspace.outtra.append();
//
//	// Forcefield config
//	//bool useBrokenBonded = true; // Flag the use of the bonded forcefield that allows bond breaks.
//	//Forcefield ff = createff(wspace,ffps,useBrokenBonded,useBrokenBonded);
//
//	//Minimisation min2( ff );
//	//min2.Steps = 401;
//	//min2.StepSize = 2E1;
//	//min2.UpdateScr = 10;
//	//min2.UpdateTra = 1;
//	//min2.UpdateMon = 10;
//	//min2.run();
//
//	//MolecularDynamics md(ff);
//	//md.Steps = 5001;
//	//md.UpdateScr = 10;
//	//md.UpdateTra = 10;
//	//md.UpdateMon = 10;
//	//md.run();
//
//	return;
//}

//void main_CraigSim()
//{
//	FFParamSet ffps;
//	ffps.readLib("amber03aa.ff");
//
//	// Load our PDB
//	PDB_In sys(ffps,"before.pdb");
//	sys.load();
//
//	WorkSpace wspace( sys );
//	wspace.info();
//	wspace.addStdTra("out");
//
//	//Forcefield ff = createff(wspace,false,true);
//	Forcefield ff = createffVac(wspace,false,true);
//
//	CartesianRestraint rest(wspace);
//	rest.k = 5;
//	ff.add(rest);
//
//	Minimisation minim(ff);
//	minim.Steps = 20001;
//	minim.UpdateNList = 50;
//	minim.UpdateMon = 0;
//	minim.UpdateScr = 100;
//	minim.UpdateTra = 0;
//	minim.StepSize = minim.StepSize * 10;
//	minim.SlopeCutoff = 0.05 * PhysicsConst::kcal2J / PhysicsConst::Na;
//	minim.run();
//
//	MolecularDynamics md( ff );
//	md.Steps = 2500;
//	md.Integrator = MolecularDynamics::Langevin;
//	md.Thermostat = MolecularDynamics::NoThermostat;
//	md.Timestep   = 1.00E-15;
//	md.UpdateScr = 100;
//	md.UpdateNList = 10;
//	md.UpdateTra = 0;
//	md.UpdateMon = 0;
//	md.TargetTemp = new Temp(600);
//	md.run();
//	md.RandVel = false; // have to run MD once before we call this
//
//	SnapShot s;
//	for( size_t i = 0; i < 1000; i++ )
//	{
//		md.run();
//		s = wspace.save();
//		minim.run();
//		wspace.Step = i;
//		wspace.outtra.append();
//		wspace.load(s);
//	}
//}
//
//int Tweaky(int argc, char** argv)
//{
//	if( argc != 4 )
//	{
//		printf("Incorect params");
//		return -1;
//	}
//	std::string fileStem(argv[1]);
//
//	int startIndex = -1;
//	if( 0 != str2int(argv[2], startIndex ) )
//	{
//		printf("Incorect param 2");
//		return -1;
//	}
//
//	int loopLength = -1;
//	if( 0 != str2int(argv[3], loopLength ) )
//	{
//		printf("Incorect param 3");
//		return -1;
//	}
//
//	int endIndex = startIndex + loopLength - 1;
//
//	// Load FFParam
//	FFParamSet ffps;
//	ffps.readLib("amber03aa.ff");
//
//	Library_Legacy::RotamerLibrary rotLib;
//	rotLib.readRotamerLibrary( "dict.pdb" );
//
//	// Load our favourite PDB molecule
//	PDB_In sys(ffps,fileStem + ".pdb");
//	sys.setFilter(Library::Polypeptide); // ONLY load the polypeptide 
//	sys.load(); // load from model 1 (any chains)
//
//	// Make a workspace
//	WorkSpace wspace = WorkSpace( sys );
//	wspace.addStdTra(fileStem + "_" + int2str(startIndex) + "_" + int2str(loopLength) );
//	wspace.info();
//
//	SegmentDef segdef( wspace, startIndex, endIndex, false );
//
//	// Create the forcefields
//	Forcefield ffs = createffs(wspace,true,false);
//	//Forcefield ff = createff(wspace,true,false);
//	Forcefield ff = createffVac( wspace, true, false );
//
//	// Our moveset
//	Manipulator::MoveSet moves( &wspace );
//	Manipulator::SegPerturbation_BBTorsional p1(segdef,0.5,Maths::DegToRad(10.0));
//	Manipulator::SegPerturbation_SCRotamer p2(segdef, &rotLib, 0.5, 0.0 );
//	Manipulator::SegPerturbation_SCTorsional p3( segdef, 1.0, Maths::DegToRad(30.0), 0.2 );
//	moves.add(p1);
//	moves.add(p2);
//	moves.add(p3);
//
//	Protocol::DualFFMinimiser eval(ff, ffs, PickResidueRange(wspace, startIndex, endIndex ) );
//	eval.!OutputLevel = true;
//
//	eval.SDPreMinSteps = 50;
//	
//	eval.StericKillFull = 4500.0; // kcal / mol
//	eval.StericMinSteps = 750;
//	eval.StericStepSize = 2E5;
//	eval.StericSlopeCutoff = 0.1 * Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na;
//	
//	eval.Steps = 750;
//	eval.StepSize = 2E5;
//
//	eval.UpdateScr = 0;
//	eval.UpdateTra = 0;
//	eval.UpdateMon = 0;
//	
//	eval.SlopeCutoff = 0.1 * Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na;
//
//	eval.run();
//	//eval.!OutputLevel = true;
//
//	// Perform a montecarlo procedure
//	Protocol::MonteCarloLowKeep mc( eval, moves );
//	mc.DoZeroGeometry = false;
//	mc.Steps = 1000;
//	mc.UpdateTraAcc = true;
//	mc.UpdateTraRej = false;
//	mc.UpdateTra = 1; // performed by low-keep accept function
//	mc.UpdateScr = 1;
//	mc.UpdateMon = 1;
//	mc.FinalState = Protocol::MonteCarlo::LowestEpot;	
//
//	BondFilter bf(wspace);
//	bf.setTollerance(0.5);
//	mc.addPostFilter(bf);
//
//	Monitors::CRMSMonitor mon(&wspace);
//	mc.addMonitor(mon);
//
//	mc.run(); // Pack sidechains using the rotamers	
//
//	return 0;
//}
//
void main_GentleHarmonicRestraintMinimisation()
{
	StatClock timer;
	timer.Begin();

	// Load our forcefield parameters
	FFParamSet ffps;
	ffps.readLib("amber03aa.ff");
	ffps.readLib("tip3.ff");

	// Load our PDB
	PDB_In sys(ffps,"TrpcageMine.pdb");
	sys.load();

	// Make a full workspace from the imported system
	WorkSpace wspace( sys );
	wspace.info();
	wspace.addStdTra("sim");

	// Write what we have interpreted from the import file	
	PDB_Writer("imported.pdb").write(wspace); // pdb

	// Forcefield config
	bool useBrokenBonded = false; // Flag the use of the bonded forcefield that allows bond breaks.
	Forcefield ff = createffVac(wspace,useBrokenBonded,true);

	FF_Restraint_Positional rest(wspace);
	rest.k = 20.0;
	ff.add(rest);

	// Harmonic
	Minimisation minim(ff);
	minim.Algorithm = Minimisation::SteepestDescent;
	minim.Steps = 50;
	minim.UpdateTra = 0;
	minim.UpdateNList = 10;
	minim.UpdateMon = 0;
	minim.UpdateScr = 10;
	minim.StepSize = minim.StepSize * 10;
	minim.run();

	rest.k = 10;
	minim.Algorithm = Minimisation::ConjugateGradients;
	minim.run();

	rest.k = 5;
	rest.setSelection( PickCoreBackbone() );
	minim.StepSize = minim.StepSize * 10;
	minim.run();

	rest.deactivate();
	minim.Steps = 100;
	minim.run();

	PDB_Writer("cooked.pdb").write(wspace); // pdb

	timer.End();
	timer.ReportMilliSeconds();

	return;
}

//void TempFunc()
//{	
//	std::string seqString = "*Y-(PCPFCFKEFTRKDNMTAHVKIIHK)-I*";
//	std::string outStem = "2drp40to65";
//	std::string loadPDB = outStem + ".pdb";
//
//	//  1) Load the forcefield;
//	//  We will be using bog standard amber03 all atom;
//	FFParamSet ffps;
//	ffps.readLib("amber03aa.ff");
//	
//	//  2) Load our template PDB file;
//	PDB_In sim(ffps,loadPDB);
//	sim.setAlignerDirect(); //  required to disable internal sequence alignment;
//	
//	//  3) Parse the sequence;
//	Sequence::BioSequence seq( ffps );
//	seq.setTo(seqString);
//	
//	//  4) Thread our structure;
//	sim.loadExplicit('A',Library::Polypeptide,seq);
//	sim.printPDB(outStem + ".pre-packed.pdb");
//	
//	//  5) Generate a workspace so that we can pack sidechains;
//	WorkSpace wspace(sim);
//	wspace.addStdTra(outStem);
//	
//	//  6) Configure dual forcefield;
//	Forcefield ffs = createffs(wspace); //  steric only forcefield for severe clash resolution;
//	ffs.printEnergySummary();
//	Forcefield ff = createff(wspace); //  normal forcefield;
//	ff.printEnergySummary();
//	
//	//  7) Call the packing procedure;
//	Library_Legacy::RotamerLibrary rotLib;
//	rotLib.readRotamerLibrary("dict.pdb");
//	
//	//  8) Our moveset;
//	MoveSet moves( &wspace );
//	SidechainRotamerLibMove rotamer( &wspace, &rotLib, 1.0, 1.0 );
//	moves.add( rotamer );
//	
//	//  9) Minimise only the sidechains;
//	DualFFMinimiser eval(ff, ffs, PickSidechains() );
//	eval.!OutputLevel = true;
//	eval.SDPreMinSteps = 51;
//	eval.StericMinSteps = 1001;
//	eval.StericSlopeCutoff = 0.1 * Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na;
//	eval.StericKillFull = 1000.0;
//	eval.Steps = 501;
//	eval.StepSize *= 10;
//	eval.SlopeCutoff = 0.05 * Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na;
//	eval.UpdateScr = 0;
//	eval.UpdateTra = 0;
//	eval.UpdateMon = 0;
//	//eval.run();
//	
//	//  10) Perform a montecarlo procedure;
//	MonteCarlo mc( eval, moves );
//	mc.Temperature = new Temp(500);
//	mc.Steps = 500;
//	mc.UpdateTra = 1;
//	mc.UpdateScr = 1;
//	mc.UpdateMon = 1;
//	mc.UpdateTraAcc = true;
//	//  mc.UpdateTraRej = True;
//	mc.FinalState = MonteCarlo::LowestEpot;
//	
//	//  11) Filtering
//	ClashFilter cFilter(wspace);
//	cFilter.setForPicker( Pick_AND(PickHeavyAtoms(), PickSidechains()) );
//	cFilter.setOLapFac(0.5); //  heavy-ish;
//	mc.addPreFilter(cFilter);
//
//	BondFilter bFilter(wspace);
//	mc.addPostFilter(bFilter);
//	
//	//  13) Print the cooked structure;
//	mc.run(); //  Pack sidechains using the rotamers;
//	
//	//  14) Print the cooked structure;
//	sim.printPDB(outStem + ".final.pdb");	
//}

void TestProximityGrid()
{
	FFParamSet ffps;
	ffps.readLib("lib/amber03aa.ff");

	System sys( ffps );
	sys.add( NewProteinHelix(ffps,"*A-(CDEFGHIKLMNPQRSTVW)-Y*") );

	WorkSpace wspace( sys );
	wspace.info();
	wspace.printPDB("rot_start.pdb");
	IO::OutTra_BTF& tra = wspace.addStdTra("rot_test");
	IO::BTF_Block_Vector* vect = new IO::BTF_Block_Vector(10000);
	tra.addOwnedBlock(vect);
	wspace.outtra.append();

	const double cutOff = 2.5;
	PickCoreBackbone pick;
	ProximityGrid grid( wspace, pick, cutOff );
	grid.refresh();

	//grid.drawGrid( *vect );
	//wspace.outtra.append();

	const ParticleStore& atoms = wspace.atom;

	Maths::dvector seekPos;
	Maths::FastRandom fr;

	seekPos.x = (fr.nextDouble() * (2.0 * grid.getXBound())) - grid.getXBound();
	seekPos.y = (fr.nextDouble() * (2.0 * grid.getYBound())) - grid.getYBound();
	seekPos.z = (fr.nextDouble() * (2.0 * grid.getZBound())) - grid.getZBound();
	seekPos.add( grid.getWorldCentre() );

	//for( double inc = 0.1; inc < 20.0; inc += 0.1 )
	//{
	//	grid.setCutoff(inc);
	//	grid.drawMatches( *vect, seekPos );
	//	grid.drawGrid( *vect );
	//	wspace.outtra.append();
	//}

	StatClock boba;
	StatClock bob;
	while( true )
	{
		seekPos.x = (fr.nextDouble() * (2.0 * grid.getXBound())) - grid.getXBound();
		seekPos.y = (fr.nextDouble() * (2.0 * grid.getYBound())) - grid.getYBound();
		seekPos.z = (fr.nextDouble() * (2.0 * grid.getZBound())) - grid.getZBound();
		seekPos.add( grid.getWorldCentre() );

		bob.Begin();

		double sumMe = 0.0;
		const int ddd = 1000;

		const GridPoint* gridList;
		int atomsClash = 0;
		for( int k = 0; k < ddd; k++ )
		{
			grid.refresh();
			for( int i = 0; i < ddd; i++ )
			{			
				grid.atomList( seekPos, gridList );

				for( int j = 0; j < gridList->size(); j++ )
				{
					double dist = seekPos.sqrdist( atoms[gridList->atomIndexes[j]].pos() );
					if( dist < 1.0 )
					{
						sumMe += dist;
						atomsClash++;
					}
					//if( dist < 1.0 )
					//{
					//	printf("Clash\n");
					//	drawPoint(*vect,seekPos,Colour::Magenta);
					//	wspace.outtra.append();
					//}
				}

				//if( grid.drawMatches( *vect, seekPos ) )
				//{
				//	grid.drawGrid( *vect );
				//	wspace.outtra.append();
				//}
				//else
				//{
				//	grid.drawGrid( *vect );
				//	wspace.outtra.append();
				//	//vect->ForceClear();
				//}
			}
		}

		bob.End();
		bob.ReportMilliSeconds();

		PosPointer poses(wspace,pick);

		boba.Begin();
		int atomsClashA = 0;
		double sumMeA = 0.0;
		for( int k = 0; k < ddd; k++ )
		{
			for( int i = 0; i < ddd; i++ )
			{
				for( int j = 0; j < poses.size(); j++ )
				{
					double dist = seekPos.sqrdist(poses.p(j));
					if( dist < 1.0 )
					{
						sumMeA += dist;
						atomsClashA++;
					}
				}				
			}
		}		

		boba.End();
		boba.ReportMilliSeconds();
		
		if( atomsClash > 0 )
		{
			printf("Atom -> %d\n",atomsClash);
			printf("%8.3f\n",sumMe);
			printf("AtomA -> %d\n",atomsClashA);
			printf("%8.3f\n\n",sumMeA);
		}
	}

	return;
}

void randomRotamerApplication( bool _testSterics, bool importForeignLib )
{	
	FFParamSet ffps;
	ffps.readLib("lib/amber03aa.ff");

	System sys( ffps );
	sys.add( NewProteinHelix(ffps,"*A-(CDEFGHIKLMNPQRSTVW)-Y*") );

	WorkSpace wspace( sys );
	wspace.info();
	wspace.printPDB("rotamer_perturb_random_start.pdb");
	wspace.addStdTra("rotamer_perturb_random_end");

	Library::RotamerLibrary rot(ffps);

	if(importForeignLib)
	{	
		rot.convertLib("lib/rotlib/shetty/scl-B30-occ1.0-rmsd1.0-prop20.0",Library::RotLibConvert_Shetty());
		rot.writeLib("lib/rotlib/shetty.rotamer");
	}
	else
	{
		rot.readLib("lib/rotlib/shetty.rotamer");
	}

	StatClock timeMe;
	timeMe.Begin();

	RandomRotamerApplicator app( wspace, rot, ApplyCartesian );
	if( _testSterics )
		app.addFilter_DefaultSteric(); // Add the linear-time clash filter for all rotamer states
	app.test(500); // 500 random rotamer applications over the entire PickedResidueList

	//SidechainRotamerLibMove move(wspace, rot, 1.0, 1.0 );
	//move.test(200);

	//SegmentDef def(wspace, 5, 15, false );
	//SegPerturbation_SCRotamer segmove( def, rot, 0.5, 1.0 );
	//segmove.test(350);

	timeMe.End();
	timeMe.ReportMilliSeconds();

	wspace.printPDB("rotamer_perturb_random_end.pdb");
}

void TheWandersOfRotamers()
{
	FFParamSet ffps;
	ffps.readLib("lib/amber03aa.ff");

	Library::RotamerLibrary rot(ffps);
	//rot.convertLib("lib/rotlib/shetty/scl-B30-occ1.0-rmsd1.0-prop20.0",Library::RotLibConvert_Shetty());
	rot.convertLib("lib/rotlib/shetty/scl-B30-occ1.0-rmsd0.5-prop20.0",Library::RotLibConvert_Shetty());
	//rot.writeLib("lib/rotlib/shetty.rotamer");
	//rot.readLib("lib/rotlib/shetty.rotamer");

	//System sys(ffps);
	//sys.add( NewProteinHelix(ffps,"*A-(CDEFGHIKLMNPQRSTVW)-Y*") );
	//sys.add( NewProteinHelix(ffps,"*A-(RRRRRRRRRR)-A*") );

	// ERROR state!
	// YECKYCTKQFPDLNMFTFH	

	Sequence::BioSequence seq(ffps);
	seq.setTo("*Y-(E)-CYM-(KY)-CYM-(RKQFPDRNMFTF)-HID-(VDSE)-HID-(P)-N*");

	//PDB_In sys( ffps, "dek/for_jon/zf_backbones/24/zinc_finger.pdb" );
	PDB_In sys( ffps, "dek/_scripts/1. Make Initial Zinc Finger Models/2drp40to65.pdb" );
	sys.setAlignerDirect();
	sys.loadExplicit('A',Library::Polypeptide,seq);

	//PDB_In sys( ffps, "_egress/1SMD_truncated.pdb" );
	//sys.loadExplicit(' ',Library::Polypeptide);

	WorkSpace wspace( sys );
	wspace.info();
	wspace.printPDB("rot_start.pdb");
	IO::OutTra_BTF& tra = wspace.addStdTra("rot_test");
	//IO::BTF_Block_Vector* vect = new IO::BTF_Block_Vector(5000);
	//tra.addOwnedBlock(vect);
	wspace.outtra.append();

	StatClock bob;
	bob.Begin();

	RotamerApplicator_SCWRL app(wspace,rot);
	//app.OutputLevel = Verbosity::Loud;
	app.apply();
	wspace.outtra.append();

	bob.End();
	bob.ReportMilliSeconds();

	wspace.printPDB("rot_end.pdb");

	int a;
	std::cout << "cooked";
	std::cin >> a;
	return;
}

void testGraphTheory()
{
	Maths::UndirectedGraph g(16);
	
	g.addEdge(0, 5);
	g.addEdge(0, 1);	
	g.addEdge(0, 6);
	g.addEdge(1, 2);
	g.addEdge(1, 3);
	g.addEdge(1, 4);
	g.addEdge(2, 3);
	g.addEdge(4, 5);
	g.addEdge(6, 8);
	g.addEdge(6, 7);	
	g.addEdge(7, 8);

	g.addEdge(2, 15);

	g.addEdge(15, 9);
	g.addEdge(9, 10);
	g.addEdge(10, 11);
	g.addEdge(11, 9);

	g.addEdge(4, 12);
	g.addEdge(13, 12);
	g.addEdge(13, 14);
	g.addEdge(5, 13);

	g.calcArticulationPoints();

	g.info();

	return;
}

void doMD()
{
	std::string traFile = "2ghf19to44_YPCPFCFKEFTRKDNMTAHVKIIHKI_100";

	StatClock bob;
	bob.Begin();

	FFParamSet ffps("lib/amber03aa.ff");
	ffps.readLib("lib/amber.zinc.ff");

	System sim(ffps);
	InTra_BTF tra( traFile + ".tra" );
	tra.loadIntoSystem( sim, -1 ); // load the final entry

	WorkSpace wspace( sim );
	wspace.info();
	wspace.addStdTra(traFile + ".sim.tra" );
	wspace.outtra.append();
	wspace.setAsReference(); // for CRMS calcs

	Monitors::PotentialEnergyMonitor mon1(wspace);
	Monitors::CRMSMonitor mon2(wspace);

	Forcefield ff = createff( wspace );

	Minimisation min2( ff );
	min2.Steps = 21;
	min2.StepSize = 1.0E7;
	min2.UpdateNList = 10;
	min2.UpdateScr = 10;
	min2.UpdateTra = 10;
	min2.UpdateMon = 0;
	min2.SlopeCutoff = 0.01 * PhysicsConst::kcal2J / PhysicsConst::Na;
	min2.run();

	MolecularDynamics md(ff);
	md.Integrator = MolecularDynamics::Langevin;
	md.Thermostat = MolecularDynamics::NoThermostat;
	md.Timestep   = 1.00E-15;
	md.UpdateScr = 100;
	md.UpdateNList = 10;
	md.UpdateTra = 100;
    md.setTargetTemp(298);
	md.Steps = 150;
	md.run();

	// Now set up for production run!
	md.Steps = 101; // 0.01ns
	md.RandVel = false; // CRITICAL !!

	for( int i = 0; i < 20; i++ )
	{
		SnapShot s = wspace.save();
		min2.run();
		mon1.measure();
		mon2.measure();
		wspace.load(s);
		md.run();
	}

	wspace.printPDB(traFile + ".post-md.pdb");
		
	printf("EPot Energy monitor:\n");
	mon1.printHeader();
	mon1.printAllData();
	mon1.printHistogram(10);
	printf("EPot average:\n");
	mon2.printRunningAverage();

	printf("CRMS monitor:\n");
	mon2.printHeader();
	mon2.printAllData();
	mon2.printHistogram(0.2);
	printf("CRMS average:\n");
	mon2.printRunningAverage();

	return;
}

void coiledCoilMake()
{
	StatClock bob;
	bob.Begin();

	FFParamSet ffps("lib/amber03aa.ff");

	Library::RotamerLibrary rot(ffps);
	rot.convertLib("lib/scl-B30-occ1.0-rmsd1.0-prop20.0",Library::RotLibConvert_Shetty());

	std::string seq1String = "*K-(VDEERYDIEAKVTKNITEIADLTQ)-K*";
	std::string seq2String = "*K-(VDEERYDIEAKVTKNITEIADLTQ)-K*";
	
	Sequence::BioSequence seq1(ffps);
	seq1.setTo(seq1String);
	Sequence::BioSequence seq2(ffps);
	seq2.setTo(seq2String);

	// load only polypeptide from the first model using the override sequences!
	PDB_In sim(ffps, "ccoil_template.pdb" );	
	sim.setAlignerDirect(); // we provide the alignment by specifying the sequence!
	sim.loadExplicit('B',Library::Polypeptide,seq1);
	sim.loadExplicit('C',Library::Polypeptide,seq2);

	WorkSpace wspace( sim );
	wspace.info();
	wspace.addStdTra("ccoil");
	wspace.outtra.append();
	wspace.printPDB("pre-packed.pdb");

	// Pack all except the coordinating residues and zinc!
	RotamerApplicator_SCWRL app(wspace,rot);
	app.OutputLevel = Verbosity::Loud;
	app.apply();
	wspace.outtra.append();

	// Create the primary forcefields
	Forcefield ffs = createffs(wspace);
	Forcefield ff = createff(wspace);

	// Now minimise all chains
	DualFFMinimiser eval(ff, ffs, PickSidechains() );
	eval.StepSize /= 10;
	eval.StericMinSteps = 10001;
	eval.StericSlopeCutoff = 0.1 * PhysicsConst::kcal2J / PhysicsConst::Na;
	eval.Steps = 10001;
	eval.StepSize *= 10;
	eval.SlopeCutoff = 0.05 * PhysicsConst::kcal2J / PhysicsConst::Na;
	eval.UpdateScr = 10;
	eval.UpdateTra = 10;
	eval.UpdateMon = 0;
	eval.run();

	// Now minimise the WORLD!!
	eval.setPickEverything();
	eval.run();

	// Cooked!
	wspace.printPDB("built.pdb");
	wspace.outtra.append();

	bob.End();
	bob.ReportMilliSeconds();

	return;
}

void zincFingerPrimaryModelBuilder()
{	
	StatClock bob;
	bob.Begin();

	FFParamSet ffps("lib/amber03aa.ff");
	ffps.readLib("lib/amber.zinc.ff");

	Library::RotamerLibrary rot(ffps);
	//rot.convertLib("lib/rotlib/shetty/scl-B30-occ1.0-rmsd1.0-prop20.0",Library::RotLibConvert_Shetty());
	//rot.convertLib("lib/rotlib/shetty/scl-B30-occ1.0-rmsd0.5-prop20.0",Library::RotLibConvert_Shetty());
	//rot.writeLib("lib/rotlib/shetty.rotamer");
	rot.readLib("lib/rotlib/shetty.rotamer");
	
	// CYS residues bonded to zinc are denoted CYM - they are -vely charged.
	Sequence::BioSequence seq(ffps);
	seq.setTo("*K-(Y)-CYM-(ST)-CYM-(DISFNYVKTYLA)-HID-(KQFY)-HID-(K)-N*");

	PDB_In sim(ffps, "dek/zinc_finger.pdb" );
	// load only polypeptide from the first model using the override sequence!
	sim.loadExplicit('A',Library::Polypeptide,seq);

	// Bond my onion!
	Molecule& mol = sim.getMolecule(0);
	mol.append( NewMolecule( ffps, "ZN" ) ); // add our zinc ion. Added arbritrarily at (0.0,0.0,0.0) (ffps default).
    int zincResIndex = mol.nResidues()-1; // we just pushed it there ...
	int zincIndex =  mol.findParticle( zincResIndex, "ZN" );
	mol.addBond( mol.findParticle( 2, "SG" ), zincIndex );
	mol.addBond( mol.findParticle( 5, "SG" ), zincIndex );
	mol.addBond( mol.findParticle( 18, "NE2" ), zincIndex );
	mol.addBond( mol.findParticle( 23, "NE2" ), zincIndex );

	WorkSpace wspace( sim );
	wspace.info();
	wspace.addStdTra("dek/zinc_finger");
	wspace.outtra.append();

	PickResidueList pickCoordinatingRes(wspace);
	pickCoordinatingRes.add( mol.res[2] ); // CYS
	pickCoordinatingRes.add( mol.res[5] ); // CYS
	pickCoordinatingRes.add( mol.res[18] ); // HIS
	pickCoordinatingRes.add( mol.res[23] ); // HIS
	pickCoordinatingRes.add( mol.res[26] ); // ZN itself

	Forcefield ffs = createffs(wspace);

	Forcefield ff(wspace);
	FF_Bonded bonds(wspace);
	ff.add( bonds );

	// SLAM that zinc in place!! Moving ONLY the zinc itself
	Minimisation zincPlacementMinimisation(ff, PickAtomIndex(zincIndex) );
	zincPlacementMinimisation.Steps = 20000;
	zincPlacementMinimisation.UpdateTra = 10;
	zincPlacementMinimisation.UpdateScr = 10;
	zincPlacementMinimisation.UpdateNList = 10;
	zincPlacementMinimisation.SlopeCutoff = 0.01 * PhysicsConst::kcal2J / PhysicsConst::Na;
	zincPlacementMinimisation.run();

	zincPlacementMinimisation.setPicking( Pick_AND( PickAtomIndex(zincIndex), pickCoordinatingRes ) );
	zincPlacementMinimisation.run();

	FF_GeneralizedBorn_Still* gbsa = new FF_GeneralizedBorn_Still( wspace ); // used to take nb
	gbsa->FastMode = true;
	ff.addWithOwnership( gbsa );

	FF_SASA_LCPO* sasa = new FF_SASA_LCPO(wspace);
	sasa->GlobalASP = 0.009;
	ff.addWithOwnership( sasa );

	// Now, NOT pickCoordinatingRes, cos we want to apply rotamers to residues not
	// involved in coordination!
	pickCoordinatingRes.invertSelection();

	// Pack all except the coordinating residues and zinc!
	RotamerApplicator_SCWRL app(wspace,rot, pickCoordinatingRes );
	app.OutputLevel = Verbosity::Loud;
	app.apply();
	wspace.outtra.append();

	// Now minimise all chains
	DualFFMinimiser eval(ff, ffs, PickSidechains() );
	eval.SDPreMinSteps = 51;
	eval.StericMinSteps = 10001;
	eval.StericSlopeCutoff = 0.1 * PhysicsConst::kcal2J / PhysicsConst::Na;
	eval.Steps = 10001;
	eval.StepSize *= 10;
	eval.SlopeCutoff = 0.05 * PhysicsConst::kcal2J / PhysicsConst::Na;
	eval.UpdateScr = 10;
	eval.UpdateTra = 10;
	eval.UpdateMon = 0;
	eval.run();

	// Now minimise the WORLD!!
	eval.setPickEverything();
	eval.run();

	// Cooked!
	wspace.printPDB("built.pdb");
	wspace.outtra.append();

	bob.End();
	bob.ReportMilliSeconds();

	return;
}
