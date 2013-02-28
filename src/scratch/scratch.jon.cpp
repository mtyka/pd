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
#include "mmlib.h"
#include "tools/io.h"
#include "../pd/accessory.h"
#include "arcus.h"
#include "scratch.jon.basics.h"
//#include "scratch.jon.more.h"

using namespace std;
using namespace Physics;
using namespace Protocol;
using namespace IO;
using namespace Manipulator;

int main_inner(int argc, char** argv);
int main(int argc, char** argv)
{
	PrintFullpdHeader();
	
#ifdef _DEBUG
	return main_inner(argc,argv);
#else
	try { return main_inner(argc, argv); }
	catch( ExceptionBase ex ) 
	{ return -1; }
#endif
}

void TestStitcher_ToStage2(int argc, char** argv);
void TestRefiner_FromStage2(int argc, char** argv);

void Test_TorsionalMinimisation()
{
	FFParamSet ffps;
	ffps.readLib("lib/amber03aa.ff");

	//System sys(ffps);
	//sys.add(NewProteinHelix(ffps,"*A-(AAAAPAAAA)-A*" ));

	PDB_In sys(ffps,"./pd_double/test.pdb");
	sys.load();

	WorkSpace wspace = WorkSpace( sys );
	RotBond replacement; // Must stay in scope for the lifetime of this wspace...
	wspace.setRotatableBondList( &replacement ); // Override the internal dummy rotbond array
	wspace.info();
	wspace.addStdTra("imported464");
	wspace.outtra.append();

	Forcefield ffs = createffts(wspace);

	SegmentDefBase sd(wspace,37,44,SegBreakCentre);

	TorsionalMinimisation min2( ffs, sd );
	min2.Steps = 3000;
	min2.InitialCapFactor = 0.05;
	min2.UpdateScr = 1;
	min2.UpdateTra = 1;
	min2.UpdateMon = -1;
	min2.run();

	return;
}

void temp()
{
	// set a standard random number seed for this test.
	Maths::FastRandom::getInstance()->reinitialise(100);

	std::string libPath = "../exec/";

	// Load our forcefield parameters
	FFParamSet ffps;
	ffps.readLib(libPath + "lib/amber03aa.ff");

	Library::RotamerLibrary rotLib( ffps );
	//rotLib.readLib( libPath + "lib/rotlib/shetty.rotamer" );
	rotLib.convertLib(libPath + "lib/rotlib/shetty/scl-B30-occ1.0-rmsd1.0-prop20.0", Library::RotLibConvert_Shetty() );
	rotLib.writeLib( libPath + "lib/rotlib/shetty.rotamer" );

	// Load our PDB
	PDB_In sys(ffps,"the_exception2.pdb");
	sys.disableRebuild(); // Not critical, but will otherwise be repeated - ArcusBase will do the work.
	sys.CentreOnBuild = true; // Not critical, moves valid particles to the COG post-build.
	sys.load();

	// Make a full workspace from the imported system
	WorkSpace wspace( sys );
	wspace.info();
	RotBond replacement; // Must stay in scope for the lifetime of this wspace...
	wspace.setRotatableBondList( &replacement ); // Override the internal dummy rotbond array
	wspace.addStdTra("parp");

	SegmentDef _region(wspace,99,106,SegBreakCentre);

	PickAtomRanges m_DynamicAtomsPicker;
	size_t breakRes = _region.getBreakResIndex();
	size_t startAtomIndex = _region.getStartAtomIndex();
	size_t centreAtomIndexA = wspace.res[breakRes-1].ilast;
	m_DynamicAtomsPicker.addRange(startAtomIndex,centreAtomIndexA-startAtomIndex+1,false);
	size_t endAtomIndex = _region.getEndAtomIndex();
	size_t centreAtomIndexB = wspace.res[breakRes].ifirst;
	ASSERT( centreAtomIndexA == centreAtomIndexB-1, CodeException, "Oddness happened" );
	m_DynamicAtomsPicker.addRange(centreAtomIndexB,endAtomIndex-centreAtomIndexB+1,true);

	PickResidueList move(wspace);
	move.add( 99 ); 
	move.add( 100 ); 
	move.add( 101 ); 
	move.add( 102 ); 
	move.add( 103 ); 
	move.add( 104 ); 
	move.add( 105 ); 
	move.add( 106 ); 

	ProximityGrid m_StaticGrid4A(wspace,Pick_NOT(m_DynamicAtomsPicker),4.0);
	Manipulator::RotamerApplicator_SCWRL rot( wspace, rotLib, 
		move
		);
	rot.OutputLevel = Verbosity::Loud;
	rot.overrideStaticGrid(m_StaticGrid4A); // rotamers can only now clash with the rigid-body

	rot.apply();

	PDB_Writer("parp.pdb").write(wspace); // pdb

	return;
}

void tempSegDistFan()
{
	// set a standard random number seed for this test.
	Maths::FastRandom::getInstance()->reinitialise(100);

	std::string libPath = "../exec/";

	// Load our forcefield parameters
	FFParamSet ffps;
	ffps.readLib(libPath + "lib/amber03aa.ff");

	// Load our anglsest
	Library::AngleSet as;
	as.loadFromFile(libPath + "lib/angleset/big.loop.angleset");
	as.info();

	// Load our PDB
	PDB_In sys(ffps,"../pdb/1NLS_.pdb");
	sys.disableRebuild(); // Not critical, but will otherwise be repeated - ArcusBase will do the work.
	sys.CentreOnBuild = true; // Not critical, moves valid particles to the COG post-build.
	sys.load();

	// Make a full workspace from the imported system
	WorkSpace wspace( sys );
	wspace.info();
	RotBond replacement; // Must stay in scope for the lifetime of this wspace...
	wspace.setRotatableBondList( &replacement ); // Override the internal dummy rotbond array
	IO::OutTra_BTF& tra = wspace.addStdTra("testSegDistFan");
	IO::BTF_Block_Vector* vect = new IO::BTF_Block_Vector(); // due to a bug in dave, BTF_Block_Vector must be added 1st...
	tra.addOwnedBlock(vect);
	IO::BTF_Block_Comment* cmnt = new IO::BTF_Block_Comment();
	tra.addOwnedBlock(cmnt);
	tra.append();

	SegmentDef region(wspace,96,103,SegBreakCentre);
	SegmentDef subRegion( wspace, region.getStartResIndex(), region.getCentreResIndex(), SegBreakEnd );
	//SegmentDef subRegion( wspace, region.getCentreResIndex()+1, region.getEndResIndex(), SegBreakStart );

	region.info();
	subRegion.info();

	bool reversed = subRegion.getStartResIndex() > region.getStartResIndex();
	SegmentDistanceFilter filter( libPath + "lib/segment/segdist.dat" );
	filter.initialise( region, reversed );

	// Assert -->
	filter.draw( *vect );
	wspace.outtra.append();
	std::cout << filter.reason();
	ASSERT( filter.passes(), CodeException, "Native filtered!!");

	ConfBuild_RandomSingle conf( as, subRegion );
	conf.PropensityWeighting = true;
	conf.NoiseSigma = Maths::DegToRad(15.0);
	conf.setRandCount(1000);
	conf.setBuildMode(Manipulator::ConfBuilderBase::Backbone);

	int filt = 0;
	while( conf.next() )
	{
		if( !filter.passes() )
		{
			filt++;
			std::string reason = filter.reason();
			cmnt->setFormat( reason );
			std::cout << reason << std::endl;
			filter.draw( *vect );
			wspace.outtra.append();
		}
		//else
		//{
		//	//filt++;
		//	//std::string reason = filter.reason();
		//	//cmnt->setFormat( reason );
		//	//std::cout << reason << std::endl;
		//	//filter.draw( *vect );
		//	//wspace.outtra.append();
		//}		
	}

	return;
}

int main_inner(int argc, char** argv)
{
	//temp();
	//tempSegDistFan();
	//Test_TorsionalMinimisation();
	TestStitcher_ToStage2(argc,argv);
	//TestRefiner_FromStage2(argc,argv);
	return 0;
}

std::string traPathToStem( const std::string& traPath, double dielec )
{
	StringBuilder stem(traPath);
	size_t at;
	while( SIZE_T_FAIL != (at = stem.FirstOf('/') ) )
	{
		stem.replace( at, '\\' );
	}
	at = stem.LastOf( '\\' );
	if( at != SIZE_T_FAIL )
		stem.TruncateLeftBy( at + 1 );
	at = stem.FirstOf( '.' );
	if( at != SIZE_T_FAIL )
		stem.TruncateRightTo( at );

	stem.append('_');
	stem.append( dielec, "%4.2f" );

	return stem.toString();
}

std::string pdbPathToStem( const std::string& pdbPath, int startRes, int length, int branchNum, int execID )
{	
	StringBuilder stem;
	stem.append(execID);
	stem.append('_');
	stem.append( pdbPath );
	size_t at;
	while( SIZE_T_FAIL != (at = stem.FirstOf('/') ) )
	{
		stem.replace( at, '\\' );
	}
	at = stem.LastOf( '\\' );
	if( at != SIZE_T_FAIL )
		stem.TruncateLeftBy( at + 1 );
	at = stem.FirstOf( '.' );
	if( at != SIZE_T_FAIL )
		stem.TruncateRightTo( at );

	stem.append('_');
	stem.append(startRes);
	stem.append('_');
	stem.append(length);
	stem.append('_');
	stem.append(branchNum);

	return stem.toString();
}

void TestRefiner_FromStage2( int argc, char** argv)
{
	ASSERT( argc == 6, ArgumentException, "Give me traname, startres, length, branchCount\n" );

	std::string traFile( argv[1] );
	int startRes = atoi( argv[2] );
	int nRes = atoi( argv[3] );
	double dielec = atof( argv[4] );
	int execID = atoi( argv[5] );
	int endRes = startRes + nRes - 1;	
	std::string stem = traPathToStem( traFile, dielec ) + "_out";

	// set a standard random number seed for this test.
	Maths::FastRandom::getInstance()->reinitialise(100);

	std::string libPath = "./";

	// Load the replay InTra
	if( !IO::fileExists( traFile ) )
	{
		throw IOException("Cant find tra!");
	}

	// Load our anglsest
	Library::AngleSet as;
	as.loadFromFile(libPath + "big.loop.angleset");
	as.info();

	// Load our forcefield parameters
	FFParamSet ffps;
	ffps.readLib(libPath + "amber03aa.ff");

	Library::RotamerLibrary rotLib( ffps );
	//rotLib.readLib( libPath + "lib/rotlib/shetty.rotamer" );
	// Use the HIGH resolution lib for refinement
	rotLib.convertLib(libPath + "scl-B30-occ1.0-rmsd0.5-prop20.0", Library::RotLibConvert_Shetty() );
	//rotLib.convertLib(libPath + "scl-B30-occ1.0-rmsd0.2-prop20.0", Library::RotLibConvert_Shetty() );
	//rotLib.convertLib(libPath + "scl-B30-occ1.0-rmsd1.0-prop20.0", Library::RotLibConvert_Shetty() );
	//rotLib.writeLib( libPath + "shetty.rotamer" );

	// Load our tra
	System sys(ffps);
	IO::InTra_BTF traInput( traFile );
	traInput.loadIntoSystem( sys, -1 );

	// Make a full workspace from the imported system
	WorkSpace wspace( sys );
	wspace.info();
	RotBond replacement; // Must stay in scope for the lifetime of this wspace...
	wspace.setRotatableBondList( &replacement ); // Override the internal dummy rotbond array
	IO::OutTra_BTF& tra = wspace.addStdTra(stem); // tra
	IO::BTF_Block_Comment* cmnt = new IO::BTF_Block_Comment();
	tra.addOwnedBlock(cmnt);

	// Write what we have interpreted from the import file	
	PDB_Writer(stem + ".imported.pdb").write(wspace); // pdb

	// Forcefield config
	bool useBrokenBonded = true; // Flag the use of the bonded forcefield that allows bond breaks.
	Forcefield ffs = createffts(wspace,useBrokenBonded,false);
	Forcefield ff = createff(wspace,useBrokenBonded,dielec,false);

	SegmentDef def(wspace,startRes,endRes,SegBreakCentre);
	def.info();

	// Loop stitcher
	Protocol::ArcusReplay loopBuild(ffs,ff,as,rotLib,traInput);
	loopBuild.UpdateNList = 5;
	loopBuild.OutputLevel = Verbosity::Eleven; // I love that flag :-D
	loopBuild.enableComments( cmnt );
	loopBuild.buildAddExplicit( def ); // turn off auto-detection of what to build and explicitly define it

	// Override to split protocol ...
	ArcusRefine_CGMin minimisingRefiner;
	minimisingRefiner.OutputLevel = Verbosity::Normal;
	minimisingRefiner.enableRotamerPack( rotLib );
	minimisingRefiner.UpdateNList = 10;
	loopBuild.setRefiner( minimisingRefiner );

	Printf("All setup complete... Rolling the primary loop. Hold onto your hat!\n");

	int result = loopBuild.run();

	ff.calcEnergies();
	cmnt->setFormat("**Final Best** Structure According to FF! Ene: %8.3f, BB-cRMS: %8.3f, Heavy-cRMS: %8.3f\n")
		(wspace.ene.epot / (Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na))
		(def.cCRMS_BB())
		(def.cCRMS_Heavy());
	cmnt->toScreen();
	wspace.outtra.append();

	Printf("Complete! :-D\n");

	return;
}

void TestStitcher_ToStage2(int argc, char** argv)
{
	ASSERT( argc == 6, ArgumentException, "Give me only: pdbname, startres, length, branchCount\n" );

	std::string pdbFile( argv[1] );
	int startRes = atoi( argv[2] );
	int nRes = atoi( argv[3] );
	int ProduceNJoinedBranches = atoi( argv[4] );
	int execID = atoi( argv[5] );
	int endRes = startRes + nRes - 1;	
	std::string stem = pdbPathToStem( pdbFile, startRes, nRes, ProduceNJoinedBranches,execID );

	// set a standard random number seed for this test.
	Maths::FastRandom::getInstance()->reinitialise(100);

	std::string libPath = "./";

	// Load our anglsest
	Library::AngleSet as;
	as.loadFromFile(libPath + "big.loop.angleset");
	as.info();

	// Load our forcefield parameters
	FFParamSet ffps;
	ffps.readLib(libPath + "amber03aa.ff");

	Library::RotamerLibrary rotLib( ffps );
	//rotLib.readLib( libPath + "lib/rotlib/shetty.rotamer" );
	rotLib.convertLib(libPath + "scl-B30-occ1.0-rmsd1.0-prop20.0", Library::RotLibConvert_Shetty() );
	//rotLib.writeLib( libPath + "shetty.rotamer" );

	// Load our PDB
	PDB_In sys(ffps,pdbFile);
	sys.disableRebuild(); // Not critical, but will otherwise be repeated - ArcusBase will do the work.
	sys.CentreOnBuild = true; // Not critical, moves valid particles to the COG post-build.
	sys.load();

	// Make a full workspace from the imported system
	WorkSpace wspace( sys );
	wspace.info();
	RotBond replacement; // Must stay in scope for the lifetime of this wspace...
	wspace.setRotatableBondList( &replacement ); // Override the internal dummy rotbond array
	IO::OutTra_BTF& tra = wspace.addStdTra(stem); // tra
	IO::BTF_Block_Comment* cmnt = new IO::BTF_Block_Comment();
	tra.addOwnedBlock(cmnt);

	// Write what we have interpreted from the import file	
	PDB_Writer(stem + ".imported.pdb").write(wspace); // pdb

	// Forcefield config
	bool useBrokenBonded = true; // Flag the use of the bonded forcefield that allows bond breaks.
	Forcefield ffs = createffts(wspace,useBrokenBonded,false);
	Forcefield ff = createff(wspace,useBrokenBonded,false);

	SegmentDef def(wspace,startRes,endRes,SegBreakCentre);
	def.info();

	// Loop stitcher
	Arcus loopBuild(ffs,ff,as,rotLib);
	loopBuild.UpdateNList = 5;
	loopBuild.OutputLevel = Verbosity::Eleven; // i love that flag :-D
	loopBuild.SegFanFilename = libPath + "segdist.dat";
	loopBuild.ProduceNJoinedBranches = ProduceNJoinedBranches;
	//loopBuild.enableComments( cmnt );
	loopBuild.buildAddExplicit( def ); // turn off auto-detection of what to build and explicitly define it

	// override to split protocol ...
	ArcusRefine_ToTra storageRefiner( stem );
	loopBuild.setRefiner( storageRefiner );

	// Assert that basic filters are not removing the native!
	loopBuild.AssertCA6_32Filter();
	loopBuild.AssertSegDistFilter();

	Printf("All setup complete... Rolling the primary loop. Hold onto your hat!\n");

	loopBuild.nameMe = stem;
	int result = loopBuild.run();

	Printf("Complete! :-D\n");

	return;
}

