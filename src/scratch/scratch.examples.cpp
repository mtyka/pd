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
#include "tools/quote.h"

#include "workspace/workspace.h"
#include "workspace/segdef.h"

#include "sequence/sequence.h"
#include "sequence/alignment.h"

#include "forcefields/ffbonded.h"
#include "forcefields/breakablebonded.h"
#include "forcefields/ffnonbonded.h"
#include "forcefields/gbff.h"
#include "forcefields/lcpo.h"
#include "forcefields/ffsoftvdw.h"

#include "fileio/tra.h"
#include "fileio/pdb.h"

#include "protocols/minimise.h"
#include "protocols/torsionalminimisation.h"
#include "protocols/md.h"

using namespace std;
using namespace Physics;
using namespace Protocol;
using namespace Tra;
using namespace PDB;

Forcefield createffs( WorkSpace& wspace, bool useBreakableFF = false, bool summary = true)
{
	Forcefield ff = Forcefield(wspace);

	if( useBreakableFF )
	{
		FF_BreakableBonded* bonds = new FF_BreakableBonded(wspace);
		ff.addWithOwnership( *bonds ) ;
	}
	else
	{
		BondedForcefield* bonds = new BondedForcefield(wspace);
		ff.addWithOwnership( *bonds ) ;
	}

	SoftVDWForcefield *sff = new SoftVDWForcefield(wspace);
	//sff->Hardness = 6;
	sff->Hardness = 12;
	ff.addWithOwnership(*sff);

	if( summary ) ff.printEnergySummary();

	return ff;
}

Forcefield createffVac(WorkSpace& wspace, bool useBreakableFF = false, bool summary = true)
{
	Forcefield ff = Forcefield(wspace);

	if( useBreakableFF )
	{
		FF_BreakableBonded* bonds = new FF_BreakableBonded(wspace);
		ff.addWithOwnership( *bonds ) ;
	}
	else
	{
		BondedForcefield* bonds = new BondedForcefield(wspace);
		ff.addWithOwnership( *bonds ) ;
	}

	FF_NonBonded* nb = new FF_NonBonded(wspace);
	nb->Cutoff = 12.0;
	nb->InnerCutoff = 6.0;
	ff.addWithOwnership( *nb );

	if( summary ) ff.printEnergySummary();

	return ff;
}

Forcefield createff(WorkSpace& wspace, bool useBreakableFF = false, double dielec = 1.0, bool summary = true)
{
	Forcefield ff = Forcefield(wspace);

	if( useBreakableFF )
	{
		FF_BreakableBonded* bonds = new FF_BreakableBonded(wspace);
		ff.addWithOwnership( *bonds ) ;
	}
	else
	{
		BondedForcefield* bonds = new BondedForcefield(wspace);
		ff.addWithOwnership( *bonds ) ;
	}

	GB_Still* gbsa = new GB_Still( wspace ); // used to take nb
	gbsa->FastMode = true;
	gbsa->DielectricSolute = dielec;
	ff.addWithOwnership( *gbsa );

	SASA_LCPO* sasa = new SASA_LCPO(wspace);
	sasa->GlobalASP = 0.009;
	ff.addWithOwnership( *sasa );

	if( summary ) ff.printEnergySummary();

	return ff;
}

void Test_Minimisation()
{
	StatClock clockMe;
	clockMe.Begin();

	FFParamSet ffps;
	ffps.readLib("lib/amber03aa.ff");

	System sys(ffps);
	sys.add(NewProteinHelix(ffps,"*A-(AAAAPAAAA)-A*" ));

	WorkSpace wspace = WorkSpace( sys );
	wspace.info();
	wspace.addStdTra("imported464");

	Forcefield ff = createff(wspace);

	for( int i = 0; i < 1; i++ ) // average over 10 minimisations
	{
		Minimisation min( ff );
		min.Steps = 1500;
		min.StepSize = 2E1;
		min.UpdateNList = 5;
		min.UpdateScr = 100;
		min.UpdateTra = 0;
		min.UpdateMon = 0;
		min.run();

		clockMe.Stamp();
	}

	clockMe.ReportMilliSeconds(10);

	return;
}

void Test_MD()
{
	StatClock clockMe;
	clockMe.Begin();

	FFParamSet ffps;
	ffps.readLib("lib/amber03aa.ff");

	//System sys(ffps);
	//sys.add(NewProteinHelix(ffps,"*P-(ACDEFGHIKLMNPQRSTVWY)-P*" ));
	//sys.add(NewProteinHelix(ffps,"*A-(A)-A*" ));

	PDB_In sys( ffps, "pep4_A_84_237_ss.pdb" );
	sys.load();

	WorkSpace wspace = WorkSpace( sys );
	PDB_Out out1("imported1",wspace);
	wspace.info();

	wspace.addStdTra("imported1");

	Forcefield ff = createff(wspace);

	Minimisation min2( ff ); // Whole Protein
	min2.Steps = 3501;
	min2.StepSize = 2E1;
	min2.UpdateNList = 5;
	min2.UpdateScr = 20;
	min2.UpdateTra = 10;
	min2.UpdateMon = 0;
	min2.run();

	MolecularDynamics md(ff); // Whole Protein
	md.Steps = 5001;
	md.Integrator = MolecularDynamics::Langevin;
	md.Thermostat = MolecularDynamics::NoThermostat;
	md.Timestep   = 1.00E-15;
	md.UpdateScr = 20;
	md.UpdateNList = 10;
	md.UpdateTra = 20;
	md.UpdateMon = 0;
	md.TargetTemp = new Temp(300);
	md.run();

	clockMe.End();

	clockMe.ReportMilliSeconds(10);

	return;
}

void Test_TorsionalMinimisation()
{
	FFParamSet ffps;
	ffps.readLib("lib/amber03aa.ff");

	System sys(ffps);
	sys.add(NewProteinHelix(ffps,"*A-(AAAAPAAAA)-A*" ));

	WorkSpace wspace = WorkSpace( sys );
	PDB_Out out1("imported464",wspace);
	wspace.info();

	wspace.addStdTra("imported464");

	Forcefield ff = createff(wspace);

	SegmentDefBase sd(wspace,0,6,false);
	TorsionalMinimisation min2( ff, sd );

	min2.Steps = 3000;
	min2.UpdateScr = 1;
	min2.UpdateTra = 1;
	min2.UpdateMon = -1;
	min2.run();
}

void Test_PDB_In()
{
	FFParamSet ffps;
	ffps.readLib("lib/amber03aa.ff");
	ffps.readLib("lib/tip3.ff");

	string stem = "mike_pdb_sample\\pdb\\";
	std::vector<string> names;

	names.push_back("bba1.pdb"); // Correctly complains about PYA
	names.push_back("water.pdb"); // Correctly throws an exception - MOL is not a molecule in amber03aa.ff
	names.push_back("trpzip1.pdb"); // Good
	names.push_back("sh3.pdb"); // Good
	names.push_back("acyltra.pdb"); // Good
	names.push_back("1fsd.pdb"); // Good
	names.push_back("villin.pdb"); // Good
	names.push_back("protg.pdb"); // Good
	names.push_back("proteinAZ.pdb");// Good
	names.push_back("ubifrag.pdb"); // Good
	names.push_back("trpzip.pdb"); // Good
	names.push_back("trpcage.pdb"); // Good
	names.push_back("alcalase.pdb"); // Good
	names.push_back("1hz6.pdb"); // Good
	names.push_back("trpzip2.pdb"); // Good

	// 1) Try loading everything from each file
	for( size_t i = 0; i < names.size(); i++ )
	{
		try
		{
			// Read in
			PDB_In sysA(ffps,stem + names[i]);
			sysA.load(); // Load everything from model 1

			// Write out
			PDB_Writer writerA(stem + names[i] + ".out.pdb");
			writerA.write(sysA);
		}
		catch(exception ex)
		{
			printf("\nCritical Load Failure!!!!\n");
		}
	}


	// 2) An example of using custom name mappings for oddball residues
	PDB_In sysB(ffps,stem + names[0]);
	sysB.addAlias("PYA","ALA"); // Define a brand new custom alias
	sysB.load();

	PDB_Writer write(stem + "bba1.pdb" + ".customalias.pdb");
	write.write(sysB);

	// 3) Content filters
	PDB_In sysC(ffps,"1dqe.pdb");
	sysC.setFilter(Library::Polypeptide);
	sysC.load(); // loads only polypeptide from model 1

	PDB_In sysD(ffps,"1dqe.pdb");
	sysD.setFilter(Library::Water);
	sysD.load(); // loads only the water molecules from model 1

	PDB_In sysE(ffps,"1dqe.pdb");
	sysE.getClass().addClassMember("custom1","ALA"); // ALA is a member of class 'custom1'
	sysE.getClass().addClassMember("custom1","GLU"); // GLU is a member of class 'custom1'
	sysE.getClass().addClassMember("custom1","ARG"); // ARG is a member of class 'custom1'
	sysE.setFilter("custom1"); // loads only from the custom gilter
	sysE.load(); // loads only the specified molecules from model 1

	// 4) Load specifics - reolace "1dqe.pdb" with your file
	PDB_In sysF(ffps,"1dqe.pdb");
	sysF.loadAll(); // load absolutely fricking everything from all models at once - this could get messy!

	PDB_In sysG(ffps,"1dqe.pdb");
	sysG.load('A'); // load chain A from model 1

	PDB_In sysH(ffps,"1dqe.pdb");
	sysH.load('A',4); // load chain A from model 4

	PDB_In sysI(ffps,"1dqe.pdb");
	sysI.loadExplicit('A',"BOM",47,' '); // Load BOM 47 with no icode from chain A of model 1

	PDB_In sysJ(ffps,"1dqe.pdb");
	Sequence::BioSequence bioseq;
	bioseq.setTo("*P-(AAAAAAAAAAAAPAAAAAAAAAA)-P*"); // Map the particles defined in the PDB over this sequence!
	sysJ.loadExplicit('C',Library::Polypeptide,bioseq); // Load polypeptide from chain C using a sequece override

	return;
}

void Test_PDB_Sim( const std::string& fileStem, bool doMD = false )
{
	// Load FFParam
	FFParamSet ffps;
	ffps.readLib("lib/amber03aa.ff");

	// Load our favourite PDB molecule
	PDB_In sys(ffps,fileStem + ".pdb");
	sys.setFilter(Library::Polypeptide); // ONLY load the polypeptide
	sys.load(); // load from model 1 (any chains)

	// Make a workspace
	WorkSpace wspace = WorkSpace( sys );
	wspace.addStdTra(fileStem);
	wspace.info();

	// Write out what we have imported to check
	PDB_Writer writer(fileStem + ".imp.pdb",false);
	writer.write(wspace);

	// Create the bond only forcefield
	Forcefield ffb = Forcefield(wspace);
	BondedForcefield* bonds = new BondedForcefield(wspace);
	ffb.addWithOwnership( *bonds );

	Minimisation minb( ffb );
	minb.Steps = 50000;
	minb.SlopeCutoff = 0.1 * PhysicsConst::kcal2J / PhysicsConst::Na;
	minb.StepSize = 1E2;
	minb.UpdateScr = 10;
	minb.UpdateTra = 10;
	minb.UpdateMon = 10;
	minb.run();

	// Create the nice steric forcefield
	Forcefield ffs = createffs(wspace);

	Minimisation mins( ffs );
	mins.Steps = 2000;
	mins.SlopeCutoff = 0.001 * PhysicsConst::kcal2J / PhysicsConst::Na;
	mins.StepSize = 0.4E1;
	mins.UpdateScr = 10;
	mins.UpdateTra = 10;
	mins.UpdateMon = 10;
	mins.run();

	// Create the full forcefield
	Forcefield ff = createff(wspace);

	// Minimise - PDB files often contain bad clashes
	Minimisation min( ff );
	min.Steps = 500;
	min.StepSize = 2E1;
	min.UpdateScr = 10;
	min.UpdateTra = 10;
	min.UpdateMon = 10;
	min.run();

	if( doMD )
	{
		// Do a nice MD simulation
		MolecularDynamics md(ff);
		md.Steps = 5000;
		md.UpdateScr = 10;
		md.UpdateTra = 10;
		md.UpdateMon = 10;
		md.run();
	}

	// Whats the workspace like now?
	wspace.info();

	PDB_Writer writerEnd(fileStem + ".min.pdb",false);
	writerEnd.write(wspace);

	return;
}

void Test_Tra_In()
{
	FFParamSet ffps;
	ffps.readLib("lib/amber03aa.ff");

	System sim2( ffps );
	Tra::Traj_In mytra( "test.tra", true );
	mytra.loadIntoSystem( sim2, -1 );

	WorkSpace wspace2( sim2 );
	PDB_Out out("imported",wspace2);

	Forcefield ff2 = createff(wspace2);
	wspace2.info();

	Minimisation min2( ff2 );
	min2.Steps = 1;
	min2.UpdateScr = 10;
	min2.UpdateTra = 10;
	min2.UpdateMon = 10;
	min2.run();

	wspace2.info();

	out.append();
}

void Test_Alignment_Code()
{
	FFParamSet ffps;
	ffps.readLib("lib/amber03aa.ff");

	// Forcefield resolution
	Sequence::BioSequence s1(ffps);
	char* seq1 = "*s-ASP-GLU-AlA-Gly-ASP-D-GLU-ALA-GLU-A*";
	printf("Resolving using 'FFParamSet':\n");
	printf("Seq1: '%s'\n",seq1);
	s1.setTo(seq1);
	printf("ResolvesTo: '");
	s1.printToScreen();
	printf("'\n\n");

	Library::NamingConventions* map2 = Library::NamingConventions::getSingleton();
	Sequence::BioSequence s2(*map2);
	char* seq2 = "GLU-AlA-G-ASP-D-AlA-GLU";
	printf("Resolving using Library::NamingConventions:\n");
	printf("Seq2: '%s'\n",seq2);
	s2.setTo(seq2);
	printf("ResolvesTo: '");
	s2.printToScreen();
	printf("'\n\n");

	Sequence::ExpSeqPair expPair(s1,s2);
	Sequence::SimpleAligner ali;
	ali.Align(expPair);
	expPair.printToScreen();

	return;
}

void Test_StatClock()
{
	int count = 20;
	TextProgressBar bar( count, 50 );
	StatClock stats("Clock Test");
	bar.Begin();
	stats.Begin();
	for( int i = 0 ; i < count; i++ )
	{
		system("sleep 0.1"); // replace with desired OutputLevel! operation to be timed - must not print to screen!
		bar.next();
		stats.Stamp();
	}
	bar.End();
	stats.ReportSeconds();
}

void Test_CharacterMap()
{
	int printed = 0;
	for( int i = 1; i < 256; i++ )
	{
		if( (i >= 7 && i <= 13) || i == 27 ) continue;
		if( i < 10 ) std::cout << ' ';
		if( i < 100 ) std::cout << ' ';
		std::cout << i << std::string(": ") << (char)i;
		printed++;
		if( (printed % 10 == 0) || i == 255 )
		{
			std::cout << std::endl;
		}
		else
		{
			std::cout << std::string(", ");
		}
	}
	std::cout << std::endl;
	std::cout << std::endl;
}

void Test_Quote()
{
	Quote* q = Quote::getSingleton();
	q->printQuote();
	q->printQuote();
	q->printQuote();
	q->printQuote();
	q->printQuote();
	q->printQuote();
	q->printQuote();
	q->printQuote();
	q->printQuote();
}

void Test_StringSystem()
{
	double d = 1.723625154;
	int i = 1;
	Printf("%d) Using 'safe Printf' with the same format strings as printf(), you can print numbers: %lf\n")(i++)(d);

	StringBuilder sb;
	sb.setFormat("%d) StringBuilder now supports format strings!")(i++);
	std::cout << sb << std::endl;

	sb.clear();
	float f = 4.33243343f;
	sb.append("Jon says: ");
	sb.appendFormat("%d %s %5.3f")(2343)(" Derek! ")(f);
	sb.endl();

	// Printing : StringBuilder supports both output methods
	std::cout << sb;
	sb.toScreen(); // Print whole thing
	sb.toScreen(16,14); // or subranges!

	std::cout << std::endl << std::endl;
}

int main_inner(int argc, char** argv);
int main( int argc, char** argv )
{
	//if( argc <= 1 ) return -1;
	PrintFullpdHeader();

#ifdef _DEBUG
	return main_inner(argc,argv);
#else
	try { return main_inner(argc,argv); }
	catch( ExceptionBase ex )
	{ return -1; }
#endif
}

int main_inner(int argc, char** argv)
{
	// Text and tools
	//Test_StringSystem();
	//Test_Quote();
	//Test_CharacterMap();
	//Test_Alignment_Code();
	//Test_StatClock();

	// Simulations
	//Test_MD();
	//Test_Minimisation();
	//Test_TorsionalMinimisation();

	// FileIO
	//Test_PDB_In();
	Test_PDB_Sim("2drp40to65");
	//Test_Tra_In();

	return 0;
}

