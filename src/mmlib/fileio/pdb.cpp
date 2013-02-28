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
#include <iomanip> // Needed for some stream manipulation functionality e.g. setw()
#include "revision.h"
#include "tools/stringbuilder.h"
#include "system/genpolymer.h"
#include "system/system.h"
#include "fileio/infile.h"
#include "pdb.h"

using namespace std;
using namespace Sequence;

namespace IO 
{
	// -----------------------------------------------
	// 'PDB_SEQRES'
	// Responsible for importing from and exporting to
	// SEQRES tags in PDB files.
	// -----------------------------------------------

	PDB_SEQRES::PDB_SEQRES() : BioSequenceCollection<char>(), m_ParserNameMapper(NULL)
	{
		resetStatics();
	}

	void PDB_SEQRES::resetStatics()
	{
	}

	void PDB_SEQRES::setAlias(const Library::AliasMapper& _ParserNameMapper)
	{
		m_ParserNameMapper = &_ParserNameMapper;
	}

	const Library::AliasMapper& PDB_SEQRES::getAlias() const
	{
		return *m_ParserNameMapper;
	}

	void PDB_SEQRES::clear()
	{
		BioSequenceCollection<char>::clear();
	}

	void PDB_SEQRES::loadFromFile( const std::string &_filename )
	{
		std::ifstream file(_filename.c_str(), ifstream::in);
		if( !file.is_open() ) throw(IOException( "PDB-SEQRES file not found: '" + _filename + "'!" ));

		std::string line;
		while(std::getline(file,line))
		{
			parseLine(line);
		}

		file.close(); // close the m_Stream to release the file_handle
	}

	void PDB_SEQRES::parseLine( const std::string &_Line )
	{
		static StringBuilder lineParser(81); // 81 from: 1 + Zero Termination During Parsing + 80 is the normal length of a string in these files...
		static StringBuilder resNameBuilder(3); // StringBuilder usable to temporarily hold the 3 residue chars

		static int serNum = 1;
		static char chainID = CHAR_MAX;
		static int foundResCount = 0;
		static int totalResCount = INT_MAX;

		// get out line into a buffer
		lineParser.setTo( _Line );
		lineParser.PadRight(80); // ensure that our string is a minimum length

		// Extract and verify the relevent line header writermation
		if( !lineParser.compare("SEQRES",0,6,0,false) ) 
		{
			throw ArgumentException("Non-'SEQRES' line sent to the PDB_SEQRES class!");
		}

		char newChainID = lineParser[11];
		if( newChainID != chainID )
		{
			serNum = 1;
			totalResCount = lineParser.parseInt(13);
			foundResCount = 0;
			chainID = newChainID;
			addSequence(chainID);
		}
		else
		{
			serNum++;
			int newSerNum = lineParser.parseInt(8);
			if( newSerNum != serNum ) printf("WARNING: 'SEQRES' line serial is incorrect! There may be sequence errors!\n");
			int newTotalResCount = lineParser.parseInt(13);
			if( newTotalResCount != totalResCount ) printf("WARNING: 'SEQRES' line residue total is incorrect! There may be sequence errors!\n");
		}

		// get the sequence we want to add residue definitions to
		BioSequence &seq = getLastSequence();

		int pos = 19; // start here
		while( foundResCount < totalResCount )
		{
			//pos is given by 19 + (4 * nameSlotIndex);
			resNameBuilder.setTo(lineParser,pos,3);

			// check that we actually have the corresponding number of amino acids here!
			// 0,1,2 as DNA can be listed as " A " "A " " A" etc....
			if( resNameBuilder.compare("   ", 0, 3, 0, false) )
			{
				THROW(ProcedureException,"Error in 'SEQRES' parsing, the parseable amino acid count is incorrect!");
			}

			std::string resName = resNameBuilder.toString();
			char resID = '?';
			if( m_ParserNameMapper != NULL )
			{
				m_ParserNameMapper->lookupAlias(resName,resName); // rename the residue if there the current name a valid alias
				m_ParserNameMapper->lookupShortName(resName,resID);
			}
			seq.append( ResidueInfo( resID, foundResCount, resName ) );

			pos += 4;
			foundResCount++;
			if( (foundResCount % 13) == 0 ) break; // 13 residues per parseLine *max*
		}
	}

	const Library::AliasMapper* PDB_SEQRES::getNameMapper() const
	{
		return m_ParserNameMapper;
	}





	// -----------------------------------------------
	// PDB-Specific Parsing Helper Function
	// -----------------------------------------------

	inline bool getLineType(const std::string &line, std::string &lineType)
	{
		if( line.length() < 6 ) return false; // too short to be valid
		lineType[0] = toupper(line[0]);
		lineType[1] = toupper(line[1]);
		lineType[2] = toupper(line[2]);
		lineType[3] = toupper(line[3]);
		lineType[4] = toupper(line[4]);
		lineType[5] = toupper(line[5]);
		return true;
	}





	// -----------------------------------------------
	// PDBAtomLine
	// -----------------------------------------------

	PDBAtomLine::PDBAtomLine()
	{
		init();
	}

	void PDBAtomLine::init()
	{
		lineType = ATOM;  // lineType (1-6)
		atomNum = 0;      // Atom Number (7-11)
		atomName = "";    // Atom Name (13-16)
		altLoc = '\0';    // AltLoc Indicator (17)
		resName = "";     // Residue Name (18-20)
		chainID = '\0';   // Chain ID (22)
		resNum = 0;       // Residue Number (23-26)
		iCode = '\0';     // Insertion Code (27)
		x = 0.0;          // X Coordinate (31-38)
		y = 0.0;          // Y Coordinate (39-46)
		z = 0.0;          // Z Coordinate (47-54)
		occupancy = 0.0;  // Occupancy (55-60)
		TempFactor = 0.0; // Temperature Factor (61-66)
		segID = "";       // Segment Identifier (73-76)
		element = "";     // Element symbol (77-78)
		charge = "";      // Atomic Charge (79-80)
	}

	bool PDBAtomLine::loadAtom(const Particle &atom, int atomIndex)
	{
		lineType = ATOM;
		atomNum = atomIndex;
		atomName = std::string(atom.pdbname);
		altLoc = ' ';
		resName = std::string(atom.parentl3name);
		chainID = ' ';
		resNum = atom.ir;
		iCode = ' ';
		x = atom.pos().x;
		y = atom.pos().y;
		z = atom.pos().z;
		occupancy =  1.0;
		TempFactor =  1.0;
		segID = "";
		element = "";
		charge = "";

		return true;
	}

	bool PDBAtomLine::loadLine(const std::string& _Line)
	{
		bool error = false;
		std::string line = _Line;
		line.resize(80,' ');

		// Reinitialise the structure.
		init();

		// Line Type check
		if( 0 == line.compare( 0,6,"ATOM  " ) ) 
		{
			lineType = ATOM;
		}
		else if( 0 == line.compare( 0,6,"HETATM" ) ) 
		{
			lineType = HETATM;
		}
		else if( 0 == line.compare( 0,6,"TER   " ) ) 
		{
			lineType = TER;
			// These calls may fail, it doesnt matter, only 'TER' is important for parsing...
			str2int( line.substr(6,5), atomNum );
			resName = line.substr(17,3);
			chainID = line[21];
			str2int( line.substr(22,4), resNum ); 
			iCode = line[26];
			return true;
		}
		else 
		{
			return false; // Not an atom record type!
		}

		if( 0 != str2int( line.substr(6,5), atomNum ) ) error = true;
		atomName = trim(line.substr(12,4));
		altLoc = line[16];
		resName = trim(line.substr(17,3));
		chainID = line[21];
		if( 0 != str2int( line.substr(22,4), resNum ) ) error = true;
		iCode = line[26];
		if( 0 != str2double( line.substr(30,8), x ) ) error = true;
		if( 0 != str2double( line.substr(38,8), y ) ) error = true;
		if( 0 != str2double( line.substr(46,8), z ) ) error = true;
		if( 0 != str2double( line.substr(54,6), occupancy ) ) error = true;
		if( 0 != str2double( line.substr(60,6), TempFactor ) ) error = true;
		segID = line.substr(72,4);
		element = line.substr(76,2);
		charge = line.substr(78,2);

		return true;
	}

	void PDBAtomLine::saveLine(std::ostream &_file )const
	{
		if( lineType == ATOM ) _file << "ATOM  ";
		else if( lineType == HETATM ) _file << "HETATM";
		else THROW(CodeException,"Unknown lineType is PDBAtomLine::saveLine()");

		int atomPrint = atomNum;
		if( atomPrint > 99999 )
		{
			std::cout << "WARNING: Atom number is too large for the PDB format. Rolling into range 0->99999 to allow output!" << std::endl;
			while( atomPrint > 99999 )
			{
				atomPrint -= 100000;
			}
		}
		_file << setw(5) << atomPrint << ' ';

		// Prepare the atom name
		std::string nameBuf = atomName;
		if( nameBuf.length() > 4 )
		{
			// truncate to 4 for the file format
			nameBuf = nameBuf.substr(0,4);
		}
		else if( nameBuf.length() <= 3)
		{
			// expand to 4 using formatted spaces
			if( !isASCIINumeric(nameBuf[0]) ) nameBuf.insert(0," ");
			for( size_t i = nameBuf.length(); i < 4; i++ )
			{
				nameBuf.push_back(' ');
			}
		}
		_file << nameBuf;
		_file << altLoc;

		// Prepare the residue name
		_file << makeStringOfLength(resName,3,' ',1); // 1: because of Nxxx and Cxxx
		_file << ' ';
		_file << chainID;

		int resPrint = resNum;
		if( resPrint > 9999 )
		{
			std::cout << "WARNING: Residue number is too large for the PDB format. Rolling into range 0->9999 to allow output!" << std::endl;
			while( resPrint > 9999 )
			{
				resPrint -= 10000;
			}
		}
		_file << setw(4) << resPrint;
		_file << iCode << "   ";

		setFloatFmt(_file,8,3);
		_file << x;
		setFloatFmt(_file,8,3);
		_file << y;
		setFloatFmt(_file,8,3);
		_file << z;
		setFloatFmt(_file,6,2);
		_file << occupancy;
		setFloatFmt(_file,6,2);
		_file << TempFactor;

		_file << makeStringOfLength(segID,4,' ');
		_file << makeStringOfLength(element,2,' ');
		_file << makeStringOfLength(charge,2,' ');

		_file << endl;
	}





	// -----------------------------------------------
	// ModResLine
	// -----------------------------------------------

	ModResLine::ModResLine() :
	chainID('\0'),
		resNum(-1),
		iCode('\0'),
		chemicalName("")
	{
	}

	bool ModResLine::parseLine( const std::string& line )
	{
		bool parseStatus = true;

		// make life easier and faster by caching the line in a mutable stringbuilder
		StringBuilder sb(line);
		if( sb.size() < 80 ) sb.PadRight(80,' ');

		if( sb.compare("MODRES",0,6,0,false) )
		{
			return false; // not a MODRES line!
		}

		// Assign our properties
		m_Alias = sb.toString(12, 3);
		chainID = sb[16];
		try 
		{ 
			resNum = sb.parseInt(18); 
		}
		catch( ParseException ex ) 
		{ 
			printf("Syntax error found in MODRES line. Index parsing failed.");
			parseStatus = false; 
		}
		iCode = sb[22];
		m_Name = sb.toString(24, 3);
		if ( 0 == m_Alias.compare("   ") )
		{
			printf("Syntax error found in MODRES line. Standard residue name is undefined.");
			parseStatus = false; // The bitch is blank - it tells us nothing useful, why, why is it there!!
		}
		chemicalName = sb.toString(29, 41);
		chemicalName = trim(chemicalName);

		return parseStatus;
	}

	bool ModResLine::matches(const FileParticle& particle) const
	{
		return ((particle.chainID == chainID) &&
			(particle.iCode == iCode) &&
			(particle.resNum == resNum));
	}

	bool ModResLine::matches(const std::string& _LookUp) const
	{
		return ( 0 == m_Alias.compare(_LookUp) );
	}

	// -----------------------------------------------
	// PDB_MODRES
	// -----------------------------------------------

	PDB_MODRES::PDB_MODRES() : Library::AliasMapper()
	{
	}

	bool PDB_MODRES::AttemptRename( FileParticle& particle ) const
	{
		if (m_Lines.size() == 0) 
			return false; // We have nothing to work with

		for( size_t i = 0; i < m_Lines.size(); i++ )
		{
			if (m_Lines[i].matches(particle))
			{                   
				char iCode = m_Lines[i].iCode;
				if (iCode == ' ') iCode = '_';
				std::cout << "Successful Name-Munge from MODRES data: Converted:'" << 
					m_Lines[i].chainID << "'s '" << particle.resName << ':' << m_Lines[i].resNum << ':' 
					<< iCode << "' to a '" << m_Lines[i].getAlias() << "'";

				// Assign the standard name to this residue.
				particle.resName = m_Lines[i].getName();

				return true; // We now have a valid particle in the eyes of the forcefield.
			}
		}

		return false; // We didn't a valid entry
	}

	void PDB_MODRES::parseLine( const std::string& line )
	{
		ModResLine modres;
		if (modres.parseLine(line))
		{
			m_Lines.push_back(modres);
		}
	}

	bool PDB_MODRES::lookupLongName( char _SourceID, std::string& _DestID ) const
	{
		return false; // MODRES doesn't provide single letter <-> multiletter mappings
	}

	bool PDB_MODRES::lookupShortName( const std::string& _SourceID, char& _DestID ) const
	{
		return false; // MODRES doesn't provide single letter <-> multiletter mappings
	}

	bool PDB_MODRES::lookupAlias( const std::string& _SourceID, std::string& _DestID ) const
	{
		if (m_Lines.size() == 0) 
			return false; // We have nothing to work with

		for( size_t i = 0; i < m_Lines.size(); i++ )
		{
			if (m_Lines[i].matches(_SourceID))
			{                   
				_DestID = m_Lines[i].getName(); // Assign the standard name
				return true;
			}
		}

		return false; // We didn't a valid entry
	}



	// -----------------------------------------------
	// PDB_Connections
	// -----------------------------------------------

	PDB_Connections::PDB_Connections() : FileCovalency()
	{
	}

	void PDB_Connections::parseCONECT( const std::string& line )
	{
		return; // We dont do these yet ...
	}

	void PDB_Connections::parseSSBOND( const std::string& line )
	{
		FileBond bond;
		StringBuilder sb(line);
		sb.PadRight(80,' ');
		if( !sb.compare("SSBOND",0,6,0) ) return; // Not a real SSBOND line
		if( !sb.compare("CYS",11,3,0) ) return; // How can you have an SSBOND not involving a CYS residue !?!
		if( !sb.compare("CYS",25,3,0) ) return;
		try
		{
			bond.r1ResNum = sb.parseInt(17);
			bond.r2ResNum = sb.parseInt(31);
		}
		catch(ParseException ex)
		{
			return;
		}
		bond.Type = FileBond::Disulphide;
		bond.r1ChainID = sb[15];
		bond.r1ICode = sb[21];
		bond.r1Name = "CYS";
		bond.r2ChainID = sb[29];
		bond.r2ICode = sb[35];
		bond.r2Name = "CYS";

		m_Lines.push_back(bond);
	}

	// -----------------------------------------------
	// PDB_Model
	// -----------------------------------------------

	PDB_Model::PDB_Model(const PDB_ImportBase& _PDBParent, const FileImportBase& _FileParent) : 
	FileMoleculeMaker(_FileParent),
		m_PDBParent(&_PDBParent)
	{
	}

	PDB_Model::~PDB_Model()
	{
	}

	const PDB_ImportBase& PDB_Model::getPDBParent() const
	{
		return *m_PDBParent;
	}

	bool PDB_Model::useSEQRES() const
	{
		return m_PDBParent->UseSEQRES;
	}

	File_Molecule PDB_Model::make( char chainID )
	{
		return make( chainID, NULL );
	}

	File_Molecule PDB_Model::make( char chainID, const Sequence::BioSequence* buildSequenceOverride )
	{
		if( 0 == pdbLines.size() ) THROW(ProcedureException,"PDB_Model contains no PDB atoms!");

		bool usingSEQRES = false;
		File_Molecule mol(chainID);

		if( buildSequenceOverride != NULL )
		{
			// If we have a sequence override, use it for the build process.
			if( buildSequenceOverride->size() == 0 )
			{
				THROW( ArgumentException, "Error creating 'File_Molecule' - a 'buildSequenceOverride' was provided, but appears to be zero-length!");
			}
			AssignBioSequence(mol,*buildSequenceOverride);
		}
		else if( useSEQRES() && isResequencableFilter() )
		{
			// Otherwise check to see if we have SEQRES data to define missing residues.
			const PDB_SEQRES& seqres = getPDBParent().getSEQRES();
			if( seqres.HasIndex(chainID) )
			{
				AssignBioSequence(mol,seqres.getSequence(chainID));
				usingSEQRES = true;
			}
		}		

		// Main chain extraction loop
		size_t added = 0; // keep track of how many we have added to our molecule 
		// - '0' returns an empty molecule now rather than later.
		for( size_t i = 0; i < pdbLines.size(); i++ )
		{
			if( (pdbLines[i].lineType != PDBAtomLine::TER) &&
				(pdbLines[i].chainID == chainID) )
			{				
				// We need to resolve the underlying filters, CLASS and ALIAS defs
				// This is done by appendFileParticle() which calls bool PassesFilter()
				if( appendFileParticle(mol,pdbLines[i]) )
				{
					added++; // Record this, we may not actually add any...
				}
			}
		}

		if( 0 == added )
		{
			// this is an empty mol, File_Molecule.isEmpty() will return true to the parent 'FileInBase'
			return mol;
		}

		// 1) Autogenerate the structural sequence from the FileParticle array
		// 2) perform a sequence alignment to BioSequence data
		// 3) Perform 'bond specific' renaming of particles e.g. CYS -> CYX
		finaliseSequence(mol); 

		if( usingSEQRES )
		{
			ASSERT(mol.isValid(),ProcedureException,"The SEQRES derived biological sequence does not correspond to available structural data - mismatches are present!!");
		}

		return mol;
	}	

	bool PDB_Model::HasChainID(char _RequestChain) const
	{
		for( size_t i = 0; i < pdbLines.size(); i++ )
		{	
			if( (pdbLines[i].chainID == _RequestChain )  
				&& wouldAppendFileParticle( pdbLines[i]) )
			{
				return true;
			}
		}
		return false;
	}

	std::vector<char> PDB_Model::findChainIDs() const
	{
		std::vector<char> chainIDs;
		for( size_t i = 0; i < pdbLines.size(); i++ )
		{
			if( wouldAppendFileParticle(pdbLines[i]) )
			{
				bool contains = false;
				for( size_t j = 0; j < chainIDs.size(); j++ )
				{
					if( pdbLines[i].chainID == chainIDs[j] )
					{
						contains = true;
						break;
					}
				}
				if( !contains ) chainIDs.push_back(pdbLines[i].chainID);
			}
		}
		return chainIDs;
	}









	// -----------------------------------------------
	// PDB_ImportBase
	// -----------------------------------------------

	const float NULL_RESLN = 999.99f;
	const char NULL_PDB_CHAR = '\0';

	PDB_ImportBase::PDB_ImportBase() 
		: m_SEQRES()
	{
		m_Resolution = NULL_RESLN;
		m_ExpMethod = Undefined;
		m_UnrecognisedLines = 0;
		MaxModelImport = 1; // This is the default number of models to import, it can be overridden
		HeaderFilterMode = DefaultHeader;
		m_Verbose = Verbosity::Normal;
		UseSEQRES = true;
	}

	void PDB_ImportBase::parseExpData(const std::string &line)
	{
		// COLUMNS       DATA TYPE      FIELD         DEFINITION
		// -------------------------------------------------------------------------------
		//  1 -  6       Record name    "EXPDTA"
		//  9 - 10       Continuation   continuation  Allows concatenation of multiple
		//                                            records.
		// 11 - 70       SList          technique     The experimental technique(s) with
		//                                            optional comment describing the
		//                                            sample or experiment.

		// We have to deal with lines that have ending Type definitions like this:
		// 'EXPDTA    NMR                                                           1HME   6'
		std::string inputLine;
		if( line.length() > 70 ) // remove the Padding past character 70.
		{
			inputLine = line.substr(0,70);
		}
		else
		{
			inputLine = line;
		}
		inputLine = trim(inputLine.substr(6,inputLine.length()-6)); // remove the 'EXPDATA' and surrounding whitespace
		makeupper(inputLine);

		// Now ... Extract what Type we have from these rather varied example types(!!!):
		//EXPDTA INFRARED SPECTROSCOPY
		//EXPDTA SINGLE-CRYSTAL ELECTRON DIFFRACTION
		//EXPDTA ELECTRON MICROSCOPY
		//EXPDTA FIBER DIFFRACTION
		//EXPDTA X-RAY DIFFRACTION - oh my god how inconsistent are these people!
		//EXPDTA X-RAY DIFFRACTION
		//EXPDTA NMR
		//EXPDTA NMR; THEORETICAL MODEL
		//EXPDTA NMR, 30 MODELS
		//EXPDTA NMR, REPRESENTATIVE STRUCTURE
		//EXPDTA NMR, RESTRAINED REGULARIZED MEAN STRUCTURE
		//EXPDTA NMR, REGULARIZED MEAN STRUCTURE
		//EXPDTA NMR, MINIMIZED STRUCTURE
		//EXPDTA NMR, MINIMIZED AVERAGE STRUCTURE
		//EXPDTA NMR, AVERAGE STRUCTURE
		// {
		//    EXPDTA NMR,  26 STRUCTURES (another variable indentation!)
		//    EXPDTA NMR, 25 STRUCTURE
		//    EXPDTA NMR, 20 STRUCTURES
		//    EXPDTA NMR, 15 STRUCTURES
		//    EXPDTA NMR, 10 STRUCTURES
		//    EXPDTA NMR, 5 STRUCTURES
		//    EXPDTA NMR, 1 STRUCTURE
		// }
		//EXPDTA SYNCHROTRON X-RAY DIFFRACTION

		m_ExpMethod = Undefined;

		if(0==inputLine.compare("X-RAY DIFFRACTION"))
		{
			m_ExpMethod = Crystalographic;
		}
		else if(0==inputLine.compare("SYNCHROTRON X-RAY DIFFRACTION"))
		{
			m_ExpMethod = Crystalographic;
		}
		else if(0==inputLine.compare("FIBER DIFFRACTION"))
		{
			m_ExpMethod = FiberDiffraction;
		}
		else if(0==inputLine.compare("ELECTRON MICROSCOPY"))
		{
			m_ExpMethod = ElectronMicroscopy;
		}
		else if(0==inputLine.compare("INFRARED SPECTROSCOPY"))
		{
			m_ExpMethod = InfraredMicroscopy;
		}
		else if(0==inputLine.compare("SINGLE-CRYSTAL ELECTRON DIFFRACTION"))
		{
			m_ExpMethod = SingleCrystalElectronMicroscopy;
		}
		else if (0==inputLine.compare(0,3,"NMR"))
		{
			// main Type
			m_ExpMethod = NMR;
			// NMR-subtype if you really want it...
			if( inputLine.length() >= 4 )
			{
				m_NMRMethod = ltrim(inputLine.substr(4,inputLine.length()-4));
			}
			else
			{
				m_NMRMethod = ""; // sometimes there isn't one...
			}
		}
		else
		{
			printf("An unknown EXPDTA Type was found in the PDB file, please writerm the application authors!");
			m_ExpMethod = UnknownMethod;
		}
	}

	void PDB_ImportBase::printExpdata()const
	{
		switch( m_ExpMethod )
		{
		case Crystalographic:
			printf( "Expermental Method: 'X-RAY DIFFRACTION'\n" );
			break;
		case ElectronMicroscopy:
			printf( "Expermental Method: 'ELECTRON MICROSCOPY'\n" );
			break;
		case FiberDiffraction:
			printf( "Expermental Method: 'FIBER DIFFRACTION'\n" );
			break;
		case InfraredMicroscopy:
			printf( "Expermental Method: 'INFRARED SPECTROSCOPY'\n" );
			break;
		case NMR:
			if( m_NMRMethod.length() > 0 )
			{
				printf( "Expermental Method: 'NMR, %s'\n", m_NMRMethod.c_str() );
			}
			else
			{
				printf( "Expermental Method: 'NMR'\n" );
			}
			break;
		case UnknownMethod:
			printf( "Expermental Method: 'UNKNOWN' (Original \"EXPDTA\" Tag Was Not Recognised)\n" );
			break;
		case Undefined:
			printf( "Expermental Method: 'UNDEFINED' (No \"EXPDTA\" Tag In File)\n" );
			break;
		default:
			THROW(CodeException,"Unknown internal m_ExpMethod!");
			break;
		}
	}

	void PDB_ImportBase::printExpdata( ofstream &_file )const
	{
		switch( m_ExpMethod )
		{
		case Crystalographic:
			_file << "EXPDTA X-RAY DIFFRACTION" << endl;
			break;
		case ElectronMicroscopy:
			_file << "EXPDTA ELECTRON MICROSCOPY" << endl;
			break;
		case FiberDiffraction:
			_file << "EXPDTA FIBER DIFFRACTION" << endl;
			break;
		case InfraredMicroscopy:
			_file << "EXPDTA INFRARED SPECTROSCOPY" << endl;
			break;
		case NMR:
			_file << "EXPDTA NMR: " << m_NMRMethod << endl;
			break;
		case UnknownMethod:
			_file << "EXPDTA UNKNOWN: 'EXPDTA' Tag In Original File Was Not Recognised" << endl;
			break;
		case Undefined:
			_file << "EXPDTA UNDEFINED: No 'EXPDTA' Tag In Original File\n" << endl;
			break;
		default:
			THROW(CodeException,"Unknown internal 'm_ExpMethod'!");
			break;
		}
	}

	void PDB_ImportBase::ScanFile( const FileImportBase& _DerivedImporter )
	{
		const std::string& _filename = _DerivedImporter.getFileName();
		printf("Perfoming Initial PDB File Scan of '%s' ...\n", _filename.c_str() );

		bool errorCondition = false;
		std::ifstream file(_filename.c_str(), ifstream::in);

		if( !file.is_open() ) throw(IOException( "PDB File not found: '" + _filename + "'!" ));

		std::string lineType(6,' '); // blank string, **MUST** be 6 in length
		std::string line;

		// member initialisation
		m_UnrecognisedLines = 0;

		// parsing helper variables
		bool atomImportStarted = false;
		bool expectingModelStatement = true;
		int modelCount = 0;
		int endModelCount = 0;

		PDBAtomLine pdbLine;
		PDB_Model modelMember(*this,_DerivedImporter);

		while(std::getline(file,line))
		{
			if( !getLineType(line,lineType) ) 
				continue; // line not long enough, ignore it!

			// ---------------------------
			// System Lines
			// ---------------------------
			if(0==lineType.compare("ATOM  "))
			{
				atomImportStarted = true;
				if( (int)m_Models.size() < MaxModelImport )
				{
					if( !pdbLine.loadLine(line) ) 
					{
						std::cout << "Failure parsing 'ATOM  ' line!" << std::endl;
						errorCondition = true;
					}
					else
					{
						modelMember.appendLine(pdbLine);
					}
				}
			}
			else if(0==lineType.compare("HETATM"))
			{
				atomImportStarted = true;
				if( (int)m_Models.size() < MaxModelImport )
				{
					if( !pdbLine.loadLine(line) ) 
					{
						std::cout << "Failure parsing 'HETATM' line!" << std::endl;
						errorCondition = true;
					}
					else
					{
						modelMember.appendLine(pdbLine);
					}
				}
			}
			else if(0==lineType.compare("TER   "))
			{
				atomImportStarted = true;
				if( (int)m_Models.size() < MaxModelImport )
				{
					if( !pdbLine.loadLine(line) ) 
					{
						std::cout << "Failure parsing 'TER   ' line!" << std::endl;
						errorCondition = true;
					}
					else
					{
						modelMember.appendLine(pdbLine);
					}
				}
			}
			else if(0==lineType.compare("MODEL "))
			{
				if( atomImportStarted && modelCount == 0 ) throw(ParseException("\nNew \"MODEL \" statements are not allowed following the import of the first atom of the PDB!\n"));
				if( !expectingModelStatement ) throw(ParseException("\nUnexpected \"MODEL \" statement!\n"));
				modelCount++;
				expectingModelStatement = false;
			}
			else if(0==lineType.compare("ENDMDL"))
			{
				if( expectingModelStatement ) throw(ParseException("\nUnexpected \"ENDMDL\" statement!\n"));
				expectingModelStatement = true;
				if( modelMember.getLineCount() > 0 )
				{
					m_Models.push_back( modelMember );
					modelMember = PDB_Model(*this,_DerivedImporter);
				}
				else if( (int)m_Models.size() < MaxModelImport )
				{
					printf("\nWARNING: 'MODEL'\\'ENDMDL' pair detected with no atom definitions in-between!\n");
				}
			}

			// -----------------------------
			//  Sequence and renaming lines
			// -----------------------------

			else if(0==lineType.compare("SEQRES") && (0 < (HeaderFilterMode & SEQRES)))
			{
				m_SEQRES.parseLine(line);
			}
			else if(0==lineType.compare("MODRES") && (0 < (HeaderFilterMode & MODRES)))
			{
				m_MODRES.parseLine(line);
			}

			// ---------------------------
			//  Bonding Lines
			// ---------------------------

			else if(0==lineType.compare("SSBOND") && (0 < (HeaderFilterMode & SSBOND)))
			{
				// just count them for the moment...
				m_Connections.parseSSBOND(line);
			}
			else if(0==lineType.compare("CONECT") && (0 < (HeaderFilterMode & CONECT)))
			{
				// just count them for the moment...
				m_Connections.parseCONECT(line);
			}

			// ----------------------------------------------------------
			//  'Pointless Lines' - Ignore them
			//  (often annoyingly found within the atom definition block)
			// ----------------------------------------------------------

			else if(0==lineType.compare("ANISOU")) { continue; }
			else if(0==lineType.compare("SIGUIJ")) { continue; }
			else if(0==lineType.compare("SIGATM")) { continue; }

			// ---------------------------
			//  Header Lines
			// ---------------------------

			else if(0==lineType.compare("HEADER") && (0 < (HeaderFilterMode & HEADER)) )
			{
				m_PDBHeader.push_back(line);
			}
			else if(0==lineType.compare("TITLE ") && (0 < (HeaderFilterMode & TITLE)) )
			{
				m_PDBHeader.push_back(line);
			}
			else if(0==lineType.compare("COMPND") && (0 < (HeaderFilterMode & COMPND)) )
			{
				m_PDBHeader.push_back(line);
			}
			else if(0==lineType.compare("SOURCE") && (0 < (HeaderFilterMode & SOURCE)) )
			{
				m_PDBHeader.push_back(line);
			}
			else if(0==lineType.compare("KEYWDS") && (0 < (HeaderFilterMode & KEYWDS)) )
			{
				m_PDBHeader.push_back(line);
			}
			else if(0==lineType.compare("EXPDTA") && (0 < (HeaderFilterMode & EXPDTA)) )
			{
				parseExpData(line);
				m_PDBHeader.push_back(line);
			}
			else if(0==lineType.compare("AUTHOR") && (0 < (HeaderFilterMode & AUTHOR)) )
			{
				m_PDBHeader.push_back(line);
			}
			else if(0==lineType.compare("REVDAT") && (0 < (HeaderFilterMode & REVDAT)) )
			{
				m_PDBHeader.push_back(line);
			}
			else if(0==lineType.compare("JRNL  ") && (0 < (HeaderFilterMode & JRNL)) )
			{
				m_PDBHeader.push_back(line);
			}
			else if(0==lineType.compare("REMARK") )
			{
				// we still want to obtain the resolution if available from "REMARK 2 RESOLUTION"
				if(0==line.compare(0,21,"REMARK   2 RESOLUTION"))
				{
					std::string resln = line.substr( 22, 5 );
					str2float( resln, m_Resolution );
				}
				else if( 0 < (HeaderFilterMode & REMARK) )
				{
					// only add the rest of the ~1,000,000 kinds of remarks if "REMARK" is flagged in the filter
					// (only happens in uber-import mode)
					m_PDBHeader.push_back(line);
				}
			}
			else if(0==lineType.compare("FORMUL") && (0 < (HeaderFilterMode & FORMUL)) )
			{
				m_PDBHeader.push_back(line);
			}
			else if(0==lineType.compare("HELIX ") && (0 < (HeaderFilterMode & HELIX)) )
			{
				m_PDBHeader.push_back(line);
			}
			else if(0==lineType.compare("SHEET ") && (0 < (HeaderFilterMode & SHEET)) )
			{
				m_PDBHeader.push_back(line);
			}
			else if(0==lineType.compare("HET   ") && (0 < (HeaderFilterMode & HET)) )
			{
				m_PDBHeader.push_back(line);
			}
			else if(0==lineType.compare("DBREF ") && (0 < (HeaderFilterMode & DBREF)) )
			{
				m_PDBHeader.push_back(line);
			}
			else
			{
				m_UnrecognisedLines++;
				if( 0 < (HeaderFilterMode & UNKNOWN_TAGS) ) m_PDBHeader.push_back(line);
			}
		}

		if( modelCount > 0 && !expectingModelStatement )
		{
			// we are in a MODEL/ENDMDL mode and there was no tailing ENDMDL statement
			throw(ParseException("Unexpected \"ENDMDL\" statement!\n"));
		}
		else
		{
			// we are in a normal atom list, import has ended, and now we must add the only 'model'
			if( modelMember.getLineCount() > 0 )
			{
				m_Models.push_back( modelMember );
			}
		}

		// -----------------------------------------
		// Data Verification and Summary Phase ...
		// -----------------------------------------

		printf(" Done!\n\n");

		// print the header if we like ...
		if( getVerbose() && m_PDBHeader.size() > 0 ) 
		{
			printf("Imported PDB_Header Lines: \n");
			printf("----------------------------\n");
			std::string headerLine;
			for( size_t i = 0; i < m_PDBHeader.size(); i++ )
			{
				headerLine = rtrim(m_PDBHeader[i]," \n\r\t");
				printf("%s\n",headerLine.c_str());
			}
			printf("----------------------------\n\n");
		}

		printf("PDB FileScan Summary:\n\t" );
		printExpdata();// print the extracted exprimental resolution method
		if( m_Resolution != NULL_RESLN )
		{
			printf("\tExperimental Resolution: %4.2f\n",m_Resolution);
		}
		else
		{
			printf("\tExperimental Resolution: Undefined\n");
		}

		// "MODEL  " and "ENDMDL" statements control '# Model Definitions'
		if( MaxModelImport >= 1 )
		{
			printf("\t%d Imported Model Definitions (Imposed Import Restriction of %d Models)\n",(int)m_Models.size(),MaxModelImport);
		}
		else
		{
			printf("\t%d Imported Model Definitions (No Imposed Import Max)\n",(int)m_Models.size());
		}
		printf("\t%d Validated Disulphide and Connect statements\n",m_Connections.getCount());
		printf("\t%d Modified Residue statements\n",m_MODRES.getCount());
		printf("\t%d Uninterpreted lines\n\n",m_UnrecognisedLines);

		file.close();
		return;
	}

	// -------------------------
	// PDB Babel 
	// -------------------------

	PDB_Babel::PDB_Babel()
	{
	}

	PDB_Babel::~PDB_Babel()
	{
	}

	void PDB_Babel::RenameFile( const std::string& _InFilename, const std::string& _OutFilename )
	{
		printf("Babeling PDB File '%s' to '%s' ...", _InFilename.c_str(), _OutFilename.c_str() );

		std::ifstream _in;
		std::ofstream _out;
		try
		{
			_in.open(_InFilename.c_str(), ifstream::in);
			_out.open(_OutFilename.c_str(), ios::out);
			if( !_in.is_open() )
			{
				throw(IOException("Could not open the file for reading!"));
			}
			if( !_out.is_open() )
			{
				throw(IOException("Could not open the file for writing!"));
			}

			std::string lineType("      "); // blank string, **MUST** be 6 in length
			std::string line;
			PDBAtomLine pdbLine;
			while(std::getline(_in,line))
			{
				if( !getLineType(line,lineType) ) 
				{
					_out << line;
					continue;
				}

				// ---------------------------
				// Atom Lines
				// ---------------------------

				if(0==lineType.compare("ATOM  "))
				{
					pdbLine.loadLine(line);
					ReinterpretName( pdbLine );
					pdbLine.saveLine(_out);
				}
				else if(0==lineType.compare("HETATM"))
				{
					pdbLine.loadLine(line);
					ReinterpretName( pdbLine );
					pdbLine.saveLine(_out);
				}
				else
				{
					_out << line;
				}
			}
		}
		catch(ExceptionBase ex)
		{
			if( _in.is_open() ) // release the file handles!
			{
				_in.close();
			}
			if( _out.is_open() ) // release the file handles!
			{
				_out.close();
			}
			throw; // rethrow our exception - we just wanted to catch it to close the fileStream...
		}

		if( _in.is_open() ) // release the file handles!
		{
			_in.close();
		}
		if( _out.is_open() ) // release the file handles!
		{
			_out.close();
		}
	}


	// -----------------------------------------------
	// PDB_Tools
	// -----------------------------------------------

	PDB_Tools::PDB_Tools( const std::string& _filename ) 
		: FileImportBase(_filename), 
		PDB_ImportBase()
	{
		addAliasMapper(getMODRES());
		m_SEQRES.setAlias(getAlias());	
		setCovalency(m_Connections);		
	}

	void PDB_Tools::ensureScan()
	{
		ScanFile(*this);
	}


	// -----------------------------------------------
	// PDB_In
	// -----------------------------------------------

	PDB_In::PDB_In(const FFParamSet &_ffps, const std::string& _filename )
		: FileInBase(_ffps,_filename), 
		PDB_ImportBase(),
		m_ScanPerformed(false)
	{		
		// addAliasMapper(_ffps) - No: _ffps is added by the FileInBase
		addAliasMapper(getMODRES()); // MODRES is added second - it is only used if there is no FF definition for the modified residue
		m_SEQRES.setAlias(getAlias());	
		setCovalency(m_Connections);
		setVerbosity(Verbosity::Normal);
	}

	void PDB_In::ensureScanned()
	{
		if( m_ScanPerformed ) return;
		ScanFile(*this);
		m_ScanPerformed = true;
	}

	void PDB_In::loadAll()
	{
		if( m_Models.size() == 0 ) ensureScanned(); // IMPORTANT
		for( size_t i = 0; i < m_Models.size(); i++ )
		{
			loadModel( (unsigned int)i );
		}
	}

	void PDB_In::loadModel( unsigned int modelIndex )
	{
		load( '\0', modelIndex );
	}

	void PDB_In::load()
	{
		load( '\0', 0 );
	}

	void PDB_In::load( char chainID )
	{
		load( chainID, 0 );
	}

	void PDB_In::load( char chainID, unsigned int modelIndex )
	{
		ensureScanned(); // IMPORTANT

		if( modelIndex >= m_Models.size() )
		{
			throw OutOfRangeException("Requested modelIndex is outside of bounds!");
		}

		if( chainID == '\0' )
		{
			std::vector<char> chains = m_Models[modelIndex].findChainIDs();
			if( chains.size() == 0 )
			{
				StringBuilder sb;
				sb.setFormat("Load() failed! No chains found in the context of the current filter: '%s'!")(getFilterText());
				throw ProcedureException(sb.toString());
			}
			else
			{
				for( size_t i = 0; i < chains.size(); i++ )
				{
					char printChain = chains[i];
					if( printChain == ' ' ) printChain = '_';
					printf("Loading from chain '%c'\n",printChain);
					baseLoad(m_Models[modelIndex],chains[i]);
				}
			}
		}
		else
		{
			if( !m_Models[modelIndex].HasChainID(chainID) )
			{
				StringBuilder sb;
				sb.setFormat("Load() failed! ChainID '%c' not found in the context of the current filter: '%s'!")(chainID)(getFilterText());
				throw ProcedureException(sb.toString());
			}
			char printChain = chainID;
			if( printChain == ' ' ) printChain = '_';
			printf("Loading from chain '%c'\n",printChain);
			baseLoad(m_Models[modelIndex],chainID);
		}
	}

	void PDB_In::loadExplicit( char chainID, Library::ResidueClass _class, unsigned int modelIndex )
	{
		loadExplicitInner( chainID, _class, NULL, modelIndex );
	}

	void PDB_In::loadExplicit( char chainID, Library::ResidueClass _class, const Sequence::BioSequence& _SequenceOverride, unsigned int modelIndex )
	{
		loadExplicitInner( chainID, _class, &_SequenceOverride, modelIndex );
	}

	void PDB_In::loadExplicitInner( char chainID, Library::ResidueClass _class, const Sequence::BioSequence* _SequenceOverride, unsigned int modelIndex )
	{
		ensureScanned(); // IMPORTANT

		// For the inner baseLoad call to work, we actually have to set the internal filter flags!
		bool storeFlag = m_FilterFlag;
		Library::ExtdResidueClasses storeClass = m_FilterE;
		m_FilterFlag = false; // we are in enum class mode!
		m_FilterE = _class;

		if( modelIndex >= m_Models.size() )
		{
			throw OutOfRangeException("Requested modelIndex is outside of bounds!");
		}
		if( !m_Models[modelIndex].HasChainID(chainID) )
		{
			throw ArgumentException("ChainID not found!");
		}

		printf("Loading explicit entity model:%d chain:%c using external sequence...\n",modelIndex,chainID);
		baseLoad(m_Models[modelIndex],chainID,_SequenceOverride);

		m_FilterFlag = storeFlag; // we are in enum class mode!
		m_FilterE = storeClass;
	}

	void PDB_In::loadExplicit( char chainID, const std::string& resName, int residueNum, char iCode, unsigned int modelIndex )
	{
		THROW(NotImplementedException,"THis function is not yet implemented !!");
		//baseLoad(m_Models[modelIndex],chainID,resName,residueNum,iCode);
	}











	// -----------------------------------------------
	// PDB_RawWriter
	// -----------------------------------------------

	PDB_RawWriter::PDB_RawWriter() 
		: m_ModelIndex(0),
		UseBioInfo(true)
	{
	}

	void PDB_RawWriter::rawRemarks() const
	{
		assertStream();
		tm *ltime;
		time_t stime;
		stime = (int) time(NULL);
		ltime = localtime(&stime);

		stream() << "REMARK File written with " << Revision::FullStamp() << endl;
		stream() << "REMARK (c) J.Rea & M.Tyka 2003-2006" << endl;
		stream() << "REMARK Compiled: " << __DATE__ << __TIME__ << endl;
		stream() << "REMARK File Written: " << asctime(ltime);
	}

	void PDB_RawWriter::rawBeginModel() const
	{		
		// MODEL Record Format:
		//
		// COLUMNS     DATA TYPE       FIELD       DEFINITION
		// -------------------------------------------------------------
		//  1 - 6      Record name     "MODEL "
		// 11 - 14     Integer         serial      Model serial number.

		assertStream();
		stream() << "MODEL " << setw(8) << right << m_ModelIndex++ << endl;
	}

	void PDB_RawWriter::rawEnd() const
	{		
		assertStream();
		stream() << "END" << endl;
	}

	void PDB_RawWriter::rawEndModel() const
	{		
		assertStream();
		stream() << "ENDMDL" << endl;
	}

	void PDB_RawWriter::rawWriteFullHeader( const MoleculeBase& _Mol ) const
	{	
		assertStream();
		rawRemarks();
		// rawSEQRES(_Mol);
		// rawSSBOND(_Mol);
	}

	void PDB_RawWriter::rawWriteFullHeader( const System& _Sys ) const
	{	
		assertStream();
		rawRemarks();
		// rawSEQRES(_Mol);
		// rawSSBOND(_Mol);
	}

	int PDB_RawWriter::rawWrite( const MoleculeBase& _Mol, int _AtomOffset ) const
	{		
		assertStream();

		// We assume >0 atoms in code below...
		if( _Mol.nAtoms() == 0 ) return 0;

		const Sequence::BioSequence& seq = _Mol.getSequence();
		bool useSeqIndexing = UseBioInfo && (seq.size() == _Mol.nResidues());

		char prevChainID = _Mol.getChainID(0);
		PDBAtomLine line;
		for( int i = 0; i < _Mol.atom.size(); i++ )
		{
			line.loadAtom(_Mol.atom[i],i);
			line.chainID = _Mol.getChainID(i);
			line.atomNum += _AtomOffset;

			// Position must come from atomxyz() virtual function call. This cannot be done in loadAtom()
			Maths::dvector pos( _Mol.atomxyz(i) );
			line.x = pos.x;
			line.y = pos.y;
			line.z = pos.z;

			if(useSeqIndexing)
			{
				// Override the residue number and insertion code using imported biological information
				// found when this data was originally loaded, if present. Important for true PDB output.
				// Behaviour can be modified using the public flag 'useSeqIndexing'.
				const Sequence::ResidueInfo& res = seq.getResidue(line.resNum);
				line.iCode = res.getBiologicalICode();
				line.resNum = res.getBiologicalIndex();				
			}

			// This is required when outputting from a 'Workspace' as they contain multiple molecules from many chains
			// When outputting 'Molecule',this will never evaluate to true.
			if( line.chainID != prevChainID )
			{
				prevChainID = line.chainID;
				stream() << "TER" << endl;
			}

			line.saveLine(stream());
		}

		// I have added this 'if' to ensure that H20 collections where each H20 is a separate molecule
		// are  printed without the TER. This may or may not be PROPER ... who knows ... but large
		// amounts of TER statements looks silly ...
		if( _Mol.res.size() > 1 )
		{
			stream() << "TER" << endl;
		}

		return (int) _Mol.atom.size();
	}

	int PDB_RawWriter::rawWrite( const System& _Sys, int _AtomOffset ) const
	{		
		assertStream();
		for( size_t i = 0; i < _Sys.nMolecules(); i++ )
		{
			_AtomOffset += rawWrite( _Sys.getMolecule(i), _AtomOffset );
		}
		return _AtomOffset;
	}

	int PDB_RawWriter::rawWrite( const std::vector<PDBAtomLine> &_lines, char chainID, int _AtomOffset ) const
	{
		assertStream();
		PDBAtomLine line;
		for( size_t i = 0; i < _lines.size(); i++ )
		{
			line = _lines[i];
			line.chainID = chainID;
			line.atomNum += _AtomOffset;
			line.saveLine(stream());
		}
		return (int) _lines.size();
	}

	int PDB_RawWriter::rawWrite( const std::vector<Maths::dvector> &_Coordinates, char chainID, int _AtomOffset ) const
	{
		assertStream();
		PDBAtomLine line;
		line.atomName = "POS";
		line.element = 'X';
		line.resName = "MOL";
		line.resNum = 0;
		line.chainID = chainID;
		for( size_t i = 0; i < _Coordinates.size(); i++ )
		{
			line.atomNum = (int)i + _AtomOffset;
			line.x = _Coordinates[i].x;
			line.y = _Coordinates[i].y;
			line.z = _Coordinates[i].z;
			line.saveLine(stream());
		}
		return (int)_Coordinates.size();
	}

	void PDB_RawWriter::rawWrite( 
		const Maths::dvector &_Coordinate, 
		const std::string& atomName, 
		const std::string& resName, 
		int& atomnum, 
		int& resnum 
		) const
	{
		assertStream();
		static PDBAtomLine line;
		line.atomName = atomName;
		if( atomName.size() >= 2 )
		{
			line.element = atomName[1];
		}
		else
		{
			line.element = 'X';
		}
		line.resName = resName;
		line.resNum = resnum++;
		line.atomNum = atomnum++;
		line.x = _Coordinate.x;
		line.y = _Coordinate.y;
		line.z = _Coordinate.z;
		line.saveLine(stream());
	}

	void PDB_RawWriter::rawWrite( 
		const std::vector<Maths::dvector> &_Coordinates, 
		int& atomnum, 
		int& resnum 
		)const
	{
		assertStream();
		PDBAtomLine line;
		line.atomName = "POS";
		line.element = 'X';
		line.resName = "MOL";
		line.resNum = resnum++;
		for( size_t i = 0; i < _Coordinates.size(); i++ )
		{
			line.atomNum = atomnum++;
			line.x = _Coordinates[i].x;
			line.y = _Coordinates[i].y;
			line.z = _Coordinates[i].z;
			line.saveLine(stream());
		}
	}


	// ---------------
	//  PDB_RawWriter
	// ---------------

	PDB_Writer::PDB_Writer( bool _MultiModel ) 
		: PDB_RawWriter(), 
		m_MultiModel(_MultiModel),
		m_DoneHeader(false)
	{
	}

	PDB_Writer::PDB_Writer( const std::string& _filename, bool _MultiModel ) 
		: PDB_RawWriter(),
		m_MultiModel(_MultiModel),
		m_DoneHeader(false)
	{
		// Initialise the base stream. We are appending if we are multi-model else just write once.
		// Flag clear to remove existing file contents
		PDB_RawWriter::streamTo(_filename,_MultiModel,true);
	}

	PDB_Writer::PDB_Writer( const char* _filename, bool _MultiModel ) 
		: PDB_RawWriter(),
		m_MultiModel(_MultiModel),
		m_DoneHeader(false)
	{
		// Initialise the base stream. We are appending if we are multi-model else just write once.
		// Flag clear to remove existing file contents
		PDB_RawWriter::streamTo(std::string(_filename),_MultiModel,true);
	}

	PDB_Writer::PDB_Writer( std::ostream& _stream, bool _MultiModel ) 
		: PDB_RawWriter(), 
		m_MultiModel(_MultiModel),
		m_DoneHeader(false)
	{
		PDB_RawWriter::streamTo(_stream);
	}

	void PDB_Writer::write( const MoleculeBase& _Mol )
	{
		try
		{
			beginStream();
			if( !m_DoneHeader ) 
			{
				rawWriteFullHeader( _Mol );
				m_DoneHeader = true;
			}
			if( m_MultiModel ) rawBeginModel();
			rawWrite(_Mol);
			if( m_MultiModel ) rawEndModel();
			rawEnd();
			stream().flush();
			endStream();
		}
		catch(ExceptionBase ex)
		{
			endStream(); // release the file handle!
			throw; // rethrow our exception - we just wanted to catch it to close the fileStream...
		}		
	}

	void PDB_Writer::write( const System& _Sys )
	{
		try
		{
			beginStream();
			if( !m_DoneHeader )
			{
				rawWriteFullHeader( _Sys );
				m_DoneHeader = true;
			}
			if( m_MultiModel ) rawBeginModel();
			rawWrite(_Sys);
			if( m_MultiModel ) rawEndModel();
			rawEnd();
			stream().flush();
			endStream();
		}
		catch(ExceptionBase ex)
		{
			endStream(); // release the file handle!
			throw; // rethrow our exception - we just wanted to catch it to close the fileStream...
		}		
	}

	// ---------
	//  OutTra_PDB
	// ---------

	OutTra_PDB::OutTra_PDB( const std::string &_filestem, WorkSpace &_wspace, bool _CallCreate )
		: OutputTrajectoryFile(_filestem, _wspace ),
		PDB_RawWriter()
	{
		// Initialise the base stream. We are appending if we are multi-model else just write once.
		// Flag clear to remove existing file contents
		PDB_RawWriter::streamTo(_filestem + ".pdb",true,true);

		if( _CallCreate ) 
		{
			create();
		}
	}

	int OutTra_PDB::create()
	{
		if( !created )
		{
			beginStream();
			rawWriteFullHeader( getWSpace() );
			endStream();
			created = true;
		}
		return 0;
	}

	int OutTra_PDB::append()
	{
		if( !created )
		{
			int result = create();
			if( result != 0 ) return result;
		}
		beginStream();
		rawBeginModel();
		rawWrite( getWSpace() );
		rawEndModel();
		endStream();
		return 0;
	}

	// -------------------
	//  OutTra_PDB_Single
	// -------------------

	int OutTra_PDB_Single::append()
	{
		if( !created )
		{
			int result = create();
			if( result != 0 ) return result;
			beginStream();
			rawWrite( getWSpace() );
			endStream();
		}
		return 0;
	}



	// -------------------
	//  OutputFile_PDB 
	// -------------------



void OutputFile_PDB::save( WorkSpace &_wspace){
	OutTra_PDB_Single outtra( filestem, _wspace );
	outtra.append( );
}

void OutputFile_PDB::save( System &_system ){
	PDB_Writer out( filestem + ".pdb" , false );
	out.write( _system );
}


}


