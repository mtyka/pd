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

// Own Header
#include "tra.h"

// MMLib Includes
#include "sequence/sequence.h"
#include "fileio/pdb.h"
#include "workspace/workspace.h"

// Namespace Includes
using namespace std;
using namespace Maths;

namespace IO 
{
	OutTra_BTF::OutTra_BTF(const std::string &_filestem, WorkSpace& _wspace)
		: OutputTrajectoryFile(_filestem, _wspace)
	{
		filename = _filestem + std::string(".tra");
	}

	int OutTra_BTF::create()
	{
		return create(DefaultIncludes);
	}

	int OutTra_BTF::create(BFT_StandardIncludes _includes)
	{
		if( created )
		{
			THROW(ProcedureException,"OutTra_BTF already created!");
		}

		// REMEMBER: for any changes to this code:
		// check if the calculation of blocksize is still correct!
		// (see the class 'BTF_Header')...
		ofstream traFile;
		try
		{
			traFile.open(filename.c_str(), ios::out | ios::binary);

			if( !traFile.is_open() )
			{
				THROW(IOException,"'OutTra_BTF' could not open the file for writing!");
			}

			// calculate additional block and header size due to the 'additional Blocks'
			int additionalHeaderSize = 0;
			for( size_t i = 0; i < m_Blocks.size(); i++ )
			{
				additionalHeaderSize += m_Blocks[i].getHeaderSize();
			}
			int additionalBlockSize = 0;
			for( size_t i = 0; i < m_Blocks.size(); i++ )
			{
				additionalBlockSize += m_Blocks[i].getBlockSize();
			}

			// write Header
			m_Header.save(traFile,&getWSpace(),_includes,additionalHeaderSize,additionalBlockSize);

			// write System Definitions
			traFile.write("SYSTDEFI", 8);
			BTF_SystemDefinitionEntry tse;
			for(int i = 0; i < getWSpace().atom.size(); i++)
			{
				tse.setFrom(i,&getWSpace().atom[i]);
				tse.save(traFile);
			}

			for( size_t i = 0; i < m_Blocks.size(); i++ )
			{
				m_Blocks[i].appendHeader( traFile );
			}

			traFile.write("TRASTART", 8);

		}
		catch(ExceptionBase ex)
		{
			if( traFile.is_open() ) traFile.close(); // release the file handle!
			throw; // rethrow our exception - we just wanted to catch it to close the fileStream...
		}

		if( traFile.is_open() ) traFile.close(); // release the file handle!
		created = true; // flag this member in the base class PD_API to signal all is well
		return 0;
	}

	int OutTra_BTF::append()
	{
		if( !created )
		{
			create();
		}

		// Look for problems
		if( m_Header.atoms != getWSpace().atom.size() )
		{
			THROW(ProcedureException,"OutTra_BTF::append() - wspace and TraHeader atom-counts do not match!");
		}
		if( m_Header.residues != getWSpace().res.size() )
		{
			THROW(ProcedureException,"OutTra_BTF::append() - wspace and TraHeader residue-counts do not match!");
		}

		// Initiate file writing procedure...
		ofstream traFile;
		try
		{
			const int maxAttempt = 1000;
			int count = 0;
			do
			{
				traFile.open(filename.c_str(), ios::app | ios::binary);
				if( count++ == maxAttempt )
				{
					THROW(IOException,"'OutTra_BTF' could not open the file for writing!");
				}
			}
			while( !traFile.is_open() );

			traFile.write("TRAE",4);

			float writeFloat;
			for(int i = 0; i < m_Header.atoms; i++)
			{
				writeFloat = (float)getWSpace().cur.atom[i].p.x;
				traFile.write((char*)&writeFloat,sizeof(float));
				writeFloat = (float)getWSpace().cur.atom[i].p.y;
				traFile.write((char*)&writeFloat,sizeof(float));
				writeFloat = (float)getWSpace().cur.atom[i].p.z;
				traFile.write((char*)&writeFloat,sizeof(float));
			}
			if(0 != (m_Header.Type & PhiPsis))
			{
				double phi, psi;
				for( int i = 0; i < m_Header.residues; i++ )
				{
					getWSpace().calcResiduePhiPsi(i,phi,psi);
					writeFloat = (float)phi;
					traFile.write((char*)&writeFloat,sizeof(float));
					writeFloat = (float)psi;
					traFile.write((char*)&writeFloat,sizeof(float));
				}
			}
			if(0 != (m_Header.Type & Energies))
			{
				m_Ene.save(traFile,&getWSpace());
			}
			if(0 != (m_Header.Type & ForceVectors))
			{
				for(int i = 0; i < m_Header.atoms; i++)
				{
					// 1E9 is a reasonable fudge factor to get a vector magnitude suitable for rendering
					writeFloat = (float)(getWSpace().cur.atom[i].p.x + (getWSpace().cur.atom[i].f.x*5E8));
					traFile.write((char*)&writeFloat,sizeof(float));
					writeFloat = (float)(getWSpace().cur.atom[i].p.y + (getWSpace().cur.atom[i].f.y*5E8));
					traFile.write((char*)&writeFloat,sizeof(float));
					writeFloat = (float)(getWSpace().cur.atom[i].p.z + (getWSpace().cur.atom[i].f.z*5E8));
					traFile.write((char*)&writeFloat,sizeof(float));
				}
			}

			// Now also add all the extension blocks to the file
			for( size_t i = 0; i < m_Blocks.size(); i++ )
			{
				m_Blocks[i].appendData( traFile );
			}
		}
		catch(ExceptionBase ex)
		{
			if( traFile.is_open() ) traFile.close(); // release the file handle!
			throw; // rethrow our exception - we just wanted to catch it to close the fileStream...
		}

		if( traFile.is_open() ) traFile.close(); // release the file handle!
		return 0;
	}

	void OutTra_BTF::addBlock( BTF_Block& _Block )
	{
		if( created ) THROW(ProcedureException,"OutTra_BTF::addBlock() cannot be called after the TrajectoryFormat file is created on disk.");
		m_Blocks.add( _Block );		
	}

	void OutTra_BTF::addOwnedBlock( BTF_Block* _Block )
	{
		if( created ) THROW(ProcedureException,"OutTra_BTF::addBlock() cannot be called after the TrajectoryFormat file is created on disk.");
		m_Blocks.addWithOwnership( _Block );		
	}







	BTF_ImportBase::BTF_ImportBase(const std::string &_fileName, bool _FullFileValidation)
	{
		reOpen(_fileName,_FullFileValidation);
	}

	void BTF_ImportBase::reOpen( const std::string &_fileName, bool _FullFileValidation)
	{
		m_Filename = _fileName;
		load();
		if( _FullFileValidation )
		{
			validateFile(); // Performs recountEntries();
		}
		else // just...
		{
			recountEntries();
		}
	}

	void BTF_ImportBase::assertTag( std::ifstream &_stream, std::string _tag )
	{
		size_t length = _tag.length();
		if( length == 0 )
		{
			THROW(CodeException,"assertTag() requires a tag of '>0' in length!");
		}
		char* check = new char[length+1];
		check[length] = 0;
		_stream.read(check, (std::streamsize) length);
		if( _stream.gcount() != length )
		{
			THROW(ParseException,"assertTag() could not read the entire tag!!");
		}
		if( 0 != _tag.compare(check))
		{
			THROW(ParseException,"assertTag() Tag verification failed!!");
		}
		delete[] check;
	}

	void BTF_ImportBase::load()
	{
		ifstream traFile;
		try
		{
			traFile.open(m_Filename.c_str(), ios::in | ios::binary);

			if( !traFile.is_open() )
			{
				THROW(IOException,"'BTF_ImportBase' could not open the file '" + m_Filename + "' for reading!");
			}

			m_Header.load(traFile);

			assertTag(traFile,"SYSTDEFI"); // Internal Sanity check printed in the TrajectoryFormat File

			sysDefs.clear(); // remove anything in the current array.
			BTF_SystemDefinitionEntry ent;
			for( int i = 0; i < m_Header.atoms; i++ )
			{
				ent.load(traFile);
				sysDefs.push_back(ent);
			}

			// ---------------------------------------------------------
			// If needed in the future we can parse trablocks here ...
			// Otherwise we skip over the whole lot below ...
			// ---------------------------------------------------------

			// check that out file seems to be valid for import.
			traFile.seekg(m_Header.trajectorystart);
			assertTag(traFile,"TRASTART"); // Internal Sanity check printed in the TrajectoryFormat File
		}
		catch( ExceptionBase )
		{
			if( traFile.is_open() ) traFile.close();
			throw; // Close the file handle and re-throw the exception.
		}

		// Close the file handle...
		if( traFile.is_open() ) traFile.close();
	}

	void BTF_ImportBase::info(bool _verbose)
	{
		m_Header.info(_verbose);
		for(size_t i = 0; i < sysDefs.size(); i++)
		{
			sysDefs[i].info(_verbose);
		}
	}

	void BTF_ImportBase::validateFile()
	{
		recountEntries();
		ifstream traFile;
		try
		{
			traFile.open(m_Filename.c_str(), ios::in | ios::binary);
			if( !traFile.is_open() )
			{
				THROW(IOException,"'BTF_ImportBase' could not open the file for validation!");
			}

			traFile.seekg(sizeof(BTF_Header));
			assertTag(traFile,"SYSTDEFI"); // Internal Sanity check printed in the TrajectoryFormat File
			traFile.seekg(m_Header.trajectorystart);
			assertTag(traFile,"TRASTART"); // Internal Sanity check printed in the TrajectoryFormat File
			long seekPos = (long)m_Header.trajectorystart + 8;
			for(int i = 0; i < m_Entries; i++)
			{
				assertTag(traFile,"TRAE");
				seekPos += (long)m_Header.blocksize;
				traFile.seekg(seekPos);
			}
		}
		catch( ExceptionBase )
		{
			if( traFile.is_open() ) traFile.close();
			throw; // Close the file handle and re-throw the exception.
		}

		// Close the file handle...
		if( traFile.is_open() ) traFile.close();
	}

	int BTF_ImportBase::recountEntries()
	{
		long size = IO::getFileSize(m_Filename);
		size = (size - (long)m_Header.trajectorystart - 8); // 8 bytes for "TRASTART"
		if( (size % (long)m_Header.blocksize) != 0 ) THROW(ParseException,"TrajectoryFormat file_size seems to be corrupted - a 'remainder of bytes' is present.");
		size = size / (long)m_Header.blocksize;
		m_Entries = (int) size;
		return m_Entries;
	}

	int BTF_ImportBase::size()const
	{
		return m_Header.residues;
	}

	int BTF_ImportBase::getEntryCount()const
	{
		return m_Entries;
	}

	const Sequence::BioSequence& BTF_ImportBase::getSequence() const
	{
		THROW(NotImplementedException,"");
	}

	int BTF_ImportBase::passesFilter( int _Index, AtomFilter _filter, int _molNum )
	{
		// check mol range
		if( _molNum > 0 && _molNum != sysDefs[_Index].parentnumber ) return 0;
		// check atom filter
		if( strcmp(&sysDefs[_Index].pdbname[0], "CA") == 0 )
		{
			return ( 0 != (_filter & CA) );
		}
		else if( strcmp(&sysDefs[_Index].pdbname[0], "N") == 0 )
		{
			return ( 0 != (_filter & N) );
		}
		else if( strcmp(&sysDefs[_Index].pdbname[0], "O") == 0 )
		{
			return ( 0 != (_filter & O) );
		}
		else if( strcmp(&sysDefs[_Index].pdbname[0], "C") == 0 )
		{
			return ( 0 != (_filter & C) );
		}
		else if( (sysDefs[_Index].pdbname[0] == 'H') ||
			(isASCIINumeric(sysDefs[_Index].pdbname[0]) && sysDefs[_Index].pdbname[1] == 'H' ) )
		{
			return ( 0 != (_filter & Hydrogens) );
		}
		else
		{
			return ( 0 != (_filter & Sidechain) );
		}
	}

	void BTF_ImportBase::seekToEntry( std::ifstream &_stream, int _entry )
	{
		long seekPos = (long)m_Header.trajectorystart + 8 + (m_Header.blocksize * _entry);
		_stream.seekg(seekPos);
		assertTag(_stream,"TRAE");
	}

	void BTF_ImportBase::loadVectors( std::ifstream &_traFile, std::vector<Maths::dvector>& _storage, AtomFilter _filter, int _molNum )
	{
		//Storage
		Tvector<float> getPos;

		// Do we have to filter?
		if( _molNum < 0 && _filter == All )
		{
			if( _storage.size() != m_Header.atoms )
			{
				_storage.resize(m_Header.atoms);
			}
			for( int i = 0; i < m_Header.atoms; i++ )
			{
				_traFile.read((char*)&getPos,sizeof(float)*3);
				_storage[i].setTo(getPos);
			}
		}
		else
		{
			int atomCount = 0;
			for( int i = 0; i < m_Header.atoms; i++ )
			{
				atomCount += passesFilter(i,_filter,_molNum);
			}
			if( _storage.size() != atomCount )
			{
				_storage.resize(atomCount);
			}
			int allocated = 0;
			for( int i = 0; i < m_Header.atoms; i++ )
			{
				if( passesFilter(i,_filter,_molNum) > 0 )
				{
					_traFile.read((char*)&getPos,sizeof(float)*3);
					_storage[allocated++].setTo(getPos);
				}
				else
				{
					_traFile.seekg(sizeof(float)*3,ios::cur);
				}
			}
		}
	}




	InTra_BTF::InTra_BTF( const std::string &_fileName, bool _ValidateFile )
		: BTF_ImportBase(_fileName,_ValidateFile), InputTrajectory_RandomAccess( _fileName )
	{
		reset();
	}

	void InTra_BTF::loadIntoSystem(System &sysspec, int position)
	{
		if( position < -1 || position >= getEntryCount() )
		{
			THROW(OutOfRangeException,"");
		}
		if( position == -1 ) // -1 encodes "last"
		{
			position = getEntryCount() - 1;
		}

		// Now lets find the atomic positions
		ifstream traFile;
		try
		{
			traFile.open(filename.c_str(), ios::in | ios::binary);

			if( !traFile.is_open() )
			{
				THROW(IOException,"'BTF_ImportBase' could not open the file '" + filename + "' for reading!");
			}

			seekToEntry(traFile,position); // make sure we import the correct frame!
			Molecule newmolecule(sysspec.ffps()); // make ourselves a molecule to contain the imported system
			Tvector<float> getpos; // A Temporary holder for float, required as the Particle uses doubles

			Sequence::BioSequence seq( sysspec.ffps() );
			int prevIR = INT_MAX;
			StringBuilder seqParse;

			for( int i = 0; i < m_Header.atoms; i++ )
			{
				std::string pdbName = std::string(&sysDefs[i].pdbname[0]);
				std::string parentName = std::string(&sysDefs[i].parentname[0]);

				int iMol = sysspec.ffps().findMoleculeType( parentName );
				ASSERT(iMol!=-1,ProcedureException,"");
				int iAt = sysspec.ffps().molecule[iMol].findAtomPDB(pdbName);
				ASSERT(iAt!=-1,ProcedureException,"");
				const AtomParameter& atParam = sysspec.ffps().molecule[iMol].atom[iAt];
				
				Particle newparticle(atParam);

				newparticle.parentname = parentName;

				if( parentName.size() > 3 )
				{
					newparticle.parentl3name = parentName.substr( parentName.size() - 3, 3 );
				}
				else
				{
					newparticle.parentl3name = parentName;
				}

				if( prevIR != sysDefs[i].parentnumber )
				{
					prevIR = sysDefs[i].parentnumber;
					seqParse.append( parentName );
					seqParse.append('-');
				}

				newparticle.parentl3name = parentName;
				// newparticle.parentletter = cant do

				// set the target position from the system header
				newparticle.posRef().setTo(sysDefs[i].targetx, sysDefs[i].targety, sysDefs[i].targetz);

				if(sysDefs[i].structureknown != 0)
				{
					newparticle.setKnownStructure(true);
					newparticle.setRebuildRequired(false);
				}
				else
				{
					newparticle.setKnownStructure(false);
					newparticle.setRebuildRequired(true);
				}

				// residue number
				newparticle.ir = sysDefs[i].parentnumber;

				// read from our tra file
				traFile.read((char*)&getpos,sizeof(float)*3);
				newparticle.pos().setTo( getpos );

				newmolecule.addParticle(newparticle);
			}

			newmolecule.beginBonding();
			for( int i = 0; i < m_Header.atoms; i++ )
			{				
				for( int k = 0; k < sysDefs[i].n_cov12atoms; k++)
				{
					newmolecule.addBond( i, sysDefs[i].cov12atom[k] );
				}
			}
			newmolecule.endBonding();

			if( seqParse[seqParse.size()-1] == '-' )
			{
				seqParse.TruncateRightBy(1);
			}
			seq.setTo( seqParse.toString() );

			newmolecule.setBioSequence( seq );

			// and finally add the molecule to the parent container...
			sysspec.add(newmolecule);
		}
		catch( ExceptionBase )
		{
			if( traFile.is_open() ) traFile.close();
			throw; // Close the file handle and re-throw the exception.
		}

		// Close the file handle
		if( traFile.is_open() ) traFile.close();
	}

	bool InTra_BTF::readNext( SnapShot &ss )
	{
		ss = makeSnapShot( currentPos++ );
		return !isEndOfFile();
	}

	bool InTra_BTF::skip()
	{
		currentPos++;
		return !isEndOfFile();
	}

	bool InTra_BTF::isEndOfFile() const
	{
		return currentPos >= nEntries();
	}

	void InTra_BTF::reset()
	{
		currentPos = 0;
	}

	void InTra_BTF::readRandomAccess( SnapShot &ss, size_t entry )
	{
		ss = makeSnapShot( entry );
	}

	size_t InTra_BTF::nEntries() const
	{
		return getEntryCount();
	}

	SnapShot InTra_BTF::makeSnapShot( size_t position )
	{
		ASSERT( position < getEntryCount(), ArgumentException, "Entry request is outside of tra range");
		
		SnapShot ss( m_Header.atoms );

		// Now lets find the atomic positions
		ifstream traFile;
		try
		{
			traFile.open(filename.c_str(), ios::in | ios::binary);

			if( !traFile.is_open() )
			{
				THROW(IOException,"'BTF_ImportBase' could not open the file '" + filename + "' for reading!");
			}

			seekToEntry(traFile,position); // make sure we import the correct frame!
			Tvector<float> getpos; // A Temporary holder for float, required as the Particle uses doubles
			for( int i = 0; i < m_Header.atoms; i++ )
			{
				// read from our tra file
				traFile.read((char*)&getpos,sizeof(float)*3);
				SnapShotAtom& atom = ss.atom[i];
				atom.p.setTo( getpos );
				atom.f.zero();
				atom.v.zero();
			}
		}
		catch( ExceptionBase )
		{
			if( traFile.is_open() ) traFile.close();
			throw; // Close the file handle and re-throw the exception.
		}

		// Close the file handle
		if( traFile.is_open() ) traFile.close();
		
		return ss;
	}


	BTF_Tools::BTF_Tools(std::string &_fileName, bool _ValidateFile )
		: BTF_ImportBase(_fileName,_ValidateFile)
	{
	}

	void BTF_Tools::getCoordinates( int _Entry )
	{
		getCoordinates( _Entry, All, -1 );
	}

	void BTF_Tools::getCoordinates( int _Entry, AtomFilter _filter )
	{
		getCoordinates( _Entry, _filter, -1 );
	}

	void BTF_Tools::getCoordinates( int _Entry, AtomFilter _filter, int _molNum )
	{
		if( _Entry == -1 )
		{
			_Entry = recountEntries() - 1;
		}
		if( _Entry >= 0 )
		{
			if( _Entry >= m_Entries ) recountEntries(); // Try this, more *might* have been added...
			if( _Entry >= m_Entries ) THROW(ArgumentException,"The getCoordinates() entry is not within range!");
		}
		else
		{
			THROW(ArgumentException,"The getCoordinates() entry is not within range!");
		}

		ifstream traFile;
		try
		{
			traFile.open(m_Filename.c_str(), ios::in | ios::binary);
			if( !traFile.is_open() )
			{
				THROW(IOException,"'BTF_ImportBase' could not open the file for validation!");
			}

			// Seek to the correct TrajectoryFormat entry
			traFile.seekg(m_Header.trajectorystart + 8 + m_Header.blocksize);
			assertTag(traFile,"TRAE");

			loadVectors( traFile, positions, _filter, _molNum );

			// skip phis and psis if present
			if( 0 <= (m_Header.Type & PhiPsis ) )
			{
				traFile.seekg(m_Header.residues*sizeof(float)*2,ios::cur); // seek past the Phi,Psi definitions
			}

			// '4' was 'rotamers' in the enumeration - this is now deprecated, but kept for backwards compatibility with old tra files
			if( 0 <= (m_Header.Type & 4 ) )
			{
				traFile.seekg(m_Header.residues*sizeof(int),ios::cur);
			}

			// read positions in, ognoring the ones we dont want...
			if( 0 <= (m_Header.Type & Energies ) )
			{
				ene.load(traFile);
			}

			// ForceVectors can also be stored
			if( 0 <= (m_Header.Type & ForceVectors ) )
			{
				loadVectors( traFile, forces, _filter, _molNum );
			}
		}
		catch( ExceptionBase )
		{
			if( traFile.is_open() ) traFile.close();
			throw; // Close the file handle and re-throw the exception.
		}

		// Close the file handle...
		if( traFile.is_open() ) traFile.close();
	}

	void BTF_Tools::savePDBFile( std::string &_FileName, int _Entry )
	{
		savePDBFile(_FileName,_Entry,All,-1);
	}

	void BTF_Tools::savePDBFile( std::string &_FileName, int _Entry, AtomFilter _filter )
	{
		savePDBFile(_FileName,_Entry,_filter,-1);
	}

	void BTF_Tools::savePDBFile( std::string &_FileStem, int _Entry, AtomFilter _filter, int _molNum )
	{
		PDBAtomLine atom;
		std::vector<PDBAtomLine> lines;

		getCoordinates(_Entry);

		int used = 0;
		for( int i = 0; i < m_Header.atoms; i++ )
		{
			if( passesFilter(i,_filter,_molNum) > 0 )
			{
				atom.altLoc = ' ';
				atom.atomName = sysDefs[i].pdbname;
				atom.atomNum = sysDefs[i].atomnumber;
				atom.chainID = ' ';
				atom.charge = " ";
				if(isASCIINumeric(sysDefs[i].pdbname[0])) atom.element = sysDefs[i].pdbname[1];
				else atom.element = sysDefs[i].pdbname[0];
				atom.iCode = ' ';
				atom.lineType = PDBAtomLine::ATOM;
				atom.occupancy = 0;
				atom.resName = sysDefs[i].parentname;
				atom.resNum = sysDefs[i].parentnumber;
				atom.segID = " ";
				atom.TempFactor = 0;

				Maths::dvector *pos = &positions[used++];
				atom.x = pos->x;
				atom.y = pos->y;
				atom.z = pos->z;

				lines.push_back(atom);
			}
		}

		PDB_RawWriter out;
		out.streamTo( _FileStem + ".pdb", false, true );
		out.beginStream();
		out.rawRemarks();
		out.rawWrite( lines );
		out.endStream();
	}




	int loadtra(System &sysspec, const std::string &filename, const std::string &entrystring)
	{
		int entry;

		// Interpret the 'entrystring' - this allows you to ask for the 'last' entry if desired...
		if(cmpstring(entrystring,"last"))
		{
			entry = -1;
		}
		else
		{
			if(str2int(entrystring,entry)!=0)
			{
				printf("ERROR: Unknown entry identifier '%s' \n",entrystring.c_str() );
				return -1;
			}
		}

		// Now load the chosen entry from the file...
		printf("Loading %s:%d \n",filename.c_str(),entry);
		InTra_BTF mytra( filename, true );
		mytra.loadIntoSystem( sysspec, entry );

		// Excellent :-D
		return 0;
	}





	// -------------------
	//  OutputFile_BTF 
	// -------------------



	void OutputFile_BTF::save( WorkSpace &_wspace){
		OutTra_BTF outtra( filestem, _wspace );
		outtra.append( );
	}





} // namespace IO


