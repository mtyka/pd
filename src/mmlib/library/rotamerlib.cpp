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

#include <fstream>
#include <algorithm>
#include <functional>

#include "tools/vector.h"
#include "sequence/sequence.h"
#include "forcefields/ffparam.h"
#include "system/molecule.h"

#include "rotamerlib.h"

namespace Library
{
	// Oddball helper function used extensively below
	std::string ensureNaming( const MoleculeDefinition& mol, bool altNames, const std::string& altName )
	{
		if( altNames ) return altName;
		int i = mol.findAtomRaw(altName);
		ASSERT( i != -1, CodeException, "ensureNaming(): FFPS internal mapping failure?!?");
		return mol.atom[i].pdbname;
	}

	RotConvention::RotConvention()
	{
		CBIsThirdDefaultAnchor = false; // Safety-first default
	}

	void RotConvention::clear()
	{
		m_Convention.clear();
	}

	int RotConvention::anchorAtomType( const std::string& _resName, const std::string& _atomName ) const
	{
		for( size_t i = 0; i < m_Convention.size(); i++ )
		{
			if( 0 == m_Convention[i].resName.compare( _resName ) )
			{
				if( 0 == m_Convention[i].a1.compare( _atomName ) ) return 1;
				if( 0 == m_Convention[i].a2.compare( _atomName ) ) return 2;
				if( 0 == m_Convention[i].a3.compare( _atomName ) ) return 3;
				return 0;
			}			
		}
		return -1;
	}

	void RotConvention::addConvention( const ConventionDef& _Conv )
	{
		for( size_t i = 0; i < m_Convention.size(); i++ )
		{
			if( 0 == m_Convention[i].resName.compare( _Conv.resName ) )
			{
				throw ArgumentException("The RotConvention container already contains a definition for resname: '" + _Conv.resName + "'");
			}
		}
		m_Convention.push_back(_Conv);
	}

	int RotConvention::findConvention( const std::string& _resName ) const
	{
		for( size_t i = 0; i < m_Convention.size(); i++ )
		{
			if( 0 == m_Convention[i].resName.compare( _resName ) )
			{
				return i;
			}
		}
		return -1;
	}

	const ConventionDef& RotConvention::getConvention( int _id ) const
	{
		ASSERT( _id >= 0 && _id < m_Convention.size(), OutOfRangeException, "RotConvention convention id is out of range");
		return m_Convention[_id];
	}

	void RotConvention::setToDefaultConvention()
	{
		clear();

		// The default anchor atom names
		ConventionDef cd;
		cd.a1 = "CA";
		cd.a2 = "N";
		// See the associated comment for `CBIsThirdDefaultAnchor'
		if( CBIsThirdDefaultAnchor )
		{
			cd.a3 = "CB";
		}
		else
		{
			cd.a3 = "C";
		}

		ConventionDef cd_gly;
		cd_gly.a1 = "CA";
		cd_gly.a2 = "N";
		cd_gly.a3 = "C";

		// Add the 20 standard residues to the container using the default anchor atoms
		cd.resName = "ALA"; m_Convention.push_back(cd);
		cd.resName = "CYS"; m_Convention.push_back(cd);
		cd.resName = "ASP"; m_Convention.push_back(cd);
		cd.resName = "GLU"; m_Convention.push_back(cd);
		cd.resName = "PHE"; m_Convention.push_back(cd);
		cd_gly.resName = "GLY"; m_Convention.push_back(cd_gly);
		cd.resName = "HIS"; m_Convention.push_back(cd);
		cd.resName = "ILE"; m_Convention.push_back(cd);
		cd.resName = "LYS"; m_Convention.push_back(cd);
		cd.resName = "LEU"; m_Convention.push_back(cd);
		cd.resName = "MET"; m_Convention.push_back(cd);
		cd.resName = "ASN"; m_Convention.push_back(cd);
		cd.resName = "PRO"; m_Convention.push_back(cd);
		cd.resName = "GLN"; m_Convention.push_back(cd);
		cd.resName = "ARG"; m_Convention.push_back(cd);
		cd.resName = "SER"; m_Convention.push_back(cd);
		cd.resName = "THR"; m_Convention.push_back(cd);
		cd.resName = "VAL"; m_Convention.push_back(cd);
		cd.resName = "TRP"; m_Convention.push_back(cd);
		cd.resName = "TYR"; m_Convention.push_back(cd);

		// Other ionisation states
		cd.resName = "HID"; m_Convention.push_back(cd);
		cd.resName = "HIE"; m_Convention.push_back(cd);
		cd.resName = "HIP"; m_Convention.push_back(cd);
		cd.resName = "CYX"; m_Convention.push_back(cd);
		cd.resName = "CYM"; m_Convention.push_back(cd);

		// Now for each of these contained residues, add their chi definitions
		// None for m_Convention[0] ALA
		m_Convention[1].chi.push_back(ChiDef("N","CA","CB","SG")); // CYS
		m_Convention[2].chi.push_back(ChiDef("N","CA","CB","CG")); // ASP
		m_Convention[2].chi.push_back(ChiDef("CA","CB","CG","OD1")); // ASP
		m_Convention[3].chi.push_back(ChiDef("N","CA","CB","CG")); // GLU
		m_Convention[3].chi.push_back(ChiDef("CA","CB","CG","CD")); // GLU
		m_Convention[3].chi.push_back(ChiDef("CB","CG","CD","OE1")); // GLU
		m_Convention[4].chi.push_back(ChiDef("N","CA","CB","CG")); // PHE
		m_Convention[4].chi.push_back(ChiDef("CA","CB","CG","CD1")); // PHE
		// None for m_Convention[5] GLY
		m_Convention[6].chi.push_back(ChiDef("N","CA","CB","CG")); // HIS
		m_Convention[6].chi.push_back(ChiDef("CA","CB","CG","ND1")); // HIS
		m_Convention[7].chi.push_back(ChiDef("N","CA","CB","CG1")); // ILE
		m_Convention[7].chi.push_back(ChiDef("CA","CB","CG1","CD1")); // ILE
		m_Convention[8].chi.push_back(ChiDef("N","CA","CB","CG")); // LYS
		m_Convention[8].chi.push_back(ChiDef("CA","CB","CG","CD")); // LYS
		m_Convention[8].chi.push_back(ChiDef("CB","CG","CD","CE")); // LYS
		m_Convention[8].chi.push_back(ChiDef("CG","CD","CE","NZ")); // LYS
		m_Convention[9].chi.push_back(ChiDef("N","CA","CB","CG")); // LEU
		m_Convention[9].chi.push_back(ChiDef("CA","CB","CG","CD1")); // LEU
		m_Convention[10].chi.push_back(ChiDef("N","CA","CB","CG")); // MET
		m_Convention[10].chi.push_back(ChiDef("CA","CB","CG","SD")); // MET
		m_Convention[10].chi.push_back(ChiDef("CB","CG","SD","CE")); // MET
		m_Convention[11].chi.push_back(ChiDef("N","CA","CB","CG")); // ASN
		m_Convention[11].chi.push_back(ChiDef("CA","CB","CG","OD1")); // ASN
		m_Convention[12].chi.push_back(ChiDef("N","CA","CB","CG")); // PRO
		m_Convention[12].chi.push_back(ChiDef("CA","CB","CG","CD")); // PRO
		m_Convention[13].chi.push_back(ChiDef("N","CA","CB","CG")); // GLN
		m_Convention[13].chi.push_back(ChiDef("CA","CB","CG","CD")); // GLN
		m_Convention[13].chi.push_back(ChiDef("CB","CG","CD","OE1")); // GLN
		m_Convention[14].chi.push_back(ChiDef("N","CA","CB","CG")); // ARG
		m_Convention[14].chi.push_back(ChiDef("CA","CB","CG","CD")); // ARG
		m_Convention[14].chi.push_back(ChiDef("CB","CG","CD","NE")); // ARG
		m_Convention[14].chi.push_back(ChiDef("CG","CD","NE","CZ")); // ARG
		// m_Convention[14].chi.push_back(ChiDef("CD","NE","CZ","NH1")); // ARG - No, its planar!
		m_Convention[15].chi.push_back(ChiDef("N","CA","CB","OG")); // SER
		m_Convention[16].chi.push_back(ChiDef("N","CA","CB","OG1")); // THR
		m_Convention[17].chi.push_back(ChiDef("N","CA","CB","CG1")); // VAL
		m_Convention[18].chi.push_back(ChiDef("N","CA","CB","CG")); // TRP
		m_Convention[18].chi.push_back(ChiDef("CA","CB","CG","CD1")); // TRP
		m_Convention[19].chi.push_back(ChiDef("N","CA","CB","CG")); // TYR
		m_Convention[19].chi.push_back(ChiDef("CA","CB","CG","CD1")); // TYR

		// Other ionisation states
		m_Convention[20].chi.push_back(ChiDef("N","CA","CB","CG")); // HID
		m_Convention[20].chi.push_back(ChiDef("CA","CB","CG","ND1")); // HID
		m_Convention[21].chi.push_back(ChiDef("N","CA","CB","CG")); // HIE
		m_Convention[21].chi.push_back(ChiDef("CA","CB","CG","ND1")); // HIE
		m_Convention[22].chi.push_back(ChiDef("N","CA","CB","CG")); // HIP
		m_Convention[22].chi.push_back(ChiDef("CA","CB","CG","ND1")); // HIP
		m_Convention[23].chi.push_back(ChiDef("N","CA","CB","SG")); // CYX
		m_Convention[24].chi.push_back(ChiDef("N","CA","CB","SG")); // CYM
	}

	RotLibConvertBase::RotLibConvertBase()
	{
		setToDefaultConvention();
	}

	void RotLibConvertBase::setUseAltnames( RotamerLibrary& _Lib, bool _Use ) const
	{
		_Lib.setUseAltnames(_Use);
	}

	RotLibConvert_PD::RotLibConvert_PD()
		: RotLibConvertBase()
	{
		setToDefaultConvention();
	}

	void RotLibConvert_OldPDFormat::readLib( const std::string& _LibraryFileName, RotamerLibrary& _RotLib )
	{
		if( _RotLib.isFinalised() ) 
			throw ProcedureException("readLib() is not allowed, the rotamer library has been finalised, no further import can occur");

		// Flag that we would like to import using the FFPS altnames! (AKA rawNames)
		setUseAltnames(_RotLib,true);

		// Force the addition of the HID HIE and HIP ionisation state mappers
		_RotLib.addIonisationAliasWorkingDefaults();

		std::ifstream file(_LibraryFileName.c_str(), std::ifstream::in);
		if( !file.is_open() ) 
			throw(IOException( "Rotamer file not found: '" + _LibraryFileName + "'!" ));

		RotamerAtom importAtom;
		std::vector<RotamerAtom> importAtoms;
		std::string line;
		StringBuilder sb;
		StringBuilder resname(' ', 4);
		StringBuilder temp(4);
		int resIndex = INT_MAX;

		while(sb << file)
		{
			sb.Trim();
			if( sb.size() == 0 ) continue; // Empty line
			if( sb[0] == '#' ) continue; // Comment line
			if( !sb.compare("ATOM  ",0,6,0,false) ) continue;

			temp.setTo(sb,12,4);
			temp.Trim(); // Turn "N   " into "N" for PD FFParam compatibility
			importAtom.Name = temp.toString();
			importAtom.IdealPos.setTo(
				sb.parseDouble(30,8),
				sb.parseDouble(38,8),
				sb.parseDouble(46,8)				
				);

			int currentResIndex = sb.parseInt(20,6);
			if( resIndex != currentResIndex || !sb.compare(resname,17,3,0) )
			{
				// on the first itteration of the cycle, we must add the first atom, so test for INT_MAX
				if( resIndex != INT_MAX ) 
				{
					// Add the new rotamer
					_RotLib.addRotamer(resname.toString(),importAtoms, ConstantProbability() );
					// Reinitialise residue
					importAtoms.clear();
				}
				importAtoms.push_back(importAtom);
				resIndex = currentResIndex;
				resname.setTo(sb,17,3);
			}
			else
			{
				importAtoms.push_back(importAtom);
			}

			continue;
		}

		if( importAtoms.size() > 0 )
		{
			// Add the final rotamer
			_RotLib.addRotamer(resname.toString(),importAtoms, ConstantProbability());
		}
	}

	RotamerAtom::RotamerAtom() 
		: Name(), 
		IsHydrogen(false), 
		IdealPos(DBL_MAX,DBL_MAX,DBL_MAX), 
		Defined(true)
	{
	}

	RotamerAtom::RotamerAtom(const std::string &_name, const Maths::dvector &_idealPos, bool _isHydrogen )
		: Name(_name), IdealPos(_idealPos), IsHydrogen(_isHydrogen), Defined(true)
	{
	}

	// Used by RotamerLibrary::shouldRemoveRotamerAtom()
	bool RotamerAtom::operator==( const std::string& resName ) const
	{ 
		return 0 == resName.compare(Name);
	}

	inline void parseScope( std::istream& file, StringBuilder& sb, CloneHolder<ProbabilityBase>& scope )
	{
		sb.TruncateLeftBy(11);
		sb.Trim();
		ASSERT( sb.size() > 0, ParseException, "Expected a mode after SCOPE_START");
		if( sb.compare("CONST",0,5,0,true) )
		{
			scope = ConstantProbability();
			scope->fromStream( file );	
		}
		else if( sb.compare("BBScope",0,7,0,true) )
		{
			scope = ProbabilityByBBScope();
			scope->fromStream( file );
		}
		else if( sb.compare("PhiPsiMap",0,9,0,true) )
		{
			sb.TruncateLeftBy(9);
			sb.Trim();
			ASSERT(sb.size()>0,ParseException,"Expected a PhiPsiMap resolution!");
			double degreeResln = sb.parseDouble();
			ASSERT( 0.0 < degreeResln && degreeResln <= 360.0, ParseException, 
				"Imported PhiPsiMap propensity bin resolution is not within the sensible range");
			double remainder =  360.0 - floor( 360.0 / degreeResln ) * degreeResln;
			ASSERT( Maths::SigFigEquality(remainder, 0.0, 8 ),
				ParseException, "Imported PhiPsiMap propensity bin resolution does not yield an integer number of bins");
			int binResln = (int)(360.0 / degreeResln);
			scope = ProbabilityByPhiPsiMap(binResln);
			scope->fromStream( file );
		}
		else
		{
			sb.insert(0,"Unknown scope type: '");
			sb.append('\'');
			throw ParseException( sb.toString() );
		}
	}

	void RotLibConvert_PD::readBlock( std::istream& file, const StringBuilder& resName, RotamerLibrary& _RotLib ) const
	{
		double chi;
		RotamerAtom importAtom;
		std::vector<RotamerAtom> importAtoms;
		std::vector<double> importChis;
		int resIndex = INT_MAX;
		StringBuilder sb;
		StringBuilder temp(4);
		bool additionOccured = false;
		CloneHolder<ProbabilityBase> scope = ConstantProbability(); // default
		int wantChi = 0;

		bool quit = false;
		bool goScope = false;
		bool atomTodo = false;
		bool chiTodo = false;
		while( sb << file && !quit )
		{
			sb.Trim();
			if( sb.size() == 0 ) continue; // Empty line
			if( sb[0] == '#' ) continue; // Comment line

			if( sb.compare("END   ",0,6,0,false) )
			{
				sb.TruncateLeftBy(3);
				sb.Trim();
				if( sb.size() == 0 )
				{
					throw ParseException("Expected a residue name following END statement");
				}
				if( !sb.compare( resName, 0, 3, 0 ) )
				{
					throw ParseException("Residue name misamatch");
				}
				if( !additionOccured )
				{
					_RotLib.addAsBlankRotamer(resName.toString());
				}

				quit = true; // End of residue definition
				// fall through to rotamer addition code
			}
			else if( sb.compare("CHI   ",0,6,0,false) )
			{
				additionOccured = true;
				if( importAtoms.size() > 0 )
					throw ParseException("CHI Definition found in a block which already contains ATOM definitions");
				int i = sb.parseInt(6,4);	
				if( i != ++wantChi )
				{
					throw ParseException("CHI definition sequenceID error");
				}
				int currentResIndex = sb.parseInt(11,4);
				chi = Maths::DegToRad(sb.parseDouble(16));

				if( resIndex == INT_MAX ) 
				{
					// on the first itteration of the cycle, we must add the first chi, so test for INT_MAX
					resIndex = currentResIndex;
				}
				if( resIndex == currentResIndex )
				{				
					importChis.push_back(chi);
					continue;
				}

				resIndex = currentResIndex;
				chiTodo = true;
				// fall through to rotamer addition code
			}
			else if( sb.compare("ATOM  ",0,6,0,false) )
			{
				additionOccured = true;
				if( importChis.size() > 0 )
					throw ParseException("ATOM definition found in a block which already contains CHI definitions");

				temp.setTo(sb,12,4);
				temp.Trim(); // Turn "N   " into "N" for PD FFParam compatibility
				importAtom.Name = temp.toString();

				if(!sb.compare(resName,17,3,0))
				{
					throw ParseException("Residue name misamatch before END statement");
				}

				int currentResIndex = sb.parseInt(20,6);

				importAtom.IdealPos.setTo(
					sb.parseDouble(30,8),
					sb.parseDouble(38,8),
					sb.parseDouble(46,8)				
					);

				if( resIndex == INT_MAX ) 
				{
					// on the first itteration of the cycle, we must add the first chi, so test for INT_MAX
					resIndex = currentResIndex;
				}
				if( resIndex == currentResIndex )
				{				
					importAtoms.push_back(importAtom);
					continue;
				}

				resIndex = currentResIndex;
				atomTodo = true;
				// fall through to rotamer addition code
			}
			else if( sb.compare("SCOPE_START",0,11,0,false))
			{
				goScope = true;
				// fall through
			}
			else
			{
				Printf("Warning, unparsed line!!: \"%s\"")(sb.toString());
				continue;
			}

			// If we have falled through to here, we need to add the rotamer!
			if( importAtoms.size() > 0 )
			{
				// Add the final rotamer
				_RotLib.addRotamer(resName.toString(),importAtoms, scope.data());
				importAtoms.clear();
			}
			else if( importChis.size() > 0 )
			{
				// Add the final rotamer
				_RotLib.addRotamer(resName.toString(),importChis, scope.data());
				importChis.clear();
				wantChi = 0;
			}

			// Final stuff to finish off
			if( goScope )
			{
				parseScope( file, sb, scope );
				goScope = false;
			}
			else if( atomTodo )
			{
				importAtoms.push_back(importAtom);
				atomTodo = false;
			}
			else if( chiTodo )
			{
				importChis.push_back(chi);
				chiTodo = false;
			}
		}
	}

	void RotLibConvert_PD::readLib( const std::string& _LibraryFileName, RotamerLibrary& _RotLib )
	{
		if( _RotLib.isFinalised() ) 
			throw ProcedureException("readLib() is not allowed, the rotamer library has been finalised, no further import can occur");

		std::ifstream file(_LibraryFileName.c_str(), std::ifstream::in);
		if( !file.is_open() ) throw(IOException( "Rotamer file not found: '" + _LibraryFileName + "'!" ));

		const char* const headerString = "## PD_ROTAMER_FORMAT v";
		const double fileParserVersion = 1.0;
		StringBuilder sb;
		sb << file;
		if(!sb.compare(headerString,0,22,0,false))
		{
			throw ParseException("Rotamer file does not define 'PD_ROTAMER_FORMAT' it is probably not a PD format file");
		}
		else
		{
			try
			{
				double version = sb.parseDouble(22);
				if( version != fileParserVersion )
					throw ParseException("PD_ROTAMER_FORMAT parser version 1.0 being asked to read a rotamer file of incompatible version!");
			}
			catch( ExceptionBase ex )
			{
				throw ParseException("PD_ROTAMER_FORMAT does not specify a version number??" );
			}
		}

		sb << file;
		if( sb.size() < 10 || !sb.compare("NAMING ",0,7,0) )
		{
			THROW(ParseException,"No NAMING directive found");
		}
		else if( sb.compare("ALT",7,3,0) )
		{
			setUseAltnames(_RotLib,true); // Flag that we would like to import using the FFPS altnames!
		}
		else if( sb.compare("PDB",7,3,0) )
		{
			setUseAltnames(_RotLib,false); // Flag that we would like to import using the PDB altnames!
		}
		else
		{
			THROW(ParseException,"Unknown NAMING type");
		}

		while(sb << file)
		{
			sb.Trim();
			if( sb.size() == 0 ) continue; // Empty line
			if( sb[0] == '#' ) continue; // Comment line
			if( sb.compare("BEGIN ",0,6,0,false) )
			{
				sb.TruncateLeftBy(6);
				sb.Trim();
				if( sb.size() == 0 )
					throw ParseException("Expected a residue name following BEGIN token");
				readBlock(file,sb,_RotLib);
			}	
		}		
	}

	void RotLibConvert_PD::writeRotamerSet( std::ofstream& stream, Library::WriteLibMode _mode, const RotamerSet& rotset ) const
	{
		stream << "BEGIN " << rotset.getName() << std::endl;
		StringBuilder sb(80);

		const std::vector<RotamerAtom> atoms = rotset.getAtoms();
		const std::vector<RotamerChi> chis = rotset.getChis();

		size_t iRot = 0;
		size_t iAt = 0;
		for( size_t i = 0; i < rotset.nRot(); i++ )
		{
			const Rotamer& rot = rotset.getRotamer(i);

			// Per-rotamer probability and scope definition
			rot.getProbMap().toStream( stream );

			if( _mode == WriteCartesian || (_mode == WriteOriginal && rot.getSource() == WriteCartesian) )
			{
				for( size_t j = 0; j < rot.nAtom(); j++ )
				{
					const RotamerAtom& atom = atoms[j];
					const Maths::dvector& pos = rot.getPos(j);
					sb.setFormat("ATOM  %5d %-4s %3s%6d    %8.3f%8.3f%8.3f")
						(iAt++)(atom.Name)(rotset.getName())(iRot)(pos.x)(pos.y)(pos.z);
					sb.PadRight(80,' ');
					stream << sb << std::endl;
				}
				iRot++;
			}
			else if( _mode == WriteTorsional || (_mode == WriteOriginal && rot.getSource() == WriteTorsional) )
			{
				for( size_t j = 0; j < rot.nChi(); j++ )
				{
					sb.setFormat("CHI   %4d %4d %8.3f")
						(j+1)(i+1)( Maths::RadToDeg(rot.getChi(j)) );
					sb.PadRight(80,' ');
					stream << sb << std::endl;
				}
				iRot++;
			}
			else
			{
				THROW(CodeException,"Something strange in usage of the WriteLibMode enumeration");
			}
		}

		stream << "END   " << rotset.getName() << std::endl << std::endl;
	}

	void RotLibConvert_PD::writeLib( const std::string& _LibraryFileName, Library::WriteLibMode _mode, const RotamerLibrary& _RotLib ) const
	{
		std::ofstream file(_LibraryFileName.c_str(), std::ifstream::out);
		if( !file.is_open() ) throw(IOException( "Rotamer lib output file could not be opened for writing: '" + _LibraryFileName + "'!" ));
		file << "## PD_ROTAMER_FORMAT v1.0" << std::endl;
		file << (_RotLib.useAltNames() ? "NAMING ALT" : "NAMING PDB") << std::endl << std::endl;
		for( size_t i = 0; i < _RotLib.nRot(); i++ )
		{
			writeRotamerSet(file,_mode,_RotLib.getRotamerSet(i));
		}
		file << "## PD_END_FILE" << std::endl;
		file.close();
	}

	RotamerChi::RotamerChi( const ChiDef& _chi, const std::vector<RotamerAtom>& _atoms ) 
		: ChiDef(_chi)
	{
		// Obtain the rotamer _atoms array indeces
		for( size_t n = 0; n < _atoms.size(); n++ )
		{
			if( 0 == _atoms[n].Name.compare(p) )
			{
				rot.i = (int)n;
				break;
			}
		}
		for( size_t n = 0; n < _atoms.size(); n++ )
		{
			if( 0 == _atoms[n].Name.compare(q) )
			{
				rot.a = (int)n;
				break;
			}
		}
		for( size_t n = 0; n < _atoms.size(); n++ )
		{
			if( 0 == _atoms[n].Name.compare(r) )
			{
				rot.b = (int)n;
				break;
			}
		}
		for( size_t n = 0; n < _atoms.size(); n++ )
		{
			if( 0 == _atoms[n].Name.compare(s) )
			{
				rot.j = (int)n;
				break;
			}
		}
	}

	BackboneScope::BackboneScope() 
	{
		set(Undefined);
	}

	BackboneScope::BackboneScope(BackboneScopeType _type) 
	{
		set(_type);
	}

	void BackboneScope::set(BackboneScopeType _type)
	{
		m_BBScope = _type;
		switch( m_BBScope )
		{
		case Specific:
			// Fall through to the default 'Undefined' range
		case Undefined:
			{
				m_Phi = TorsionalTolleranceRange(0.0,180.0);
				m_Psi = TorsionalTolleranceRange(0.0,180.0);
				return;
			}
		case Alpha:
			{
				m_Phi = PhiAplpha;
				m_Psi = PsiAplpha;
				return;
			}
		case Beta:
			{
				m_Phi = PhiBeta;
				m_Psi = PsiBeta;
				return;
			}
		case LeftAlpha:
			{
				m_Phi = PhiLAplpha;
				m_Psi = PsiLAplpha;
				return;
			}
		default:
			{
				THROW(CodeException,"Unknown BackboneScopeType");
			}
		}
	}

	void BackboneScope::setSpecific( TorsionalTolleranceRange _Phi, TorsionalTolleranceRange _Psi )
	{
		m_BBScope = Specific;
		m_Phi = _Phi;
		m_Psi = _Psi;
	}

	const TorsionalTolleranceRange BackboneScope::PhiAplpha =  TorsionalTolleranceRange( Maths::DegToRad( -65.0), Maths::DegToRad(35.0), Maths::DegToRad(-85.0) );
	const TorsionalTolleranceRange BackboneScope::PsiAplpha =  TorsionalTolleranceRange( Maths::DegToRad( -45.0), Maths::DegToRad(95.0), Maths::DegToRad(-25.0) );
	const TorsionalTolleranceRange BackboneScope::PhiBeta =    TorsionalTolleranceRange( Maths::DegToRad(-110.0), Maths::DegToRad(80.0), Maths::DegToRad(-80.0) );
	const TorsionalTolleranceRange BackboneScope::PsiBeta =    TorsionalTolleranceRange( Maths::DegToRad( 135.0), Maths::DegToRad(55.0), Maths::DegToRad(-85.0) );
	const TorsionalTolleranceRange BackboneScope::PhiLAplpha = TorsionalTolleranceRange( Maths::DegToRad(  87.0), Maths::DegToRad(63.0), Maths::DegToRad(-27.0) );
	const TorsionalTolleranceRange BackboneScope::PsiLAplpha = TorsionalTolleranceRange( Maths::DegToRad(  47.0), Maths::DegToRad(28.0), Maths::DegToRad(-77.0) );

	TorsionalTolleranceRange BackboneScope::getPhi() const
	{
		return m_Phi;
	}

	TorsionalTolleranceRange BackboneScope::getPsi() const
	{		
		return m_Psi;
	}

	std::string BackboneScope::toString() const
	{
		switch(m_BBScope)
		{
		case Undefined:
			return "U";
		case Alpha:
			return "A";
		case Beta:
			return "B";
		case LeftAlpha:
			return "L";
		case Specific:
			{
				StringBuilder sb;
				sb.setFormat("S %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f")
					(m_Phi.getValue())(m_Phi.getPlus())(m_Phi.getMinus())(m_Psi.getValue())(m_Psi.getPlus())(m_Psi.getMinus());
				return sb.toString();
			}
		default:
			THROW(CodeException,"Unknown BackboneScopeType encountered"); 
		}
	}

	void BackboneScope::parseString( const std::string& s )
	{
		ASSERT(s.size()!=0,ParseException,"RotamerSet::parseBBScopeString() String cannot be zero length");
		switch(s[0])
		{
		case 'U':
			set(Undefined);
			return;
		case 'A':
			set(Alpha);
			return;
		case 'B':
			set(Beta);
			return;
		case 'L':
			set(LeftAlpha);
			return;
		case 'S':
			{
				StringBuilder sb(s);
				setSpecific(
					TorsionalTolleranceRange(
					sb.parseDouble(2,7),
					sb.parseDouble(10,7),
					sb.parseDouble(18,7)),
					TorsionalTolleranceRange(
					sb.parseDouble(26,7),
					sb.parseDouble(34,7),
					sb.parseDouble(42,7) ) );
			}
		default:
			THROW(CodeException,"Unknown BackboneScopeType encountered"); 
		}
	}

	ConstantProbability::ConstantProbability()
	{
		m_Prob = 1.0; // default to enabled, probability 1.0 (they will all be this, so equal probability, unless changed)
	}

	ConstantProbability::ConstantProbability( double _prob )
	{
		m_Prob = _prob;
	}

	ConstantProbability* ConstantProbability::clone() const
	{
		return new ConstantProbability(*this);
	}

	void ConstantProbability::setProbability( double value )
	{
		m_Prob = value;
	}

	void ConstantProbability::clear()
	{
		m_Prob = -1.0;
	}

	void ConstantProbability::toStream( std::ostream& _Stream ) const ///< Output this structure to a stream
	{
		_Stream << "SCOPE_START CONST" << std::endl;
		_Stream << "PROB  " << double2str(m_Prob,"%8.5lf") << std::endl;
	}

	void ConstantProbability::fromStream( std::istream& _Stream ) ///< Set this structure using data in a stream
	{
		StringBuilder sb;
		sb << _Stream;
		ASSERT(sb.size() >=7 && isxdigit(sb[7]) && sb.compare("PROB  ",0,6,0 ),ParseException,"Probability string not formatted correctly");
		m_Prob = sb.parseDouble(7);
	}

	double ConstantProbability::probability( const MoleculeBase& _mol, int ir ) const ///< Returns a set probability
	{
		return m_Prob;
	}

	ProbabilityByBBScope::ProbabilityByBBScope():ConstantProbability()
	{
		m_BBScope.set( BackboneScope::Undefined );
	}

	ProbabilityByBBScope::ProbabilityByBBScope( double _prob ) : ConstantProbability(_prob)
	{
		m_BBScope.set( BackboneScope::Undefined );
	}

	ProbabilityByBBScope::ProbabilityByBBScope( BackboneScope::BackboneScopeType _scope, double _prob ) : ConstantProbability(_prob)
	{
		m_BBScope.set( _scope );
	}

	ProbabilityByBBScope* ProbabilityByBBScope::clone() const
	{
		return new ProbabilityByBBScope(*this);
	}

	void ProbabilityByBBScope::toStream( std::ostream& _Stream ) const ///< Output this structure to a stream
	{
		_Stream << "SCOPE_START BBScope" << std::endl;
		ConstantProbability::toStream(_Stream); // call the base class to get the probability
		_Stream << m_BBScope.toString() << std::endl;
	}

	void ProbabilityByBBScope::fromStream( std::istream& _Stream ) ///< Set this structure using data in a stream
	{
		ConstantProbability::fromStream(_Stream); // call the base class to set the probability
		StringBuilder sb;
		sb << _Stream;
		ASSERT( sb.compare("BBSCOPE ",0,8,0,false), ParseException, "BBSPOPE expected");
		sb.TruncateLeftBy(8);
		m_BBScope.parseString(sb.toString());
	}

	double ProbabilityByBBScope::probability( const MoleculeBase& _mol, int ir ) const ///< Returns a set probability if in-scope, or -1.0 if out of scope
	{
		double phi, psi;
		_mol.calcResiduePhiPsi(ir, phi, psi );
		if( m_BBScope.isInRange( phi, psi ) )
		{
			return m_Prob;
		}
		else
		{
			return -1.0;
		}
	}

	ProbabilityByPhiPsiMap::ProbabilityByPhiPsiMap( size_t binCount ) 
		: m_BBMap(binCount,binCount)
	{	
		clear(); // sets m_Assignments and clears the BBMap
		m_RadResln = Maths::MathConst::TwoPI / (double)binCount;
	}

	void ProbabilityByPhiPsiMap::clear()
	{
		m_Assignments = 0;
		m_BBMap = -1.0; // all states in the 2D grid are now -1.0; default to disabled
	}

	ProbabilityByPhiPsiMap* ProbabilityByPhiPsiMap::clone() const
	{
		return new ProbabilityByPhiPsiMap(*this);
	}

	inline double angleFor( const double binResln, size_t _index )
	{
		// Place in the range -PI -> +PI
		return Maths::RadToDeg((binResln * (double)_index) - Maths::MathConst::PI);
	}

	void ProbabilityByPhiPsiMap::toStream( std::ostream& _Stream ) const ///< Output this structure to a stream
	{
		_Stream << "SCOPE_START PhiPsiMap " << double2str(Maths::RadToDeg(m_RadResln),"%6.2lf") << std::endl;
		for( int i = 0; i < m_BBMap.dim1(); i++ )
		{
			for( int j = 0; j < m_BBMap.dim2(); j++ )
			{
				// if -1, it was never defined as a valid state, but will obviously be in a 2D grid.
				if( m_BBMap[i][j] != -1.0 )
				{
					_Stream << "SCOPE "
						<< " " << double2str(angleFor(m_RadResln,i),"%7.2lf") // phi
						<< " " << double2str(angleFor(m_RadResln,j),"%7.2lf") // psi
						<< double2str(m_BBMap[i][j],"%10.7lf") // propensity
						<< std::endl;
				}
			}
		}
		_Stream << "SCOPE_STOP" << std::endl;
	}

	void ProbabilityByPhiPsiMap::fromStream( std::istream& _Stream ) ///< Set this structure using data in a stream
	{
		StringBuilder sb;
		while( sb << _Stream )
		{
			if( sb.size() < 8 )
			{
				// fall through
			}
			else if( sb.compare("SCOPE_STOP",0,10,0) )
			{
				// Our work here is done!
				return;
			}
			else if( sb.compare("SCOPE ",0,6,0) ) 
			{
				sb.TruncateLeftBy(6);
				sb.Trim();
				ASSERT(sb.size()>0,ParseException,"Expected data following SCOPE directive");
				assignFrom( sb );
				continue;
			}
			sb.insert(0,"Warning!! Line not parsed: ");
			std::cout << sb;
		}
		THROW( ParseException, "Unexpected end of Sream whilst parsing ProbabilityByPhiPsiMap!");
	}

	double ProbabilityByPhiPsiMap::probability( const MoleculeBase& _mol, int ir ) const ///< Returns a probability based upon the current phi / psi angles
	{
		double phi, psi;
		_mol.calcResiduePhiPsi(ir, phi, psi );
		phi = torsionalRangeRadians(phi)+Maths::MathConst::PI; // +PI to place in the range 0->2PI as opposed to +-PI.
		psi = torsionalRangeRadians(psi)+Maths::MathConst::PI;
		return m_BBMap[(int)floor(phi / m_RadResln)][(int)floor(psi / m_RadResln)];			
	}

	void ProbabilityByPhiPsiMap::setTo( const ProbabilityByPhiPsiMap& _clone )
	{
		m_Assignments = _clone.m_Assignments;
		m_BBMap = _clone.m_BBMap;
		m_RadResln = _clone.m_RadResln;
	}

	void ProbabilityByPhiPsiMap::assignFrom( StringBuilder& _sb ) ///< Read in phi psi and probability from a whitespace delimited string, consuming the string in the process
	{
		m_Assignments++;

		// A phi and psi value must be converted into an integer grid location
		_sb.Trim(); // ensure oddball whitespace doesnt mess up below

		size_t cull = _sb.FirstOf(TOKEN_WHITESPACE);
		ASSERT(cull != SIZE_T_FAIL,ParseException,"SCOPE directive incorrect format");
		double phi = Maths::DegToRad(_sb.parseDouble(0,cull));
		_sb.TruncateLeftBy(cull);
		_sb.Trim();
		ASSERT(_sb.size() > 0,ParseException,"SCOPE directive incorrect format");

		cull = _sb.FirstOf(TOKEN_WHITESPACE);
		ASSERT(cull != SIZE_T_FAIL,ParseException,"SCOPE directive incorrect format");
		double psi = Maths::DegToRad(_sb.parseDouble(0,cull));
		_sb.TruncateLeftBy(cull);
		_sb.Trim();
		ASSERT(_sb.size() > 0,ParseException,"SCOPE directive incorrect format");

		double probability = _sb.parseDouble();

		assignFrom( phi, psi, probability );

		_sb.clear();
	}

	void ProbabilityByPhiPsiMap::assignFrom( double phi, double psi, double probability )
	{
		m_Assignments++;

		double indexXDouble = ((phi+Maths::MathConst::PI) / m_RadResln) + 0.5; // 0.5 assures that 0.9999 ends up in the 1.0 bin!!
		double indexDelta = indexXDouble - floor(indexXDouble);
		ASSERT( Maths::SigFigEquality( indexDelta, 0.5, 3), ParseException, "SCOPE directive phi doesn't sit on the grid");
		int indexX = (int)indexXDouble;

		double indexYDouble = ((psi+Maths::MathConst::PI) / m_RadResln) + 0.5; // 0.5 assures that 0.9999 ends up in the 1.0 bin!!
		indexDelta = (indexYDouble) - floor(indexYDouble); 
		ASSERT( Maths::SigFigEquality( indexDelta, 0.5, 3), ParseException, "SCOPE directive psi doesn't sit on the grid");
		int indexY = (int)indexYDouble;

		m_BBMap[indexX][indexY] = probability;
	}

	void Rotamer::info() const
	{
		THROW(NotImplementedException,"void Rotamer::info() const");
	}

	RotamerSet::RotamerSet( const std::string& _ResName, bool _AltNames, const RotConvention& _Conv, const MoleculeDefinition& _FFMol )
		: name(_ResName), 
		m_AltNames(_AltNames),
		m_CartesianAnchorRot1(-1), 
		m_CartesianAnchorRot2(-1), 
		m_CartesianAnchorRot3(-1), 
		m_FFMol(&_FFMol),
		m_Conv(&_Conv)
	{
	}

	int RotamerSet::closestRotamer( MoleculeBase& _molBase, int ir, double& cRMS ) const
	{
		THROW(NotImplementedException,"RotamerSet::closestRotamer()");
	}

	void RotamerSet::calcAtomIndexMap( const MoleculeBase& _molBase, int ir, std::vector<size_t>& _map ) const
	{
		const ParticleStore& atom = _molBase.atom;
		const Residue& res = _molBase.res[ir];
		size_t len = nAtom();
		size_t from = res.ifirst;
		size_t to = res.ilast;
		_map.resize(len,SIZE_T_FAIL);

		for( size_t i = 0; i < len; i++ )
		{
			for( size_t j = from; j <= to; j++ )
			{
				if( 0 == (m_AltNames ? atom[j].rawname : atom[j].pdbname).compare(m_Atom[i].Name) )
				{
					_map[i] = j;
					break;
				}
			}
			ASSERT(_map[i] != SIZE_T_FAIL, CodeException, "Atom lookup failure");
		}
	}

	void RotamerSet::applyRotamerCartesian( MoleculeBase& _molBase, int ir, size_t _rotIndex ) const
	{
		// static to ensure we aren't reallocating all the time
		static std::vector<size_t> _map; 
		calcAtomIndexMap(_molBase, ir, _map);
		applyRotamerCartesian(_molBase, ir, _rotIndex, _map );
	}

	void RotamerSet::applyRotamerCartesian( MoleculeBase& _molBase, int ir, size_t _rotIndex, const std::vector<size_t>& _indexMap ) const
	{
		D_ASSERT(_molBase.getSequence().getResName(ir).compare( name ) == 0, 
			CodeException, 
			"Refered residue does not correspond in type to this RotamerSet");

		const Residue& res = _molBase.res[ir];
		const Rotamer& rot = m_RotRes[_rotIndex];

		Maths::dvector r1 = rot.getPos(m_CartesianAnchorRot1);
		Maths::dvector r2 = rot.getPos(m_CartesianAnchorRot2);
		Maths::dvector r3 = rot.getPos(m_CartesianAnchorRot3);
		Maths::dvector a1 = _molBase.atom[_indexMap[m_CartesianAnchorRot1]].pos();
		Maths::dvector a2 = _molBase.atom[_indexMap[m_CartesianAnchorRot2]].pos();
		Maths::dvector a3 = _molBase.atom[_indexMap[m_CartesianAnchorRot3]].pos();

		// **Note**: Rotamer anchor atom 1 of the rotamer definition was ensured to be 
		// sitting on the origin during setup
		D_ASSERT( 
			Maths::SigFigEquality(r1.x, 0.0, 8) && 
			Maths::SigFigEquality(r1.y, 0.0, 8) && 
			Maths::SigFigEquality(r1.z, 0.0, 8), 
			CodeException, "Rotamer anchor atom should be sitting on the origin, it is not!");

		// superimpose the anchor sets onto each other
		Maths::matrix3x3 rmat;
		Maths::superimpose(
			a1, 
			a2, // The anchor atoms within the MoleculeBase
			a3,
			r1, 
			r2, // The anchor atoms within the rotamer definition
			r3,
			rmat);

		for( size_t i = 0; i < rot.nAtom(); i++ )
		{
			if( m_CartesianAnchorRot1 == i ||
				m_CartesianAnchorRot2 == i ||
				m_CartesianAnchorRot3 == i )
			{
				// We dont move the anchors
				continue;
			}

			// The index of the atom will be the rotamer atom index offset 
			// + the first atom index in the current residue
			const int moleculeAtomIndex = _indexMap[i];

			// Assert that the atom names match when in debug mode!
			D_ASSERT(
				(m_AltNames 
				? _molBase.atom[moleculeAtomIndex].rawname 
				: _molBase.atom[moleculeAtomIndex].pdbname)
				.compare(m_Atom[i].Name) == 0,
				CodeException, "Atom name mapping failure");

			// Perform the transform on the atomic coordinate
			Maths::dvector& pos = _molBase.atom[moleculeAtomIndex].pos();
			pos.setTo(rot.m_Pos[i]);
			//pos.sub(r1); its on 0,0,0
			pos.mulmat(rmat);
			pos.add(a1);
		}
	}

	void RotamerSet::calcCartesian(Rotamer& r) const
	{
		// Ensure a sensible starting conformation!

		// 1) Set r.m_pos to the same length as the RotamerSet atom collection
		r.m_Pos.resize( m_Atom.size() );

		// 2) Assign those atom coordinates from the corresponding 'IdealPos' (derived from the forcefield posGeom())
		for( size_t i = 0; i < m_Atom.size(); i++ )
		{
			r.m_Pos[i].setTo( m_Atom[i].IdealPos );
		}

		// 3) Calculate what we need to move
		ASSERT(m_FFMol!=NULL,CodeException,"RotamerSet::m_FFMol is NULL");
		std::vector<IndexHexet> scopeMaps;
		calcTorsionalScope_Rot(*m_FFMol,scopeMaps);

		// 4) Apply torsions ... Simple hu ;-)
		for( size_t i = 0; i < m_Chis.size(); i++ )
		{
			IndexHexet scopeMap = scopeMaps[i];
			Maths::matrix3x3 rmat;
			Maths::dvector axis;
			Maths::dvector origin( r.m_Pos[scopeMap.b] );

			double rotateCurrentAngle = 
				calcTorsionAngle( 
				r.m_Pos[scopeMap.a],
				r.m_Pos[scopeMap.b],
				r.m_Pos[scopeMap.c],
				r.m_Pos[scopeMap.d]
			);
			rotateCurrentAngle = r.m_Chis[i] - rotateCurrentAngle; // and now ... what do we need to rotate by?

			// assign the rotation matrix
			axis.setTo(r.m_Pos[scopeMap.c]);
			axis.sub(r.m_Pos[scopeMap.b]);
			rmat.setToAxisRot( axis, rotateCurrentAngle );

			// rotate all the atoms that are affetced by the rotation within the Scope of the loop definition
			for( int k = scopeMap.i; k <= scopeMap.j; k++)
			{
				r.m_Pos[k].sub(origin); // translate to axis origin
				r.m_Pos[k].mulmat(rmat); // rotate position
				r.m_Pos[k].add(origin); // translate back to original position
			}
		}

		// 5) When this function exits, the rotamer will be positioned on the origin by a second function.
	}

	IndexQuartet calcFFQuartet( const ChiDef& chi, const MoleculeDefinition& _mol, bool _altNaming )
	{
		IndexQuartet ffIndex;
		// Now the ff indexes
		const std::vector<AtomParameter> ffatoms = _mol.atom;
		for( size_t n = 0; n < ffatoms.size(); n++ )
		{
			if( 0 == chi.p.compare( _altNaming ? ffatoms[n].rawname : ffatoms[n].pdbname ) )
			{
				ffIndex.i = (int)n;
				break;
			}
		}
		for( size_t n = 0; n < ffatoms.size(); n++ )
		{
			if( 0 == chi.q.compare( _altNaming ? ffatoms[n].rawname : ffatoms[n].pdbname ) )
			{
				ffIndex.a = (int)n;
				break;
			}
		}
		for( size_t n = 0; n < ffatoms.size(); n++ )
		{
			if( 0 == chi.r.compare( _altNaming ? ffatoms[n].rawname : ffatoms[n].pdbname ) )
			{
				ffIndex.b = (int)n;
				break;
			}
		}
		for( size_t n = 0; n < ffatoms.size(); n++ )
		{
			if( 0 == chi.s.compare( _altNaming ? ffatoms[n].rawname : ffatoms[n].pdbname ) )
			{
				ffIndex.j = (int)n;
				break;
			}
		}
		return ffIndex;
	}

	void RotamerSet::calcTorsionalScope_Rot( const MoleculeDefinition& mol, std::vector<IndexHexet>& _scopeMap ) const
	{
		// Call the FF version
		calcTorsionalScope_FF( mol, _scopeMap );

		// Adapt to the rotamer index system
		const std::vector<AtomParameter>& atoms = mol.atom;
		for( size_t i = 0; i < _scopeMap.size(); i++ )
		{
			IndexHexet& set = _scopeMap[i];

			// Convert from the ffindex to the rot index	
			size_t testDeltaFF = set.j - set.i;
			set.i = FindFirstInVector(m_Atom, ensureNaming( mol, m_AltNames, atoms[set.i].rawname ) );
			set.j = FindFirstInVector(m_Atom, ensureNaming( mol, m_AltNames, atoms[set.j].rawname ) );
			size_t testDeltaRot = set.j - set.i;
			ASSERT( set.j >=0 && set.i >=0, CodeException, "Assumption failure");
			ASSERT( testDeltaFF == testDeltaRot, CodeException, "Assumption failure");

			// Lookup the indexes of the FFPS atoms within the rotamer definition
			set.a = FindFirstInVector(m_Atom, ensureNaming( mol, m_AltNames, atoms[set.a].rawname ) );
			set.b = FindFirstInVector(m_Atom, ensureNaming( mol, m_AltNames, atoms[set.b].rawname ) );
			set.c = FindFirstInVector(m_Atom, ensureNaming( mol, m_AltNames, atoms[set.c].rawname ) );
			set.d = FindFirstInVector(m_Atom, ensureNaming( mol, m_AltNames, atoms[set.d].rawname ) );
			ASSERT( 
				set.a != SIZE_T_FAIL &&
				set.b != SIZE_T_FAIL &&
				set.c != SIZE_T_FAIL &&
				set.d != SIZE_T_FAIL,
				CodeException, "Atom name mapping failure");
		}
	}

	void RotamerSet::calcTorsionalScope_FF( const MoleculeDefinition& mol, std::vector<IndexHexet>& _scopeMap ) const
	{			
		const std::vector<AtomParameter>& atoms = mol.atom;
		_scopeMap.clear();	
		for( size_t i = 0; i < m_Chis.size(); i++ )
		{
			IndexQuartet ffIndex = calcFFQuartet( m_Chis[i], mol, m_AltNames );
			IndexHexet set; // fill this bitch

			size_t seekPos = 1;
			std::vector<int> ffIndexes;
			ffIndexes.push_back(ffIndex.i);
			ffIndexes.push_back(ffIndex.a);
			ffIndexes.push_back(ffIndex.b);
			ffIndexes.push_back(ffIndex.j);

			// Search for anchors - we cant pre-store this easily, as the N C and intra-polymer versions have different internal indexes
			// (because the N terminal one has H1, H2 and H3.
			const ConventionDef& cDef = m_Conv->getConvention(m_Conv->findConvention(name));
			int cartesianAnchorFF1 = -1, cartesianAnchorFF2 = -1, cartesianAnchorFF3 = -1;
			for( size_t i = 0; i < atoms.size(); i++ )
			{
				if( 0 == cDef.a1.compare( m_AltNames ? atoms[i].rawname : atoms[i].pdbname ) )
				{
					if( cartesianAnchorFF1 != -1 ) 
						throw ProcedureException("Internal atom array contains duplicates for FF anchor-atom 1: '" + cDef.a1 + "'");
					cartesianAnchorFF1 = i;
				}
				else if( 0 == cDef.a2.compare( m_AltNames ? atoms[i].rawname : atoms[i].pdbname ) )
				{
					if( cartesianAnchorFF2 != -1 ) 
						throw ProcedureException("Internal atom array contains duplicates for FF anchor-atom 2: '" + cDef.a2 + "'");
					cartesianAnchorFF2 = i;
				}
				else if( 0 == cDef.a3.compare( m_AltNames ? atoms[i].rawname : atoms[i].pdbname ) )
				{
					if( cartesianAnchorFF3 != -1 ) 
						throw ProcedureException("Internal atom array contains duplicates for FF anchor-atom 3: '" + cDef.a3 + "'");
					cartesianAnchorFF3 = i;
				}
			}

			while(++seekPos < ffIndexes.size())
			{
				const std::vector<CovalentLink>& cov = atoms[ffIndexes[seekPos]].r_covalent;

				bool bomb = false;
				for( size_t k = 0; k < cov.size(); k++ )
				{		
					int comp = cov[k].i;
					if( comp == ffIndex.i || (comp == ffIndex.a && ffIndexes[seekPos] != ffIndex.b) )
					{
						// We just touched our own tail :-D
						// OMG!! We are in a cyclic-structure like proline!!
						// Add the remaining atoms, then exit the loop.
						bomb = true;
					}
					else if( !VectorContains(ffIndexes,comp) && // Test for back-links
						comp != cartesianAnchorFF1 && 
						comp != cartesianAnchorFF2 && 
						comp != cartesianAnchorFF3
						)
					{
						int index = cov[k].i;
						if( index < 0 ) 
							continue;
						ffIndexes.push_back(index);
					}
				}
				if( bomb )
					break;
			}

			// Test that our side-chain atom-order assumption holds true.
			// i.e. that the side-chain is a continuous stretch of atoms with increasing atom index.
			ffIndexes.erase(ffIndexes.begin(),ffIndexes.begin()+3); // remove the 1st three elements
			std::sort(ffIndexes.begin(),ffIndexes.end()); // there may be mild order anomalies, so sort...
			for( int k = 1; k < ffIndexes.size(); k++ )
			{
				ASSERT( ffIndexes[k-1] == ffIndexes[k]-1, CodeException, "Assumption failure");
			}

			set.i = (size_t)ffIndexes[0]; // start range
			set.j = (size_t)ffIndexes[ffIndexes.size()-1]; // end range
			size_t testDeltaFF = set.j - set.i;
			ASSERT( set.i >=0 && set.j >=0, CodeException, "Assumption failure");

			// The torsion atoms
			set.a = ffIndex.i;
			set.b = ffIndex.a;
			set.c = ffIndex.b;
			set.d = ffIndex.j;

			_scopeMap.push_back(set);
		}
	}

	void RotamerSet::applyRotamerTorsional( MoleculeBase& _molBase, int ir, size_t _rotIndex ) const
	{
		std::vector<IndexHexet> scopeList;
		calcTorsionalScope_FF( *_molBase.res[ir].param, scopeList );
		applyRotamerTorsional( _molBase, ir, _rotIndex, scopeList );
	}

	void RotamerSet::applyRotamerTorsional( MoleculeBase& _molBase, int ir, size_t _rotIndex, const std::vector<IndexHexet>& _scopeMap ) const
	{
		D_ASSERT(_molBase.getSequence().getResName(ir).compare( name ) == 0, 
			CodeException, 
			"Refered residue does not correspond in type to this RotamerSet");

		const Residue& res = _molBase.res[ir];
		ParticleStore& atoms = _molBase.atom;
		size_t iFirst = (size_t)res.ifirst;

		// Apply torsions using _scopeMap ... Simple hu ;-)
		Maths::matrix3x3 rmat;
		Maths::dvector axis;
		Maths::dvector changePos;
		Maths::dvector origin;
		for( size_t i = 0; i < m_Chis.size(); i++ )
		{
			// Local copy (easier and computationally faster, verified by benchmark)
			IndexHexet scopeMap = _scopeMap[i];

			// Add offset for the current residue
			scopeMap.a += iFirst;
			scopeMap.b += iFirst;
			scopeMap.c += iFirst;
			scopeMap.d += iFirst;
			scopeMap.i += iFirst;
			scopeMap.j += iFirst;

			origin.setTo( atoms[scopeMap.b].pos() );

			double rotateCurrentAngle = 
				calcTorsionAngle( 
				atoms[scopeMap.a].pos(),
				atoms[scopeMap.b].pos(),
				atoms[scopeMap.c].pos(),
				atoms[scopeMap.d].pos()
				);
			rotateCurrentAngle = m_RotRes[_rotIndex].m_Chis[i] - rotateCurrentAngle; // and now ... what do we need to rotate by?

			// assign the rotation matrix
			axis.setTo(atoms[scopeMap.c].pos());
			axis.sub(origin); // origin is atoms[scopeMap.b].pos();
			rmat.setToAxisRot( axis, rotateCurrentAngle );

			// rotate all the atoms that are affetced by the rotation within the Scope of the loop definition
			for( int n = scopeMap.i; n <= scopeMap.j; n++)
			{
				changePos.setTo(atoms[n].pos());
				changePos.sub(origin); // translate to axis origin
				changePos.mulmat(rmat); // rotate position
				changePos.add(origin); // translate back to original position
				atoms[n].pos().setTo(changePos);
			}
		}
	}

	void RotamerSet::applyFFIdealised( MoleculeBase& _molBase, int ir, const std::vector<size_t>& _indexMap ) const
	{
		const Residue& res = _molBase.res[ir];

		Maths::dvector r1 = _molBase.atom[_indexMap[m_CartesianAnchorRot1]].posGeom();
		Maths::dvector r2 = _molBase.atom[_indexMap[m_CartesianAnchorRot2]].posGeom();
		Maths::dvector r3 = _molBase.atom[_indexMap[m_CartesianAnchorRot3]].posGeom();
		Maths::dvector a1 = _molBase.atom[_indexMap[m_CartesianAnchorRot1]].pos();
		Maths::dvector a2 = _molBase.atom[_indexMap[m_CartesianAnchorRot2]].pos();
		Maths::dvector a3 = _molBase.atom[_indexMap[m_CartesianAnchorRot3]].pos();

		// superimpose the anchor sets onto each other
		Maths::matrix3x3 rmat;
		Maths::superimpose(
			a1, 
			a2, // The anchor atoms within the residue
			a3,
			r1, 
			r2, // The anchor atoms within the FFParam geometric definition of the residue
			r3,
			rmat);

		for( size_t i = 0; i < _indexMap.size(); i++ )
		{
			// Perform the transform on the atomic coordinate
			Maths::dvector& pos = _molBase.atom[_indexMap[i]].pos();
			pos.setTo(_molBase.atom[_indexMap[i]].posGeom());
			pos.sub(r1);
			pos.mulmat(rmat);
			pos.add(a1);
		}
	}

	void RotamerSet::placeOrigin(Rotamer& _r) const
	{
		Maths::dvector offset = _r.m_Pos[m_CartesianAnchorRot1];
		if( offset.x == 0.0 && offset.y == 0.0 && offset.z == 0.0 ) 
			return;
		for( size_t i = 0; i < _r.m_Pos.size(); i++ )
		{
			_r.m_Pos[i].sub(offset);
		}
	}

	void RotamerSet::calcTorsions(Rotamer& r) const
	{
		// Proxies
		std::vector<double>& chi = r.m_Chis;
		std::vector<Maths::dvector>& atom = r.m_Pos;

		chi.clear();

		for( size_t i = 0; i < m_Chis.size(); i++ )
		{
			D_ASSERT( m_Chis[i].rot.i != -1 && m_Chis[i].rot.a != -1 && m_Chis[i].rot.b != -1 && m_Chis[i].rot.j != -1,
				CodeException, "Mapping error, chi definition should not have '-1' atom indexes");
			chi.push_back( Maths::calcTorsionAngle(
				atom[m_Chis[i].rot.i],
				atom[m_Chis[i].rot.a],
				atom[m_Chis[i].rot.b],
				atom[m_Chis[i].rot.j]) );
		}
	}

	void RotamerSet::assignTorsionDefs()
	{
		ASSERT(m_Conv!=NULL,CodeException,"RotamerSet::assignTorsionDefs() m_Conv is null");
		int id = m_Conv->findConvention(name);
		if( id == -1 ) throw ProcedureException("Internal rotamer conventions do not define the residue type: '" + name + "'");
		const ConventionDef& cDef = m_Conv->getConvention(id);
		m_Chis.clear();		
		for( size_t i = 0; i < cDef.chi.size(); i++ )
		{
			m_Chis.push_back(RotamerChi(cDef.chi[i],m_Atom));
		}
	}

	void RotamerSet::assignAnchorIndexes()
	{
		ASSERT(m_Conv!=NULL,CodeException,"RotamerSet::assignAnchorIndexes() m_Conv is null");
		int id = m_Conv->findConvention(name);
		if( id == -1 ) throw ProcedureException("Internal rotamer conventions do not define the residue type: '" + name + "'");
		const ConventionDef& cDef = m_Conv->getConvention(id);

		for( size_t i = 0; i < m_Atom.size(); i++ )
		{
			if( 0 == m_Atom[i].Name.compare(cDef.a1) )
			{
				if( m_CartesianAnchorRot1 != -1 ) 
					throw ProcedureException("Internal atom array contains duplicates for anchor-atom 1: '" + cDef.a1 + "'");
				m_CartesianAnchorRot1 = i;
			}
			else if( 0 == m_Atom[i].Name.compare(cDef.a2) )
			{
				if( m_CartesianAnchorRot2 != -1 ) 
					throw ProcedureException("Internal atom array contains duplicates for anchor-atom 2: '" + cDef.a2 + "'");
				m_CartesianAnchorRot2 = i;
			}
			else if( 0 == m_Atom[i].Name.compare(cDef.a3) )
			{
				if( m_CartesianAnchorRot3 != -1 ) 
					throw ProcedureException("Internal atom array contains duplicates for anchor-atom 3: '" + cDef.a3 + "'");
				m_CartesianAnchorRot3 = i;
			}
		}

		if( m_CartesianAnchorRot1 == -1 ) 
			throw ProcedureException("Internal atom array does not define anchor-atom 1: '" + cDef.a1 + "'");
		if( m_CartesianAnchorRot2 == -1 ) 
			throw ProcedureException("Internal atom array does not define anchor-atom 2: '" + cDef.a2 + "'");
		if( m_CartesianAnchorRot3 == -1 ) 
			throw ProcedureException("Internal atom array does not define anchor-atom 3: '" + cDef.a3 + "'");
	}

	int RotamerSet::searchAtomName( const std::string& atomName ) const
	{
		for( size_t i = 0; i < m_Atom.size(); i++ )
		{
			if( 0 == m_Atom[i].Name.compare(atomName) )
			{
				return i;
			}
		}
		return -1;
	}

	void RotamerSet::add( const std::vector<double>& _Chis, const ProbabilityBase& _probMap )
	{
		// Assign the internal position array (checking that the atom mappings are correct)
		ASSERT(_Chis.size() == m_Chis.size(), CodeException, "RotamerSet::add() internal chi count mismatch, the RotamerLibrary should have assured this was true.");

		Rotamer r;
		r.m_ProbMap = _probMap;
		r.m_SourceRepresentation = WriteTorsional;
		r.m_Chis = _Chis;

		// Calculate the cartesian coordinares using the internal '_conv' chi conventions
		calcCartesian(r);

		// For optimisation of the apply() routine later, anchor-atom-1 **must** be placed on the origin!!
		placeOrigin(r);

		// Add the completed rotamer :-D
		m_RotRes.push_back(r);
	}

	void RotamerSet::add( const std::vector<RotamerAtom>& _atom, const ProbabilityBase& _probMap )
	{
		// Assign the internal position array (checking that the atom mappings are correct)
		ASSERT(_atom.size() == m_Atom.size(), CodeException, "RotamerSet::add() internal atom count mismatch, the RotamerLibrary should have assured this was true.");

		Rotamer r;
		r.m_ProbMap = _probMap;
		r.m_SourceRepresentation = WriteCartesian;

		for( size_t i = 0; i < _atom.size(); i++ )
		{
			D_ASSERT(0 == _atom[i].Name.compare(m_Atom[i].Name), CodeException, "RotamerSet::add(): Internal atom count mismatch, the RotamerLibrary should have assured this was true.");
			D_ASSERT(_atom[i].Defined && _atom[i].IdealPos.isReal(), CodeException, "RotamerSet::add(): RotamerLibrary should have ensured that all atoms have sensible defined coordinates!");
			r.m_Pos.push_back( _atom[i].IdealPos );			
		}

		// For optimisation of the apply() routine later, anchor-atom-1 **must** be placed on the origin!!
		placeOrigin(r);

		// Calculate the torsions using the internal '_conv' atom conventions
		calcTorsions(r);

		// Add the completed rotamer :-D
		m_RotRes.push_back(r);
	}

	inline bool findAtomNames( const MoleculeDefinition& mol, int ffBuildIndex, int rotBuildIndex, 
		bool altNames, const std::vector<RotamerAtom>& rotamerAtoms, 
		int& ffpsIndex1, int& ffpsIndex2, int& ffpsIndex3, int& rotIndex1, int& rotIndex2, int& rotIndex3)
	{
		// Proxies
		const std::vector<AtomParameter> atoms = mol.atom;

		size_t seekPos = 0;
		std::vector<int> ffIndexes;
		std::vector<int> rotIndexes;
		std::vector<int> levels; // Covalent levels - the atom to be rebuilt is level 0

		// Set the initial covalentLink pointer to that of the atom we are building
		const std::vector<CovalentLink>* p_cov = &atoms[ffBuildIndex].r_covalent;
		levels.push_back(0);
		ffIndexes.push_back(mol.findAtomRaw(rotamerAtoms[rotBuildIndex].Name));
		ASSERT(ffIndexes[0] != -1, ProcedureException, "RotamerLibrary::positionUndefinedAtoms(): Error in FFPS atom lookup");

		size_t rotIndex = FindFirstInVector( rotamerAtoms, 
			altNames ? atoms[ffBuildIndex].rawname : atoms[ffBuildIndex].pdbname // ani is an altName, it may not match the rotamer names
			);
		ASSERT( rotIndex != SIZE_T_FAIL, CodeException, "" );
		rotIndexes.push_back( rotIndex );

		while(levels[seekPos] < 3)
		{
			const std::vector<CovalentLink>& cov = *p_cov;
			for( size_t i = 0; i < cov.size(); i++ )
			{					
				if( !VectorContains(ffIndexes,cov[i].i) ) // test for back-links
				{
					// Lookup atom concerned in our rotamer atoms
					// ani is an altName, it may not match the rotamer names if it is using PDB convention
					// therefore use ensureNaming to do the interconversion using FFPS.
					size_t rotIndex = FindFirstInVector( rotamerAtoms, 
						ensureNaming(mol, altNames, cov[i].ani) 
						);
					if( rotIndex == SIZE_T_FAIL )
					{
						// Cant find it in the rotamer
						continue;
					}
					if( rotamerAtoms[rotIndex].Defined )
					{
						ffIndexes.push_back(cov[i].i);
						rotIndexes.push_back(rotIndex);
						levels.push_back(levels[seekPos]+1);

						if( ffIndexes.size() == 4 )
						{
							ffpsIndex1 = ffIndexes[1];
							ffpsIndex2 = ffIndexes[2];
							ffpsIndex3 = ffIndexes[3];
							rotIndex1 = rotIndexes[1];
							rotIndex2 = rotIndexes[2];
							rotIndex3 = rotIndexes[3];
							return true;
						}
					}
				}
			}

			seekPos++;
			if( seekPos == ffIndexes.size() )
			{
				break;
			}

			// Increment
			p_cov = &atoms[ffIndexes[seekPos]].r_covalent;
		}

		// Whoops, we have run out of atoms. Fail.
		return false;
	}

	void RotamerLibrary::positionUndefinedAtoms( std::vector<RotamerAtom>& rotamerAtoms ) const
	{
		// Rebuild any missing atoms using standard geometry, provided by FFParam.
		// **NOTE** - undefined atoms are garunteed to be hydrogens, the check was 
		// performed earlier; although this routine can rebuild
		// any atom from the FFParam definition...
		const FFParamSet& ffps = *m_ffps;

		// The rotation matrix we will use for superimposition
		Maths::matrix3x3 rmat;

		int prevUndefined = -1;
		while(true)
		{
			int undefined = 0;
			for( size_t i = 0; i < rotamerAtoms.size(); i++ )
			{
				if( rotamerAtoms[i].Defined ) 
					continue;

				RotamerAtom& buildThisAtom = rotamerAtoms[i];
				Maths::dvector& buildPos = buildThisAtom.IdealPos;

				const MoleculeDefinition& mol = ffps.molecule[m_CurrentFFPSLink];
				int ffpsTemplateAtomIndex = mol.findAtomRaw(buildThisAtom.Name);
				ASSERT(ffpsTemplateAtomIndex != -1, ProcedureException, "RotamerLibrary::positionUndefinedAtoms(): Error in FFPS atom lookup");

				// We need to find 3 other covalently bound atoms
				// Each will act as a reference state for the orientation of the missing atom
				int ffpsIndex1, ffpsIndex2, ffpsIndex3, rotIndex1, rotIndex2, rotIndex3;
				bool findResult = findAtomNames( 
					mol, 
					ffpsTemplateAtomIndex, 
					i,
					m_AltNames, 
					rotamerAtoms, 
					ffpsIndex1, 
					ffpsIndex2, 
					ffpsIndex3,
					rotIndex1, 
					rotIndex2, 
					rotIndex3 
					);
				if( !findResult )
				{
					undefined++;
					continue; // we may have multi-level rebuilds, cycle, we will build later
				}

				superimpose(
					rotamerAtoms[rotIndex1].IdealPos,
					rotamerAtoms[rotIndex2].IdealPos,
					rotamerAtoms[rotIndex3].IdealPos,
					mol.atom[ffpsIndex1].posGeom(),
					mol.atom[ffpsIndex2].posGeom(),
					mol.atom[ffpsIndex3].posGeom(),
					rmat );

				buildPos.setTo(mol.atom[ffpsTemplateAtomIndex].posGeom());
				buildPos.sub(mol.atom[ffpsIndex1].posGeom());
				buildPos.mulmat(rmat);
				buildPos.add(rotamerAtoms[rotIndex1].IdealPos);

				rotamerAtoms[i].Defined = true;
			}
			if( undefined == 0 )
			{
				break;
			}
			else if( undefined == prevUndefined)
			{
				THROW( ProcedureException, "Cannot find suitable reference atoms for Rotamer rebuild process" );
			}
			else
			{
				prevUndefined = undefined;
				continue;
			}
		}
	}

	RotamerLibrary::RotamerLibrary(const FFParamSet& _ffps) : m_ffps(&_ffps)
	{
		commonInit();
	}

	RotamerLibrary::RotamerLibrary(const FFParamSet& _ffps, const std::string _FileName ) : m_ffps(&_ffps)
	{
		commonInit();
		readLib( _FileName );
	}

	void RotamerLibrary::commonInit()
	{
		// Setup internal conventions :-D
		m_DefaultConvention.setToDefaultConvention();
		setDefaultConventions(); // sets m_UserOverrideConventions = false

		// Internal-workings flags
		m_AutoFinalise = true;
		m_Finalised = false;
		RebuildTollerance = Hydrogens; // hydrogens do not need to be defined by cartesian rotamers, they will be rebuild using ff definitions!
		m_AltNames = false;

		clear();
	}

	void RotamerLibrary::overrideConventions( const RotConvention& _ConvPermanentReference )
	{
		m_UserOverrideConventions = true;
		m_Convention = &_ConvPermanentReference;
	}

	void RotamerLibrary::setDefaultConventions()
	{
		m_UserOverrideConventions = false;
		m_Convention = &m_DefaultConvention;
	}

	void RotamerLibrary::readLib( const std::string& _FileName )
	{
		// Meerly trigger a "conversion" using the standard internal PD format converter
		convertLib( _FileName, m_PDConverter );
	}

	void RotamerLibrary::writeLib( const std::string& _FileName, Library::WriteLibMode _mode )
	{
		if( m_AutoFinalise ) 
			finalise(); // ensure the lib is completed.
		else
			throw ProcedureException("RotamerLibrary::writeLib required library finalisation. AutoFinalise has been disabled and the library is not finalised!");
		// For the moment, I can only see a point in writing PD format libraries...
		// Therefore always use the PD converter to format the output file.
		m_PDConverter.writeLib( _FileName, _mode, *this ); 
	}

	void RotamerLibrary::convertLib( const std::string& _FileName, RotLibConvertBase& _Converter )
	{
		clear(); // Reinitialise, we dont import multiple libraries - that doesnt really make sense.
		_Converter.readLib( _FileName, *this );
		if( m_AutoFinalise ) 
			finalise(); // ensure the lib is completed.
	}

	void RotamerLibrary::writePDB( const std::string& _FileName ) const
	{
		THROW(NotImplementedException,"RotamerLibrary::writePDB is not currently available");
	}

	void RotamerLibrary::writePDB( const std::string& _FileName, size_t _residueID ) const
	{
		THROW(NotImplementedException,"RotamerLibrary::writePDB is not currently available");
	}

	void RotamerLibrary::info()
	{
		Printf("RotamerLibrary:\n");
		Printf("Is%sFinalised\n")(m_Finalised ? " " : " Not Yet ");
		Printf("Contains %d residue definitions\n")(m_RotRes.size());
	}

	void RotamerLibrary::finalise()
	{
		if( m_Finalised ) return; 
		// So final setup to do... but we need to "seal" the library after RotamerLibraryApplicators 
		// are created to keep internal cross-references constant! Otherwise badness will break out :-D
		m_Finalised = true;
	}

	void RotamerLibrary::clear()
	{
		// Reinitialise the containers
		m_RotRes.clear();

		// Completely reset m_WorkingIonMaps to only what the user has defined.
		// The RotLibConvertBase can add its own to the working set as desired...
		m_WorkingIonMaps = m_UserIonMaps;

		// This will be re-overridden by an import converter, but needs to be reset, as import converters 
		// dont *have* to override, and the previous converter might have not used the default.
		if( !m_UserOverrideConventions )
		{
			setDefaultConventions();
		}

		// Sort-helper variables
		m_CurrentFFPSLink = -1;
		m_CurrentRotLink = SIZE_T_FAIL;

		// We are no longer finalised
		m_Finalised = false;
	}

	// Used by 'std::sort' in the function 'RotamerLibrary::standardiseAtoms()'
	bool RotamerLibrary::operator()( const RotamerAtom& a, const RotamerAtom& b ) const
	{
		const MoleculeDefinition& mol = m_ffps->molecule[m_CurrentFFPSLink];
		const std::vector<AtomParameter>& atom = mol.atom;

		// Find the index of a and b in the ffps
		size_t ai = SIZE_T_FAIL;
		size_t bi = SIZE_T_FAIL;

		for( size_t i = 0; i < atom.size(); i++ )
		{
			if( 0 == a.Name.compare(m_AltNames ? atom[i].rawname : atom[i].pdbname) )
			{
				ai = i;
				break;
			}
		}
		if(ai == SIZE_T_FAIL)
		{
			return false; // put the crap at the end ...
		}

		for( size_t i = 0; i < atom.size(); i++ )
		{
			if( 0 == b.Name.compare(m_AltNames ? atom[i].rawname : atom[i].pdbname) )
			{
				bi = i;
				break;
			}
		}
		if(bi == SIZE_T_FAIL ) 
		{
			return false; // put the crap at the end ...
		}

		return ai < bi;
	}

	bool RotamerLibrary::shouldRemoveRotamerAtom( const RotamerAtom rot ) const
	{
		const std::vector<RotamerAtom>& _atom = m_RotRes[m_CurrentRotLink].m_Atom;
		bool result = _atom.end() == std::find(_atom.begin(),_atom.end(),rot.Name);
		return result;
	}

	std::vector<RotamerAtom> RotamerLibrary::standardiseAtoms( const std::vector<RotamerAtom>& _Pos ) const
	{
		const RotamerSet& theSet = m_RotRes[m_CurrentRotLink];
		const std::vector<RotamerAtom>& _atom = theSet.m_Atom;

		// Setup
		std::vector<RotamerAtom> pos = _Pos; // Copy const array		
		// Remove any atoms that are present in the pos array, but not in the internal container
		// (Things like 'HA' will have been removed from the parent)
		pos.erase(std::remove_if(pos.begin(),pos.end(), // erase elements in the vector using std::erase and remove_if
			std::bind1st(std::mem_fun(&RotamerLibrary::shouldRemoveRotamerAtom),this)),pos.end()); // bind the member functionshouldRemoveRotamerAtom() to decide what to remove!
		std::sort(pos.begin(),pos.end(),(*this)); // Makes the atom search vs ffps more efficient...

		if( RebuildTollerance == NoRebuild )
		{
			// Assert that the atom names match!
			if( _atom.size() > pos.size() ) 
				throw ProcedureException( "Not enough atoms are present in the rotamer library to fulfil ffps definition! Residue: " + theSet.name );
			for( int i = 0; i < _atom.size(); i++ )
			{
				if( 0 != _atom[i].Name.compare(pos[i].Name) )
					throw ProcedureException( "Atom name mismatch! Residue: " + theSet.name + " FFParam: " + _atom[i].Name + "Rotamer: " + pos[i].Name);
			}
		}
		else if ( RebuildTollerance == Hydrogens )
		{
			RotamerAtom rebuildAtom;
			rebuildAtom.IdealPos.x = DBL_MAX; // initialise to an arbritrary position
			rebuildAtom.IdealPos.y = DBL_MAX;
			rebuildAtom.IdealPos.z = DBL_MAX;
			rebuildAtom.Defined = false; // we need to rebuild later

			for( int i = 0; i < _atom.size(); i++ )
			{
				if( i >= pos.size() || 0 != _atom[i].Name.compare(pos[i].Name) )
				{
					if( _atom[i].IsHydrogen )
					{
						rebuildAtom.Name = _atom[i].Name; // assign the name						
						pos.insert(pos.begin()+i,RotamerAtom(rebuildAtom));
					}
					else
					{
						throw ProcedureException( "Non-hydrogen atom definition missing! Residue: " + theSet.name + " FFParam: '" + _atom[i].Name + "'");
					}
				}
			}
		}
		else
		{
			THROW(CodeException,"Unknown 'RebuildTollerance' encountered in RotamerLibrary::assertAtoms");
		}

		// Ensure that all atoms have sensible defined positions using the ffps!
		positionUndefinedAtoms(pos);

		return pos;
	}

	inline void addRotamerAssertFinalised( bool _is )
	{
		if(_is) throw CodeException("RotamerLibrary::addRotamer() cannot be called following finalisation.\nDid you forget to disable AutoFinalise on the library?\nFinalisation occurs following readLib() or convertLib() if AutoFinalise is enabled (default).");
	}

	void RotamerLibrary::addIonisationAlias( const std::string& _rotSource, const std::string& _ionisedAlias )
	{
		addRotamerAssertFinalised(isFinalised());
		StringPair alias(_rotSource,_ionisedAlias);
		if( !VectorContains(m_UserIonMaps, alias ) )
		{
			m_UserIonMaps.push_back(alias);
		}
	}

	void RotamerLibrary::addIonisationAlias( const StringPair& alias )
	{
		addRotamerAssertFinalised(isFinalised());
		if( !VectorContains(m_UserIonMaps, alias )  )
		{
			m_UserIonMaps.push_back(alias);
		}
	}

	void RotamerLibrary::clearIonisationAliases() 
	{ 
		m_UserIonMaps.clear(); 
	}

	void RotamerLibrary::addIonisationAlias( const std::vector<StringPair>& aliases )
	{
		addRotamerAssertFinalised(isFinalised());
		for( size_t i = 0; i < aliases.size(); i++ )
		{
			if( !VectorContains(m_UserIonMaps, aliases[i] )  )
			{
				m_UserIonMaps.push_back(aliases[i]);
			}
		}
	}

	void RotamerLibrary::addIonisationAliasWorkingDefaults()
	{
		// Force the addition of the HID HIE and HIP ionisation state mappers
		StringPair alias = StringPair("HIS","HIE");
		if( !VectorContains(m_WorkingIonMaps, alias ) ) 
			m_WorkingIonMaps.push_back(alias);

		alias = StringPair("HIS","HID");
		if( !VectorContains(m_WorkingIonMaps, alias ) ) 
			m_WorkingIonMaps.push_back(alias);

		alias = StringPair("HIS","HIP");
		if( !VectorContains(m_WorkingIonMaps, alias ) ) 
			m_WorkingIonMaps.push_back(alias);

		alias = StringPair("CYS","CYX");
		if( !VectorContains(m_WorkingIonMaps, alias ) ) 
			m_WorkingIonMaps.push_back(alias);

		alias = StringPair("CYS","CYM");
		if( !VectorContains(m_WorkingIonMaps, alias ) ) 
			m_WorkingIonMaps.push_back(alias);
	}

	void RotamerLibrary::addAsBlankRotamer( const std::string& resName)
	{
		addRotamerAssertFinalised(isFinalised());		
		ensureRotamerSetIsBuilt(resName);
	}

	void RotamerLibrary::addRotamer( const std::string& resName, const std::vector<RotamerAtom>& _Positions, const ProbabilityBase& _probMap )
	{
		// Always add the primary rotamer
		addRotamer_Core( resName, _Positions, _probMap );

		// ----------------------------------------------------------------------------
		// NOW; we need to see if there are any explicit ionisation state duplications
		// for this rotamer. These may or not require a certain amount
		// of hydrogen-rebuilding. This is performed automatically...
		// ----------------------------------------------------------------------------
		for(std::vector<StringPair>::iterator i = m_WorkingIonMaps.begin(); i < m_WorkingIonMaps.end(); i++ )
		{
			if(
				0 == (*i).p.compare(resName) && // find the parent name
				-1 != m_ffps->findMoleculeType((*i).q)) // only if its defined by the forcefield
			{
				// We have a match - duplicate and rename
				addRotamer_Core( (*i).q, _Positions, _probMap );
			}
		}
	}

	void RotamerLibrary::ensureRotamerSetIsBuilt(const std::string& resName)
	{
		// Unfortunatly I cant see a way round the use of this 'm_CurrentFFPSLink' hack... 
		// it's required to transmit extra info to the std::sort predicate within its functional form
		m_CurrentFFPSLink = m_ffps->findMoleculeType(resName); 
		if( -1 == m_CurrentFFPSLink )
		{
			Printf("WARNING: Forcefield definition for residue '%s' does not exist - Rotamer definition has been discarded!!\n")(resName);
			return;
		}
		// Find or add the library definition
		m_CurrentRotLink = getIDForResidue(resName);

		// Check for the existence of a current set!!
		if( m_CurrentRotLink == SIZE_T_FAIL )
		{
			// We need to make a new RotamerSet. This will involve checking we have all the FFPS atoms defined

			// Generate the molecule definition in the internal atom array of the new RotamerSet...
			const MoleculeDefinition& mol = m_ffps->molecule[m_CurrentFFPSLink];
			const std::vector<AtomParameter>& atom = mol.atom;

			RotamerSet rotSet(resName,m_AltNames,*m_Convention,mol);
			for( size_t i = 0; i < atom.size(); i++ )
			{
				std::string name = m_AltNames ? atom[i].rawname : atom[i].pdbname;
				int anchorType = m_Convention->anchorAtomType(resName, name);
				if(anchorType==-1) 
					ProcedureException("RotamerLibrary: Internal convention for anchor atoms does not represent residue '" + resName + "'");
				if( 0 < anchorType || !atom[i].isBackbone() )
				{
					rotSet.m_Atom.push_back( RotamerAtom( name, atom[i].posGeom(), atom[i].isHydrogen() ) );
				}
			}

			// Trigger assignment of its relevent internal conventions
			rotSet.assignAnchorIndexes();
			rotSet.assignTorsionDefs();

			m_RotRes.push_back(rotSet); // It is cooked, add it to the parent container

			// Assign 'resID' for use below
			m_CurrentRotLink = getIDForResidue(resName);
			D_ASSERT(m_CurrentRotLink!=SIZE_T_FAIL,CodeException,"RotamerLibrary::ensureRotamerSetIsBuilt(): Something smells distinctly fishy!");
		}
	}

	void RotamerLibrary::addRotamer_Core( const std::string& resName, const std::vector<RotamerAtom>& _Positions, const ProbabilityBase& _probMap )
	{
		addRotamerAssertFinalised(isFinalised());		

		// There may or may not be an existing rotamer set to append to
		ensureRotamerSetIsBuilt(resName);
		if( m_CurrentFFPSLink == -1 ) 
			return; // the FFPS does not define this rotamer class!

		RotamerSet& set = m_RotRes[m_CurrentRotLink];
		set.add(
			standardiseAtoms(_Positions), // Ensure that all atoms are present :-D
			_probMap);
	}

	void RotamerLibrary::addRotamer( const std::string& resName, const std::vector<double>& _Chis, const ProbabilityBase& _probMap )
	{
		addRotamer_Core( resName, _Chis, _probMap );

		// ----------------------------------------------------------------------------
		// NOW, we need to see if there are any explicit ionisation state duplications
		// for this rotamer.
		// ----------------------------------------------------------------------------
		for(std::vector<StringPair>::iterator i = m_WorkingIonMaps.begin(); i < m_WorkingIonMaps.end(); i++ )
		{
			if( 0 == (*i).p.compare(resName) && // find the parent name
				-1 != m_ffps->findMoleculeType((*i).q)) // only if its defined by the forcefield
			{
				// We have a match - duplicate and rename
				addRotamer_Core( (*i).q, _Chis, _probMap );
			}
		}
	}

	void RotamerLibrary::addRotamer_Core( const std::string& resName, const std::vector<double>& _Chis, const ProbabilityBase& _probMap )
	{
		addRotamerAssertFinalised(isFinalised());		

		// There may or may not be an existing rotamer set to append to
		ensureRotamerSetIsBuilt(resName);

		RotamerSet& set = m_RotRes[m_CurrentRotLink];
		set.add( _Chis, _probMap );
	}

	void RotamerLibrary::addRotamer( const MoleculeBase& _Mol, size_t _IRes, const ProbabilityBase& _probMap )
	{
		addRotamerAssertFinalised(isFinalised());

		// Proxies :-D
		//std::string resName = _Mol.getResName(_IRes);
		//const Residue& res = _Mol.getRes(_IRes);
		//int start = res.ifirst;
		//size_t length = res.size();
		//const MoleculeDefinition& _ffDef = *res.param;
		//size_t nFFAtom = _ffDef.atom.size();

		THROW(NotImplementedException,"");
	}

	size_t RotamerLibrary::getIDForResidue( const std::string& resName ) const
	{
		// Use the std::find algorithm to find our residue name in the main collection. 
		// If it fails, it will return m_RotRes.end(),
		// but we wish to use SIZE_T_FAIL as the fail flag, so do a test...
		for( size_t i = 0; i < m_RotRes.size(); i++ )
		{
			if( m_RotRes[i].operator==(resName) ) 
				return i;
		}
		return SIZE_T_FAIL;
	}

	size_t RotamerLibrary::nRot() const
	{
		return m_RotRes.size();
	}

	size_t RotamerLibrary::nRotIn( size_t _RotamerID ) const
	{
		return m_RotRes[_RotamerID].nRot();
	}

	RotamerSet& RotamerLibrary::getRotamerSet( size_t _id )
	{
		ASSERT( _id < m_RotRes.size(), OutOfRangeException, "RotamerLibrary RotamerSet request is outside bounds.");
		return m_RotRes[_id];
	}

	const RotamerSet& RotamerLibrary::getRotamerSet( size_t _id ) const
	{
		ASSERT( _id < m_RotRes.size(), OutOfRangeException, "RotamerLibrary RotamerSet request is outside bounds.");
		return m_RotRes[_id];
	}
}

