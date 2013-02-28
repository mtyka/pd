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
#include <deque>
#include "tools/vector.h"
#include "library/rotamer_dunbrack.h"

namespace Library
{
	RotLibConvert_Dunbrack_BBInd::RotLibConvert_Dunbrack_BBInd()
		: RotLibConvertBase()
	{
	}

	void RotLibConvert_Dunbrack_BBInd::readDefinition( StringBuilder& sb, RotamerLibrary& _RotLib ) const
	{
		ASSERT( sb.size() >= 67, ParseException, "String too short to parse!!");

		StringBuilder resname(3);
		resname.setTo( sb, 0, 3 );
		sb.TruncateLeftBy(18);

		int pdbCount = sb.parseInt(0,7);
		if( pdbCount == 0 )
			return;
		double probability = sb.parseDouble(7,8);
		sb.TruncateLeftBy(37);
		
		std::vector<double> importChis;
		std::vector<std::string> parts = sb.tokenise();
		ASSERT( (parts.size() % 2 == 0), ParseException, "Unexpected element count");

		for( size_t i = 0; i < parts.size(); i+=2 )
		{
			double chi;
			ASSERT( 0 == str2double( parts[i], chi ), ParseException, "Error parsing chi def");
			importChis.push_back(Maths::DegToRad(chi));
		}

		ASSERT( importChis.size() != 0, ParseException, "No chis found!!");
		_RotLib.addRotamer( resname.toString(), importChis, ConstantProbability( probability ) );
	}

	void RotLibConvert_Dunbrack_BBInd::readLib( const std::string& _LibraryFilename, RotamerLibrary& _RotLib )
	{
		StringBuilder sb;
		std::ifstream* p_torsionFile;
		try
		{
			if( _RotLib.isFinalised() ) 
				throw ProcedureException("readLib() is not allowed, the rotamer library has been finalised, no further import can occur");

			// Mappings for things like HIS -> HIE HID HIP
			_RotLib.addIonisationAliasWorkingDefaults();

			// Make sure Alanine and Glycine are defined non-rotamer residues (the library itself ignores them)
			_RotLib.addAsBlankRotamer("ALA");
			_RotLib.addAsBlankRotamer("GLY");

			p_torsionFile = new std::ifstream(_LibraryFilename.c_str(), std::ifstream::in);
			std::ifstream& torsionFile = *p_torsionFile;
			if( !torsionFile.is_open() ) throw(IOException( "Dunbrack BB-independent torsional definition file not found: '" + _LibraryFilename + "'!" ));
		
			sb << torsionFile;
			ASSERT( sb.size() >= 36 && sb.compare("Backbone-independent rotamer library",0,36,0,false),
				ParseException,
				"The input coordinate file does not appear to be in Dunbrack BB-independent format");

			const char* initLine = "Res Rotamer   n(r1) n(r1234) p(r1234) sig p(r234|r1) sig  chi1 sig      chi2 sig      chi3 sig      chi4  sig";
			int initLineLength = strlen(initLine);

			bool begun = false;
			while( sb << torsionFile )
			{
				if( sb.size() >= initLineLength && 
					sb.compare( initLine, 0, initLineLength, 0 ) )
				{
					begun = true;
					break;
				}
			}
			ASSERT(begun,ParseException,"Unexpected EOF whilst parsing the Dunbrack BB-independent format library");
			sb << torsionFile; // Two blanking lines are present here
			sb << torsionFile; // So eradicate them - mwa ha ha ha!

			while( sb << torsionFile )
			{
				sb.Trim();
				if( sb.size() == 0 ) continue;
				readDefinition( sb, _RotLib );
			}

			// Ensure file-handle cleanup
			p_torsionFile->close();
			delete p_torsionFile;
		}
		catch( ExceptionBase ex )
		{
			// Ensure file-handle cleanup
			p_torsionFile->close();
			delete p_torsionFile;
			throw ex;
		}
	}

	RotLibConvert_Dunbrack_BBDep::RotLibConvert_Dunbrack_BBDep()
		: RotLibConvertBase(), switchImportMode(false)
	{
	}

	struct ParseData
	{
		ParseData();
		bool parse( StringBuilder& sb );

		StringBuilder rotID;
		short phi;
		short psi;	
		float probability;
		std::vector<double> chis;	
	};

	typedef std::vector<ParseData> ContainerType;

	ParseData::ParseData() 
		: rotID(7)
	{
		chis.reserve(4);
	}

	bool ParseData::parse(StringBuilder& sb)
	{
		if( sb.size() < 107 ) 
			return false;
		phi = (short)sb.parseInt(5,4);
		psi = (short)sb.parseInt(10,4);
		ASSERT( 
			180 >= phi &&
			-180 <= phi &&
			180 >= psi &&
			-180 <= psi, 
			ParseException, "Phi/Psi range error");
		rotID.setTo(sb, 24, 7 );
		rotID.removeAll(' ');
		ASSERT(rotID.size() == 4, ParseException, "RotID is bad");
		probability = (float)sb.parseDouble(33, 8);
		float chi;
		chis.clear();
		chi = (float)sb.parseDouble( 42, 7 );
		if( chi != 0.0 ) // As a coder, this line makes me vomit, but unfortunatly their library is written this way. Eugh!
			chis.push_back(chi);
		chi = (float)sb.parseDouble( 50, 7 );
		if( chi != 0.0 ) // Why could they not have a gap, or at least use NULL or something ?!?
			chis.push_back(chi);
		chi = (float)sb.parseDouble( 58, 7 );
		if( chi != 0.0 ) 
			chis.push_back(chi);
		chi = (float)sb.parseDouble( 66, 7 );
		if( chi != 0.0 ) 
			chis.push_back(chi);
		return true;		
	}

	inline bool removeParseData( const ParseData& pd )
	{
		return pd.phi == 180.0 || pd.psi == 180.0;
	}

	inline bool removeLowPropensity( const ParseData& pd )
	{
		return Maths::SigFigEquality(0.0,(double)pd.probability,6);
	}

	void processData( const std::string& resName, ContainerType& data, RotamerLibrary& _RotLib, double tollerance )
	{
		std::vector<double> averagedChis = data[0].chis;
		if( averagedChis[0] < 0.0 ) 
			averagedChis[0] += 360.0;

		for( size_t j = 1; j < data.size(); j++ )
		{	
			while( data[0].chis.size() > data[j].chis.size() )
			{
				// Chi count disagreement: Unfortunatly some, due to being 0.0 are not read in. Just pad it...
				data[j].chis.push_back(0.0);
			}
			for( size_t k = 0; k < data[0].chis.size(); k++ )
			{
				double d1 = data[0].chis[k];				
				double d2 = data[j].chis[k];

				double delta = d1 - d2;
				if( delta < 0.0 ) 
					delta = -delta;
				if( delta > 180.0 ) 
					delta = 360.0 - delta;
				if( delta != 0.0 )
				{
					// We have to allow a +-tollerance wobble, as some are *a little* different depending on the exact phi/psi!
					ASSERT( tollerance > delta, ParseException, "Internal data correlation mismatch (Chi)");
				}
				if( d2 < 0.0 ) d2 += 360.0; // assure all are +ve

				averagedChis[k] += d2;
			}
		}

		for( size_t i = 0; i < averagedChis.size(); i++ )
		{
			averagedChis[i] = averagedChis[i] / (double)data.size();
			averagedChis[i] = torsionalRangeDegrees( averagedChis[i] );
			averagedChis[i] = Maths::DegToRad(averagedChis[i]);
		}		

		ProbabilityByPhiPsiMap map( 36 ); // 10 degree

		for( size_t j = 0; j < data.size(); j++ )
		{
			map.assignFrom( 
				Maths::DegToRad(data[j].phi), 
				Maths::DegToRad(data[j].psi), 
				Maths::DegToRad(data[j].probability) 
				);
		}

		_RotLib.addRotamer( resName, averagedChis, map );

		return;
	}

	inline bool sortPhi( const ParseData& _a, const ParseData& _b )
	{
		return _a.phi < _b.phi;
	}

	inline bool sortPsi( const ParseData& _a, const ParseData& _b )
	{
		return _a.psi < _b.psi;
	}

	inline bool sortRotID( const ParseData& _a, const ParseData& _b )
	{
		return _a.rotID.compareValue(_b.rotID,0,4,0) >= 0;
	}

	inline double ensure( const std::vector<double>& d, size_t i )
	{
		if( d.size() > i ) return d[i];
		return 0.0;
	}

	inline bool sortChi1( const ParseData& _a, const ParseData& _b )
	{
		return ensure(_a.chis,0) < ensure(_b.chis,0);
	}

	inline bool sortChi2( const ParseData& _a, const ParseData& _b )
	{
		return ensure(_a.chis,1) < ensure(_b.chis,1);
	}

	inline bool sortChi3( const ParseData& _a, const ParseData& _b )
	{
		return ensure(_a.chis,2) < ensure(_b.chis,2);
	}

	inline bool sortChi4( const ParseData& _a, const ParseData& _b )
	{
		return ensure(_a.chis,3) < ensure(_b.chis,3);
	}

	bool equalChis( ParseData& _a, ParseData& _b, double tollerance )
	{
		for( size_t i = 0; i < _a.chis.size(); i++ )
		{
			double a = _a.chis[i];
			double b = _b.chis[i];
			double delta = fabs(a - b);
			if( delta > 180.0 )
				delta = 360.0 - delta;
			if( delta >= tollerance )
				return false;
		}

		return true;
	}

	void processResData_ByChis( const std::string& resName, ContainerType& data, RotamerLibrary& _RotLib )
	{
		const double tollerance = 12.0;

		// Ensure full chi arrays
		size_t largestSoFar = data[0].chis.size();
		for( size_t i = 1; i < data.size(); i++ )
		{
			if( data[i].chis.size() > largestSoFar )
			{
				largestSoFar = data[i].chis.size();
				i = 0; // reset
			}
			while( data[i].chis.size() < largestSoFar )
			{
				data[i].chis.push_back(0.0);
			}
		}

		size_t done = 0;
		std::vector<bool> taken(data.size());
		while( data.size() > done )
		{
			ContainerType subdata;
			bool on = false;
			for( int i = (int)data.size()-1; i >= 0; i-- )
			{
				if( taken[i] != true && (!on || equalChis( subdata[0], data[i], tollerance ) ) )
				{
					subdata.push_back( data[i] );
					taken[i] = true;
					on = true;
					done++;
				}
			}
			processData( resName, subdata, _RotLib, tollerance );
			subdata.clear();
		}
		data.clear();
	}

	void processResData_ByRotID( const std::string& resName, ContainerType& data, RotamerLibrary& _RotLib )
	{
		// Sort the data to increase efficiency below
		std::stable_sort( data.begin(), data.end(), sortPsi);
		std::stable_sort( data.begin(), data.end(), sortPhi);		
		std::stable_sort( data.begin(), data.end(), sortRotID);

		ContainerType subdata;
		
		while( true )
		{
			if( data.size() > 0 && (subdata.size() == 0 || subdata[0].rotID.compare( data.back().rotID,0,4,0 )) )
			{
				subdata.push_back( data.back() );
				data.erase(data.end()-1); // rob the end, not the start (ContainerType<> perfomance reasons!)
			}
			else
			{
				ASSERT( data.size() == (36*36), ParseException, "Assumption error!");
				ASSERT( data[0].phi == -180.0 && data[0].psi == -180.0, ParseException, "Wrong phi/psi");
				for( size_t j = 1; j < data.size(); j++ )
				{	
					ASSERT( data[j].phi == (-180.0+(j/36)*10) && data[j].psi == (-180.0+(j%36)*10), ParseException, "Wrong phi/psi");			
				}
				processData( resName, subdata, _RotLib, 40.0 ); // 40 degree tollerance delta from data[0]
				subdata.clear();
				if( data.size() == 0 )
				{
					return;
				}
			}
		}		
	}

	// This function is given an entire block of the same residue. The rotamerIDs will of course be different.
	void processResData( const std::string& resName, ContainerType& data, RotamerLibrary& _RotLib, bool switchImportMode )
	{
		// Either there is something fundamental I don't understand, or ...
		// The dunbrack library has this really annoying thing whereby the 180 and -180 phi/psi torsions,
		// which are by definition the same thing, are both defined! 
		// We only want one, as both represent the same 10 degree bin!
		data.erase(std::remove_if(data.begin(),data.end(),removeParseData),data.end());

		if( switchImportMode )
		{
			processResData_ByRotID(resName,data,_RotLib);
		}
		else
		{
			data.erase(std::remove_if(data.begin(),data.end(),removeLowPropensity),data.end());
			processResData_ByChis(resName,data,_RotLib);
		}
	}

	void RotLibConvert_Dunbrack_BBDep::readLib( const std::string& _LibraryFilename, RotamerLibrary& _RotLib )
	{
		StringBuilder sb;
		std::ifstream* p_torsionFile;
		try
		{
			if( _RotLib.isFinalised() ) 
				throw ProcedureException("readLib() is not allowed, the rotamer library has been finalised, no further import can occur");

			// Mappings for things like HIS -> HIE HID HIP
			_RotLib.addIonisationAliasWorkingDefaults();

			// Make sure Alanine and Glycine are defined non-rotamer residues (the library itself ignores them)
			_RotLib.addAsBlankRotamer("ALA");
			_RotLib.addAsBlankRotamer("GLY");

			p_torsionFile = new std::ifstream(_LibraryFilename.c_str(), std::ifstream::in);
			std::ifstream& torsionFile = *p_torsionFile;
			if( !torsionFile.is_open() ) 
				throw(IOException( "Dunbrack BB-independent torsional definition file not found: '" + _LibraryFilename + "'!" ));
		
			ParseData pd;
			ContainerType data; 
	
			StringBuilder prevResName(3);
			StringBuilder currResName(4);
			StringBuilder cacheLine;

			while( true )
			{
				if( cacheLine.size() > 0 )
				{
					sb.setTo(cacheLine);
					cacheLine.clear();
				}
				else if( !(sb << torsionFile) )
				{
					break;
				}

				sb.Trim();
				if( sb.size() == 0 ) 
					continue; // blank line

				currResName.setTo(sb,0,3);

				if( prevResName.size() == 0 )
				{
					prevResName.setTo( currResName );
					ASSERT( pd.parse( sb ), ParseException, "Malformed line");
					data.push_back( pd );
				}
				else if( prevResName.compare( currResName, 0, 3, 0 ) )
				{
					// Woooo - we found one :-D
					ASSERT( pd.parse( sb ), ParseException, "Malformed line");
					data.push_back( pd );
				}
				else
				{					
					if( data.size() > 0 )
					{
						processResData( prevResName.toString(), data, _RotLib, switchImportMode );
						ASSERT( data.size() == 0, CodeException, "CodeFail");
					}
					prevResName.setTo( currResName );
					cacheLine.setTo( sb ); // we havent actually processed the current line! Store it.
				}
			}
			if( data.size() > 0 )
			{
				processResData( prevResName.toString(), data, _RotLib, switchImportMode );
				ASSERT( data.size() == 0, CodeException, "CodeFail");
			}
			
			// Ensure file-handle cleanup
			p_torsionFile->close();
			delete p_torsionFile;
		}
		catch( ExceptionBase ex )
		{
			// Ensure file-handle cleanup
			p_torsionFile->close();
			delete p_torsionFile;
			throw ex;
		}
	}
}

