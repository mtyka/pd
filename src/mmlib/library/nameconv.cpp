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

#include "tools/stringtool.h"
#include "forcefields/ffparam.h"

#include "nameconv.h"

using namespace std;

namespace Library
{
	// --------------------------------------------------------------------------------
	// Handles naming conversions
	// --------------------------------------------------------------------------------

	const char* NameSet_PDB[NAME_SET_COUNT][NAME_SET_DEPTH] =
	{
		{ "PDB", "The standard PDB naming convention." },
		{ "polypeptide", "NTE", "n", "N-Terminus", "1H", "2H", "3H", "end"},
		{ "polypeptide", "CTE", "c", "C-Terminus", "O", "OXT", "end"},
		{ "polypeptide", "ALA", "A", "Alanine", "H", "HA", "1HB", "2HB", "3HB", "C", "CA", "CB", "N", "O", "end"},
		{ "polypeptide", "ARG", "R", "Arginine", "H", "HA", "1HB", "2HB", "1HG", "2HG", "1HD", "2HD", "HE", "1HH1", "2HH1", "1HH2", "2HH2", "C", "CA", "CB", "CG", "CD", "CZ", "N", "NE", "NH1", "NH2", "O", "end"},
		{ "polypeptide", "ASP", "D", "Aspartic Acid", "H", "HA", "1HB", "2HB", "C", "CA", "CB", "CG", "N", "O", "OD1", "OD2", "end"},
		{ "polypeptide", "ASN", "N", "Aspargine", "H", "HA", "1HB", "2HB", "2HD2", "1HD2", "C", "CA", "CB", "CG", "N", "ND2", "O", "OD1", "end"},
		{ "polypeptide", "CYS", "C", "Cysteine", "H", "HA", "1HB", "2HB", "HG", "C", "CA", "CB", "N", "O", "SG", "end"},
		{ "polypeptide", "GLU", "E", "Glutamic Acid", "H", "HA", "1HB", "2HB", "1HG", "2HG", "C", "CA", "CB", "CG", "CD", "N", "O", "OE1", "OE2", "end"},
		{ "polypeptide", "GLN", "Q", "Glutamine", "H", "HA", "1HB", "2HB", "1HG", "2HG", "2HE2", "1HE2", "C", "CA", "CB", "CG", "CD", "N", "NE2", "O", "OE1", "end"},
		{ "polypeptide", "GLY", "G", "Glycine", "H", "1HA", "2HA", "C", "CA", "N", "O", "end"},
		{ "polypeptide", "HIS", "H", "Histidine", "H", "HA", "1HB", "2HB", "HD1", "HD2", "HE1", "C", "CA", "CB", "CG", "CD2", "CE1", "N", "ND1", "NE2", "O", "end"},
		{ "polypeptide", "ILE", "I", "isoleucine", "H", "HA", "HB", "1HG1", "2HG1", "1HG2", "2HG2", "3HG2", "1HD1", "2HD1", "3HD1", "C", "CA", "CB", "CG1", "CG2", "CD1", "N", "O", "end"},
		{ "polypeptide", "LEU", "L", "Leucine", "H", "HA", "1HB", "2HB", "HG", "1HD1", "2HD1", "3HD1", "1HD2", "2HD2", "3HD2", "C", "CA", "CB", "CG", "CD1", "CD2", "N", "O", "end"},
		{ "polypeptide", "LYS", "K", "Lysine", "H", "HA", "1HB", "2HB", "1HG", "2HG", "1HD", "2HD", "1HE", "2HE", "1HZ", "2HZ", "3HZ", "C", "CA", "CB", "CG", "CD", "CE", "N", "NZ", "O", "end"},
		{ "polypeptide", "MET", "M", "Methionine", "H", "HA", "1HB", "2HB", "1HG", "2HG", "1HE", "2HE", "3HE", "C", "CA", "CB", "CG", "CE", "N", "O", "SD", "end"},
		{ "polypeptide", "PHE", "F", "Phenyl Alanine", "H", "HA", "1HB", "2HB", "HD1", "HD2", "HE1", "HE2", "HZ", "C", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "N", "O", "end"},
		{ "polypeptide", "PRO", "P", "Proline", "H2", "H1", "HA", "1HB", "2HB", "1HG", "2HG", "1HD", "2HD", "C", "CA", "CB", "CG", "CD", "N", "O", "end"},
		{ "polypeptide", "SER", "S", "Serine", "H", "HA", "1HB", "2HB", "HG", "C", "CA", "CB", "N", "O", "OG", "end"},
		{ "polypeptide", "THR", "T", "Threonine", "H", "HA", "HB", "HG1", "1HG2", "2HG2", "3HG2", "C", "CA", "CB", "CG2", "N", "O", "OG1", "end"},
		{ "polypeptide", "TRP", "W", "Tryptophan", "H", "HA", "1HB", "2HB", "HD1", "HE1", "HE3", "HZ2", "HZ3", "HH2", "C", "CA", "CB", "CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2", "N", "NE1", "O", "end"},
		{ "polypeptide", "TYR", "Y", "Tyrosine", "H", "HA", "1HB", "2HB", "HD1", "HD2", "HE1", "HE2", "HH", "C", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "N", "O", "OH", "end"},
		{ "polypeptide", "VAL", "V", "Valine", "H", "HA", "HB", "1HG1", "2HG1", "3HG1", "1HG2", "2HG2", "3HG2", "C", "CA", "CB", "CG1", "CG2", "N", "O", "end"}
	};

	const char* NameSet_XPLOR[NAME_SET_COUNT][NAME_SET_DEPTH] =
	{
		{ "XPLOR", "The standard XPLOR naming convention." },
		{ "polypeptide", "NTE", "n", "N terminus", "HT1", "HT2", "HT3", "end"},
		{ "polypeptide", "CTE", "c", "C terminus", "OT1", "OT2", "end"},
		{ "polypeptide", "ALA", "A", "Alanine", "HN", "HA", "HB1", "HB2", "HB3", "C", "CA", "CB", "N", "O", "end"},
		{ "polypeptide", "ARG", "R", "Arginine", "HN", "HA", "HB2", "HB1", "HG2", "HG1", "HD2", "HD1", "HE", "HH11", "HH12", "HH21", "HH22", "C", "CA", "CB", "CG", "CD", "CZ", "N", "NE", "NH1", "NH2", "O", "end"},
		{ "polypeptide", "ASP", "D", "Aspartic Acid", "HN", "HA", "HB2", "HB1", "C", "CA", "CB", "CG", "N", "O", "OD1", "OD2", "end"},
		{ "polypeptide", "ASN", "N", "Aspargine", "HN", "HA", "HB2", "HB1", "HD21", "HD22", "C", "CA", "CB", "CG", "N", "ND2", "O", "OD1", "end"},
		{ "polypeptide", "CYS", "C", "Cysteine", "HN", "HA", "HB2", "HB1", "HG", "C", "CA", "CB", "N", "O", "SG", "end"},
		{ "polypeptide", "GLU", "E", "Glutamic Acid", "HN", "HA", "HB2", "HB1", "HG2", "HG1", "C", "CA", "CB", "CG", "CD", "N", "O", "OE1", "OE2", "end"},
		{ "polypeptide", "GLN", "Q", "Glutamine", "HN", "HA", "HB2", "HB1", "HG2", "HG1", "HE21", "HE22", "C", "CA", "CB", "CG", "CD", "N", "NE2", "O", "OE1", "end"},
		{ "polypeptide", "GLY", "G", "Glycine", "HN", "HA2", "HA1", "C", "CA", "N", "O", "end"},
		{ "polypeptide", "HIS", "H", "Histidine", "HN", "HA", "HB2", "HB1", "HD1", "HD2", "HE1", "C", "CA", "CB", "CG", "CD2", "CE1", "N", "ND1", "NE2", "O", "end"},
		{ "polypeptide", "ILE", "I", "isoleucine", "HN", "HA", "HB", "HG12", "HG11", "HG21", "HG22", "HG23", "HD11", "HD12", "HD13", "C", "CA", "CB", "CG1", "CG2", "CD1", "N", "O", "end"},
		{ "polypeptide", "LEU", "L", "Leucine", "HN", "HA", "HB2", "HB1", "HG", "HD11", "HD12", "HD13", "HD21", "HD22", "HD23", "C", "CA", "CB", "CG", "CD1", "CD2", "N", "O", "end"},
		{ "polypeptide", "LYS", "K", "Lysine", "HN", "HA", "HB2", "HB1", "HG2", "HG1", "HD2", "HD1", "HE2", "HE1", "HZ1", "HZ2", "HZ3", "C", "CA", "CB", "CG", "CD", "CE", "N", "NZ", "O", "end"},
		{ "polypeptide", "MET", "M", "Methionine", "HN", "HA", "HB2", "HB1", "HG2", "HG1", "HE1", "HE2", "HE3", "C", "CA", "CB", "CG", "CE", "N", "O", "SD", "end"},
		{ "polypeptide", "PHE", "F", "Phenyl Alanine", "HN", "HA", "HB2", "HB1", "HD1", "HD2", "HE1", "HE2", "HZ", "C", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "N", "O", "end"},
		{ "polypeptide", "PRO", "P", "Proline", "HT2", "HT1", "HA", "HB2", "HB1", "HG2", "HG1", "HD2", "HD1", "C", "CA", "CB", "CG", "CD", "N", "O", "end"},
		{ "polypeptide", "SER", "S", "Serine", "HN", "HA", "HB2", "HB1", "HG", "C", "CA", "CB", "N", "O", "OG", "end"},
		{ "polypeptide", "THR", "T", "Threonine", "HN", "HA", "HB", "HG1", "HG21", "HG22", "HG23", "C", "CA", "CB", "CG2", "N", "O", "OG1", "end"},
		{ "polypeptide", "TRP", "W", "Tryptophan", "HN", "HA", "HB2", "HB1", "HD1", "HE1", "HE3", "HZ2", "HZ3", "HH2", "C", "CA", "CB", "CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2", "N", "NE1", "O", "end"},
		{ "polypeptide", "TYR", "Y", "Tyrosine", "HN", "HA", "HB2", "HB1", "HD1", "HD2", "HE1", "HE2", "HH", "C", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "N", "O", "OH", "end"},
		{ "polypeptide", "VAL", "V", "Valine", "HN", "HA", "HB", "HG11", "HG12", "HG13", "HG21", "HG22", "HG23", "C", "CA", "CB", "CG1", "CG2", "N", "O", "end"}
	};

	const char* NameSet_AMBER[NAME_SET_COUNT][NAME_SET_DEPTH] =
	{
		{ "AMBER", "The standard AMBER naming convention." },
		{ "polypeptide", "NTE", "n", "N terminus", "HT1", "HT2", "HT3", "end"},
		{ "polypeptide", "CTE", "c", "C terminus", "OT1", "OT2", "end"},
		{ "polypeptide", "ALA", "A", "Alanine", "H", "HA", "HB1", "HB2", "HB3", "C", "CA", "CB", "N", "O", "end"},
		{ "polypeptide", "ARG", "R", "Arginine", "H", "HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE", "HH11", "HH12", "HH21", "HH22", "C", "CA", "CB", "CG", "CD", "CZ", "N", "NE", "NH1", "NH2", "O", "end"},
		{ "polypeptide", "ASP", "D", "Aspartic Acid", "H", "HA", "HB2", "HB3", "C", "CA", "CB", "CG", "N", "O", "OD1", "OD2", "end"},
		{ "polypeptide", "ASN", "N", "Aspargine", "H", "HA", "HB2", "HB3", "HD21", "HD22", "C", "CA", "CB", "CG", "N", "ND2", "O", "OD1", "end"},
		{ "polypeptide", "CYS", "C", "Cysteine", "H", "HA", "HB2", "HB3", "HG", "C", "CA", "CB", "N", "O", "SG", "end"},
		{ "polypeptide", "GLU", "E", "Glutamic Acid", "H", "HA", "HB2", "HB3", "HG2", "HG3", "C", "CA", "CB", "CG", "CD", "N", "O", "OE1", "OE2", "end"},
		{ "polypeptide", "GLN", "Q", "Glutamine", "H", "HA", "HB2", "HB3", "HG2", "HG3", "HE21", "HE22", "C", "CA", "CB", "CG", "CD", "N", "NE2", "O", "OE1", "end"},
		{ "polypeptide", "GLY", "G", "Glycine", "H", "HA2", "HA3", "C", "CA", "N", "O", "end"},
		{ "polypeptide", "HIS", "H", "Histidine", "H", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "C", "CA", "CB", "CG", "CD2", "CE1", "N", "ND1", "NE2", "O", "end"},
		{ "polypeptide", "ILE", "I", "isoleucine", "H", "HA", "HB", "HG12", "HG13", "HG21", "HG22", "HG23", "HD11", "HD12", "HD13", "C", "CA", "CB", "CG1", "CG2", "CD1", "N", "O", "end"},
		{ "polypeptide", "LEU", "L", "Leucine", "H", "HA", "HB2", "HB3", "HG", "HD11", "HD12", "HD13", "HD21", "HD22", "HD23", "C", "CA", "CB", "CG", "CD1", "CD2", "N", "O", "end"},
		{ "polypeptide", "LYS", "K", "Lysine", "H", "HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE2", "HE3", "HZ1", "HZ2", "HZ3", "C", "CA", "CB", "CG", "CD", "CE", "N", "NZ", "O", "end"},
		{ "polypeptide", "MET", "M", "Methionine", "H", "HA", "HB2", "HB3", "HG2", "HG3", "HE1", "HE2", "HE3", "C", "CA", "CB", "CG", "CE", "N", "O", "SD", "end"},
		{ "polypeptide", "PHE", "F", "Phenyl Alanine", "H", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "HE2", "HZ", "C", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "N", "O", "end"},
		{ "polypeptide", "PRO", "P", "Proline", "HN1", "HN2", "HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "C", "CA", "CB", "CG", "CD", "N", "O", "end"},
		{ "polypeptide", "SER", "S", "Serine", "H", "HA", "HB2", "HB3", "HG", "C", "CA", "CB", "N", "O", "OG", "end"},
		{ "polypeptide", "THR", "T", "Threonine", "H", "HA", "HB", "HG1", "HG21", "HG22", "HG23", "C", "CA", "CB", "CG2", "N", "O", "OG1", "end"},
		{ "polypeptide", "TRP", "W", "Tryptophan", "H", "HA", "HB2", "HB3", "HD1", "HE1", "HE3", "HZ2", "HZ3", "HH2", "C", "CA", "CB", "CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2", "N", "NE1", "O", "end"},
		{ "polypeptide", "TYR", "Y", "Tyrosine", "H", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "HE2", "HH", "C", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "N", "O", "OH", "end"},
		{ "polypeptide", "VAL", "V", "Valine", "H", "HA", "HB", "HG11", "HG12", "HG13", "HG21", "HG22", "HG23", "C", "CA", "CB", "CG1", "CG2", "N", "O", "end"}
	};




	AliasMapper::AliasMapper() 
	{
	}





	AliasMapperCollection::AliasMapperCollection()
	{
	}

	void AliasMapperCollection::addAliasMapper(const Library::AliasMapper& _mapper)
	{
		m_Mappers.push_back(&_mapper);
	}

	bool AliasMapperCollection::lookupLongName( char _SourceID, std::string& _DestID ) const
	{
		for( size_t i = 0; i < m_Mappers.size(); i++ )
		{
			if( m_Mappers[i]->lookupLongName(_SourceID,_DestID) ) 
				return true;
		}
		return false;
	}

	bool AliasMapperCollection::lookupShortName( const std::string& _SourceID, char& _DestID ) const
	{
		for( size_t i = 0; i < m_Mappers.size(); i++ )
		{
			if( m_Mappers[i]->lookupShortName(_SourceID,_DestID) ) 
				return true;
		}
		return false;
	}

	bool AliasMapperCollection::lookupAlias( const std::string& _SourceID, std::string& _DestID ) const
	{
		// Return true if any substitution is made, but keep interogating the list in case we have
		// any recursive aliases to aliases.
		bool foundOne = false;
		for( size_t i = 0; i < m_Mappers.size(); i++ )
		{
			if( m_Mappers[i]->lookupAlias(_SourceID,_DestID) ) 
			{
				foundOne = true;
				i = 0;
				continue;
			}
		}
		return foundOne;
	}



	NamingConventions* NamingConventions::getSingleton()
	{
		static NamingConventions inst;
		return &inst;
	}

	NamingConventions::NamingConventions()
	{
		// Assign the NameSet used as the internal naming convention for MMLib: PDB! <-- 'who hoo', my fave!
		m_InternalNameSet.parseSet( NameSet_PDB );

		// Push the standard namesets into the internal array.
		m_NameSets.push_back(NameSet());
		m_NameSets[0].parseSet( NameSet_PDB );
		m_NameSets.push_back(NameSet());
		m_NameSets[1].parseSet( NameSet_XPLOR );
		m_NameSets.push_back(NameSet());
		m_NameSets[2].parseSet( NameSet_AMBER );
		// More namesets can be added from file if needed.

		// Obtain our class mappings
		ObtainClasses( NameSet_PDB );
	}

	NamingConventions::~NamingConventions()
	{
		// Nothing to clean up...
	}

	void NamingConventions::printNaming()const
	{
		cout << "The internal naming convention is: " << m_InternalNameSet.getName() << endl << endl;
		cout << "Available Namesets:" << endl;
		for( size_t i = 0; i < m_NameSets.size(); i++ )
		{
			cout << m_NameSets[i].getName() << ": " << m_NameSets[i].getDescription() << endl;
		}
		cout << endl;

		cout << "Data print:" << endl << endl;
		for( size_t i = 0; i < m_InternalNameSet.size(); i++ )
		{
			const ResidueNameSet& resNameSet = m_InternalNameSet.findResidueDef(i);
			cout << "Summary for: " << resNameSet.getName() << " '" << resNameSet.getSingleLetter() << "' " << resNameSet.getLongName() << endl;
			cout << "Intnl.\t";
			for( size_t k = 0; k < m_NameSets.size(); k++ )
			{
				cout << m_NameSets[k].getName();
				if( k != m_NameSets.size() - 1 ) cout << "\t";
				else cout << endl;
			}
			for( size_t j = 0; j < resNameSet.size(); j++ )
			{
				cout << resNameSet.findAtomName(j) << '\t';
				for( size_t k = 0; k < m_NameSets.size(); k++ )
				{
					cout << m_NameSets[k].findAtomName(i,j);
					if( k != m_NameSets.size() - 1 ) cout << "\t";
					else cout << endl;
				}
			}
			cout << endl;
		}
	}

	void NamingConventions::extendNameSets( const std::string &_DefinitionsFile )
	{
		// See comments by the extendNameSets() function declaration
		THROW(NotImplementedException,"");
	}

	const NameSet& NamingConventions::findNameSet(const std::string &_lookupSetName)const
	{
		string setName = _lookupSetName;
		makeupper(setName); // ensure that our string comparisons are case-insensitive
		for( size_t i = 0; i < m_NameSets.size(); i++ )
		{
			if( 0 == setName.compare( m_NameSets[i].getName() ) )
			{
				return m_NameSets[i];
			}
		}
		THROW(ProcedureException,"Namset not found!");
	}

	const NameSet& NamingConventions::findNameSet( StandardNameSets _lookupSet )const
	{
		int lookupIndex = (int)_lookupSet;
		return m_NameSets[lookupIndex];
	}

	const std::string &NamingConventions::findAtomName( StandardNameSets _lookupSet, const std::string &_ResidueLookup, const std::string &_AtomLookup )const
	{
		const NameSet &ns = findNameSet( _lookupSet );
		return findAtomName( ns, _ResidueLookup, _AtomLookup );
	}

	const std::string &NamingConventions::findAtomName( const std::string &_lookupSetName, const std::string &_ResidueLookup, const std::string &_AtomLookup )const
	{
		const NameSet &ns = findNameSet( _lookupSetName );
		return findAtomName( ns, _ResidueLookup, _AtomLookup );
	}

	const std::string &NamingConventions::findAtomName( const NameSet &_ns, const std::string &_ResidueLookup, const std::string &_AtomLookup )const
	{
		int resIndex = _ns.findResidueIndex(_ResidueLookup);
		int atomIndex = _ns.findAtomIndex( resIndex, _AtomLookup );
		return m_InternalNameSet.findAtomName( resIndex, atomIndex );
	}

	bool NamingConventions::lookupLongName( char _SourceID, std::string& _DestID ) const
	{
		for( size_t i = 0; i < m_InternalNameSet.size(); i++ )
		{
			const ResidueNameSet & name = m_InternalNameSet.findResidueDef(i);
			if( _SourceID == name.getSingleLetter() )
			{
				_DestID = name.getName();
				return true;
			}
		}
		return false;
	}

	bool NamingConventions::lookupShortName( const std::string& _SourceID, char& _DestID ) const
	{
		for( size_t i = 0; i < m_InternalNameSet.size(); i++ )
		{
			const ResidueNameSet & name = m_InternalNameSet.findResidueDef(i);
			if( _SourceID.compare(name.getName()) == 0 )
			{
				_DestID = name.getSingleLetter();
				return true;
			}
		}
		return false;
	}

	bool NamingConventions::lookupAlias( const std::string& _SourceID, std::string& _DestID ) const
	{
		return false; // NamingConventions does not supply 'residue : residue' aliases
	}

	void NamingConventions::ObtainClasses( const char* internalDefinition[NAME_SET_COUNT][NAME_SET_DEPTH] )
	{
		// Obtain the residue and corresponding residue class mappings
		for( size_t i = 1; i < NAME_SET_COUNT; i++ )
		{
			addClassMember(internalDefinition[i][0],internalDefinition[i][1]);
		}

		setClassMapperTermini("Polypeptide",
			Library::PrePostFix("N",false),
			Library::PrePostFix("C",false));

		setClassMapperTermini("DNA",
			Library::PrePostFix("5",true),
			Library::PrePostFix("3",true));

		setClassMapperTermini("RNA",
			Library::PrePostFix("5",true),
			Library::PrePostFix("3",true));
	}



	const ResidueNameSet& NameSet::findResidueDef( size_t _Index )const
	{
		return m_ResDefs[_Index];
	}

	const ResidueNameSet& NameSet::findResidueDef( const std::string &_ResName )const
	{
		string resName = string(_ResName);
		makeupper(resName);
		for( size_t i = 0; i < m_ResDefs.size(); i++ )
		{
			if( resName.compare( m_ResDefs[i].getName() ) ) return m_ResDefs[i];
		}
		THROW(ProcedureException,"Cannot find a residue definition with that name!");
	}

	const std::string& NameSet::findResidueName( size_t _Index )const
	{
		return m_ResDefs[_Index].getName();
	}

	int NameSet::findResidueIndex( const std::string &_ResName )const
	{
		string resName = string(_ResName);
		makeupper(resName);
		for( size_t i = 0; i < m_ResDefs.size(); i++ )
		{
			if( resName.compare( m_ResDefs[i].getName() ) ) return i;
		}
		THROW(ProcedureException,"Cannot find a residue definition with that name!");
	}

	int NameSet::findAtomIndex( int _ResIndex, const std::string &_AtomName )const
	{
		return m_ResDefs[_ResIndex].findAtomIndex( _AtomName );
	}

	const std::string &NameSet::findAtomName( int _ResIndex, int _AtomIndex )const
	{
		return m_ResDefs[_ResIndex].findAtomName(_AtomIndex);
	}

	void NameSet::parseSet( const char* internalDefinition[NAME_SET_COUNT][NAME_SET_DEPTH] )
	{
		// clear the curent array to support reinitialisation
		m_ResDefs.clear();

		// Obtain the name and description of the NameSet
		m_Name = string(internalDefinition[0][0]);
		m_Description = string(internalDefinition[0][1]);

		// Obtain the residue and corresponding atom names
		for( size_t i = 1; i < NAME_SET_COUNT; i++ )
		{
			m_ResDefs.push_back( ResidueNameSet(
				string(internalDefinition[i][1]),
				internalDefinition[i][2][0],
				string(internalDefinition[i][3]) ) );
			bool endFound = false;
			for( size_t j = 4; j < NAME_SET_DEPTH; j++ )
			{
				if( 0 == strcmp(internalDefinition[i][j], "end") )
				{
					endFound = true;
					break;
				}
				m_ResDefs[i-1].addName(string(internalDefinition[i][j]));
			}
			if( !endFound )
			{
				THROW(ProcedureException,"End statement not found during residue-definition loading!!");
			}
		}
	}

	void NameSet::parseSet( const std::string &filename )
	{
		m_ResDefs.clear(); // clear the curent array to support reinitialisation
		THROW(NotImplementedException,"");
	}

	void NameSet::parseSet( std::ifstream &fileStream )
	{
		m_ResDefs.clear(); // clear the curent array to support reinitialisation
		THROW(NotImplementedException,"");
	}




	ResidueNameSet::ResidueNameSet( const std::string _ShortName, char _SingleLetter, const std::string _LongName )
	{
		m_Name = _ShortName;
		m_SingleLetter = _SingleLetter;
		m_LongName = _LongName;
	}

	ResidueNameSet::~ResidueNameSet()
	{
	}

	const std::string &ResidueNameSet::findAtomName( int index )const
	{
		return m_AtomNames[index];
	}

	int ResidueNameSet::findAtomIndex( const std::string &_AtomName )const
	{
		string atomName = string(_AtomName);
		makeupper(atomName);
		for( size_t i = 0; i < m_AtomNames.size(); i++ )
		{
			if( 0 == atomName.compare( m_AtomNames[i] ) ) return i;
		}
		THROW(ProcedureException,"Cannot find a atom definition with that name!");
	}

	void ResidueNameSet::addName(const std::string &_Name)
	{
		m_AtomNames.push_back(_Name);
	}
}

