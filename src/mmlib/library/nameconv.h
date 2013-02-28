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

#ifndef __MOL_NAMING_H
#define __MOL_NAMING_H

#include <string>
#include <vector>

#include "residues.h"
#include "mapper.h"

class PD_API FFParamSet;

namespace Library
{






//-------------------------------------------------
//
/// \brief  Defines residue alias name 
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
/// 
/// 1) Sometimes we use short names to represent the same residue as a longer name
///
/// 2) In structural filetypes, one often finds modified residues, i.e. Standard 
/// residue types with additional chemical groups added post-translation. Aliases
/// provide the mappings between the names for the modified residues and one of the 
/// 20 Standard Residues.
/// 
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API ResidueAliasDefinition
	{
	public:
		ResidueAliasDefinition() 
		{
		}

		bool operator==( const ResidueAliasDefinition& Compare )
		{
			return ((0==m_Alias.compare(Compare.m_Name)) && (0==m_Name.compare(Compare.m_Name)));
		}

		void set(const std::string &_Alias,const std::string &_Name)
		{
			m_Alias = _Alias;
			m_Name = _Name;
		}

		bool cmpAlias(const std::string &_queryAlias,std::string &_returnAlias) const
		{
			if( m_Alias.compare(_queryAlias)!=0) return false;
			_returnAlias = m_Name;
			return true;
		}

		const std::string& getAlias() const
		{
			return m_Alias;
		}

		const std::string& getName() const
		{
			return m_Name;
		}

	protected:
		std::string m_Alias;
		std::string m_Name;
	};

	// --------------------------------------------------------------------------------
	// Below is 'NamingConventions' and associated classes:
	// Handles knowledge-based naming mappings 
	// 1) From Residue name to single letter residue name
	// 2) From Atom names in one name-set to atom names in another
	// --------------------------------------------------------------------------------

	const int NAME_SET_COUNT = 23;
	const int NAME_SET_DEPTH = 31;

	enum StandardNameSets
	{
		// NOTE!! - These indexes of these names **MUST** corespond to the
		// indexes used in the NamingConventions class PD_API constructor!
		PDB = 0,
		XPLOR = 1,
		AMBER = 2,
	};






//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API ResidueNameSet
	{
		friend class PD_API NameSet;
	public:
		ResidueNameSet( const std::string _ShortName, char _SingleLetter, const std::string _LongName );
		~ResidueNameSet();
		const std::string &findAtomName( int index ) const;
		int findAtomIndex( const std::string &_LookupName ) const;
		inline size_t size() const { return m_AtomNames.size(); }
		inline const std::string &getName() const { return m_Name; }
		inline char getSingleLetter() const { return m_SingleLetter; }
		inline const std::string getLongName() const { return m_LongName; }
	protected:
		// Functions
		void addName(const std::string &name);
		// Data
		std::string m_Name;
		char m_SingleLetter;
		std::string m_LongName;
		std::vector<std::string> m_AtomNames;
	};






//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API NameSet
	{
		friend class PD_API NamingConventions;
	public:
		NameSet(){}
		inline const std::string &getName() const { return m_Name; }
		inline const std::string &getDescription() const { return m_Description; }
		const ResidueNameSet& findResidueDef( const std::string &_ResName ) const;
		const ResidueNameSet& findResidueDef( size_t _Index ) const;
		const std::string& findResidueName( size_t _Index ) const;
		int findResidueIndex( const std::string &_ResName ) const;
		int findAtomIndex( int _ResIndex, const std::string &_AtomName ) const;
		const std::string &findAtomName( int _ResIndex, int _AtomIndex ) const;
		inline size_t size() const { return m_ResDefs.size(); }
	protected:
		// Functions
		void parseSet( const char* internalDefinition[NAME_SET_COUNT][NAME_SET_DEPTH] );
		void parseSet( const std::string &filename );
		void parseSet( std::ifstream &fileStream );
		// Data
		std::string m_Name;
		std::string m_Description;
		std::vector<ResidueNameSet> m_ResDefs;
	};






//-------------------------------------------------
//
/// \brief  Map by naming convention
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API NamingConventions : public AliasMapper, public ClassMapper
	{
	public:
		static NamingConventions* getSingleton(); /// 'Singleton Pattern' Implementation
		NamingConventions();
		virtual ~NamingConventions();

		void printNaming() const;
		inline size_t size() const { return m_NameSets.size(); };

		// extendNameSets() could be utilised in the future to extend the internal naming definitions for
		// non-standard residue types.
		void extendNameSets( const std::string &_DefinitionsFile );
		// Implementation:
		// 1) Each residue 'nameset extension' should define both an internal naming convention AND
		// all the corresponding names in additional defined/named NameSet(s).
		//
		// a) If there is no internal definition for an internal residue name, it should be added
		// to the internal definition; otherwise the definition in the file should be checked
		// so that it corresponds exactly (atom name for atom name) to the current definiton
		// in memory - if it doesnt an error should be thrown.
		// b) Each additional nameset residue definition should be added to the corresponding nameset,
		// each atom must be defined to correspond to the internal definiton by numerical **index**
		// (how else would you tell what goes with what).
		// c) For each additional set where there is no defined naming convension, a single 'null'
		// keyword should beused in the first entry slot. This 'define' indicates that no definition
		// corresponding to the internal set is present.
		// d) Each new definition should be checked against existing definitions to prevent
		// naming conflicts for either long or singleChar names.
		//
		// 2) Exact external file format is still to be defined.
		//
		// 3) This should overall hopefully provide an extendible and versitile resdiue re-naming
		// system for MMLib...

		const NameSet& findNameSet(const std::string &_lookupSetName) const;
		const NameSet& findNameSet(StandardNameSets _lookupSet) const;

		const std::string &findAtomName( StandardNameSets _lookupSet, const std::string &_ResidueLookup, const std::string &_AtomLookup ) const;
		const std::string &findAtomName( const std::string &_lookupSetName, const std::string &_ResidueLookup, const std::string &_AtomLookup ) const;
		const std::string &findAtomName( const NameSet &ns, const std::string &_ResidueLookup, const std::string &_AtomLookup ) const;

		virtual bool lookupLongName( char _SourceID, std::string& _DestID ) const; ///< Look up the long name for a given single letter name
		virtual bool lookupShortName( const std::string& _SourceID, char& _DestID ) const; ///< Look up the single letter for a long residue name
		virtual bool lookupAlias( const std::string& _SourceID, std::string& _DestID ) const; ///< Look to see if the source name is an alias, and return the true name

	protected:
		void ObtainClasses( const char* internalDefinition[NAME_SET_COUNT][NAME_SET_DEPTH] );

		std::vector<NameSet> m_NameSets;
		NameSet m_InternalNameSet;
	};
}

#endif

