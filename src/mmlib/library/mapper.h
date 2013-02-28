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

#ifndef __MAPPER
#define __MAPPER

#include <vector>
#include "tools/enum.h"

namespace Library
{
	// Filter what we actually require from our PDB file. Often, we just want the polypeptide.
	const int ResidueClassCount = 7; // ***IMPORTANT*** - Set this, or break: FileInBase::baseLoad()
	enum ResidueClass
	{
		// Basic Classes
		Polypeptide = 1, // enum MUST start with 1 - the first bit
		DNA = 2,
		RNA = 4,
		SmallMolecule = 8,
		Water = 16,
		Ion = 32,
		Carbohydrate = 64
	};

	bool isPolymerClass( ResidueClass _class );
	std::string getResidueClassString( ResidueClass _Class );

	enum ResidueClassGroup
	{
		// Supersets
		Solvent = Water | Ion,
		Nucleotide = DNA | RNA,
		AllExceptSolvent = Polypeptide | SmallMolecule | Nucleotide | Carbohydrate,
		AllClasses = Polypeptide | SmallMolecule | Solvent | Nucleotide | Ion | Carbohydrate
	};

	// Use enum tool to allow enum inheritance
	typedef InheritEnum< ResidueClass, ResidueClassGroup > ExtdResidueClasses;



	//-------------------------------------------------
	//
	/// \brief  Residue Class Mapping
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
	class PrePostFix
	{
	public:
		PrePostFix();
		PrePostFix(const std::string& _Fix, bool _IsPost );
		static const PrePostFix NoFix;
		std::string build( const std::string& _from ) const;
		bool isDefined() const;
		std::string fix;
		bool isPost;
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
	class ClassMapper; // Pre-declaration
	class ResidueClassDef
	{
		friend class ClassMapper;
	public:
		ResidueClassDef( const std::string& _Name );
		const std::string& getName() const;
		bool contains( const std::string& _Member, bool _IncludeCappingResidues ) const;
		bool containsCap( const std::string& _Member ) const;
	private:
		void addMember( const std::string& _NewMember );
		void addCap( const std::string& _NewMember );
		void AssignTerminalCodes( const PrePostFix& _StartCode, const PrePostFix& _EndCode );
		std::vector<std::string> m_Members;
		std::vector<std::string> m_Caps;
		std::string m_Name;
		PrePostFix m_PolymerStart;
		PrePostFix m_PolymerEnd;
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
	class PD_API ClassMapper
	{
	public:
		ClassMapper();
		bool isOfClass( ExtdResidueClasses _Class, const std::string& _ResName ) const;
		bool isOfClass( const std::string& _Class, const std::string& _ResName ) const;
		bool getClassString( const std::string& _ResName, std::string& _ReturnClassName ) const;
		bool getStartTerminusFromResidue( const std::string& _StartResName, PrePostFix& _Fix ) const;
		bool getEndTerminusFromResidue( const std::string& _EndResName, PrePostFix& _Fix ) const;
		bool getStartTerminusFromClass( const std::string& _StartResName, PrePostFix& _Fix ) const;
		bool getEndTerminusFromClass( const std::string& _EndResName, PrePostFix& _Fix ) const;
		void detectSplitRequirement( std::string& _resName, bool& isStart, bool& isEnd ) const; ///< Deal with NARG as a name --> ARG + PrePostFix("N",false)

		void setClassMapperTermini( const std::string& _Class, const PrePostFix& _StartCode, const PrePostFix& _EndCode );
		void addClassCap( const std::string& _Class, const std::string& _NewMember );
		void addClassMember( const std::string& _Class, const std::string& _NewMember );
	protected:

	private:
		std::vector<ResidueClassDef> m_Classes;
	};









	//-------------------------------------------------
	//
	/// \brief  Residue Alias and Name Mapping
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
	class PD_API AliasMapper
	{
	public:
		AliasMapper();
		virtual ~AliasMapper(){};
		virtual bool lookupLongName( char _SourceID, std::string& _DestID ) const = 0; ///< Look up the long name for a given single letter name
		virtual bool lookupShortName( const std::string& _SourceID, char& _DestID ) const = 0; ///< Look up the single letter for a long residue name
		virtual bool lookupAlias( const std::string& _SourceID, std::string& _DestID ) const = 0; ///< Look to see if the source name is an alias, and return the true name
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
	class PD_API AliasMapperCollection : public AliasMapper
	{
	public:
		AliasMapperCollection();
		virtual ~AliasMapperCollection(){};

		void addAliasMapper(const AliasMapper& _mapper);

		virtual bool lookupLongName( char _SourceID, std::string& _DestID ) const; ///< Look up the long name for a given single letter name
		virtual bool lookupShortName( const std::string& _SourceID, char& _DestID ) const; ///< Look up the single letter for a long residue name
		virtual bool lookupAlias( const std::string& _SourceID, std::string& _DestID ) const; ///< Look to see if the source name is an alias, and return the true name

	private:
		std::vector<const AliasMapper*> m_Mappers;
	};






	//-------------------------------------------------
	//
	/// \brief  Map residue naming
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
	class PD_API CustomAliasMapper : public AliasMapper
	{
	private:
		struct NameSet
		{
			char letter;
			std::string name;
			std::string l3name;
		};
		struct AliasSet
		{
			std::string name;
			std::string alias;
		};
	public:
		CustomAliasMapper();
		virtual ~CustomAliasMapper(){};
		void addNameSet(char letter, const std::string& l3name, const std::string& name );
		void addAlias( const std::string& alias, const std::string& name );
		virtual bool lookupLongName( char _SourceID, std::string& _DestID ) const; ///< Look up the long name for a given single letter name
		virtual bool lookupShortName( const std::string& _SourceID, char& _DestID ) const; ///< Look up the single letter for a long residue name
		virtual bool lookupAlias( const std::string& _SourceID, std::string& _DestID ) const; ///< Look to see if the source name is an alias, and return the true name
	protected:
		std::vector<NameSet> m_Names;
		std::vector<AliasSet> m_Alias;
	};
}

#endif

