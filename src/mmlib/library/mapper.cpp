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
#include "mapper.h"

namespace Library
{
	bool isPolymerClass( ResidueClass _class )
	{
		return
			((_class & Polypeptide ) > 0) ||
			((_class & DNA ) > 0) || 
			((_class & RNA ) > 0) || 
			((_class & Carbohydrate ) > 0); 
	}

	std::string getResidueClassString( ResidueClass _Class )
	{
		StringBuilder name;
		if( 0 < ( _Class & Polypeptide ) ) name.append("polypeptide");
		if( 0 < ( _Class & DNA ) ) { if( name.size() > 0 ) name.append(" | "); name.append("dna"); }
		if( 0 < ( _Class & RNA ) ) { if( name.size() > 0 ) name.append(" | "); name.append("rna"); }
		if( 0 < ( _Class & SmallMolecule ) ) { if( name.size() > 0 ) name.append(" | "); name.append("smallmolecule"); }
		if( 0 < ( _Class & Water ) ) { if( name.size() > 0 ) name.append(" | "); name.append("water"); }
		if( 0 < ( _Class & Carbohydrate ) ) { if( name.size() > 0 ) name.append(" | "); name.append("carbohydrate"); }
		if( 0 < ( _Class & Ion ) ) { if( name.size() > 0 ) name.append(" | "); name.append("ion"); }
		return name.toString();
	}

	// --------------------------------------------------------------------------------
	// Handles residue class mappings
	// --------------------------------------------------------------------------------

	const PrePostFix PrePostFix::NoFix = PrePostFix();

	PrePostFix::PrePostFix() 
		: fix(""), isPost(false)
	{
	}
	
	PrePostFix::PrePostFix(const std::string& _Fix, bool _IsPost )
		: fix(_Fix), isPost(_IsPost)
	{
	}

	std::string PrePostFix::build( const std::string& _from ) const
	{
		if( isPost )
		{
			return _from + fix;
		}
		else
		{
			return fix + _from;
		}
	}

	bool PrePostFix::isDefined() const
	{
		return fix.size() > 0;
	}

	ResidueClassDef::ResidueClassDef( const std::string& _Name )
	{
		m_Name = _Name;
	}

	const std::string& ResidueClassDef::getName() const
	{ 
		return m_Name; 
	}

	bool ResidueClassDef::contains( const std::string& _Member, bool _IncludeCappingResidues ) const
	{
		std::string name = _Member;
		makeupper(name);
		for( int i = 0; i < m_Members.size(); i++ )
		{
			if( 0 == name.compare(m_Members[i]) )
			{
				return true;
			}
		}
		if( _IncludeCappingResidues )
		{
			return containsCap(_Member);
		}
		return false;
	}

	bool ResidueClassDef::containsCap( const std::string& _Member ) const
	{
		for( int i = 0; i < m_Caps.size(); i++ )
		{
			if( 0 == _Member.compare(m_Caps[i]) )
			{
				return true;
			}
		}
		return false;
	}

	void ResidueClassDef::addMember( const std::string& _NewMember ) 
	{ 
		std::string name = _NewMember;
		makeupper(name);
		m_Members.push_back(name); 
	}

	void ResidueClassDef::addCap( const std::string& _NewMember )
	{
		std::string name = _NewMember;
		makeupper(name);
		m_Caps.push_back(name); 
	}

	void ResidueClassDef::AssignTerminalCodes( const PrePostFix& _StartCode, const PrePostFix& _EndCode )
	{
		m_PolymerStart = _StartCode;
		m_PolymerEnd = _EndCode;
	}

	ClassMapper::ClassMapper()
	{		
	}

	bool ClassMapper::isOfClass( ExtdResidueClasses _Class, const std::string& _ResName ) const
	{
		if( _Class == AllClasses ) return true;
		if( 0 < ( _Class & Polypeptide ) && isOfClass("polypeptide",_ResName) ) return true;
		if( 0 < ( _Class & DNA ) && isOfClass("dna",_ResName) ) return true;
		if( 0 < ( _Class & RNA ) && isOfClass("rna",_ResName) ) return true;
		if( 0 < ( _Class & SmallMolecule ) && isOfClass("smallmolecule",_ResName) ) return true;
		if( 0 < ( _Class & Water ) && isOfClass("water",_ResName) ) return true;
		if( 0 < ( _Class & Carbohydrate ) && isOfClass("carbohydrate",_ResName) ) return true;
		if( 0 < ( _Class & Ion ) && isOfClass("ion",_ResName) ) return true;
		return false; 
	}

	bool ClassMapper::isOfClass( const std::string& _Class, const std::string& _ResName ) const
	{
		std::string classname = _Class;
		makelower(classname);
		for( size_t i = 0; i < m_Classes.size(); i++ )
		{
			if( 0 == classname.compare(m_Classes[i].getName()) && 
				m_Classes[i].contains(_ResName,true))
			{
				return true;
			}
		}
		return false; // we have no class definition by that name, so we cannot find it...
	}

	bool ClassMapper::getClassString( const std::string& _ResName, std::string& _ReturnClassName ) const 
	{
		std::string resName = _ResName;
		makeupper(resName);
		for( size_t i = 0; i < m_Classes.size(); i++ )
		{
			if( m_Classes[i].contains(resName,true) )
			{
				_ReturnClassName = m_Classes[i].getName();
				return true;
			}
		}
		return false; // We have no class definition by that name, so we cannot find it...
	}

	bool ClassMapper::getStartTerminusFromClass( const std::string& _Class, PrePostFix& _Fix ) const
	{
		std::string classname = _Class;
		makelower(classname);
		for( size_t i = 0; i < m_Classes.size(); i++ )
		{
			if( 0 == classname.compare(m_Classes[i].getName() ) )
			{
				_Fix = m_Classes[i].m_PolymerStart;
				return true;
			}
		}
		return false; // we have no class definition by that name, so we cannot find it...
	}

	bool ClassMapper::getEndTerminusFromClass( const std::string& _Class, PrePostFix& _Fix ) const
	{
		std::string classname = _Class;
		makelower(classname);
		for( size_t i = 0; i < m_Classes.size(); i++ )
		{
			if( 0 == classname.compare(m_Classes[i].getName() ) )
			{
				_Fix = m_Classes[i].m_PolymerEnd;
				return true;
			}
		}
		return false; // we have no class definition by that name, so we cannot find it...
	}

	bool ClassMapper::getStartTerminusFromResidue( const std::string& _StartResName, PrePostFix& _Fix ) const
	{
		// What class is this residue?
		std::string resName = _StartResName;
		for( size_t i = 0; i < m_Classes.size(); i++ )
		{
			if( m_Classes[i].contains(resName,false) )
			{
				_Fix = m_Classes[i].m_PolymerStart;
				return true;
			}
		}
		for( size_t i = 0; i < m_Classes.size(); i++ )
		{
			if( m_Classes[i].containsCap(resName) )
			{
				_Fix = PrePostFix::NoFix; // capping residues have no prefix
				return true;
			}
		}
		return false;
	}

	void ClassMapper::detectSplitRequirement( std::string& resName, bool& isStart, bool& isEnd ) const
	{
		std::string returnName;

		StringBuilder sbFix;
		const StringBuilder sbRes(resName);
		int starts = 0;
		int ends = 0;

		// go through each class:
		// 1) Test to see wether the pre/post fix of the class matches the name
		// 2) If it does, test to see if the class contains the newly truncated residue name
		// 3) Test that we have only found one match
		// e.g. NARG will match Polypeptide with a *pre*-fix of N

		for( size_t i = 0; i < m_Classes.size(); i++ )
		{
			// Obtain the PrePostFix of this class
			PrePostFix fix = m_Classes[i].m_PolymerStart;
			sbFix.setTo(fix.fix);

			// The fix has to fit within the length of the _resname allowing there to be a resname of length at least one
			if( sbFix.size() > 0 && sbFix.size() < sbRes.size() ) 
			{
				if( fix.isPost && sbFix.compare( sbRes, 0, sbFix.size(), sbRes.size() - 1 - sbFix.size(), true ) )
				{
					std::string testName = sbRes.toString(0,sbRes.size()-sbFix.size());
					if( m_Classes[i].contains(testName,false))
					{
						returnName = testName;
						starts++;
					}
				}
				else if( sbFix.compare( sbRes, 0, sbFix.size(), 0, true ))
				{
					size_t length = sbRes.size()-sbFix.size();
					std::string testName = sbRes.toString(sbRes.size()-length,length);
					if( m_Classes[i].contains(testName,false))
					{
						returnName = testName;
						starts++;
					}
				}
			}

			// NOTE - below is an exact repeat of above, but for the m_PolymerEnd
			// Obtain the PrePostFix of this class
			fix = m_Classes[i].m_PolymerEnd;
			sbFix.setTo(fix.fix);

			// The fix has to fit within the length of the _resname allowing there to be a resname of length at least one
			if( sbFix.size() > 0 && sbFix.size() < sbRes.size() ) 
			{
				if( fix.isPost && sbFix.compare( sbRes, 0, sbFix.size(), sbRes.size() - 1 - sbFix.size(), true ) )
				{
					std::string testName = sbRes.toString(0,sbRes.size()-sbFix.size());
					if( m_Classes[i].contains(testName,false))
					{
						returnName = testName;
						ends++;
					}
				}
				else if( sbFix.compare( sbRes, 0, sbFix.size(), 0, true ) )
				{
					size_t length = sbRes.size()-sbFix.size();
					std::string testName = sbRes.toString(sbRes.size()-length,length);
					if( m_Classes[i].contains(testName,false))
					{
						returnName = testName;
						ends++;
					}
				}
			}
		}

		// If we found more than one match, we have trouble. Ambiguity is an error condition.
		ASSERT( (!(starts>0 && ends>0) && starts<=1 && ends<=1), CodeException, "ClassMapper::detectSplitRequirement() detection failure. Ambiguous detection mapping!!");
		if( starts>0 )
		{
			isStart = true;
			resName = returnName;
			return;
		}
		if( ends > 0 )
		{
			isEnd = true;
			resName = returnName;
			return;
		}
		return;
	}

	bool ClassMapper::getEndTerminusFromResidue( const std::string& _EndResName, PrePostFix& _Fix ) const
	{
		// What class is this residue?
		std::string resName = _EndResName;
		for( size_t i = 0; i < m_Classes.size(); i++ )
		{
			if( m_Classes[i].contains(resName,false) )
			{
				_Fix = m_Classes[i].m_PolymerEnd;
				return true;
			}
		}
		for( size_t i = 0; i < m_Classes.size(); i++ )
		{
			if( m_Classes[i].contains(resName,true) )
			{
				_Fix = PrePostFix::NoFix; // capping residues have no prefix
				return true;
			}
		}
		return false;
	}

	void ClassMapper::setClassMapperTermini( const std::string& _KnownClass, const PrePostFix& _StartCode, const PrePostFix& _EndCode )
	{
		std::string classname = _KnownClass;
		makelower(classname);
		for( int i = 0; i < m_Classes.size(); i++ )
		{
			if( 0 == classname.compare(m_Classes[i].getName()) )
			{
				m_Classes[i].AssignTerminalCodes(_StartCode,_EndCode);
				return;
			}
		}
		// make a new class
		ResidueClassDef def(_KnownClass);
		def.AssignTerminalCodes(_StartCode,_EndCode);
		m_Classes.push_back(def); 
		return;
	}

	void ClassMapper::addClassCap( const std::string& _Class, const std::string& _NewMember )
	{
		std::string classname = _Class;
		makelower(classname);
		for( int i = 0; i < m_Classes.size(); i++ )
		{
			if( 0 == classname.compare(m_Classes[i].getName()) )
			{
				m_Classes[i].addCap(_NewMember);
				return;
			}
		}
		// make a new class
		ResidueClassDef def(_Class);
		def.addCap(_NewMember);
		m_Classes.push_back(def); 			
		return;
	}

	void ClassMapper::addClassMember( const std::string& _Class, const std::string& _NewMember )
	{
		std::string classname = _Class;
		makelower(classname);
		for( int i = 0; i < m_Classes.size(); i++ )
		{
			if( 0 == classname.compare(m_Classes[i].getName()) )
			{
				m_Classes[i].addMember(_NewMember);
				return;
			}
		}
		// make a new class
		ResidueClassDef def(classname);
		def.addMember(_NewMember);
		m_Classes.push_back(def); 			
		return;
	}

	CustomAliasMapper::CustomAliasMapper()
	{
	}

	void CustomAliasMapper::addNameSet(char letter, const std::string& l3name, const std::string& name )
	{
		if( l3name.size() != 3 ) THROW( ArgumentException, "'l3name' must be 3 in length!");
		NameSet ns;
		ns.letter = letter;
		ns.name = name;
		ns.l3name = l3name;
		m_Names.push_back( ns );
	}

	void CustomAliasMapper::addAlias( const std::string& alias, const std::string& name )
	{
		AliasSet as;
		as.name = name;
		as.alias = alias;
		m_Alias.push_back( as );
	}

	bool CustomAliasMapper::lookupLongName( char _SourceID, std::string& _DestID ) const
	{
		for( size_t i = 0; i < m_Names.size(); i++ )
		{
			if( _SourceID == m_Names[i].letter )
			{
				_DestID = m_Names[i].name;
				return true;
			}
		}
		return false;
	}

	bool CustomAliasMapper::lookupShortName( const std::string& _SourceID, char& _DestID ) const
	{
		for( size_t i = 0; i < m_Names.size(); i++ )
		{
			if( 0 == _SourceID.compare(m_Names[i].name) )
			{
				_DestID = m_Names[i].letter;
				return true;
			}
		}
		return false;
	}

	bool CustomAliasMapper::lookupAlias( const std::string& _SourceID, std::string& _DestID ) const
	{
		for( size_t i = 0; i < m_Alias.size(); i++ )
		{
			if( 0 == _SourceID.compare(m_Alias[i].alias) )
			{
				_DestID = m_Alias[i].name;
				return true;
			}
		}
		return false;
	}
}

