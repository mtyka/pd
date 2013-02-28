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

#ifndef __INFILE_H
#define __INFILE_H

#include <string>
#include <vector>

#include "sequence/alignment.h"
#include "system/system.h"
#include "system/genpolymer.h"
#include "system/rebuilder.h"

class MoleculeBase;
class System;

namespace IO
{

	class PD_API FileMoleculeMaker;

	//-------------------------------------------------
	//
	/// \brief 
	/// Particle data members common to all import file types.
	///
	/// \details 
	/// NOTE: This is a much simpler form of particle than those present in fundamentals.h which
	/// requires other information not usually present in a coordinate text file.
	///
	/// \author Mike Tyka & Jon Rea 
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class PD_API FileParticle
	{
	public:
		int atomNum;
		std::string atomName;
		std::string resName;
		char chainID;
		char iCode;
		int resNum;
		double x;
		double y;
		double z;

		void info() const;
	};

	void setPolymerPositions( const FFParamSet& ffps, MoleculeBase& mol, const Sequence::AlignmentDef& sequencePair, const std::vector<FileParticle>& particles, Verbosity::Type verbose );
	void setPolymerPositions( const FFParamSet& ffps, MoleculeBase& mol, const std::vector<FileParticle>& particles, Verbosity::Type verbose = Verbosity::Normal );


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
	class PD_API File_Molecule : public Sequence::ExpSeqPair
	{
	public:
		friend class  ::IO::FileMoleculeMaker;

		File_Molecule( char _ChainID );
		File_Molecule( const File_Molecule& _clone );

		// Public Modifiers
		void setPrefix(Library::PrePostFix _First, Library::PrePostFix _Last);

		// Public printing functions
		void printToScreen() const;

		// Public readonly accessors
		const Sequence::BioSequence& getBiologicalSeq() const { return bioSequence; }
		const Sequence::BioSequence& getStructuralSeq() const { return structuralSequence; }
		size_t getBioLength() const { return bioSequence.size(); }
		const std::vector<FileParticle>& getAtomLines() const { return atomLines; }
		char getChainID() { return chainID; }

		bool isEmpty() const { return atomLines.size() == 0; }

	protected:
		// Structure member data		
		Sequence::BioSequence bioSequence;
		Sequence::BioSequence structuralSequence;
		std::vector<FileParticle> atomLines;
		char chainID;
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
	class FileBond
	{
	public:
		enum BondType
		{
			Invalid,
			Disulphide,
			Covalent,
			HBond
		};

		FileBond();

		void clear();

		bool matches( const Sequence::BioSequence& _BioSeq, size_t _Index ) const;
		bool matches( const FileParticle& atom ) const;

		BondType Type;      ///< Bond Type

		std::string r1Name; ///< 1. Residue name.
		char r1ChainID;     ///< 1. Chain identifier.
		int r1ResNum;       ///< 1. Residue sequence number.
		char r1ICode;       ///< 1. Insertion code.

		std::string r2Name; ///< 2. Residue name.
		char r2ChainID;     ///< 2. Chain identifier.
		int r2ResNum;       ///< 2. Residue sequence number.
		char r2ICode;       ///< 2. Insertion code.
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
	class FileCovalency
	{
	public:
		FileCovalency(){};
		~FileCovalency(){};

	#ifndef SWIG
		FileBond& operator[](size_t _Index) 
		{ 
			if( _Index >= m_Lines.size() ) 
				throw OutOfRangeException("FileBond[] out of range!");
			return m_Lines[_Index]; 
		}

		const FileBond& operator[](size_t _Index) const
		{ 
			if( _Index >= m_Lines.size() ) 
				throw OutOfRangeException("FileBond[] out of range!");
			return m_Lines[_Index]; 
		}
	#endif

		size_t getCount() const { return m_Lines.size(); }
		void addBond(const FileBond& fileBond) { m_Lines.push_back(fileBond); }
			
	protected:
		std::vector<FileBond> m_Lines;
	};



	//-------------------------------------------------
	//
	/// \brief  The base class for any System importing file Type
	///
	/// \details The 'FileInBase' class also manages name mappings
	/// All things that affect the naming of the atoms and residues in MMLibs resultant 'internal definition'
	/// are managed by this class ...
	///
	/// Atom and Residue naming in file formats like the PDB are rather complex beasts:
	/// 1) Different atom names are used depending on which convension is followed. The only way to manage this is to 
	///    have knowledge based mappings from one nameset to another. This is provided by class 'NamingConventions' and 
	///    is fully extendable to add new namesets.
	/// 2) in files like PDB MODRES lines can define residues that are modified in some way and therefore have names that stray away
	///    from the standard 20 residues. The mapping back to these residues is provided by this MODRES line. This can be 
	///    added as a feature to 'FileInBase' by adding an additional AliasMapper.
	/// 3) In a really simple file that just contains coordinates, there are occations when you will encounter something
	///    like a CSE - seleno cysteine - and you want to simulate it. Our standard amber forcefield does not contain
	///    a definition for CSE, but the residue is effectively a CYS. It is therefore possible to use this knowledge
	///    to hack-rename the CSE to a SYS. This is not flexible and therefore disabled by default. Knowledge is also only
	///    used if none of the other avenues are open i.e. A ff definition is added, or a MODRES line is present.
	/// An 'atomic babel-fish' ;-D
	/// http://en.wikipedia.org/wiki/Tower_of_Babel
	///
	/// \author Jon Rea 
	///
	/// \todo
	///
	/// \bug
	///
	class FileBabelBase 
	{
	public:
		FileBabelBase();

		virtual bool ReinterpretName( FileParticle& atom ) const;

		// Atom name conversion on import of the PDB file using the naming StandardNameSets in 'sequence.h'
		void setAtomNameConverterEnabled( bool _Enabled ) { m_EnableAtomConversion = _Enabled; } ///< Allows the user to enable single atom name conversion via a knowlege-base
		void setAtomNameConverter( Library::StandardNameSets _ConvertFrom ) { m_AtomNameset = _ConvertFrom; }
		Library::StandardNameSets getAtomNameConverter() const { return m_AtomNameset; } ///< Which name-set are we currently converting to

		// Generalised Knowledge-based hacking
		void enableKnownedgeBasedRenaming( bool _enable ) { m_UseKnowledge = _enable; } ///< Enable allan-style if() knowledge-based renaming function (hacky and disabled by default)

	private:
		// Knowledge-based:
		bool m_UseKnowledge;

		// Perform individual ATOM name import conversion?
		bool m_EnableAtomConversion; // Enable babel-fish style naming convention conversion (disabled by default)
		Library::StandardNameSets m_AtomNameset;
		Library::NamingConventions* m_AtomNamer;
	};



	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
	/// 'FileImportBase' houses the relevent naming mappers
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
	class FileImportBase : public FileBabelBase
	{
	public:
		FileImportBase( const std::string& _filename, bool allocateClass = true );

		const std::string& getFileName() const;

		// Sequence Management
		const Sequence::AlignerBase& getAligner() const;
		void setAligner( Sequence::AlignerBase& _AlignerOverride );
		void setAlignerDefault();
		void setAlignerDirect();
		
		// Public so that the end-user can add alias mapper overrides to this class
		void addAliasMapper(const Library::AliasMapper& _mapper);
		void setCovalency( FileCovalency& connections );
		void setClassMapper(const Library::ClassMapper& _mapper);

		const Library::ClassMapper& getClass() const;
		Library::ClassMapper& getClass();

		const Library::AliasMapperCollection& getAlias() const;
		const FileCovalency& getCovalency() const;

		void setFilter( const std::string& _filter );
		void setFilter( Library::ResidueClass _filter );
		//void setFilter( Library::ResidueClassGroup _filter );
		void setFilter( Library::ExtdResidueClasses _filter );
		bool isResequencableFilter() const; ///< We are only allowed to resequence if the molecule in question is a defined filter class and a polymer
		bool PassesFilter( const std::string& _ResName ) const;
		std::string getFilterText() const;

		virtual bool ReinterpretName( FileParticle& atom ) const; ///< Use various interbal converters to rename this atom

	protected:
		// Derived class controls access to these functions (Important for PDB)
		Verbosity::Type getVerbose() const;
		void setVerbosity( Verbosity::Type _BeVerbose );

		Library::AliasMapperCollection& getAlias();
		FileCovalency& getCovalency();

		// Protected filter control
		bool m_FilterFlag; ///< True if we are using m_FilterS, false if we are using m_FilterE.
		std::string m_FilterS; ///< A string-based class filter
		Library::ExtdResidueClasses m_FilterE; ///< Enumeration based class filter

		// Helper functions for 'ReinterpretName()'
		void ReinterpretFromCovalency( FileParticle& atom ) const; ///< Some residue AND bonding specific renaming needed in proper PDB import
		bool ReinterpretCore(FileParticle& atom ) const; ///< Reinterpet a residue name using the internal 'AliasMapperCollection'

	private:
		std::string m_Filename;
		Verbosity::Type m_Verbose;
		Library::AliasMapperCollection m_Alias;
		FileCovalency* m_Covalency;
		FileCovalency m_CovalencyDefault;
		Library::ClassMapper m_Class;

		// Sequence alignment control
		Sequence::AlignerBase* m_Aligner;
		Sequence::SimpleAligner m_AlignerDefault;
		Sequence::DirectAligner m_AlignerDirect;
	};


	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
	/// The base classes for any 'System' importing file Type
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
	class FileInBase : public FileImportBase, 
                       public System, 
                       public RebuilderOwner, 
                       public Library::CustomAliasMapper
	{
	public:
		FileInBase( const FFParamSet &_ffps, const std::string& _filename );

		// Public flags
		bool TreatMolsAsGroups; ///< See FileInBase::build() - should we treat small molecules like water as a collective molecule?
		bool CentreOnBuild; ///< Allow recentre the molecules after they have been built?

	protected:
		void baseLoad( ::IO::FileMoleculeMaker& modelMaker, char chainID, const Sequence::BioSequence* _SequenceOverride = NULL);
		void SupplementaryBondScan( Molecule& sysMol );

	private:
		virtual bool ReinterpretName( FileParticle& atom ) const;
		void baseLoadCore( Library::ResidueClass _Class, FileMoleculeMaker& modelMaker, char chainID, const Sequence::BioSequence* _SequenceOverride = NULL);
		void build( File_Molecule& mol, bool _IsPolymer ); ///< Perform 'File_xxx' related procedures and invoke use m_Rebuilder to position undefined atoms
		void buildFinalise( Molecule& _NewMol, bool _IsPolymer ); ///< The final tweaks for each new molecule generated
		void doCentring( Molecule& _NewMol ); ///< Call the post-build recentre procedure
	};

	/// The base classes for any 'System' importing file Type





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
	class FileToolBase : public FileImportBase
	{
	public:
		FileToolBase( const std::string& _filename );
	};




	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
	/// Derive from this class to be able to play with the internal members of FileMolecule in a structured way.
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
	class PD_API FileMoleculeMaker
	{
	public:
		FileMoleculeMaker( const FileImportBase& _parent );
		virtual ~FileMoleculeMaker();
		FileMoleculeMaker( const FileMoleculeMaker &_clone ) :
			m_Parent(&_clone.getFileParent())
		{
		}

	#ifndef SWIG
		const FileMoleculeMaker& operator= (const FileMoleculeMaker &_clone)
		{
			m_Parent = &_clone.getFileParent();
			return (*this);
		}
	#endif

		virtual File_Molecule make( char chainID ) = 0;
		virtual File_Molecule make( char chainID, const Sequence::BioSequence* buildSequenceOverride ) = 0; /// 'buildSequenceOverride' must be a pointer because we need to have a possible unallocated state
		
		virtual std::vector<char> findChainIDs() const = 0;
		virtual bool HasChainID(char _RequestChain ) const = 0;

	protected:
		// Multiple protected helper function calls for the derived classes virtual functions...
		void AssignBioSequence( File_Molecule& _Mol, const Sequence::BioSequence& _BioSequence ) const;
		bool appendFileParticle( File_Molecule& _Mol, const FileParticle& _Particle ) const; ///< Returns true if it was decided that the particle should be appended
		bool wouldAppendFileParticle(const FileParticle& _Particle ) const; ///< Would appendFileParticle() return true for this particle?
		void DetectPrefix( File_Molecule& mol ) const;
		void finaliseSequence( File_Molecule& mol ) const; ///< Autogenerate the structural sequence from the FileParticles and perform an alignment to BioSequence data
		const FileImportBase& getFileParent() const;
		bool isResequencableFilter() const; ///< If we are filtering for solvent, het molecules or ions, then reassigning the sequence doesnt make sense (e.g. a load of water molecules) - This flag tells us if resequenceing is allowed.
		std::string getFilterText() const;
	private:
		/// Called by the internal appendFileParticle()
		bool PassesFilter( const std::string& _ResName) const;
	private:
		/// An internal reference for access to any available sequence data etc.
		const FileImportBase* m_Parent;
	};





	////////////// TEMPORARY LOCATION


	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
	/// The base classes for any 'System' importing file Type
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

	/*
		18640 !NATOM
					 1 SEG0 0    GLY  N    CAV3   0.000000   20.0000          0
					 2 SEG0 0    GLY  CA   CAV3   0.000000   20.0000          0
					 3 SEG0 0    GLY  C    CAV3   0.000000   20.0000          0

	*/

	struct PSF_AtomLine
	{
		int index;
		std::string segname;
		int resnumber;
		std::string resname;
		std::string atomname;
		std::string ffname;
		double charge;
		double mass;
		int    extratoken;
	};

	class Load_PSF: public System 
	{
	public:
		Load_PSF( const FFParamSet &_ffps, const std::string& _filename );
	private:
	};

} //namespace IO

#endif

