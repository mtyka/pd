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

#ifndef __ROTAMERLIB_H
#define __ROTAMERLIB_H

#include <vector>
#include <string>
#include "tools/enum.h" // Required enum tool
#include "tools/cloneholder.h" // Required cloning helper
#include "primitives.h" // StringQuartet class
#include "manipulators/movebase.h" // Provides the MoveBase class
#include "maths/maths.fwd.h"
#include "maths/tntjama/tnt_array2d.h"
#include "workspace/workspace.fwd.h"

namespace Library
{
	class RotamerLibrary;
}

class PD_API FFParamSet; 
class PickResidueBase;
class PickResidueRange;
class PD_API MoleculeBase;
class PD_API MoleculeDefinition;

namespace Library
{
////	/// \brief Defines the represnetation of a given rotamer
////	enum WriteLibMode
////	{
////		Cartesian = 0,
////		Torsional = 1
////	};
////
////	/// \brief Adding "Original" to the WriteLibMode
////	enum WriteLibMode
////	{
////		Original = 2
////	};
////
	/// Inherit an enum
	//typedef InheritEnum< WriteLibMode, WriteLibMode > WriteLibMode;

	/// \brief Defines the represnetation of a given rotamer
	enum WriteLibMode 
	{
		WriteCartesian = 0,
		WriteTorsional = 1,
		WriteOriginal = 2
	};

	typedef PD_API StringQuartet ChiDef; /// A chi angle definition is meerly four atom names. Use a StringQuartet.

	//-------------------------------------------------
	/// \brief  Used by RotConvention
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	/// \author Jon Rea 
	/// \todo Development is underway, class condition may not be finalised.
	/// \bug No known bugs
	class PD_API ConventionDef
	{
	public:
		std::string resName;
		std::string a1;
		std::string a2;
		std::string a3;
		std::vector<ChiDef> chi;
	};


	//-------------------------------------------------
	/// \brief  Holds cartesian anchor and chi angle definitions
	/// \details Defines mappings between Chi angles and atom names and defines which atoms of a given residue should 
	/// be used as the anchor to apply the rotamer. One of these objects is owned by RotLib which simply initialises it by
	/// calling setToDefaultConvention(). The user can also create of of these classes and give it to the RotLib to
	/// override its internal definitions. This class can be fine-tuned by the user by adding new conventions or
	/// derived from to load conventions from a file (Not implemented, but should be trtivial if it is ever required).
	/// \author Jon Rea 
	class PD_API RotConvention
	{
	public:
		RotConvention();
		int anchorAtomType( const std::string& _resName, const std::string& _atomName ) const;
		void addConvention( const ConventionDef& _Conv );
		int findConvention( const std::string& _resName ) const; ///< Get the convention ID for a given residue name
		const ConventionDef& getConvention( int _id ) const; ///< Retrieve a convention based upon the internal ID obtained from findConvention()
		void clear(); ///< Clear all internal data
		void setToDefaultConvention(); ///< Use "standard" (simply the 20 standard amino acids) rotamer-anchor-atom and chi definitions
		
		/// If true, then efficiency can be gained. However, if one is building rotamers onto a backbone with a
		/// non-sensical CB position, then the rotamer application will fail. Therefore, the default is false to make
		/// the application more robust, as N, CA and C positions are usually well defined in all models. 
		// Set to true to gain efficiency when CB is known.
		bool CBIsThirdDefaultAnchor; 

	protected:		
		std::vector<ConventionDef> m_Convention;
	};


	//--------------------------------------------------------------------------------------------------
	/// \brief  Represents a rotamer library import converter
	/// \details This abstract base class defines the interrelationship between an arbritrary text file 
	/// containing rotamer definitions and the RotLib.addRotamer() calls. It is intended for use ONLY
	/// through the RotLib.convert( Filename, RotLibConvertBase ) call.
	/// \author Jon Rea 
	class PD_API RotLibConvertBase : public RotConvention
	{
	public: // Very limited publc interface - you can just make one...
		friend class RotamerLibrary;
		RotLibConvertBase();
	protected:
		/// The core of this class is used to call addRotamer() on the library
		virtual void readLib( const std::string& _LibraryFileName, RotamerLibrary& _RotLib ) = 0;
		void setUseAltnames( RotamerLibrary& _Lib, bool use ) const;
	};


	//-------------------------------------------------
	/// \brief   PD-format standard rotamer library importer
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	/// \author Jon Rea 
	/// \todo Development is underway, class condition may not be finalised.
	/// \bug No known bugs
	class RotamerSet; // Pre-definition
	class PD_API RotLibConvert_PD : public RotLibConvertBase
	{
	public:
		friend class RotamerLibrary; // only required as this is the default converter
		RotLibConvert_PD();		
	protected:
		virtual void readLib( const std::string& _LibraryFileName, RotamerLibrary& _RotLib );	
		void readBlock( std::istream& file, const StringBuilder& resName, RotamerLibrary& _RotLib ) const;
		void writeLib( const std::string& _LibraryFileName, Library::WriteLibMode _mode, const RotamerLibrary& _RotLib ) const;
		void writeRotamerSet( std::ofstream& stream, Library::WriteLibMode _mode, const RotamerSet& rotamer ) const;
	};


	//-------------------------------------------------
	/// \brief  Old PD-format standard rotamer library importer
	/// \details Reads in legacy rotamer files
	/// \author Jon Rea 
	/// \todo Development is underway, class condition may not be finalised.
	/// \bug No known bugs
	class PD_API RotLibConvert_OldPDFormat : public RotLibConvertBase
	{
	public:
		RotLibConvert_OldPDFormat(){};		
	protected:
		virtual void readLib( const std::string& _LibraryFileName, RotamerLibrary& _RotLib );	
	};


	//-------------------------------------------------
	/// \brief  BRIEF DESCRIPTION
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	/// \author Jon Rea 
	/// \todo Development is underway, class condition may not be finalised.
	/// \bug No known bugs
	class PD_API RotamerAtom
	{
	public:
		RotamerAtom();
		RotamerAtom(const std::string &_name, const Maths::dvector &_idealPos, bool _isHydrogen );

#ifndef SWIG
		// Used by RotamerLibrary::shouldRemoveRotamerAtom()
		bool operator==( const std::string& resName ) const;
#endif

		std::string Name;
		Maths::dvector IdealPos;
		bool Defined;
		bool IsHydrogen;
	};

	//-------------------------------------------------
	/// \brief  Extends the atom names stored in a ChiDef with the indexes of those atoms in a given RotamerSet
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	/// \author Jon Rea 
	/// \todo Development is underway, class condition may not be finalised.
	/// \bug No known bugs
	class PD_API RotamerChi : public ChiDef
	{
	public:
		RotamerChi( const ChiDef& _chi, const std::vector<RotamerAtom>& _atoms );
		IndexQuartet rot;
	};


	//-------------------------------------------------
	/// \brief  BRIEF DESCRIPTION
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	/// \author Jon Rea 
	/// \todo Development is underway, class condition may not be finalised.
	/// \bug No known bugs
	class PD_API BackboneScope
	{
	public:	
		enum BackboneScopeType
		{
			Undefined,
			Alpha,
			Beta,
			LeftAlpha,
			Specific
		};

	private:
		BackboneScopeType m_BBScope;
		TorsionalTolleranceRange m_Phi; ///< If 'Specific' backbone scope is selected, what is our scope?
		TorsionalTolleranceRange m_Psi; ///< If 'Specific' backbone scope is selected, what is our scope?

	public:
		BackboneScope();
		BackboneScope(BackboneScopeType _type);

		void set(BackboneScopeType _type);
		void setSpecific( TorsionalTolleranceRange _Phi, TorsionalTolleranceRange _Psi );

		bool isInRange( double phi, double psi  ) const { return m_Phi.isInRange(phi) && m_Psi.isInRange(psi); }

		TorsionalTolleranceRange getPhi() const;
		TorsionalTolleranceRange getPsi() const;
		BackboneScopeType getScopeID() const { return m_BBScope; }

		static const TorsionalTolleranceRange PhiAplpha;
		static const TorsionalTolleranceRange PsiAplpha;
		static const TorsionalTolleranceRange PhiBeta;
		static const TorsionalTolleranceRange PsiBeta;
		static const TorsionalTolleranceRange PhiLAplpha;
		static const TorsionalTolleranceRange PsiLAplpha;

		std::string toString() const;
		void parseString(const std::string& s);
	};

	//-------------------------------------------------
	/// \brief  ProbabilityBase - maps a probability to a rotamer state
	/// \details Base class for a probability mapped to a residues structural state
	/// This can be based upon any given property
	/// On of these classes is owned by each rotamer
	/// \author Jon Rea 
	class ProbabilityBase
	{
	public:
		ProbabilityBase(){}
		virtual ~ProbabilityBase(){}
		virtual ProbabilityBase* clone() const = 0;

		virtual void clear() = 0;
		virtual void toStream( std::ostream& _Stream ) const = 0; ///< Output this structure to a stream
		virtual void fromStream( std::istream& _Stream ) = 0; ///< Set this structure using data in a stream
		void toScreen() const { toStream(std::cout); }

		/// returns the probability of 'rotamer i' in this molecular context
		/// return value is a relative probability if valid, or <0.0 (-1.0) if not
		/// The relative probability relates this rotamer to others in the RotamerSet
		virtual double probability( const MoleculeBase& _mol, int ir ) const = 0;
	};

	class ConstantProbability : public ProbabilityBase
	{
	public:
		ConstantProbability();
		ConstantProbability( double _prob );
		virtual ConstantProbability* clone() const;

		void setProbability( double value );
		virtual void clear();
		virtual void toStream( std::ostream& _Stream ) const; ///< Output this structure to a stream
		virtual void fromStream( std::istream& _Stream ); ///< Set this structure using data in a stream
		virtual double probability( const MoleculeBase& _mol, int ir ) const; ///< Returns a set probability

	protected:
		double m_Prob;
	};

	class ProbabilityByBBScope : public ConstantProbability
	{
	public:
		ProbabilityByBBScope();
		ProbabilityByBBScope( double _prob );
		ProbabilityByBBScope( BackboneScope::BackboneScopeType _scope, double _prob );
		virtual ProbabilityByBBScope* clone() const;

		virtual void toStream( std::ostream& _Stream ) const; ///< Output this structure to a stream
		virtual void fromStream( std::istream& _Stream ); ///< Set this structure using data in a stream
		virtual double probability( const MoleculeBase& _mol, int ir ) const; ///< Returns a set probability if in-scope, or -1.0 if out of scope

	private:
		BackboneScope m_BBScope;		
	};

	class ProbabilityByPhiPsiMap : public ProbabilityBase
	{
	public:
		ProbabilityByPhiPsiMap( size_t binCount );
		virtual void clear();
		virtual ProbabilityByPhiPsiMap* clone() const;

		virtual void toStream( std::ostream& _Stream ) const; ///< Output this structure to a stream
		virtual void fromStream( std::istream& _Stream ); ///< Set this structure using data in a stream
		virtual double probability( const MoleculeBase& _mol, int ir ) const; ///< Returns a probability based upon the current phi / psi angles

		void setTo( const ProbabilityByPhiPsiMap& _clone );

		void assignFrom( StringBuilder& _sb ); ///< Read in phi psi and probability from a whitespace delimited string, consuming the string in the process
		void assignFrom( double phi, double psi, double probability );

		size_t numAssignments() const { return m_Assignments; }

	private:
		size_t m_Assignments;
		TNT::Array2D<double> m_BBMap;
		double m_RadResln;
	};

	//-------------------------------------------------
	/// \brief  BRIEF DESCRIPTION
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	/// \author Jon Rea 
	/// \todo Development is underway, class condition may not be finalised.
	/// \bug No known bugs
	class PD_API Rotamer
	{
	public:
		friend class PD_API RotamerLibrary;
		friend class PD_API RotamerSet;

		Rotamer() 
		{
			m_SourceRepresentation = WriteCartesian; // we have to initialise to something ...
			m_ProbMap = ConstantProbability();
		}

		inline const Maths::dvector& getPos( size_t i ) const { return m_Pos[i]; } ///< Get a const hook to a given atom position
		inline double getChi( size_t i ) const { return m_Chis[i]; } ///< Get a const hook to a given Chi definition
		inline size_t nAtom() const { return m_Pos.size(); } ///< How many atoms do we contain?
		inline size_t nChi() const { return m_Chis.size(); } ///< How many chi definitions do we contain?
		inline ProbabilityBase& getProbMap() { return m_ProbMap.data(); }
		inline const ProbabilityBase& getProbMap() const { return m_ProbMap.data(); }		
		inline WriteLibMode getSource() const { return m_SourceRepresentation; } ///< Where the 'ding-dang-doodly' did this rotamer come from?

		void info() const;

		/// Externally override the probability decision for this rotamer
		void setProbMap( const ProbabilityBase& _prob ) { m_ProbMap = _prob; }

	private:
		CloneHolder<ProbabilityBase> m_ProbMap;
		std::vector<double> m_Chis;
		std::vector<Maths::dvector> m_Pos; ///< This is dvector and not RotamerAtom as the parent rotamer set holds these. Here we only store positions which map to that of the parent container
		WriteLibMode m_SourceRepresentation; ///< Was this rotamer originally defined as a cartesian or torsional state?
	};


	//-------------------------------------------------
	/// \brief   A set of rotamers for a single residue type
	/// \details RotamerSet can only be created via the RotamerLibrary
	/// \author Jon Rea 
	/// \todo Development is essentially complete
	/// \bug No known bugs
	class PD_API RotamerSet
	{
		// *NOTE*: No public constructor logic
		// RotamerSet();
	public:
		const std::string getName() const { return name; }
		size_t nAtom() const { return m_Atom.size(); }
		size_t nRot() const { return m_RotRes.size(); } 
		const Rotamer& getRotamer(size_t i) const { return m_RotRes[i]; }
		Rotamer& getRotamer(size_t i) { return m_RotRes[i]; }
		int getCartesianAnchorRot1() const { return m_CartesianAnchorRot1; }
		int getCartesianAnchorRot2() const { return m_CartesianAnchorRot2; }
		int getCartesianAnchorRot3() const { return m_CartesianAnchorRot3; }
		const std::vector<RotamerAtom>& getAtoms() const { return m_Atom; }
		const std::vector<RotamerChi>& getChis() const { return m_Chis; }
		bool usesAltnames() const { return m_AltNames; }
		int searchAtomName( const std::string& atomName ) const;

		// **CODE-DESIGN-NOTE**, these need to be here. RotamerApplicatorBase *must* derive from MoveBase to be of most use,
		// which means that it can only work on WorkSpace. It would be nice if Rotamers could be applied to any
		// MoleculeBase as well as WorkSpace, and therefore the logical place for the apply function is within the
		// RotamerSet, with wrapper functions in RotamerLibrary and speed-optimised wrappers in RotamerApplicatorBase.

		int closestRotamer( MoleculeBase& _molBase, int ir, double& cRMS ) const;

		// Cartesian
		void calcAtomIndexMap( const MoleculeBase& _molBase, int ir, std::vector<size_t>& _map ) const; // Fill the map for applyRotamer()
		void applyRotamerCartesian( MoleculeBase& _molBase, int ir, size_t _rotIndex ) const; ///< Apply the cartesian rotamer, internally calculating atom correspondence (slow)
		void applyRotamerCartesian( MoleculeBase& _molBase, int ir, size_t _rotIndex, const std::vector<size_t>& _indexMap ) const; ///< Apply the cartesian rotamer, using pre-calculated atom correspondence (speedy, but needs to be correct!)

		// Torsional
		void calcTorsionalScope_FF( const MoleculeDefinition& _molBase, std::vector<IndexHexet>& _scopeMap ) const;
		void calcTorsionalScope_Rot( const MoleculeDefinition& _molBase, std::vector<IndexHexet>& _scopeMap ) const; ///< Fills the IndexHexet with index information
		void applyRotamerTorsional( MoleculeBase& _molBase, int ir, size_t _rotIndex ) const; ///< SLOW! - calls calcTorsionalScope_FF() 
		void applyRotamerTorsional( MoleculeBase& _molBase, int ir, size_t _rotIndex, const std::vector<IndexHexet>& _scopeMap ) const; ///< FAST. Apply the torsional rotamer using pre-calculated calcTorsionalScope().
		void applyFFIdealised( MoleculeBase& _molBase, int ir, const std::vector<size_t>& _indexMap ) const; ///< Used to ensure that the torsional representation doesnt distort over time

	protected:
		// Protected-ness: You are only allowed to make a RotamerSet through RotamerLibrary
		friend class RotamerLibrary;
		RotamerSet( const std::string& _ResName, bool _AltNames, const RotConvention& _Conv, const MoleculeDefinition& _FFMol );

		// Setup functions	
		void add( const std::vector<RotamerAtom>& _atom, const ProbabilityBase& _probMap ); ///< Add rotamer information
		void add( const std::vector<double>& _Chis, const ProbabilityBase& _probMap ); ///< Add rotamer information

#ifndef SWIG
		bool operator==( const std::string& resName ) const { return resName == name; }
#endif

	private:
		/// The actual rotamer information
		std::vector<Rotamer> m_RotRes;

		// Naming data
		const RotConvention* m_Conv;
		bool m_AltNames; ///< Are the names stored in m_Atom AltNames?
		std::string name; ///< The name of the residue which this 'RotamerSet' represents
		const MoleculeDefinition* m_FFMol;
		std::vector<RotamerAtom> m_Atom; ///< What atoms does each 'Rotamer' contain?
		std::vector<RotamerChi> m_Chis; ///< Definitions as to what the chi torsions actually are

		int m_CartesianAnchorRot1; ///< The index of anchor atom 1 in the m_Atom list and in all dvector lists in m_RotRes
		int m_CartesianAnchorRot2; ///< The index of anchor atom 2
		int m_CartesianAnchorRot3; ///< The index of anchor atom 3

		void placeOrigin(Rotamer& _r) const; ///< Ensure that the rotamers cartesian system centres on the origin
		void calcTorsions(Rotamer& _r) const; ///< From the cartesian system
		void calcCartesian(Rotamer& _r) const; ///< From the torsional system

		// Convention mappings, called by the parent RotamerLibrary
		void assignTorsionDefs(); ///< Torsion definitions are assigned using these convention mappings
		void assignAnchorIndexes(); ///< Assign each m_CartesianAnchorX
	};


	//-------------------------------------------------
	/// \brief  BRIEF DESCRIPTION
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	/// \author Jon Rea 
	/// \todo Development is underway, class condition may not be finalised.
	/// \bug No known bugs
	class PD_API RotamerLibrary
	{
		friend class RotLibConvertBase; // allows setUseAltNaming() to be called by this class
	public:
		enum RebuildTolleranceType
		{
			Hydrogens, ///< Allow hydrogens to be rebuilt for cartesian rotamer addition.
			NoRebuild  ///< Dont perform any rebuild, and if atoms are not present, complain loudly.
		} RebuildTollerance;

		RotamerLibrary( const FFParamSet& _ffps );
		RotamerLibrary( const FFParamSet& _ffps, const std::string _FileName); ///< Constructor which calls readLib() for the file specified.

#ifndef SWIG
		bool operator()( const RotamerAtom& a, const RotamerAtom& b ) const; ///< Internal sorting helper, requires public.
		bool operator()( const RotamerAtom& a, const std::string& b ) const; ///< Internal sorting helper, requires public.
#endif

		// Library interaction
		void readLib( const std::string& _FileName ); ///< Reads a PD format library
		void convertLib( const std::string& _FileName, RotLibConvertBase& _Converter ); ///< Import a foreign format library
		void writeLib( const std::string& _FileName, WriteLibMode _mode = WriteOriginal ); ///< Write a PD format library
		void writePDB( const std::string& _FileName ) const; ///< Write a PDB format file containing our rotamers
		void writePDB( const std::string& _FileName, size_t _residueID ) const; ///< Write a PDB format file containing our rotamers for a specific residue

		void info(); ///< print rotamer library info	

		size_t size() const { return m_RotRes.size(); }

		bool isFinalised() const { return m_Finalised; } ///< On first use, the library is finalised and various internals set. This prevents the addition of further rotamers...

		// Why-oh-why can everyone not agree on one set of names?
		bool useAltNames() const { return m_AltNames; } ///< Should we import and export using the PDB or RAW name for rotamer import?

		void autoFinalise( bool _do ) { m_AutoFinalise = _do; } ///< After a readLib() or convertLib() operation, should the library be finalised?

		void finalise(); ///< Scale all probabilities and prevent the addition of new rotamers.
		void clear(); ///< Clear all definitions from the library

		/// When we look at HIS, it has HID HIE and HIP ionisation states which differ by the position
		/// and number and names of hydrogen atoms only. As most rotamer definitions do not define the
		/// hydrogen atoms, and by definition any based on chi angles; the best way of dealing with this
		/// seems to be adding a mapping from the ionised name (e.g. HID) to the name in the rotamer lib.
		/// Of course this then assumes that the rotameric distributions are identical with equal propensities,
		/// which is clearly wrong, but then rotamer states are PDB-derived, which does not define ionisation,
		/// information, and therefore rotamer libraries are averaged over all ionisation states anyway.
		/// If one wanted to have different distributions, then dont add ionisation name mappings, but instead,
		/// specify individual definitions for HIP HIE and HID in the rotamer file.
		void addIonisationAlias( const std::string& _rotSource, const std::string& _ionisedAlias );
		void addIonisationAlias( const StringPair& alias );
		void addIonisationAlias( const std::vector<StringPair>& aliases );
		void clearIonisationAliases();
		void addIonisationAliasWorkingDefaults();

		void addAsBlankRotamer( const std::string& resName ); ///< Things like ALA and GLY have no rotamers, but must be "defined" as such.
		void addRotamer( const std::string& resName, const std::vector<RotamerAtom>& _Positions, const ProbabilityBase& _probMap );
		void addRotamer( const std::string& resName, const std::vector<double>& _Chis, const ProbabilityBase& _probMap );
		void addRotamer( const MoleculeBase& _Mol, size_t _IRes, const ProbabilityBase& _probMap );

		const std::vector<StringPair>& getIonMaps() const { return m_WorkingIonMaps; }
		size_t getIDForResidue( const std::string& resName ) const; ///< Return the internal ID for a residue name. This can be used later to access the Rotamer without incuring a string look-up operation.
		size_t nRotIn( size_t _RotamerID ) const; ///< How many rotamers in this set?
		size_t nRot() const; ///< The number of rotamer definitions
		const RotamerSet& getRotamerSet( size_t _id ) const;
		RotamerSet& getRotamerSet( size_t _id );

		void overrideConventions( const RotConvention& _ConvPermanentReference );
		void setDefaultConventions();

	private:
		void commonInit(); ///< Constructor helper
		void setUseAltnames( bool _do ) { m_AltNames = _do; } ///< Should we import and export using the PDB or RAW name for rotamer import?

		void ensureRotamerSetIsBuilt(const std::string& resName); ///< If no rotamerset already exists for this resName, build one and append.
		void addRotamer_Core( const std::string& resName, const std::vector<RotamerAtom>& _Positions, const ProbabilityBase& _probMap );
		void addRotamer_Core( const std::string& resName, const std::vector<double>& _Chis, const ProbabilityBase& _probMap );

		/// Make sure each imported set of atoms has the same contents, add missing ones and flag defined as false
		std::vector<RotamerAtom> standardiseAtoms( const std::vector<RotamerAtom>& _Positions ) const; ///< Standardise '_Positions' using the info in the current RotamerSet, referenced by 'm_CurrentRotLink'
		bool shouldRemoveRotamerAtom( const RotamerAtom rot ) const;
		void positionUndefinedAtoms( std::vector<RotamerAtom>& _rot ) const;

		std::vector<StringPair> m_UserIonMaps;    ///< (User added) Holds mappings like HIE -> HIS and HID -> HIS
		std::vector<StringPair> m_WorkingIonMaps; ///< (Working-set) Holds mappings like HIE -> HIS and HID -> HIS
		std::vector<RotamerSet> m_RotRes;
		RotLibConvert_PD m_PDConverter; ///< This is used by readLib, which calls convertLib with this default converter
		bool m_Finalised;
		bool m_AutoFinalise;
		bool m_AltNames; ///< Do rotamer names correspond to altnames; if not then they are PDB names.
		const FFParamSet* m_ffps;

		// It's all convention sweetie ;-)
		RotConvention m_DefaultConvention; ///< The default anchor and chi definitions
		const RotConvention* m_Convention; ///< A pointer to the current rotamer anchor and chi definitions
		bool m_UserOverrideConventions; ///< The user has balled overrideConventions(). When the user calls setDefaultConventions(), this is unflagged.

		// Helpers
		mutable int m_CurrentFFPSLink; ///< A mild simplifying mini-hack (used in an internal sorting operation)
		mutable size_t m_CurrentRotLink; ///< A mild simplifying mini-hack (used in an internal sorting operation)
	};
}

#endif

