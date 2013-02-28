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

#ifndef __ROTAMER_APPLICATOR_BASE
#define __ROTAMER_APPLICATOR_BASE

#include "verbosity.h"
#include "movebase.h"
#include "system/proximitygrid.h"
#include "system/molecule.h"
#include "workspace/pospointer.h"
#include "maths/fastrandom.fwd.h"

namespace Library
{
	class RotamerLibrary;
	class RotamerSet;
}
class PickResidueBase;
class SnapShot;

namespace Manipulator
{
	enum PD_API RotamerMode
	{
		ApplyTorsional, ///< When a rotamer is applied, just apply a set of torsions, disregarding deviations from idealised geometry.
		ApplyTorsional_Idealise, ///< Apply the standard idealised conformation present in the ffps, and then apply torsions.
		ApplyTorsional_PeriodicIdealise, ///< Apply the standard idealised conformation present in the ffps at a given interval, and then apply torsions. Useful in a protocol which can shaft idealised geometry, but it happens gradually.
		ApplyCartesian ///< Apply the rotamer by cartesian superimposition of stem residues.
	};

	//-------------------------------------------------
	/// \brief A base class which defines a class which states if a given rotamer should be active.
	class PD_API RotamerLink;
	class PD_API RotamerApplicatorBase; // forward declatation
	class PD_API RotamerFilterBase
	{
	public:
		RotamerFilterBase( const RotamerApplicatorBase& _app );
		virtual ~RotamerFilterBase();
		virtual RotamerFilterBase* clone() const = 0;
		virtual bool passes( const RotamerLink& _ir, size_t _rotID ) = 0;
		/// Does this filter instance require that the rotamer be applied to make its decision? Default, no.
		virtual bool requiresCoordinates() const { return false; }
	protected:
		RotamerFilterBase();
		const RotamerApplicatorBase* m_App;
	};
}

#ifdef SWIG
%template(ObjectContainer_RotamerFilter) ObjectContainer<Manipulator::RotamerFilterBase>;
#endif

namespace Manipulator
{
	class PD_API RotamerFilterContainer: public RotamerFilterBase, public ObjectContainer<RotamerFilterBase> 
	{
	public:
		RotamerFilterContainer();
		virtual ~RotamerFilterContainer();
		virtual RotamerFilterContainer* clone() const;
		/// Returns true if all child filters pass. The test *must* be independent of 
		/// other active side-chain orientations. i.e. you can only use information from 
		/// static atoms and rotamer properties!!
		virtual bool passes( const RotamerLink& _ir, size_t _rotID ); 
		virtual bool requiresCoordinates() const;
	};

	/// \brief Deactivates rotamers if they heavily clash with any static atoms.
	class PD_API RotamerStericFilter : public RotamerFilterBase
	{
	public:
		RotamerStericFilter( const RotamerApplicatorBase& _app );
		virtual ~RotamerStericFilter();
		virtual RotamerStericFilter* clone() const;
		virtual bool passes( const RotamerLink& _ir, size_t _rotID );
		virtual bool requiresCoordinates() const { return true; }
	private:
		const ProximityGrid* m_PtrStaticGrid; ///< A pointer to the parent applicators static grid
		double m_OLapFactor;
	};

	class RotamerApplicatorBase; // forward delcaration for RotamerLink
	/// \brief RotamerLink is a performance optimisation, storing pre-calculated info for RotamerApplicatorBase.
	/// \details 
	/// RotamerApplicatorBase owns an array of these of the same length as number of residues in the
	/// Molecule on which it works which match the internal residue picker.
	/// This class defines the link between a residue in a Molecule and the Library representing the possible states
	/// for that residue.
	/// \author Jon Rea 
	/// \todo Development is underway, class condition may not be finalised.
	/// \bug No known bugs
	class PD_API RotamerLink
	{
	private:
		friend class RotamerApplicatorBase;
		RotamerLink(){}
	public:
		size_t serial; ///< Its index in the parent!
		size_t resIndex; ///< The residue index of this rotamer link within the parent WorkSpace.	
		size_t libraryRotID; ///< The index of this residues type in the rotamer library
		std::vector<size_t> indexMap; ///< Store the atom indexes for rapid application
		std::vector<IndexHexet> scopeMap; ///< Store the chi scopes for rapid application
		const Library::RotamerSet* rot;
		WorkSpace* wspace;
		RotamerMode mode;
		int PeriodicIdealiseInterval; ///< A pointer to the periodic interval assigned in the parent RotamerApplicatorBase
		double maxInteractionDistance; ///< The distance from the CA atom to the furthest other heavy atom, plus that atoms radius. The maximal value of this quantity in any rotamer is stored.
		mutable PosStore sidechainPosCacheA; ///< A temporary store for a given conformation
		mutable PosStore sidechainPosCacheB; ///< A temporary store for a given conformation
		double maxProbability; ///< The highest single probability of any one rotamer in this set of rotamers, wether active or not.

		void reinit( WorkSpace* wspace, const Library::RotamerLibrary* _Lib, RotamerMode _mode, size_t _serial, int _resIndex, int _rotLibID, int _PeriodicIdealiseInterval );

		/// Function, using the above information to trigger the assignment of a requested Rotameric state
		/// Only applys if the rotIndex does not match the previous appyed value
		void apply( size_t rotIndex ) const;///< Always applies the rotamer
		void applyIfChangedID( size_t rotIndex ) const; ///< Applies the rotamer only if the interlal ID has changed (i.e. if you call apply() for rotamer 3 twice, the second one wont transform coordinates)
		void applyFFIdeal() const; /// Apply the idealised version from the forcefield.
		double probability( size_t rotIndex ) const;

		void setActive( size_t rotIndex, bool active );
		bool getActive( size_t rotIndex ) const;
		size_t countActive() const;

		inline size_t nRot() const { return m_Active.size(); } // always same length as rot->nRot()

		size_t lastRotamerApplied() const { return m_PreviousRotApply; }

	private:
		void calcMaxProbability(); ///< The result is stored in maxProbability.
		void calcMaxInteractionDistance(); ///< Assign the maxInteractionDistance member. Called once during setup. The answer doesnt change unless the rotamer lib does.
		mutable size_t m_PreviousRotApply;
		std::vector<bool> m_Active; ///< Each rotamer can be flagged as (de)commisioned from active service at any given time
		size_t m_ActiveCount;
	};

	typedef float RotEneStoreType;
	extern const RotEneStoreType RotEneStoreType_MAX;

	//----------------------------------------------------------------------
	/// \brief A base class for classes which apply rotamer states
	/// \details
	/// \author Jon Rea 
	/// \todo Development is underway, class condition may not be finalised.
	/// \bug No known bugs
	class PD_API RotamerApplicatorBase : public MoveBase
	{
	public:

		RotamerApplicatorBase( WorkSpace& _wspace, const Library::RotamerLibrary& _Lib, RotamerMode _mode = ApplyCartesian );
		RotamerApplicatorBase( WorkSpace& _wspace, const Library::RotamerLibrary& _Lib, const PickResidueBase& _Picker, RotamerMode _mode = ApplyCartesian );
		virtual ~RotamerApplicatorBase(){}
		virtual RotamerApplicatorBase* clone() const = 0;

		void printRotLinkActivity() const;

		Verbosity::Type OutputLevel;

		void overrideStaticGrid( const ProximityGrid& _StaticGrid ) { m_StaticGrid = _StaticGrid; }

		void reInitialise( RotamerMode _mode, const PickResidueBase& _Picker ); ///< Setup the internal RotamerLink array
		void reInitialise( const PickResidueBase& _Picker ); ///< Setup the internal RotamerLink array

		virtual int apply() = 0; ///< This behaviour must be set by deriving classes
		void applyRotamer( size_t _ResIndex, size_t _RotamerID ) const; ///< A mechanism to arbitrarily apply a single rotamer
		
		// Specialised apply functions
		void applyMostLikelyNonClashingRotamers(double cutoff) const;
		size_t applyMostLikelyNonClashingRotamer( size_t rotamerID, double cutoff ) const; // returns the number of active rotamers

		const PickResidueBase& getPicker() const { return m_ResPicker.data(); }
		void clearFilters(); ///< Remove all rotamer filters
		void enableFilters( bool on ) { m_FiltersEnabled = on; } ///< Temporarily enable or disable the rotamer filtering system
		void addFilter( RotamerFilterBase& _filter ) { m_ActivationFilters.add( _filter ); } ///< Add a new rotamer filter
		void addFilter_DefaultSteric(); ///< Add the default steric filter
		void deActivateAll(); ///< Deactivate all rotamers - only really does something permanent if filters are not enabled.
		void activateAll(); ///< Activate all rotamers - only really does something permanent if filters are not enabled.

		size_t nRot( size_t ir ) const;
		size_t size() const { return m_RotamerLinks.size(); }
		const ProximityGrid& getStaticGrid() const { return m_StaticGrid; }
		const ProximityGrid& getDynamicGrid() const { return m_DynamicGrid; }
		const Library::RotamerLibrary& getLibrary() const { return *m_Lib; }

		void setPeriodicIdealiseInterval( size_t interval ); ///< If ApplyMode == Torsional_PeriodicIdealise, this is our interval.

		const std::vector< std::vector<RotEneStoreType> >& getStoredPairwiseEne( size_t _rotLinkI, size_t _rotLinkK ) const;
		RotEneStoreType getStoredPairwiseEne( size_t _rotLinkI, size_t _rotJ, size_t _rotLinkK, size_t _rotM ) const;
		std::pair<RotEneStoreType, size_t> getStoredStaticEne( size_t _rotLinkIndex, size_t _slot ) const;
		size_t getStoredStaticEneIndexer( size_t _rotLinkIndex, size_t _slot ) const;

	protected:
		void recalcRotlinkMaxProbability();

		RotEneStoreType calcEPairwise( size_t _rotLinkI, size_t _rotLinkJ ) const;
		RotEneStoreType calcESteric_Static( size_t _rotLinkIndex ) const;

		void preCalculateSterics(); ///< Fills both the arrays below for efficiency
		void preCalculateSterics_PostCall( size_t i, size_t j ); ///< Add steric interactions AFTER reactivation
		mutable std::vector< std::vector< std::pair< RotEneStoreType, size_t > > > m_StoreSelfEnergy; ///< 2D jagged array. Stores self energies. Mutable as this is an internal cache that needs to be modified by public const functions.
		mutable std::vector< std::vector< std::vector< std::vector< RotEneStoreType > > > > m_StorePairEnergy; ///< 4D jagged array. Stores pair energies. Mutable as this is an internal cache that needs to be modified by public const functions.

		RotEneStoreType staticStericEnergy_Core( size_t _rotLinkIndex ) const; ///< Steric energy based upon the **static** grid, using the atomInteractionEnergy()
		virtual RotEneStoreType atomInteractionEnergy( const ParticleStore& atom, const SnapShot& snap, size_t i, size_t j ) const;

		ProximityGrid m_DynamicGrid;
		ProximityGrid m_StaticGrid;
		CloneHolder<PickResidueBase> m_ResPicker;
		RotamerMode m_Mode;
		void deactivateViaFilters();
		void deactivateViaFilters(RotamerFilterContainer& filters);
		bool m_FiltersEnabled;
		RotamerFilterContainer m_ActivationFilters;
		std::vector<RotamerLink> m_RotamerLinks;
		size_t m_PeriodicIdealiseInterval;
		const Library::RotamerLibrary* m_Lib;

	private:
		void preCalculateSterics_inner( size_t i );
		void preCalculateSterics_inner( size_t i, size_t j ); ///< called by preCalculateSterics
	};

	//-------------------------------------------------
	/// \brief Applies rotamers at random to a workspace.
	/// \details A simple rotamer applicator. Speedy because its rather simple. 
	/// Useful in outside apps which perform structural selection themselves. 
	/// Naive use of this class is almost certainly rather wasteful. 
	/// Try something like RotamerApplicator_SCWRL!
	/// \author Jon Rea 
	/// \bug No known bugs
	class PD_API RandomRotamerApplicator : public RotamerApplicatorBase
	{
	public:
		RandomRotamerApplicator( WorkSpace& _wspace, const Library::RotamerLibrary& _Lib, RotamerMode _mode = ApplyCartesian, bool _ProbabilityWeighting = false );
		RandomRotamerApplicator( WorkSpace& _wspace, const Library::RotamerLibrary& _Lib, const PickResidueBase& _Picker, RotamerMode _mode = ApplyCartesian, bool _ProbabilityWeighting = false ); 
		virtual RandomRotamerApplicator* clone() const { return new RandomRotamerApplicator(*this); }

		/// The default, applying random rotamers to the picked residues in the molecule
		virtual int apply();

		bool ProbabilityWeighting; ///< Should more likely rotamers be applied in preference by the random generator?

	private:
		void applyRandomRotamer( size_t _PickedIndex );
		Maths::FastRandom* rand;
	};
}

#endif


