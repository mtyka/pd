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

#include <algorithm>

#include "primitives.h"
#include "pickers/basicpickers.h"
#include "library/rotamerlib.h"
#include "workspace/workspace.h"
#include "workspace/bondorder.h"
#include "workspace/pospointer.h"

#include "rotamer_applicatorbase.h"

namespace Manipulator
{

	// -------------------------
	// Helpers
	// -------------------------

	/// Compare the ratio of the distance to the sum of the VDW radii. If under _OlapFac, we have a clash.
	bool simpleClash( 
		const Particle& _atom1, 
		const Particle& _atom2, 
		double _OlapFac = 0.75 
	){
		double sqrdist = _atom1.pos().sqrdist( _atom2.pos() );
		double radSum = _atom1.radius + _atom2.radius;
		radSum *= radSum;
		return (sqrdist / radSum) < _OlapFac;
	}

	const RotEneStoreType RotEneStoreType_MAX = FLT_MAX;

	// -----------------------------
	// Classes: RotamerFilters
	// -----------------------------

	RotamerFilterBase::RotamerFilterBase()
		: m_App(NULL)
	{
	}

	RotamerFilterBase::RotamerFilterBase( const RotamerApplicatorBase& _app )
		: m_App(&_app)
	{
	}

	RotamerFilterBase::~RotamerFilterBase()
	{
	}

	RotamerFilterContainer::RotamerFilterContainer()
	{
	}

	RotamerFilterContainer::~RotamerFilterContainer()
	{
	}

	RotamerFilterContainer* RotamerFilterContainer::clone() const
	{
		return new RotamerFilterContainer(*this);
	}

	bool RotamerFilterContainer::passes( const RotamerLink& _ir, size_t _rotID )
	{
		for( size_t i = 0; i < size(); i++ )
		{
			if( !element(i).passes( _ir, _rotID ) )
				return false;
		}
		return true;
	}

	bool RotamerFilterContainer::requiresCoordinates() const
	{
		for( size_t i = 0; i < size(); i++ )
		{
			if( element(i).requiresCoordinates() )
				return true;
		}
		return false;
	}

	RotamerStericFilter::RotamerStericFilter( const RotamerApplicatorBase& _app )
		: RotamerFilterBase(_app), 
		m_OLapFactor(0.7),
		m_PtrStaticGrid(&_app.getStaticGrid())
	{
	}

	RotamerStericFilter::~RotamerStericFilter()
	{
		// Nothing!
	}

	RotamerStericFilter* RotamerStericFilter::clone() const
	{
		return new RotamerStericFilter(*this);
	}

	bool RotamerStericFilter::passes( const RotamerLink& _ir, size_t _rotID )
	{
		const WorkSpace& wspace = m_App->getWspace();		
		const ParticleStore& atom = wspace.atom;
		const BondOrder& bondOrder = wspace.bondorder();
		const Residue& res = wspace.res[_ir.resIndex];

		const GridPoint* list = NULL;
		for( int iat = res.ifirst; iat <= res.ilast; iat++ )
		{
			if( atom[iat].isSideChain() )
			{
				if( m_PtrStaticGrid->atomList( atom[iat].pos(), list ) )
				{
					for( int j = 0; j < list->size(); j++ )
					{
						int k = list->atomIndexes[j];
						if(	k < res.ifirst || k > res.ilast ) // We assume that all rotamers are not self-clashing (else they wouldnt be very good rotamers would they :-D)
						{
							if( bondOrder.getBondOrder( iat, k ) > BondOrder_1_4_Pair ) // Only include 1:5 bonded pairs or above
							{
								if( simpleClash(atom[iat],atom[k], m_OLapFactor) ) // Do we clash?
								{
									return false;
								}
							}
						}
					}
				}
			}
		}
		return true;
	}


	// -----------------------------
	// Class: RotamerApplicatorBase
	// -----------------------------

	const double GRID_SEARCH = 5.0;

	RotamerApplicatorBase::RotamerApplicatorBase( WorkSpace& _wspace, const Library::RotamerLibrary& _Lib, RotamerMode _mode )
		: MoveBase( _wspace ), 
		OutputLevel( Verbosity::Normal ),
		m_Lib(&_Lib), 
		m_PeriodicIdealiseInterval(1), 
		m_FiltersEnabled(true), 
		m_Mode(_mode), 
		m_ResPicker(PickEverything()),
		m_StaticGrid( _wspace, Pick_NOT(PickSidechains()), GRID_SEARCH ),
		m_DynamicGrid( _wspace, PickSidechains(), GRID_SEARCH )
	{
		ASSERT(_Lib.isFinalised(), CodeException, "RotamerApplicatorBase: Library::RotamerLibrary::finalise() has not been called on passed library");
		reInitialise( m_ResPicker.data() );
	}

	RotamerApplicatorBase::RotamerApplicatorBase( WorkSpace& _wspace, const Library::RotamerLibrary& _Lib, const PickResidueBase& _Picker, RotamerMode _mode )
		: MoveBase( _wspace ),
		OutputLevel(Verbosity::Normal),
		m_Lib(&_Lib), 
		m_PeriodicIdealiseInterval(1), 
		m_FiltersEnabled(true), 
		m_Mode(_mode), 
		m_ResPicker(_Picker),
		m_StaticGrid( _wspace, Pick_NOT(Pick_AND(_Picker,PickSidechains())), GRID_SEARCH ),
		m_DynamicGrid( _wspace, Pick_AND(_Picker,PickSidechains()), GRID_SEARCH )
	{
		ASSERT(_Lib.isFinalised(), CodeException, "RotamerApplicatorBase: Library::RotamerLibrary::finalise() has not been called on passed library");
		reInitialise( m_ResPicker.data() );
	}

	void RotamerApplicatorBase::clearFilters() 
	{ 
		m_ActivationFilters.clear(); 
	}	

	void RotamerApplicatorBase::deActivateAll()
	{
		for( size_t i = 0; i < m_RotamerLinks.size(); i++ )
		{
			for( size_t j = 0; j < m_RotamerLinks[i].nRot(); j++ )
			{
				m_RotamerLinks[i].setActive(j,false);
			}
		}
	}

	void RotamerApplicatorBase::activateAll()
	{
		for( size_t i = 0; i < m_RotamerLinks.size(); i++ )
		{
			for( size_t j = 0; j < m_RotamerLinks[i].nRot(); j++ )
			{
				m_RotamerLinks[i].setActive(j,true);
			}
		}
	}

	void RotamerApplicatorBase::addFilter_DefaultSteric()
	{
		m_ActivationFilters.addWithOwnership( new RotamerStericFilter(*this) );
	}

	void RotamerApplicatorBase::deactivateViaFilters()
	{
		deactivateViaFilters(m_ActivationFilters);
	}

	void RotamerApplicatorBase::deactivateViaFilters(RotamerFilterContainer& filters)
	{
		if( m_FiltersEnabled && filters.size() > 0 )
		{
			// some filters require that the rotamer be applied, others dont. Therefore if none of the
			// current filter list requires application, there is no point in doing so! Test for this...
			bool applyRotamers = filters.requiresCoordinates();

			// Now, deactivate if we don't like them!			
			for( size_t i = 0; i < m_RotamerLinks.size(); i++ )
			{
				RotamerLink& rotI = m_RotamerLinks[i];
				rotI.sidechainPosCacheA.store();
				for( size_t j = 0; j < rotI.nRot(); j++ )
				{
					if( !rotI.getActive(j) )
						continue;
					if( applyRotamers )
					{
						rotI.apply( j );					
					}
					if( !filters.passes( rotI, j ) )
					{
						rotI.setActive(j, false);
					}
				}
				rotI.sidechainPosCacheA.revert();
			}			
		}
	}

	size_t RotamerApplicatorBase::nRot( size_t ir ) const
	{ 
		return m_RotamerLinks[ir].rot->nRot(); 
	}

	void RotamerApplicatorBase::setPeriodicIdealiseInterval( size_t interval )
	{
		m_PeriodicIdealiseInterval = interval;
		for( size_t i = 0; i < m_RotamerLinks.size(); i++ )
		{
			m_RotamerLinks[i].PeriodicIdealiseInterval = m_PeriodicIdealiseInterval;
		}
	}

	void RotamerApplicatorBase::recalcRotlinkMaxProbability()	
	{
		// Recalculate the max probability, its backbone-dependent!!
		for( size_t i = 0; i < m_RotamerLinks.size(); i++ )
		{
			m_RotamerLinks[i].calcMaxProbability();
		}
	}

	RotEneStoreType RotamerApplicatorBase::atomInteractionEnergy( const ParticleStore& atom, const SnapShot& snap, size_t i, size_t j ) const
	{
		// Trivially count clashes, deriving classes can do a better job.
		// This function is required to be in the base class as stericEnergy() uses it.
		return simpleClash( atom[i], atom[j], 0.7 ) ? 1.0f : 0.0f;
	}

	size_t RotamerApplicatorBase::getStoredStaticEneIndexer( size_t _rotLinkIndex, size_t _slot ) const
	{
		return m_StoreSelfEnergy[_rotLinkIndex][_slot].second;
	}

	const std::vector< std::vector<RotEneStoreType> >& RotamerApplicatorBase::getStoredPairwiseEne( size_t _rotLinkI, size_t _rotLinkK ) const
	{
		ASSERT( _rotLinkI != _rotLinkK, CodeException, "No self pair-interactions allowed!" );		
		if( _rotLinkK < _rotLinkI )
		{
			std::swap( _rotLinkK, _rotLinkI );
		}
		size_t kLookup = _rotLinkK - _rotLinkI - 1;
		return m_StorePairEnergy[ _rotLinkI ][ kLookup ];
	}

	RotEneStoreType RotamerApplicatorBase::getStoredPairwiseEne( size_t _rotLinkI, size_t _rotJ, size_t _rotLinkK, size_t _rotM ) const
	{
		ASSERT( _rotLinkI != _rotLinkK, CodeException, "No self pair-interactions allowed!" );
		if( _rotLinkK < _rotLinkI )
		{
			std::swap( _rotLinkK, _rotLinkI );
			std::swap( _rotJ, _rotM );
		}
		size_t kLookup = _rotLinkK - _rotLinkI - 1;
		return m_StorePairEnergy[ _rotLinkI ][ kLookup ][ _rotJ ][ _rotM ];
	}

	std::pair<RotEneStoreType, size_t> RotamerApplicatorBase::getStoredStaticEne( size_t _rotLinkIndex, size_t _slot ) const 
	{ 
		return m_StoreSelfEnergy[_rotLinkIndex][_slot];
	}

	RotEneStoreType RotamerApplicatorBase::calcEPairwise( size_t _rotLinkI, size_t _rotLinkJ ) const
	{
		ASSERT( _rotLinkI != _rotLinkJ, CodeException, "Rotamer cannot be pair-wise conpared with itself");

		const ParticleStore& atoms = wspace->atom;
		const SnapShot& cur = wspace->cur;
		const Residue& resI = wspace->res[ m_RotamerLinks[_rotLinkI].resIndex ];
		const Residue& resJ = wspace->res[ m_RotamerLinks[_rotLinkJ].resIndex ];

		// Only sidechains make pairwise interactions for the purposes of this...
		// any interaction with the static grid is already included in this residues self energy
	
		RotEneStoreType stericEne = 0.0;
		for( size_t iat = resI.ifirst; iat <= resI.ilast; iat++ )
		{
			if( atoms[iat].isSideChain() ) 
			{
				for( size_t jat = resJ.ifirst; jat <= resJ.ilast; jat++ )
				{
					if( atoms[jat].isSideChain() )
					{
						RotEneStoreType ene = atomInteractionEnergy(atoms, cur, iat, jat );
						stericEne += ene;
					}
				}
			}
		}
		return stericEne;
	}

	RotEneStoreType RotamerApplicatorBase::calcESteric_Static( size_t _rotLinkIndex ) const
	{
		const ParticleStore& atoms = wspace->atom;
		int resIIndex = (int)m_RotamerLinks[_rotLinkIndex].resIndex;
		const Residue& res = wspace->res[ resIIndex ];
		const SnapShot& cur = wspace->cur;

		RotEneStoreType stericEnergy = 0.0;
		const GridPoint* _list;
		for( size_t iat = res.ifirst; iat <= res.ilast; iat++ )
		{
			if( m_StaticGrid.atomList(cur.atom[iat].p,_list) )
			{
				for( size_t j = 0; j < _list->size(); j++ )
				{
					size_t jat = _list->atomIndexes[j];
					int resJIndex = atoms[jat].ir;

					// We assert that i does not interact with i+-1!
					// We are therefore assuming there is always a rotamer which is able not 
					// to clash with directly adjacent resiudes (ref SCWRL paper)
					if( abs(resJIndex-resIIndex) >= 2 )
					{
						RotEneStoreType stericEne = atomInteractionEnergy( atoms, cur, iat, jat );
						stericEnergy += stericEne;
					}
				}
			}
		}

		return stericEnergy;
	}

	void RotamerApplicatorBase::reInitialise( const PickResidueBase& _Picker )
	{
		reInitialise( m_Mode, _Picker );
	}

	void RotamerApplicatorBase::printRotLinkActivity() const
	{
		Printf( "Names: " );
		for( size_t i = 0; i < m_RotamerLinks.size(); i++ )
		{
			Printf("%s ")(m_RotamerLinks[i].rot->getName() );
		}
		Printf( "\n" );
		Printf( "Names: " );
		for( size_t i = 0; i < m_RotamerLinks.size(); i++ )
		{
			Printf("%3d ")(m_RotamerLinks[i].countActive() );
		}
		Printf( "\n\n" );
	}

	void RotamerApplicatorBase::reInitialise( RotamerMode _mode, const PickResidueBase& _Picker )
	{
		ASSERT( wspace != NULL, CodeException, "Internal NULL wspace reference");
		ASSERT( m_Lib != NULL, CodeException, "Internal NULL _Lib reference");

		// Proxies
		ResidueStore& res = wspace->res;
		const Sequence::BioSequence& seq = wspace->getSequence();

		size_t nres = res.size();
		ASSERT( nres == seq.size(), CodeException, "Internal wspace->seq mismatch");

		// Reassign our grids from the new picker
		m_StaticGrid = ProximityGrid( *wspace, Pick_NOT(Pick_AND(_Picker,PickSidechains())), GRID_SEARCH );
		m_DynamicGrid = ProximityGrid( *wspace, Pick_AND(_Picker,PickSidechains()), GRID_SEARCH );

		// Setup our rotamer-link array for speedy access
		m_RotamerLinks.clear(); // Wash any existing underpants
		for( size_t i = 0; i < res.size(); i++ )
		{
			if( !_Picker.matches( wspace->res[i] ) )
			{
				continue;
			}

			std::string name = seq.getResName(i);
			size_t libraryID = m_Lib->getIDForResidue(name);
			if( libraryID == SIZE_T_FAIL )
			{
				// Um, we have no definition!
				throw ProcedureException( "RotamerLibrary does not define residue: '" + name + "'");
			}

			RotamerLink link;
			link.reinit( wspace, m_Lib, m_Mode, m_RotamerLinks.size(), i, libraryID, m_PeriodicIdealiseInterval );
			m_RotamerLinks.push_back(link);
		}
	}

	void RotamerApplicatorBase::applyRotamer( size_t _ResIndex, size_t _RotamerID ) const
	{
		m_RotamerLinks[_ResIndex].apply(_RotamerID);
	}

	void RotamerApplicatorBase::applyMostLikelyNonClashingRotamers(double cutoff) const
	{
		for( size_t i = 0; i < m_RotamerLinks.size(); i++ )
		{
			applyMostLikelyNonClashingRotamer(i,cutoff);
		}
	}

	size_t RotamerApplicatorBase::applyMostLikelyNonClashingRotamer( size_t rotIndex, double cutoff ) const
	{
		D_ASSERT( rotIndex < m_RotamerLinks.size(), OutOfRangeException, "applyMostLikelyNonClashingRotemer(), out of bounds request" );

		const RotamerLink& rotamerLink = m_RotamerLinks[rotIndex];

		size_t nrot = rotamerLink.rot->nRot();
		std::vector< std::pair<double,size_t> > probs;
		probs.reserve(nrot);		
		for( size_t i = 0; i < nrot; i++ )
		{
			if( rotamerLink.getActive(i) )
			{
				probs.push_back( std::pair<double,size_t>( rotamerLink.probability(i), i) );
			}
		}

		if( probs.size() == 0 )
		{
			// There are no active rotamers. 
			// Best way?: Simply use the idealised-geometry forcefield definition.
			rotamerLink.applyFFIdeal();
			return 0;
		}

		if( probs.size() == 1 )
		{
			// We only have one, slam it it, its all we can do!
			rotamerLink.apply( probs[0].second );
			return 1;
		}

		// Sort by probability, **backwards**, hence rbegin and rend. (i.e. reverse itterators)
		std::sort( probs.rbegin(), probs.rend() );

		for( size_t i = 0; i < probs.size(); i++ )
		{
			rotamerLink.apply( probs[i].second );
			if( cutoff >= calcESteric_Static( rotIndex ) )
			{
				return probs.size();
			}
		}

		// They all fecking clash!! 
		// Therefore just shove the most likely in there :-S
		rotamerLink.apply( probs[0].second );
		return probs.size();
	}

	void RotamerLink::reinit( WorkSpace* _wspace, const Library::RotamerLibrary* _Lib, RotamerMode _mode, size_t _serial, int _resIndex, int _rotLibID, int _PeriodicIdealiseInterval )
	{
		m_PreviousRotApply = SIZE_T_FAIL;
		wspace = _wspace;
		mode = _mode;
		serial = _serial;
		PeriodicIdealiseInterval = _PeriodicIdealiseInterval;
		rot = &_Lib->getRotamerSet( _rotLibID );
		resIndex = _resIndex; // Store the residue index
		libraryRotID = _rotLibID; // Store the internal library ID
		Pick_AND picker = Pick_AND( PickResidue( _resIndex ), PickSidechains() );
		sidechainPosCacheA.setPicking( *wspace, picker );
		sidechainPosCacheB.setPicking( *wspace, picker );	

		m_Active.clear();
		m_Active.resize(_Lib->nRotIn( _rotLibID ), true); // true, all rotamers begin active
		m_ActiveCount = m_Active.size();

		rot->calcAtomIndexMap(*wspace, _resIndex, indexMap ); // Get the atom indeces!
		rot->calcTorsionalScope_FF(*wspace->res[_resIndex].param, scopeMap );

		calcMaxInteractionDistance(); // How tall is this residue when its stretched out like a cheap whore?					
		calcMaxProbability();
	}

	void RotamerLink::calcMaxProbability()
	{
		maxProbability = -1.0;
		const size_t size = nRot();
		for( size_t i = 0; i < size; i++ )
		{
			double currentProbability = probability(i);
			if( currentProbability > maxProbability )
			{
				maxProbability = currentProbability;
			}
		}
	}

	void RotamerApplicatorBase::preCalculateSterics_inner( size_t i, size_t j )
	{
		// Calc pair-energy to the current APPLIED rotamer of picked residue i	

		// Calculate ESelf
		m_StoreSelfEnergy[i][j] = std::pair<RotEneStoreType,size_t>(calcESteric_Static(i),j);

		// Proxies
		const SnapShot& cur = wspace->cur;
		RotamerLink& rotI = m_RotamerLinks[i];
		size_t nrotI = rotI.nRot();
		int ancI = rotI.indexMap[rotI.rot->getCartesianAnchorRot1()];

		size_t pairedCount = m_RotamerLinks.size() - i - 1;
		for( size_t k = 0; k < pairedCount; k++ )
		{		
			size_t rotKIndex = k+i+1;
			RotamerLink& rotK = m_RotamerLinks[rotKIndex];
			size_t nrotK = rotK.nRot();
#ifdef _DEBUG
			bool checker = false;
#endif
			double interactionSum = rotI.maxInteractionDistance + rotK.maxInteractionDistance;
			interactionSum *= interactionSum;
			int ancK = rotK.indexMap[rotK.rot->getCartesianAnchorRot1()];
			double sqrDist = cur.atom[ancI].p.sqrdist( cur.atom[ancK].p );

			m_StorePairEnergy[i][k][j].resize(nrotK,RotEneStoreType_MAX);

			if( interactionSum < sqrDist )
			{		
				for( size_t m = 0; m < nrotK; m++ )
				{							
					if( rotK.getActive(m) )
					{
						m_StorePairEnergy[i][k][j][m] = 0.0;
					}
				}						
#ifndef _DEBUG
				continue;							
#else
				checker = true;
#endif
			}
			else
			{						
				for( size_t m = 0; m < nrotK; m++ )
				{							
					if( rotK.getActive(m) )
					{
						rotK.apply( m );
						RotEneStoreType e = calcEPairwise( i, rotKIndex );
#ifdef _DEBUG
						if( checker )
							ASSERT( e == 0.0, CodeException, "Assumption Failure");
#endif
						m_StorePairEnergy[i][k][j][m] = e;
					}
				}
			}
		}
	}

	void RotamerApplicatorBase::preCalculateSterics_inner( size_t i )
	{
		// Proxies
		RotamerLink& rotI = m_RotamerLinks[i];
		size_t nrotI = rotI.nRot();

		// Self
		m_StoreSelfEnergy[i].resize( nrotI, std::pair<RotEneStoreType,size_t>(RotEneStoreType_MAX,SIZE_T_FAIL) );

		// Pair
		size_t pairedCount = m_RotamerLinks.size() - i - 1;
		m_StorePairEnergy[i].resize(pairedCount); // diagonal of square matrix		

		// Alloc arrays
		for( size_t k = 0; k < pairedCount; k++ )
		{
			m_StorePairEnergy[i][k].resize(nrotI);
		}

		for( size_t j = 0; j < nrotI; j++ )
		{
			if( !m_RotamerLinks[i].getActive(j) )
			{
				m_StoreSelfEnergy[i][j] = std::pair<RotEneStoreType,size_t>(RotEneStoreType_MAX,j);
			}
			else
			{
				m_RotamerLinks[i].applyIfChangedID(j);
				preCalculateSterics_inner( i, j );	
			}
		}
	}

	// A copy of preCalculateSterics_inner( size_t i ) with a different 'rotKIndex' loop and indexer range checks
	void RotamerApplicatorBase::preCalculateSterics_PostCall( size_t roti, size_t j )
	{
		// Calc pair-energy to the current APPLIED rotamer of picked residue i	

		// Calculate ESelf
		m_StoreSelfEnergy[roti][j] = std::pair<RotEneStoreType,size_t>(calcESteric_Static(roti),j);

		// Proxies
		const SnapShot& cur = wspace->cur;
		RotamerLink& rotI = m_RotamerLinks[roti];
		size_t nrotI = rotI.nRot();
		int ancI = rotI.indexMap[rotI.rot->getCartesianAnchorRot1()];

		size_t pairedCount = m_RotamerLinks.size() - roti - 1;
		for( size_t rotKIndex = 0; rotKIndex < m_RotamerLinks.size(); rotKIndex++ )
		{		
			if( rotKIndex == roti )
				continue;

			RotamerLink& rotK = m_RotamerLinks[rotKIndex];
			size_t nrotK = rotK.nRot();

			// array indexers
			if( roti < rotKIndex )
			{
				ASSERT( m_StorePairEnergy[roti][rotKIndex-roti-1][j].size() == 0, CodeException, "Assumption Fail");
				m_StorePairEnergy[roti][rotKIndex-roti-1][j].resize(nrotK,RotEneStoreType_MAX);
			}

#ifdef _DEBUG
			bool checker = false;
#endif
			double interactionSum = rotI.maxInteractionDistance + rotK.maxInteractionDistance;
			interactionSum *= interactionSum;
			int ancK = rotK.indexMap[rotK.rot->getCartesianAnchorRot1()];
			double sqrDist = cur.atom[ancI].p.sqrdist( cur.atom[ancK].p );	

			if( interactionSum < sqrDist )
			{		
				for( size_t m = 0; m < nrotK; m++ )
				{							
					if( rotK.getActive(m) )
					{
						if( roti < rotKIndex )
						{
							m_StorePairEnergy[roti][rotKIndex-roti-1][j][m] = 0.0;
						}
						else
						{
							ASSERT( m_StorePairEnergy[rotKIndex][roti-rotKIndex-1][m].size() != 0, CodeException, "Assumption failure" );
							//m_StorePairEnergy[rotKIndex][roti-rotKIndex][m].resize(nrotI,RotEneStoreType_MAX);
							m_StorePairEnergy[rotKIndex][roti-rotKIndex-1][m][j] = 0.0;
						}
					}
				}						
#ifndef _DEBUG
				continue;							
#else
				checker = true;
#endif
			}
			else
			{						
				for( size_t m = 0; m < nrotK; m++ )
				{							
					if( rotK.getActive(m) )
					{
						rotK.apply( m );
						RotEneStoreType e = calcEPairwise( roti, rotKIndex );
#ifdef _DEBUG
						if( checker )
							ASSERT( e == 0.0, CodeException, "Assumption Failure");
#endif
						if( roti < rotKIndex )
						{
							m_StorePairEnergy[roti][rotKIndex-roti-1][j][m] = e;
						}
						else
						{
							ASSERT( m_StorePairEnergy[rotKIndex][roti-rotKIndex-1][m].size() != 0, CodeException, "Assumption failure" );
							//m_StorePairEnergy[rotKIndex][roti-rotKIndex][m].resize(nrotI,RotEneStoreType_MAX);
							m_StorePairEnergy[rotKIndex][roti-rotKIndex-1][m][j] = e;							
						}
					}
				}
			}
		}
	}

	void RotamerApplicatorBase::preCalculateSterics()
	{
		// Init arrays
		m_StoreSelfEnergy.clear();
		m_StoreSelfEnergy.resize( m_RotamerLinks.size() );
		m_StorePairEnergy.clear();
		m_StorePairEnergy.resize( m_RotamerLinks.size() );

		for( size_t i = 0; i < m_RotamerLinks.size(); i++ )
		{
			preCalculateSterics_inner(i);
		}
	}

	void RotamerLink::applyFFIdeal() const
	{
		rot->applyFFIdealised(*wspace,resIndex,indexMap);
	}

	void RotamerLink::applyIfChangedID( size_t rotIndex ) const
	{
		if( rotIndex != m_PreviousRotApply )
		{
			apply( rotIndex );
		}
	}

	void RotamerLink::apply( size_t rotIndex ) const
	{
		ASSERT( getActive(rotIndex), ArgumentException, "Inactive rotamer requested for application");

		// Flag this!
		m_PreviousRotApply = rotIndex;

		switch( mode )
		{
		case ApplyCartesian:
			{
				rot->applyRotamerCartesian(*wspace,resIndex,rotIndex,indexMap);
				return;
			}
		case ApplyTorsional:
			{
				rot->applyRotamerTorsional(*wspace,resIndex,rotIndex,scopeMap);
				return;
			}
		case ApplyTorsional_Idealise:
			{
				rot->applyFFIdealised(*wspace,resIndex,indexMap);
				rot->applyRotamerTorsional(*wspace,resIndex,rotIndex,scopeMap);
				return;
			}
		case ApplyTorsional_PeriodicIdealise:
			{
				static int switcher = 0;
				if( switcher % PeriodicIdealiseInterval == 0 )
					rot->applyFFIdealised(*wspace,resIndex,indexMap);
				rot->applyRotamerTorsional(*wspace,resIndex,rotIndex,scopeMap);
				switcher++;
				return;
			}
		}
	}

	double RotamerLink::probability( size_t rotIndex ) const
	{
		return rot->getRotamer(rotIndex).getProbMap().probability(*wspace,resIndex);
	}

	void RotamerLink::setActive( size_t rotIndex, bool activityChange )
	{
		if( m_Active[rotIndex] != activityChange )
		{
			// its changed, flag the change
			m_Active[rotIndex] = activityChange;
			if( activityChange )
			{
				m_ActiveCount++; // one more active
			}
			else
			{
				m_ActiveCount--; // on less is active
			}
		}
	}

	bool RotamerLink::getActive( size_t rotIndex ) const
	{
		return m_Active[rotIndex];
	}

	size_t RotamerLink::countActive() const
	{
		return m_ActiveCount;
	}

	void RotamerLink::calcMaxInteractionDistance()
	{
		size_t natom = rot->nAtom();
		size_t nrot = rot->nRot();

		maxInteractionDistance = 0.0;
		for( size_t i = 0; i < nrot; i++ )
		{
			const Library::Rotamer& rotamer = rot->getRotamer(i);
			Maths::dvector anchorZeroPos( rotamer.getPos( rot->getCartesianAnchorRot1() ) );
			for( size_t j = 0; j < natom; j++ )
			{				
				maxInteractionDistance = std::max( maxInteractionDistance,
					rotamer.getPos(j).sqrdist( anchorZeroPos ) );
			}
		}
		maxInteractionDistance = sqrt( maxInteractionDistance ); // sqrt(): we stored a square distance, and want the distance.
	}

	// -------------------------------
	// Class: RandomRotamerApplicator
	// -------------------------------

	RandomRotamerApplicator::RandomRotamerApplicator( WorkSpace& _wspace, const Library::RotamerLibrary& _Lib, RotamerMode _mode, bool _ProbabilityWeighting )
		: RotamerApplicatorBase(_wspace,_Lib,_mode)
	{
		ProbabilityWeighting = _ProbabilityWeighting;
		rand = Maths::FastRandom::getInstance();
	}

	RandomRotamerApplicator::RandomRotamerApplicator( WorkSpace& _wspace, const Library::RotamerLibrary& _Lib, const PickResidueBase& _Picker, RotamerMode _mode, bool _ProbabilityWeighting )
		: RotamerApplicatorBase(_wspace,_Lib,_Picker,_mode)
	{
		ProbabilityWeighting = _ProbabilityWeighting;
		rand = Maths::FastRandom::getInstance();
	}

	int RandomRotamerApplicator::apply()
	{ 
		activateAll(); // All default to active
		recalcRotlinkMaxProbability(); // This needs to be called - its backbone dependent!!
		deactivateViaFilters(); // set which rotamer states are active
		for( size_t i = 0; i < m_RotamerLinks.size(); i++ )
		{
			applyRandomRotamer(i);
		}
		return 1;
	} 

	void RandomRotamerApplicator::applyRandomRotamer( size_t _PickedIndex )
	{
		RotamerLink& link = m_RotamerLinks[_PickedIndex];
		const Library::RotamerSet& rotSet = *link.rot;
		const double rotFrac = 1.0 / (double)rotSet.nRot();

		// First find the maximum probability based on active rotamer states
		double activeProbSum = 0.0;
		for( size_t i = 0; i < rotSet.nRot(); i++ )
		{
			if( link.getActive(i) ) // Only include active rotamers
			{
				double prob = ProbabilityWeighting ? link.probability(i) : rotFrac;
				activeProbSum += prob;
			}
		}
		if( activeProbSum == 0.0 ) // If there are no active states, then there is no rotamer to set...
			return;

		double probSum = 0.0;
		// nextDouble(): 'Values returned are from 0.0 up to but not including 1.0' - which is what we want here :-D
		double probCap = rand->nextDouble() * activeProbSum; // The maximum probability sum given the current active states
		if( ProbabilityWeighting )
		{
			for( size_t i = 0; i < rotSet.nRot(); i++ )
			{
				if( link.getActive(i) ) 
				{
					probSum += link.probability(i);
					if( probSum > probCap )
					{
						// We have our rotamer
						link.apply(i);
						return;
					}
				}
			}
			THROW(CodeException,"applyRandomRotamer() probability weighting code failed?!?");
		}
		else
		{
			// We are going to apply the Xth active rotamer:
			size_t Xth = (size_t) floor(probCap / rotFrac);			
			for( size_t i = 0; i < rotSet.nRot(); i++ )
			{
				if( link.getActive(i) ) 
				{
					if( Xth <= 0 ) 
					{
						// We have our rotamer
						link.apply(i);
						return;
					}
					Xth--;
				}
			}
			THROW(CodeException,"applyRandomRotamer() non-probability weighting code failed?!?");
		}
	}
}


