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
#include "tools/vector.h"
#include "system/fundamentals.h"
#include "system/molecule.h"
#include "workspace/workspace.h"
#include "library/rotamerlib.h"
#include "rotamer_filter.h"
#include "rotamer_scwrl.h"

namespace Manipulator
{
	const double DEFAULT_EBBMAX = 5.0; // 50.0 is the SCWRL original paper EBBMAX cutoff for interaction
	const double DEFAULT_PROB_CAP = 0.9; // 90% is the SCWRL original paper probability cutoff for rotamer inclusion

	size_t NthActiveBackLookup( RotamerLink& rotLink, size_t compRotamerSerial )
	{
		ASSERT( rotLink.getActive( compRotamerSerial ) && compRotamerSerial < rotLink.nRot(), CodeException, "Assumption failure" );
		size_t index = 0;
		for( size_t i = 0; i < compRotamerSerial; i++ )
		{
			if( rotLink.getActive( i ) )
			{
				index++;
			}
		}
		return index;
	}

	RotamerApplicator_SCWRL::RotamerApplicator_SCWRL( WorkSpace& _wspace, const Library::RotamerLibrary& _Lib, RotamerMode _mode )
		: RotamerApplicatorBase(_wspace,_Lib,_mode)
	{
		EBBMAX = DEFAULT_EBBMAX;
		IgnoreHydrogenForSterics = false;
		PROB_CAP = DEFAULT_PROB_CAP;
		SupressStaticGridRefreshOnApply = false;
	}

	RotamerApplicator_SCWRL::RotamerApplicator_SCWRL( WorkSpace& _wspace, const Library::RotamerLibrary& _Lib, const PickResidueBase& _Picker, RotamerMode _mode )
		: RotamerApplicatorBase(_wspace,_Lib,_Picker,_mode)
	{
		EBBMAX = DEFAULT_EBBMAX;
		IgnoreHydrogenForSterics = false;
		PROB_CAP = DEFAULT_PROB_CAP;
		SupressStaticGridRefreshOnApply = false;
	}

	RotamerApplicator_SCWRL* RotamerApplicator_SCWRL::clone() const
	{
		return new RotamerApplicator_SCWRL(*this);
	}

	size_t RotamerApplicator_SCWRL::detectRotamerContacts()
	{
		// NOTE!!
		// This function makes the assumption that all active and inactive picked resdues are
		// currently using their most likely rotamers.

		const ParticleStore& atom = wspace->atom;
		const ResidueStore& res = wspace->res;
		const SnapShot& cur = wspace->cur;

		// Initialise m_FlagActiveResidue to flag all picked residues as inactive
		m_FlagActiveResidue.clear();
		m_FlagActiveResidue.resize( m_RotamerLinks.size(), false );

		size_t countActive = 0;

		// Detect which should begin as active
		for( size_t i = 0; i < m_RotamerLinks.size() - 1; i++ )
		{
			for( size_t j = i+1; j < m_RotamerLinks.size(); j++ )
			{
				double ePair = calcEPairwise( i, j );
				if( ePair > 0.0 )
				{
					if( m_RotamerLinks[i].countActive() > 1 && m_RotamerLinks[j].countActive() > 1 )
					{
						if( !m_FlagActiveResidue[i] )
						{
							m_FlagActiveResidue[i] = true;
							countActive++;
						}
						if( !m_FlagActiveResidue[j] )
						{
							m_FlagActiveResidue[j] = true;
							countActive++;
						}
						m_Graph.addEdge( i, j );
					}
					else if( m_RotamerLinks[i].countActive() > 1 )
					{
						if( !m_FlagActiveResidue[i] )
						{
							m_FlagActiveResidue[i] = true;
							countActive++;
						}
					}
					else if( m_RotamerLinks[j].countActive() > 1 )
					{
						if( !m_FlagActiveResidue[j] )
						{
							m_FlagActiveResidue[j] = true;
							countActive++;
						}
					}
					else
					{
						// There is 1 or less rotamers for both, but we have a clash.
						// There is nothing we can do about this with the packer...
						// You would just need a more fine-grained library, or the
						// backbone structure iteself is flawed.
						if( OutputLevel > Verbosity::Silent )
						{
							Printf("Warning! (Picked indexes %d vs. %d) Clash detected (%6.2f), but there are <= 1 active rotamers for both residues! (i:%d)(j:%d)\n")
								(i)(j)(ePair)(m_RotamerLinks[i].countActive())(m_RotamerLinks[j].countActive());
						}
					}
				}
			}
		}

		if( countActive == 0 )
		{
			// None are active, all are in their most probable non-clashing
			// states, so the system must by definition already be in its minimum energy conformation.
			// We got 'nout to do!
			return 0;
		}

		std::vector<Maths::Edge> testedPossibleEdges;

		// We now need to see if any additional residues should be marked as active!
		// This occurs when any one rotamer of one residue interacts with a rotamer from another residue.
		bool activations;
		do
		{
			activations = false;
			for( size_t i = 0; i < m_RotamerLinks.size(); i++ )
			{
				if( m_FlagActiveResidue[i] )
				{
					const RotamerLink& rotI = m_RotamerLinks[i];
					D_ASSERT( rotI.countActive() > 1, CodeException, "Residue should not be active here if it has <= 1 active rotamers!");
					int ancI = rotI.indexMap[rotI.rot->getCartesianAnchorRot1()];

					for( size_t j = 0; j < m_RotamerLinks.size(); j++ )
					{
						if( i == j )
							continue; // you cant interact with yourself
						const RotamerLink& rotJ = m_RotamerLinks[j];

						const std::vector< std::vector<RotEneStoreType> >& ePairs = getStoredPairwiseEne( i, j );

						Maths::Edge e(i, j);
						// Computationally cheap, so test first, prior to interactionSum which requires
						// slightly more expensive pointer lookups for evaluation.
						if( VectorContains(testedPossibleEdges, e ) )
							continue; // we know they interact, why check again?
						testedPossibleEdges.push_back( e );
						std::swap( e.fromVertex, e.toVertex );
						testedPossibleEdges.push_back( e ); // also add the reverse

						// Is it actually possible for these two residues to interact?
						double interactionSum = rotI.maxInteractionDistance + rotJ.maxInteractionDistance;
						interactionSum *= interactionSum;
						int ancJ = rotJ.indexMap[rotJ.rot->getCartesianAnchorRot1()];
						double sqrDist = cur.atom[ancI].p.sqrdist( cur.atom[ancJ].p );
						if( interactionSum < sqrDist )
						{
							// Assert that there is indeed no interaction when in debug!
							D_ASSERT( calcEPairwise( i, j ) == 0.0, CodeException, "Critical code assumption failure!");
							continue; // They cant possibly interact, so continue :-D
						}

						for( size_t k = 0; k < rotI.rot->nRot(); k++ )
						{
							if( !rotI.getActive(k) )
								continue; // this rotamer or residue I is not active

							if( rotJ.countActive() > 0 ) // otherwise, it cant be activated, cos it only has one!
							{
								for( size_t m = 0; m < rotJ.rot->nRot(); m++ )
								{
									if( !rotJ.getActive(m) )
										continue; // this rotamer or residue I is not active

									RotEneStoreType ePair;
									if( i < j )
										ePair = ePairs[k][m];
									else
										ePair = ePairs[m][k];

									//DEBUG_PAIR_ENE( i, j, k, m, ePair );

									ASSERT( ePair != RotEneStoreType_MAX, CodeException, "Cached value has not been stored yet!");
									if( ePair > 0.0 )
									{
										// you can call addEdge() even if we have added this pair before, there is an internal check.
										m_Graph.addEdge( i, j);

										// here we know, by definition, that m_FlagActiveResidue[i] is true
										if( !m_FlagActiveResidue[j] && m_RotamerLinks[j].countActive() > 1 )
										{
											m_FlagActiveResidue[j] = true;
											activations = true; // we found a new one! Flag this.
											countActive++; // mark that a further residue is active
										}
										// no point in checking any more rotamer pairs, we now know these residues interact.
										goto BREAK_OUTER;
									}
								}
							}
						}
BREAK_OUTER:
						continue;
					}
				}
			}
		}
		while( activations );

		return countActive;
	}

	void RotamerApplicator_SCWRL::deadEndElimination()
	{
		// Against all other residues, we are going to identify rotamers in residue '_resID'
		// that are higher in energy than all other rotamers in residue '_resID' for all other
		// residues and their rotamers.
		// To do this we will use the Goldstein Criterion!! Wooooo, that sounds fancy!
		// See the SCQRL paper or google
		// 'An extension of dead-end elimination for protein side-chain conformation using merge-decoupling'
		// for an explanation of the simplest 'deadEndElimination' (DEE) which is the one used here
		// Note, below we are using sorted lists to increase the efficiency of the algorithm
		// Hence the backward enumeration below.

		double eliminated = 0.0;
		double totalRotamers = 0.0;

		for( size_t i = 0; i < m_RotamerLinks.size(); i++ )
		{
			RotamerLink& rot = m_RotamerLinks[i];

			size_t activeCnt = rot.countActive();
			totalRotamers += activeCnt;
			if( activeCnt <= 1 )
			{
				// There cant by definition be any to eliminate
				continue;
			}

			for( int j = ((int)rot.nRot())-1; j >= 0; j-- ) // from highest to lowest
			{
				std::pair<RotEneStoreType, size_t> selfE_j = getStoredStaticEne( i, j );
				if( rot.getActive( selfE_j.second ) )
				{
					ASSERT( selfE_j.first != RotEneStoreType_MAX, CodeException, "Unassigned ESelf");

					for( int k = ((int)rot.nRot())-1; k >= 0; k-- ) // from highest to lowest
					{
						if( j == k )
							continue;
						std::pair<RotEneStoreType, size_t> selfE_k = getStoredStaticEne( i, k );

						if( rot.getActive( selfE_k.second ) )
						{
							ASSERT( selfE_k.first != RotEneStoreType_MAX, CodeException, "Unassigned ESelf");

							double decisionSum = deadEndElimination_Core( rot, i, selfE_j, selfE_k );
							if( decisionSum > 0.0 )
							{
								if( OutputLevel >= Verbosity::Loud )
									Printf("ELIMINATION!!, %d: %d vs %d (%8.3lf vs. %8.3lf) (%8.3lf)\n")(i)(selfE_j.second)(selfE_k.second)(selfE_j.first)(selfE_k.first)(decisionSum);
								// Woooo; One more eliminated from the overall search!
								rot.setActive( selfE_j.second, false );
								eliminated += 1.0;
								break;
							}
							else
							{
								if( OutputLevel >= Verbosity::Loud )
									Printf("No Elimination, %d: %d vs %d (%8.3lf vs. %8.3lf) (%8.3lf)\n")(i)(selfE_j.second)(selfE_k.second)(selfE_j.first)(selfE_k.first)(decisionSum);
							}
						}
					}
				}
			}
		}

		if( OutputLevel )
			Printf("Dead-end Eliminated %5.2lf%% of rotamers! (%d/%d)\n\n")((eliminated/totalRotamers)*100.0)((int)eliminated)((int)totalRotamers);
	}

	// void RotamerApplicator_SCWRL::DEBUG_PAIR_ENE( size_t _resI, size_t _resJ, size_t _rotK, size_t _rotM, double theEne )
	// {
	//#ifdef _DEBUG
//		m_RotamerLinks[_resI].forceApply( _rotK );
//		m_RotamerLinks[_resJ].forceApply( _rotM );
//		double realENE = calcEPairwise( _resI, _resJ );
//		ASSERT( realENE == theEne, CodeException, "DEBUG_PAIR_ENE");
	//#endif
	// }

	double RotamerApplicator_SCWRL::deadEndElimination_Core(RotamerLink& rotI, size_t _rotID, std::pair<RotEneStoreType, size_t> Si, std::pair<RotEneStoreType, size_t> Ri )
	{
		RotEneStoreType pairMinimisedSum = 0.0f;
		for( size_t j = 0; j < m_RotamerLinks.size(); j++ )
		{
			if( j != _rotID )
			{
				RotamerLink& rotJ = m_RotamerLinks[j];

				RotEneStoreType minWRTSi = RotEneStoreType_MAX;
				RotEneStoreType minWRTRi = RotEneStoreType_MAX;
				const std::vector< std::vector<RotEneStoreType> >& ePairwise = getStoredPairwiseEne( _rotID, j );

				if( rotJ.countActive() == 0 )
				{
					// Then finalisePreCalculatedSterics() will have set cache position k=0 to be the pairwise
					// interactions of the rotamers with the FFParam conformation...

					RotEneStoreType eneSiJu;
					RotEneStoreType eneRiJu;

					// The call to getStoredPairwiseEne above only stores 1/2 the matrix!
					if( _rotID < j )
					{
						eneSiJu = ePairwise[Si.second][0];
						eneRiJu = ePairwise[Ri.second][0];
					}
					else
					{
						eneSiJu = ePairwise[0][Si.second];
						eneRiJu = ePairwise[0][Ri.second];
					}

					ASSERT( eneSiJu != RotEneStoreType_MAX, CodeException, "Cached value has not been calculated");
					minWRTSi = std::min( minWRTSi, eneSiJu );
					ASSERT( eneRiJu != RotEneStoreType_MAX, CodeException, "Cached value has not been calculated");
					minWRTRi = std::min( minWRTRi, eneRiJu );
				}
				else
				{
					for( int k = 0; k < rotJ.nRot(); k++ )
					{
						if( rotJ.getActive(k) )
						{
							RotEneStoreType eneSiJu;
							RotEneStoreType eneRiJu;

							// The call to getStoredPairwiseEne above only stores 1/2 the matrix!
							if( _rotID < j )
							{
								eneSiJu = ePairwise[Si.second][k];
								eneRiJu = ePairwise[Ri.second][k];
							}
							else
							{
								eneSiJu = ePairwise[k][Si.second];
								eneRiJu = ePairwise[k][Ri.second];
							}

							ASSERT( eneSiJu != RotEneStoreType_MAX, CodeException, "Cached value has not been calculated");
							minWRTSi = std::min( minWRTSi, eneSiJu );
							ASSERT( eneRiJu != RotEneStoreType_MAX, CodeException, "Cached value has not been calculated");
							minWRTRi = std::min( minWRTRi, eneRiJu );

							//DEBUG_PAIR_ENE( _rotID, j, Si.second, k, eneSiJu );
							//DEBUG_PAIR_ENE( _rotID, j, Ri.second, k, eneRiJu );
						}
					}
				}

				RotEneStoreType minimisedDelta = minWRTSi - minWRTRi;

				pairMinimisedSum += minimisedDelta;
			}
		}

		RotEneStoreType selfDelta = Si.first - Ri.first;
		RotEneStoreType decisionSum = selfDelta + pairMinimisedSum;

		return decisionSum; // If its greater than zero, Si cannot be part of the global energy minimim and can therefore be culled!!
	}

	void RotamerApplicator_SCWRL::finalisePreCalculatedSterics()
	{
		// Add the probabilities to the steric scores, to create ESelf arrays
		for( size_t i = 0; i < m_RotamerLinks.size(); i++ )
		{
			if( m_RotamerLinks[i].countActive() > 0 )
			{
				size_t nrot = m_RotamerLinks[i].rot->nRot();
				for( size_t j = 0; j < nrot; j++ )
				{
					if( m_RotamerLinks[i].getActive(j) )
					{
						D_ASSERT( m_StoreSelfEnergy[i][j].second == j, CodeException, "Why is the list not in order at this point?");
						m_StoreSelfEnergy[i][j].first += calcEProbability( i, j );
					}
				}
			}
		}

		ASSERT( m_RotamerLinks.size() == m_StoreSelfEnergy.size(), CodeException, "" );

		// Apply the best ones. They dont clash with the static grid (unless there are no non-clashers)
		// But they may, of course, interact with each other!
		for( size_t i = 0; i < m_RotamerLinks.size(); i++ )
		{
			if( m_RotamerLinks[i].countActive() == 0 )
			{
				m_RotamerLinks[i].applyFFIdeal();
				preCalculateSterics_PostCall( i, 0 );
			}
			else
			{
				// Sort forwards by stored energy - the smaller the better
				std::sort( m_StoreSelfEnergy[i].begin(), m_StoreSelfEnergy[i].end() );
				for( size_t j = 0; j < m_StoreSelfEnergy.size(); j++ )
				{
					size_t rotIndex = m_StoreSelfEnergy[i][j].second;
					if( m_RotamerLinks[i].getActive(rotIndex) )
					{
						m_RotamerLinks[i].apply( rotIndex );
						break;
					}
					continue;
				}
			}
		}
	}

	class SCWRLInitRotamerFilter : public RotamerFilterBase
	{
	public:
		SCWRLInitRotamerFilter( const RotamerApplicator_SCWRL& parent, double eneCutoff )
			: RotamerFilterBase(parent), m_EneCutoff(eneCutoff), m_SCWRL(&parent) {}
		virtual ~SCWRLInitRotamerFilter(){}
		virtual SCWRLInitRotamerFilter* clone() const { return new SCWRLInitRotamerFilter(*this); }
		virtual bool passes( const RotamerLink& _ir, size_t _rotID )
		{
			std::pair<double, size_t> ene = m_SCWRL->getStoredStaticEne( _ir.serial, _rotID );
			D_ASSERT( ene.second == _rotID, CodeException, "Code Assumption Error");
			return ene.first <= m_EneCutoff;
		}
		virtual bool requiresCoordinates() const { return true; }
	private:
		const RotamerApplicator_SCWRL* m_SCWRL;
		double m_EneCutoff;
	};

	bool RotamerApplicator_SCWRL::reactivateSingleBestRotamerIfNoneAreValid( RotamerCumulativeProbabilityDensityFilter& probFilter )
	{
		ASSERT(!probFilter.requiresCoordinates(), CodeException,
			"Last time I checked RotamerCumulativeProbabilityDensityFilter didn't require rotamer application. What happened!?!");
		bool needApply = m_ActivationFilters.requiresCoordinates();

		std::vector< std::pair<double,size_t> > probs;

		bool hadUndefined = false;
		size_t undefined = 0;
		for( size_t pass = 0; pass <= 1; pass ++ )
		{
			if( pass == 1 )
				needApply = true;

			for( size_t i = 0; i < m_RotamerLinks.size(); i++ )
			{
				RotamerLink& rotamerLink = m_RotamerLinks[i];

				if( rotamerLink.countActive() == 0 )
				{
					undefined++;
					hadUndefined = true;

					size_t nrot = rotamerLink.rot->nRot();
					probs.clear();
					probs.reserve(nrot);
					for( size_t j = 0; j < nrot; j++ )
					{
						std::pair<double,size_t> selfEne = getStoredStaticEne( i, j );

						if( needApply )
						{
							ASSERT( !rotamerLink.getActive(j), CodeException, "Unexpected active rotamer !?!");
							rotamerLink.setActive(selfEne.second,true);
							rotamerLink.apply(selfEne.second);
							rotamerLink.setActive(selfEne.second,false);
						}

						if( pass==1 )
						{
							if( m_ActivationFilters.passes( rotamerLink, selfEne.second ) )
							{
								if( selfEne.first == RotEneStoreType_MAX )
								{
									// This rotamer was previously pre-filtered, and therefore is not in the cache!
									// We therefore need to calc it here!
									ASSERT( selfEne.first == RotEneStoreType_MAX, CodeException,
										"Steric pre-calculation pass 0 should have worked!" );
									// The rotamer i has been applied, needApply is set true for pass-1
									preCalculateSterics_PostCall( i, selfEne.second );
									selfEne = getStoredStaticEne( i, selfEne.second );
									ASSERT( selfEne.first != RotEneStoreType_MAX, CodeException,
										"Steric pre-calculation pass 1 failure!" );
									ASSERT(rotamerLink.probability(selfEne.second) <= 0.0, CodeException, "Probability assumption failure");
									selfEne.first += calcEProbability( i, j ); // we know the probability is <= 0.0
								}
								probs.push_back( selfEne );
							}
						}
						else
						{
							if( probFilter.passes( rotamerLink, selfEne.second ) &&
								m_ActivationFilters.passes( rotamerLink, selfEne.second ) )
							{
								ASSERT( selfEne.first != RotEneStoreType_MAX, CodeException,
									"Expected that if the filters passed this rotamer, it would have been analysed by the steric pre-calc!! Something is wrong!");
								selfEne.first += calcEProbability( i, selfEne.second ); // take the probability into account
								probs.push_back( selfEne );
							}
						}
					}

					if( probs.size() == 0 )
					{
						// There are no rotamers which pass even the basic filter list!!!
						if( OutputLevel >= Verbosity::Quiet )
							Printf("WARNING: No valid rotamer states for residue %d (pick index %d) (pass %d)\n")(rotamerLink.resIndex)(i)(pass);
						if( pass == 1 )
						{
							if( OutputLevel >= Verbosity::Quiet )
								Printf("CRITICAL WARNING: Even disregarding probability AND sterics, after user-filtering there are NO rotamers. What are you upto?!?!\n");
						}
						else
						{
							if( OutputLevel >= Verbosity::Normal )
								Printf("This is after re-filtering by the *basic* non-steric filters! Library culling criteria are too great!\n");
							if( OutputLevel >= Verbosity::Normal )
								Printf("On the other hand, something may be dodgy... (one example is if Proline is in a unusual backbone conformation with no associated-probability info)\n");
						}
					}
					else
					{
						// Sort by probability
						std::sort( probs.begin(), probs.end() );

						if( OutputLevel > Verbosity::Quiet )
						{
							Printf("WARNING: Forcing rotamer re-activation due to ");
							if( pass == 0 ) Printf("poor steric environment");
							else Printf("no valid probability states");
							Printf(" for residue %d (pick index %d)!\n")(rotamerLink.resIndex)(i);
						}

						// DEBUG
						//wspace->outtra.append();
						// END DEBUG

						// Apply the best
						rotamerLink.setActive( probs[0].second, true ); // a necessary special reactivation!
						rotamerLink.apply( probs[0].second );
						undefined--;

						// DEBUG
						//wspace->outtra.append();
						//for( size_t j = 0; j < nrot; j++ )
						//{
						//	rotamerLink.setActive( j, true ); // a necessary special reactivation!
						//	rotamerLink.apply( j );
						//	wspace->outtra.append();
						//}
						// END DEBUG
					}
				}
			}
			if( undefined == 0 )
				break;
		}

		if( undefined > 0 )
		{
			if( OutputLevel > Verbosity::Quiet )
			{
				Printf("CRITICAL WARNING: User filters result in one or more residue has no rotamers at all!\n");
			}
		}

		return hadUndefined;
	}

	VertexLink::VertexLink( SQWRL_Vertex& _parent )
	{
		parent = &_parent;
		rotLink = _parent.rotLink;
	}

	bool VertexLink::operator<( const VertexLink& v ) const
	{
		return nActiveRot() < v.nActiveRot();
	}

	SQWRL_Vertex::SQWRL_Vertex( RotamerLink& _rotLink, Maths::UndirectedGraph& _Graph, std::vector<SQWRL_Vertex>& vectors, size_t vertexIndex )
	{
		const Maths::Vertex& vert = _Graph.getVertices()[vertexIndex];
		adjacency = vert.adjacency;
		articulationOrder = vert.articulationOrder;
		ASSERT( vertexIndex == vert.serial, CodeException, "Assumption failure");
		serial = vertexIndex;
		m_VertexContainer = &vectors;
		rotLink = &_rotLink;
	}

	void SQWRL_Vertex::apply( size_t _rotIndex, Verbosity::Type _verb )
	{
		if( _verb > Verbosity::Silent )
		{
			Printf("Graph resolution for residue %d: Chosen rotamer %d\n")(serial)(_rotIndex);
		}
		rotLink->apply(_rotIndex);
		if( bestRotamers.size() > 0 )
		{
			// We have Chhheeeellldren - wha ha ha ha!
			size_t seek = NthActiveBackLookup( *rotLink, _rotIndex );
			const std::vector<size_t> bestIndexes = bestRotamers[seek];
			ASSERT( childIndeces.size() == bestIndexes.size(), CodeException, "Assumption failure");
			for( size_t i = 0; i < childIndeces.size(); i++ )
			{
				SQWRL_Vertex& childVert = (*m_VertexContainer)[childIndeces[i]];
				childVert.apply( bestIndexes[i], _verb );
			}
		}
	}

	inline bool hasVertex( std::vector<Maths::Edge> edges, const SQWRL_Vertex& v )
	{
		for( size_t i = 0; i < edges.size(); i++ )
		{
			if( v.serial == edges[i].fromVertex || v.serial == edges[i].toVertex )
			{
				return true;
			}
		}
		return false;
	}

	SQWRL_SolveUnit::SQWRL_SolveUnit( Maths::UndirectedGraph& _Graph, size_t biConnectedIndex, std::vector<SQWRL_Vertex>& vectors )
	{
		const BiconnectedComponent& bi = _Graph.getBiComponents()[biConnectedIndex];
		ASSERT( bi.serial == biConnectedIndex, CodeException, "Assumption failure");
		serial = biConnectedIndex;
		edges = bi.edges;
		for( size_t i = 0; i < vectors.size(); i++ )
		{
			ASSERT( i == vectors[i].serial, CodeException, "Assumption failure");
			if( hasVertex( edges, vectors[i] ) )
			{
				m_Vertices.push_back( i );
			}
		}
		m_VertexContainer = &vectors;
		used = false;
	}

	size_t SQWRL_SolveUnit::minArticulationOrder() const
	{
		size_t min = SIZE_T_FAIL;
		for( size_t i = 0; i < m_Vertices.size(); i++ )
		{
			if( (*m_VertexContainer)[m_Vertices[i]].articulationOrder < min )
			{
				min = (*m_VertexContainer)[m_Vertices[i]].articulationOrder;
			}
		}
		return min;
	}

	size_t SQWRL_SolveUnit::countArticulationPoints() const
	{
		size_t count = 0;
		for( size_t i = 0; i < m_Vertices.size(); i++ )
		{
			if( (*m_VertexContainer)[m_Vertices[i]].articulationOrder > 0 )
			{
				count++;
			}
		}
		return count;
	}

	size_t SQWRL_SolveUnit::firstArticulationPoint() const
	{
		for( size_t i = 0; i < m_Vertices.size(); i++ )
		{
			if( (*m_VertexContainer)[m_Vertices[i]].articulationOrder > 0 )
			{
				return i;
			}
		}
		return SIZE_T_FAIL;
	}

	size_t RotamerApplicator_SCWRL::traverseRotStack( std::vector< VertexLink >& vertexLinks )
	{
		// ---------------------------------
		// Phase 1, configure primary vertex
		// and obtain proxies
		// ---------------------------------

		ASSERT( vertexLinks.size() > 1, CodeException, "Assumption Failure");

		SQWRL_Vertex& primary = *vertexLinks[0].parent;
		if( primary.articulationOrder > 0 )
			primary.articulationOrder--;
		std::vector< std::vector< size_t > >& bestRotamers2 = primary.bestRotamers;
		std::vector< double >& superRotamerEnergies2 = primary.superRotamerEnergies;
		RotamerLink& primaryRotLink = *primary.rotLink;

		std::vector< size_t >& childIndeces2 = primary.childIndeces;
		for( size_t i = 1; i < vertexLinks.size(); i++ )
		{
			childIndeces2.push_back( vertexLinks[i].serial() );
		}

		// can be 0 or n
		ASSERT( bestRotamers2.size() == superRotamerEnergies2.size(), CodeException, "Assumption Failure");
		if( bestRotamers2.size() == 0 ) bestRotamers2.resize( primaryRotLink.countActive() );

		bool creatingSuperEneArray = superRotamerEnergies2.size() == 0;

		// ---------------------------------
		// Phase 2, Enumerate the primary vertexes rotamers, finding the best solution
		// amongst all children for each one. Record only the best, and save it for this
		// vertex - this vertex is now a super-vertex.
		// ---------------------------------

		double bestGlobalEne = DBL_MAX;
		size_t bestPrimaryRotamer = SIZE_T_FAIL;
		size_t superRotamerEnergies2_Use = 0;
		for( size_t rotamerSlot = 0; rotamerSlot < primaryRotLink.nRot(); rotamerSlot++ )
		{
			std::vector<size_t> bestRotStack;

			std::pair<RotEneStoreType, size_t> ene = getStoredStaticEne( primary.serial, rotamerSlot );
			size_t rotamerID = ene.second;

			if( primaryRotLink.getActive( rotamerID ) )
			{
				ASSERT( ene.first != RotEneStoreType_MAX, CodeException, "Unassigned ESelf");

				// ---------------------------------
				// Phase 2a: Only active rotatmerpstates for
				// the primary vertex are assessed
				// ---------------------------------

				double lowestEne = DBL_MAX;
				std::vector<size_t> appliedRotStack;
				std::vector<double> eneRotStack;

				// The energy of the 1st level is, by definition, either the self energy of a
				// simple vertex, or the super-energy of a super-vertex.
				// The first energy in the eneRotStack holds this value ...
				if( creatingSuperEneArray )
				{
					// A Simple Vertex, becoming a super-vertex - its the primary here, and has children...
					superRotamerEnergies2.push_back( 0.0 );
					eneRotStack.push_back( ene.first );
				}
				else
				{
					// Super-vertex
					eneRotStack.push_back( superRotamerEnergies2[superRotamerEnergies2_Use] );
				}

				// ----------------------------------
				// Phase 2b: Find the best childrens
				// rotamers
				// ----------------------------------

				std::stack<size_t> indexStack;
				size_t currentIndexer = 0;
				while( true )
				{
					size_t pos = appliedRotStack.size()+1;

					if( pos == vertexLinks.size() || currentIndexer == vertexLinks[pos].rotLink->nRot() )
					{
						if( appliedRotStack.size() == 0 )
						{
							break;
						}
						appliedRotStack.pop_back();
						eneRotStack.pop_back();
						currentIndexer = indexStack.top()+1;
						indexStack.pop();
						continue;
					}

					ASSERT( vertexLinks[pos].rotLink->countActive() > 0, CodeException, "If we are here, all should have active states!");

					// See if this rot is active!
					std::pair<RotEneStoreType, size_t> levelSelf =
						getStoredStaticEne( vertexLinks[pos].parent->serial, currentIndexer );
					size_t compRotamerSerial = levelSelf.second;
					if( !vertexLinks[pos].rotLink->getActive( compRotamerSerial ) )
					{
						currentIndexer++;
						continue;
					}
					ASSERT( levelSelf.first != RotEneStoreType_MAX, CodeException, "Lookup Failure!");

					// NOW, we need to add the self-energy of this levels vertex / super vertex
					// lookup the self-ene ...
					double levelEne = eneRotStack[eneRotStack.size()-1];
					if( vertexLinks[pos].parent->superRotamerEnergies.size() > 0 )
					{
						// Its a super-vertex. Get the stored ESelf + children ESelf + all pair enes!
						double superEne = vertexLinks[pos].parent->superRotamerEnergies[ NthActiveBackLookup( *vertexLinks[pos].rotLink, compRotamerSerial ) ];
						levelEne += superEne;
					}
					else
					{
						// Its just a normal vertex - get the single self energy
						levelEne += levelSelf.first;
					}

					// Add the pair interaction with the primary-vertex
					double pairEne = getStoredPairwiseEne(
						vertexLinks[pos].serial(),
						compRotamerSerial,
						vertexLinks[0].serial(),
						rotamerID );
					ASSERT( pairEne != DBL_MAX, CodeException, "Pair energy non precalculated!?!");
					levelEne += pairEne;

					// And then the pair-wise energy with all lower levels (not including the primary-vertex!)
					for( size_t i = 1; i < pos; i++ )
					{
						if( levelEne >= lowestEne )
						{
							goto DIE_PAIR_ENE_COMP_OVER_BEST;
						}
						pairEne = getStoredPairwiseEne(
							vertexLinks[pos].serial(),
							compRotamerSerial,
							vertexLinks[i].serial(),
							appliedRotStack[i-1] );
						ASSERT( pairEne != DBL_MAX, CodeException, "Pair energy non precalculated!?!");
						levelEne += pairEne;
					}

					if( levelEne >= lowestEne )
					{
DIE_PAIR_ENE_COMP_OVER_BEST:
						currentIndexer++;
						continue;
					}

					appliedRotStack.push_back( compRotamerSerial );
					eneRotStack.push_back( levelEne );
					indexStack.push(currentIndexer);
					currentIndexer = 0;

					if( appliedRotStack.size()+1 == vertexLinks.size() )
					{
						// If we have got here, then its the best so far...
						ASSERT( lowestEne > eneRotStack[eneRotStack.size()-1], CodeException, "Assumption failure");
						ASSERT( eneRotStack.size() > 0, CodeException, "Assumption failure");
						lowestEne = eneRotStack[eneRotStack.size()-1];
						bestRotStack = appliedRotStack;
						if( lowestEne < bestGlobalEne )
						{
							bestGlobalEne = lowestEne;
							bestPrimaryRotamer = rotamerID;
						}
					}
				}

				if( vertexLinks.size()-1 != bestRotStack.size() )
					THROW( CodeException, "Stack-traversal logic failure");

				// We have identified the best rotamers for this primary rotamer :-D
				//1) add the lowest sum-energy to the super-rotamer
				superRotamerEnergies2[superRotamerEnergies2_Use] += lowestEne ; // store the energy of this super-rotamer!

				//2) Add the rotamers of the children to the bestRotamer set for this parent rotamer
				bestRotamers2[superRotamerEnergies2_Use].insert(
					bestRotamers2[superRotamerEnergies2_Use].end(), bestRotStack.begin(), bestRotStack.end() );
				if( childIndeces2.size() != bestRotamers2[superRotamerEnergies2_Use].size() )
					THROW( CodeException, "Stack-traversal logic failure");

				// 3) Increment the parent current rotamer indexer to point to the next slot
				superRotamerEnergies2_Use++;
			}
		}

		// We should have set these in the loop above!! Otherwise something most dubious has happened!
		ASSERT( bestGlobalEne != DBL_MAX && bestPrimaryRotamer != SIZE_T_FAIL, CodeException, "Assumption Failure");

		return bestPrimaryRotamer;
	}

	void RotamerApplicator_SCWRL::resolveGraph( std::vector<SQWRL_SolveUnit>& units, size_t i )
	{
		SQWRL_SolveUnit& unit = units[i];

		size_t primaryVertex = unit.firstArticulationPoint();
		if( primaryVertex == SIZE_T_FAIL )
		{
			// If there is no articulation point, then our work is almost complete.
			// We have condensed to a single articulation unit. Resolve for that using an
			// arbritrary vertex as the primary...
			ASSERT( unit.countArticulationPoints() == 0, CodeException, "Resolution failure");
		}

		// Sort the application residues by the number of rotamers they have!
		std::vector< VertexLink > vertexLinks;
		for( size_t i = 0; i < unit.m_Vertices.size(); i++ )
		{
			if( i != primaryVertex )
			{
				VertexLink vert = VertexLink((*unit.m_VertexContainer)[unit.m_Vertices[i]]);
				vertexLinks.push_back( vert );
			}
		}
		std::sort(vertexLinks.begin(),vertexLinks.end());
		if( primaryVertex != SIZE_T_FAIL )
		{
			vertexLinks.insert(vertexLinks.begin(), VertexLink((*unit.m_VertexContainer)[unit.m_Vertices[primaryVertex]]));
		}

		// With respect to the first vertex, fill its best rotamer arrays!
		size_t primariesBestRot = traverseRotStack(vertexLinks);

		if( primaryVertex == SIZE_T_FAIL )
		{
			// Apply the rotamer for the primary vertex. If its is a super-vertex,
			// it will also apply its children...
			vertexLinks[0].parent->apply(primariesBestRot, OutputLevel);
		}
	}

	void RotamerApplicator_SCWRL::resolveGraph()
	{
		std::vector<SQWRL_Vertex> vectors;
		size_t count = m_Graph.getVertices().size();
		for( size_t i = 0; i < count; i++ )
		{
			vectors.push_back( SQWRL_Vertex(
				m_RotamerLinks[ m_Graph.getVertices()[i].serial ],
				m_Graph, vectors, i) );
		}

		std::vector<SQWRL_SolveUnit> components;
		count = m_Graph.getBiComponents().size();
		for( size_t i = 0; i < count; i++ )
		{
			components.push_back( SQWRL_SolveUnit(m_Graph, i, vectors) );
		}

		while( true )
		{
			size_t find = SIZE_T_FAIL;
			for( size_t i = 0; i < components.size(); i++ )
			{
				if( !components[i].used )
				{
					size_t numArticulationPoints = components[i].countArticulationPoints();
					if( numArticulationPoints <= 1 )
					{
						find = i;
						components[i].used = true;
						break;
					}
				}
			}
			if( find == SIZE_T_FAIL )
				break;
			else
				resolveGraph( components, find );
		}
		for( size_t i = 0; i < components.size(); i++ )
		{
			ASSERT( components[i].used == true && components[i].countArticulationPoints() == 0, CodeException, "Resolution failure");
		}
	}

	int RotamerApplicator_SCWRL::apply()
	{
		if( !SupressStaticGridRefreshOnApply )
		{
			m_StaticGrid.refresh();
		}
		else
		{
			D_ASSERT( !m_StaticGrid.requiresRefresh(), CodeException, "SupressStaticGridRefreshOnApply: Flagged true for a system that is moving!!!");
		}

		activateAll(); // All default to active

		if( OutputLevel )
		{
			Printf("All states active:\n");
			printRotLinkActivity();
		}

		recalcRotlinkMaxProbability();

		// Keep only the first >90% most probable rotamers
		RotamerFilterContainer sqrwlFilterA;
		RotamerCumulativeProbabilityDensityFilter probFilter( *this, PROB_CAP );
		sqrwlFilterA.add(probFilter);
		deactivateViaFilters( sqrwlFilterA );

		if( OutputLevel > Verbosity::Silent )
		{
			Printf("Following %5.1lf%% probability filtering:\n")(probFilter.getMinDensity()*100.0);
			printRotLinkActivity();
		}

		// User-defined filters to deactivate further rotamers
		deactivateViaFilters(); // Set which rotamers are active per-residue based on user criteria

		if( OutputLevel > Verbosity::Silent )
		{
			Printf("Following user rotamer filters:\n");
			printRotLinkActivity();
		}

		if( OutputLevel > Verbosity::Quiet )
		{
			Printf("Pre-calculating self and pair-wise energetic terms...\n");
		}

		// **!! CODE NOTE !!**
		// All SCWRL functions after preCalculateSterics(); below have the in-built assumption that
		// NO ROTAMERS WILL EVER BE RE-ACTIVATED!!!
		// Once a rotamer is deactivated, it **must** stay that way, otherwise the code below **will**
		// fail. There are no ASSERT statements for this, as it's hard to test for.

		// These will then be filled in **the current structural context** by SCWRLInitRotamerFilter
		preCalculateSterics(); // Base class call. Il est tres critique! Used by SCWRLInitRotamerFilter below...

		if( OutputLevel > Verbosity::Quiet )
		{
			Printf("Filtering by static-grid energy...\n");
		}

		// Activate rotamers using the standard filter set, and a custom filter to remove
		// any static grid clashing rotamers via the standard SCWRL criterion
		RotamerFilterContainer sqrwlFilterB;
		SCWRLInitRotamerFilter innerFilter(*this, EBBMAX);
		sqrwlFilterB.add( innerFilter );
		deactivateViaFilters( sqrwlFilterB );

		if( OutputLevel > Verbosity::Silent )
		{
			Printf("Following internal steric rotamer filters:\n");
			printRotLinkActivity();
		}

		// Sometimes at this stage, we have NO valid rotamers for a given residue.
		// Most commonly because none are under the 'SCWRLInitRotamerFilter' BBMAX filter...
		// In which case, we have two options.
		// 1) apply the forcefield defaults which could clash horibly ... ew!
		// 2) The second is to apply the rotamer with the best ESelf with no regard for other rotamers,
		//    BUT importantly with regard to the user and CumulativeProbabilityDensity filters!
		//    The remainder of the code here then effectively treats this rotamer with the same logic as other
		//    stuff thats not allowed to move. The user will be warned that reactivation has occured!
		// Following the re-activation attempt, if there are STILL no active rotamers we have aproblem.
		// We must ALWAYS apply the user filters...
		// Therefore either:
		// 1) Repeat the above but ignorning the probability filter - what else can we do!
		// 2) If there are NO rotamers at all after user filtering, then apply the FFPS conformation.
		//    if 2) then the preCalculated pairwise steric arrays are filled for the FFPS conformation as though
		//    it was rotamer 0 during finalisePreCalculatedSterics()
		if( reactivateSingleBestRotamerIfNoneAreValid(probFilter)
			&& OutputLevel > Verbosity::Silent )
		{
			Printf("Following dynamic-reactivation rotamer filters:\n");
			printRotLinkActivity();
		}

		// Add the probabilities to the steric scores, to convert the internal cache to an 'ESelf' array
		// Then sort by energy; then apply.
		// This makes sure each rotamer is in the most likely valid state following this call...
		finalisePreCalculatedSterics();

		if( OutputLevel > Verbosity::Quiet )
		{
			Printf("Performing Dead-end Elimination (DEE)...\n");
		}

		// Remove any rotamers that we are sure are not part of the global energy minimum using DEE
		deadEndElimination();
		// We still need to apply it! Slam in the most likely after the DEE...
		applyMostLikelyNonClashingRotamers(0.5); // allowing micro-clashes

		if( OutputLevel > Verbosity::Silent )
		{
			Printf("Following DEE:\n");
			printRotLinkActivity();
		}

		if( OutputLevel > Verbosity::Quiet )
		{
			Printf("Detecting rotamer contacts...\n");
		}

		// Initialise the internal graph
		m_Graph.clear();
		m_Graph.setVertexCount( m_RotamerLinks.size() );
		// detectRotamerContacts() -> **filling m_Graph** and m_FlagActiveResidue with the contact residue data!
		if( 0 == detectRotamerContacts() )
		{
			if( OutputLevel > Verbosity::Quiet )
			{
				Printf("No rotamer interactions, most likely singletons have been selected!\n");
			}
			// There are no contacts - we are already in the minimum-"energy" conformation!
			return 1; // return 1 as we have applied the most likely non-clashers..
		}

		if( OutputLevel > Verbosity::Quiet )
		{
			Printf("Processing interaction graph...\n");
		}

		// We know at this point that every residue is in its lowest energy conformation, except those that
		// have the potential to interact. Call resolveGraph() to allocate these...
		m_Graph.calcArticulationPoints();
		m_Graph.info(OutputLevel);
		resolveGraph(); // We have contacts, they must be resolved!

		if( OutputLevel > Verbosity::Quiet )
		{
			Printf("Graph resolved!\n");
		}

		return 1;
	}

	double getRadius( int atomicNumber )
	{
		switch( atomicNumber )
		{
		case 6: // 'C' - Carbon
			return 1.6;
		case 8: // 'O' - Oxygen
			return 1.3;
		case 7: // 'N' - Nitrogen
			return 1.3;
		case 16: // 'S' - Sulphur
			return 1.7;
		case 1: // 'H' - *not* SQWRL derived - my definition, **soft** hydrogen
			return 0.45;
		case 30: // 'Zn' - *not* SQWRL derived - Zinc, my definition
			return 1.4;
		default:
			THROW(ProcedureException,std::string("SCWRL does not define radius for element atomic number: ")
				+ int2str(atomicNumber) );
		}
	}

	RotEneStoreType RotamerApplicator_SCWRL::atomInteractionEnergy( const ParticleStore& atom, const SnapShot& snap, size_t i, size_t j ) const
	{
		if( IgnoreHydrogenForSterics )
		{
			if( atom[i].isHydrogen() || atom[j].isHydrogen() )
			{
				return 0.0;
			}
		}

		// The 'random' decimal values below are very much from the SQRWL paper.
		// I have no idea where the absolute values are derived from;
		// but, what this funciton does is obvious.
		double sqrDist = snap.atom[i].p.sqrdist( snap.atom[j].p );
		double rad1 = getRadius( atom[i].Z );
		double rad2 = getRadius( atom[j].Z );
		double radSum = rad1 + rad2;
		double sqrRadSum = radSum * radSum;
		if( sqrRadSum < sqrDist )
		{
			return 0.0f;
		}
		double dist = sqrt(sqrDist);
		if ( dist < 0.8254 * radSum )
		{
			return 10.0f; // max_cap
		}
		else
		{
			return (RotEneStoreType)(57.273 * ( 1.0 - ( dist / radSum ) ));
		}
	}

	RotEneStoreType RotamerApplicator_SCWRL::calcEProbability( size_t _residue, size_t rotamer ) const
	{
		const RotamerLink& rotLink = m_RotamerLinks[_residue];
		// The energyetic term derived from the probability of this rotamer
		double probability = rotLink.probability(rotamer);
		const double minusK = -3.0; // Set by SCWRL in the original paper
		if( probability <= 0.0 )
		{
			return 1000000.0; // flag as crap
		}
		ASSERT(
			rotLink.maxProbability > 0.0f &&
			probability <= rotLink.maxProbability, CodeException, "Probabilities are non-sensical");
		return (RotEneStoreType)(minusK * log( probability / rotLink.maxProbability ));
	}

	RotEneStoreType RotamerApplicator_SCWRL::calcESelf( size_t _rotLinkIndex, size_t rotamer ) const
	{
		m_RotamerLinks[_rotLinkIndex].apply( rotamer );
		RotEneStoreType probabilityEnergy = calcEProbability( _rotLinkIndex, rotamer );
		RotEneStoreType stericEnergy = calcESteric_Static( _rotLinkIndex );
		return probabilityEnergy + stericEnergy; // The final energy is the sum of those components
	}
}

