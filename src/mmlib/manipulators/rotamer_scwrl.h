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

#ifndef __SCWRL_H
#define __SCWRL_H

#include "rotamer_applicatorbase.h" // Base class
#include "maths/graphtheory.h"

class RotamerCumulativeProbabilityDensityFilter;

namespace Manipulator
{	
	class RotamerApplicator_SCWRL; // Forward declaration for SQWRL_SolveUnit
	struct SQWRL_SolveUnit;

	struct SQWRL_Vertex : public Maths::Vertex
	{
		SQWRL_Vertex( RotamerLink& _rotLink, Maths::UndirectedGraph& _Graph, std::vector<SQWRL_Vertex>& vectors, size_t vertexIndex );
		std::vector< size_t > childIndeces; ///< Only articulation points actually have these
		std::vector< std::vector< size_t > > bestRotamers; ///< Store the best rotamer for each child for each rotamer of this residue
		std::vector< double > superRotamerEnergies;
		RotamerLink* rotLink;
		std::vector<SQWRL_Vertex>* m_VertexContainer;
		void apply( size_t _rotIndex, Verbosity::Type _verb );
	};

	struct SQWRL_SolveUnit : public Maths::BiconnectedComponent
	{
		SQWRL_SolveUnit( Maths::UndirectedGraph& _Graph, size_t biConnectedIndex, std::vector<SQWRL_Vertex>& vectors );		
		size_t minArticulationOrder() const;
		size_t countArticulationPoints() const;
		size_t firstArticulationPoint() const;
		std::vector<size_t> m_Vertices; // indexes in m_VertexContainer
		std::vector<SQWRL_Vertex>* m_VertexContainer;
		bool used;
	};

	struct VertexLink
	{
		VertexLink( SQWRL_Vertex& _parent );
		bool operator<( const VertexLink& v ) const;
		inline size_t nActiveRot() const { return rotLink->countActive(); }
		inline size_t serial() const { return parent->serial; }
		RotamerLink* rotLink;
		SQWRL_Vertex* parent;
	};

	class PD_API RotamerApplicator_SCWRL : public RotamerApplicatorBase
	{
	public:
		RotamerApplicator_SCWRL( WorkSpace& _wspace, const Library::RotamerLibrary& _Lib, const PickResidueBase& _Picker, RotamerMode _mode = ApplyCartesian );
		RotamerApplicator_SCWRL( WorkSpace& _wspace, const Library::RotamerLibrary& _Lib, RotamerMode _mode = ApplyCartesian );
		virtual ~RotamerApplicator_SCWRL(){}
		virtual RotamerApplicator_SCWRL* clone() const;

		virtual int apply(); ///< Pack rotamers

		bool SupressStaticGridRefreshOnApply; ///< Do this **ONLY** if you **KNOW** the static grid isnt moving!!

		bool IgnoreHydrogenForSterics; ///< Ignore hydrogens in atom interaction calcs (SCWRL paper, yes, us, no)
		double EBBMAX; ///< 50.0 is the SCWRL EBBMAX cutoff for interaction
		double PROB_CAP; ///< 90% probability filter

	protected:

		//void DEBUG_PAIR_ENE( size_t _resI, size_t _resJ, size_t _rotK, size_t _rotM, double theEne );

		void finalisePreCalculatedSterics();
		void deadEndElimination(); /// Remove rotamers which cannot by definition be part of the global SCRWL energy minimim
		size_t detectRotamerContacts(); /// Fill the undirected graph from the current residue contacts
		bool reactivateSingleBestRotamerIfNoneAreValid( RotamerCumulativeProbabilityDensityFilter& sqrwlFilter ); ///< If all valid rotamer states are sterically screened, reactivate the best failing possibility... its as good as we are going to get!
		void resolveGraph(); ///< If residue interactions are found, resolve the Undirected Graph! :-D
		void resolveGraph( std::vector<SQWRL_SolveUnit>& units, size_t i );

		RotEneStoreType calcESelf( size_t _residue, size_t rotamer ) const;
		RotEneStoreType calcEProbability( size_t _residue, size_t rotamer ) const;
		virtual RotEneStoreType atomInteractionEnergy( const ParticleStore& atom, const SnapShot& snap, size_t i, size_t j ) const; // override base class function

		Maths::UndirectedGraph m_Graph;
		std::vector<bool> m_FlagActiveResidue; ///< Used to store which picked residues are "active"

	private:
		double deadEndElimination_Core(RotamerLink& rot, size_t _rotID, std::pair<RotEneStoreType, size_t> Si, std::pair<RotEneStoreType, size_t> Ri ); /// inner call
		size_t traverseRotStack( std::vector< VertexLink >& vertexLinks );
	};
}

#endif

