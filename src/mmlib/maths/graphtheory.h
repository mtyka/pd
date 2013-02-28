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

#ifndef _GRAPH_THEORY_H
#define _GRAPH_THEORY_H

// Implementation (c) Jon Rea 2007
/// \file maths/graphtheory.h 
/// \brief Contains classes related to graph theory
/// \details Graph theory algorithms
/// \author Jon Rea

#include <vector>
#include <stack>

namespace Maths
{
	struct PD_API Vertex
	{
		size_t serial; ///< Index in parent Vertex array
		size_t articulationOrder; ///< The number of connections - 1
		size_t number; ///< Part of the biconnect algorithm
		size_t lowpt; ///< low point - Part of the biconnect algorithm
		std::vector<size_t> adjacency; ///< The indexes of the vertexes connected to this one
		void info( Verbosity::Type verbose ) const;
	};

	struct PD_API Edge
	{
		Edge();
		Edge( size_t from, size_t to );
		bool operator==( const Edge& e ) const { return e.fromVertex == fromVertex && e.toVertex == toVertex; }
		size_t serial;
		size_t fromVertex;
		size_t toVertex;
		void info( Verbosity::Type verbose ) const;
	};

	struct PD_API BiconnectedComponent
	{
		size_t serial;
		std::vector<Edge> edges;
		void info( Verbosity::Type verbose ) const;
	};

	class PD_API GraphBase
	{
	public:
		GraphBase();
		GraphBase(size_t _vertexCount);

		void addEdge( Edge _e );
		void addEdge( size_t from, size_t to );

		void setVertexCount(size_t _vertexCount);
		void clear();

		const std::vector<Edge>& getEdges() const { return m_Edges; }
		const std::vector<Vertex>& getVertices() const { return m_Vertices; }

	protected:

		std::vector<Vertex> m_Vertices;
		std::vector<Edge> m_Edges;
	};

	/// \brief A connected undirected graph
	class PD_API UndirectedGraph : public GraphBase
	{
	public:
		UndirectedGraph();
		UndirectedGraph(size_t _vertexCount);

		bool hasEdge( size_t i, size_t j ) const;
		bool hasEdge( Edge e ) const;

		/// \brief Algorithm for finding the articulation points for an undirected graph.
		/// \details Algorithm by Tarjan et al:
		/// \reference
		/// article{Tarjan72,
  	/// author    = {Robert Endre Tarjan},
  	/// title     = {Depth-First Search and Linear Graph Algorithms.},
  	/// journal   = {SIAM J. Comput.},
  	/// volume    = {1},
  	/// number    = {2},
  	/// year      = {1972},
  	/// pages     = {146-160},
		/// }
		void calcArticulationPoints();

		void info( Verbosity::Type verbose = Verbosity::Normal ) const;

		const std::vector<BiconnectedComponent>& getBiComponents() const { return m_BiComponents; }

	protected:
		
		void calcAdjacency();

	private:
		void calcArticulationOrders();
		void biConnectInit();
		void biConnect( size_t v, size_t u );
		size_t m_I;
		std::vector<Edge> m_EdgeStack;
		std::vector<BiconnectedComponent> m_BiComponents;
	};
}

#endif

