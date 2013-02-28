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
#include "tools/vector.h"
#include "maths/graphtheory.h"

namespace Maths
{
	void Vertex::info( Verbosity::Type verbose ) const
	{
		if( verbose > Verbosity::Silent )
		{
			Printf("%3d %c %d {")(serial)( articulationOrder ? 'A' : '-' )(articulationOrder);
			for( size_t i = 0; i < adjacency.size(); i++ )
			{
				Printf("%d")(adjacency[i]);
				if( i != adjacency.size()-1 )
					Printf(",");
			}
			Printf("}");
			if( verbose > Verbosity::Normal )
				Printf(" (Num:%3d, Lpt:%3d)")(number)(lowpt);
		}
	}

	Edge::Edge()
	{
	}

	Edge::Edge( size_t _fromVertex, size_t _toVertex )
		: fromVertex(_fromVertex), toVertex(_toVertex)
	{
	}

	void Edge::info( Verbosity::Type verbose ) const
	{
		if( verbose > Verbosity::Silent )
			Printf("{%-3d,%-3d}")(fromVertex)(toVertex);
	}

	void BiconnectedComponent::info( Verbosity::Type verbose ) const
	{
		if( verbose > Verbosity::Silent )
		{
			Printf("%d: ")(serial);
			for( size_t i = 0; i < edges.size(); i++ )
			{
				edges[i].info(verbose);
				if( i != edges.size()-1 )
					Printf("%,");
			}
		}
	}

	GraphBase::GraphBase()
	{
	}

	GraphBase::GraphBase(size_t _vertexCount)
	{
		setVertexCount( _vertexCount );
	}

	void GraphBase::setVertexCount(size_t _vertexCount)
	{
		clear();
		m_Vertices.resize(_vertexCount);
		for( size_t i = 0; i < m_Vertices.size(); i++ )
		{
			m_Vertices[i].serial = i;
		}
	}

	void GraphBase::addEdge( Edge _e )
	{
		ASSERT( m_Vertices.size() > _e.fromVertex && m_Vertices.size() > _e.toVertex,
			ArgumentException, "Edge is out of vertex range");
		if( !VectorContains( m_Edges, _e ) ) // It is important toVertex ensure that all edges are unique
		{
			_e.serial = m_Edges.size();
			m_Edges.push_back(_e);
		}
	}

	void GraphBase::addEdge( size_t fromVertex, size_t toVertex )
	{
		addEdge( Edge( fromVertex, toVertex ) );
	}

	void GraphBase::clear()
	{
		m_Edges.clear();
		m_Vertices.clear();
	}

	void UndirectedGraph::calcAdjacency()
	{
		for( size_t i = 0; i < m_Vertices.size(); i++ )
		{
			m_Vertices[i].adjacency.clear();
		}
		for( size_t i = 0; i < m_Edges.size(); i++ )
		{
			m_Vertices[m_Edges[i].fromVertex].adjacency.push_back( m_Edges[i].toVertex );
			m_Vertices[m_Edges[i].toVertex].adjacency.push_back( m_Edges[i].fromVertex );
		}
	}

	UndirectedGraph::UndirectedGraph()
	{
	}

	UndirectedGraph::UndirectedGraph(size_t _vertexCount)
		: GraphBase(_vertexCount )
	{
	}

	void UndirectedGraph::info(Verbosity::Type verbose) const
	{
		if( verbose >= Verbosity::Normal )
		{
			Printf("Vertices:\n");
			for( size_t i = 0; i < m_Vertices.size(); i++ )
			{
				m_Vertices[i].info(verbose);
				Printf("\n");
			}
			Printf("\nEdges:\n");
			int sep = 0;
			for( size_t i = 0; i < m_Edges.size(); i++ )
			{
				sep++;
				m_Edges[i].info(verbose);
				if( sep % 6 == 0 )
					Printf("\n");
				else if( i != m_Edges.size()-1 )
					Printf(" ");
			}
			Printf("\n\nBiconnected Components:\n");
			for( size_t i = 0; i < m_BiComponents.size(); i++ )
			{
				m_BiComponents[i].info(verbose);
				Printf("\n");
			}
			Printf("\n");
		}
	}

	bool UndirectedGraph::hasEdge( size_t i, size_t j ) const
	{
		return hasEdge( Edge( i, j ) );
	}

	bool UndirectedGraph::hasEdge( Edge e ) const
	{
		// In undirected graphs, edges are bi-directional
		if( VectorContains( m_Edges, e ) )
			return true;
		std::swap( e.fromVertex, e.toVertex );
		return VectorContains( m_Edges, e );
	}

	void UndirectedGraph::calcArticulationPoints()
	{
		biConnectInit();
		m_I = 0;
		m_EdgeStack.clear();
		for( size_t i = 0; i < m_Vertices.size(); i++ )
		{
			if( m_Vertices[i].number == SIZE_T_FAIL )
			{
				biConnect( i, 0 );
			}
		}
		calcArticulationOrders();
	}

	void UndirectedGraph::calcArticulationOrders()
	{
		for( size_t i = 0; i < m_Vertices.size(); i++ )
		{
			m_Vertices[i].articulationOrder = 0; // not an articulation point, so zero order
			for( size_t j = 0; j < m_BiComponents.size(); j++ )
			{
				for( size_t k = 0; k < m_BiComponents[j].edges.size(); k++ )
				{
					if( m_BiComponents[j].edges[k].fromVertex == i ||
						m_BiComponents[j].edges[k].toVertex   == i )
					{
						m_Vertices[i].articulationOrder++;
						break;
					}
				}
			}
			if( m_Vertices[i].articulationOrder > 0 )
				m_Vertices[i].articulationOrder--; // by definition of order
		}
	}

	void UndirectedGraph::biConnectInit()
	{
		for( size_t i = 0; i < m_Vertices.size(); i++ )
		{
			m_Vertices[i].number = SIZE_T_FAIL; // Uninitialised
			m_Vertices[i].lowpt = SIZE_T_FAIL; // Uninitialised
		}
		calcAdjacency();
		m_BiComponents.clear();
	}

	void UndirectedGraph::biConnect( size_t _v, size_t _u )
	{
		Vertex& v = m_Vertices[_v];
		v.number = ++m_I;
		v.lowpt = v.number;

		for( size_t i = 0; i < v.adjacency.size(); i++ )
		{
			size_t _w = v.adjacency[i];
			Vertex& w = m_Vertices[_w];

			if( w.number == SIZE_T_FAIL ) // vertex w is not yet assigned
			{
				m_EdgeStack.push_back( Edge( _v, _w ) ); // add (v,w) to the stack of edges

				biConnect( _w, _v );
				v.lowpt = std::min( v.lowpt, w.lowpt );

				if( w.lowpt >= v.number )
				{
					BiconnectedComponent bc;
					bc.serial = m_BiComponents.size();
					while( true )
					{
						Edge& e = m_EdgeStack.back();
						if( m_Vertices[e.fromVertex].number >= w.number )
						{
							bc.edges.push_back( e );
							m_EdgeStack.pop_back();
						}
						else
						{
							break;
						}
					}

					Edge& e = m_EdgeStack.back();
					bc.edges.push_back( e );
					m_EdgeStack.pop_back();

					m_BiComponents.push_back( bc );
				}

			}
			else if( w.number < v.number && _w != _u )
			{
				m_EdgeStack.push_back( Edge( _v, _w ) );
				v.lowpt = std::min( v.lowpt, w.number );
			}
		}
	}
}

